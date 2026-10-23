"""Cutting the viral rows out of IntAct's full release.

WHAT THIS IS FOR. data/raw/vhppi held seven vhn_<taxid>.tab27 files, one per
coronavirus, and intact.zip at 1.35 GB which nothing had ever read --
ingest_vhppi globbed `vhn_*.tab27` and the pattern did not match it. The
flavivirus PPI gap was recorded as a missing download while a broader source
sat in the same directory.
"""
from __future__ import annotations

import io
import sys
import zipfile
from importlib import import_module
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

fi = import_module("filter_intact")

HEADER = b"#ID(s) interactor A\tID(s) interactor B\n"


def _row(a="P11111", b="Q22222", tax_a="9606", tax_b="64320", cols=42):
    f = [""] * cols
    f[0], f[1] = f"uniprotkb:{a}", f"uniprotkb:{b}"
    f[6] = 'psi-mi:"MI:0018"(two hybrid)'
    f[8] = "pubmed:12345678"
    f[9], f[10] = f"taxid:{tax_a}(human)", f"taxid:{tax_b}(virus)"
    f[11] = 'psi-mi:"MI:0915"(physical association)'
    f[14] = "intact-miscore:0.56"
    return ("\t".join(f) + "\n").encode()


def _zip(tmp: Path, lines: list[bytes], name="intact.txt") -> Path:
    tmp.mkdir(parents=True, exist_ok=True)
    p = tmp / "intact.zip"
    with zipfile.ZipFile(p, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(name, b"".join(lines))
    return p


def _run(tmp: Path, lines: list[bytes], *extra, name="intact.txt"):
    z = _zip(tmp, lines, name=name)
    out = tmp / "intact_viral.tab27"
    code = fi.main(["--zip", str(z), "--out", str(out), *extra])
    return code, out


# ------------------------------------------------------------ the filter
def test_a_declared_viral_taxid_is_kept(tmp_path, capsys):
    code, out = _run(tmp_path, [HEADER, _row(tax_b="64320")])
    assert code == 0 and out.exists()
    assert out.read_bytes().count(b"\n") == 1


def test_an_undeclared_taxid_is_not_kept_for_looking_viral(tmp_path, capsys):
    """Matching by name, or by "not 9606", would pull in every virus IntAct
    carries. The register's own map is the filter."""
    code, out = _run(tmp_path, [HEADER, _row(tax_b="11676")])   # HIV-1
    assert code == 1, "nothing matched, so nothing is written"
    assert not out.exists()


def test_an_isolate_taxid_is_kept_because_the_register_declares_it(tmp_path):
    """declared_taxa folds in isolate_taxa and chain_source.taxon. DENV-2's
    chains live on strain 31634, and IntAct files Shah et al.'s dengue
    interactions under that strain, not the species."""
    cfg = yaml.safe_load((ROOT / "config/viruses.yaml").read_text())
    from kgav.virus_register import declared_taxa
    assert declared_taxa(cfg).get("31634") == "11060", \
        "if this mapping goes, the dengue rows stop resolving"
    code, out = _run(tmp_path, [HEADER, _row(tax_b="31634")])
    assert code == 0 and out.read_bytes().count(b"\n") == 1


def test_the_header_line_is_counted_as_a_header_not_a_row(tmp_path, capsys):
    """The taxid filter would drop IntAct's header anyway -- it carries no
    `taxid:N` -- so an earlier version of this test passed with the explicit
    skip removed. It now asserts the thing only the skip produces: the
    header counted as a header, so `rows read` and `kept` describe data.
    A producer that comments an example row, or a header that happens to
    contain the pattern, is then a reported line rather than a silent one.
    """
    code, out = _run(tmp_path, [HEADER, _row()])
    assert code == 0
    printed = capsys.readouterr().out
    assert "2 rows read, 1 header, 1 kept" in printed
    body = out.read_bytes()
    assert not body.startswith(b"#") and body.count(b"\n") == 1


def test_a_short_row_cannot_match(tmp_path):
    """MITAB 2.5 has no column 15. Reading column 9 out of a 2.5 row would
    take whatever is there as a taxid."""
    short = b"\t".join([b"uniprotkb:P1", b"uniprotkb:Q2"]) + b"\n"
    code, _ = _run(tmp_path, [HEADER, short])
    assert code == 1


def test_either_side_may_carry_the_viral_taxon(tmp_path):
    for a, b in (("9606", "64320"), ("64320", "9606")):
        code, _out = _run(tmp_path / f"d{a}", [HEADER, _row(tax_a=a, tax_b=b)])
        assert code == 0, (a, b)


def test_nothing_is_written_when_nothing_matches(tmp_path, capsys):
    """An empty file would parse to zero edges and read as a source that was
    consulted and had nothing to say."""
    code, out = _run(tmp_path, [HEADER, _row(tax_b="11676")])
    assert code == 1 and not out.exists()
    assert "NO ROW names any declared viral taxon" in capsys.readouterr().out


def test_inspect_writes_nothing(tmp_path):
    code, out = _run(tmp_path, [HEADER, _row()], "--inspect")
    assert code == 0 and not out.exists()


def test_the_output_carries_provenance(tmp_path):
    """A file derived from an archive, with no record of which archive or
    which filter, is the state data/raw was in before fetch_raw."""
    import json
    code, out = _run(tmp_path, [HEADER, _row()])
    assert code == 0
    side = out.with_suffix(out.suffix + ".source.json")
    assert side.exists()
    rec = json.loads(side.read_text())
    assert rec["derived_from"] == "intact.zip"
    assert rec["rows_kept"] == 1
    assert "64320" in rec["filtered_on"]


def test_an_ambiguous_archive_is_refused_rather_than_guessed(tmp_path):
    p = tmp_path / "intact.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("a.txt", b"")
        z.writestr("b.txt", b"")
    with pytest.raises(SystemExit, match="exactly one"):
        fi.main(["--zip", str(p), "--out", str(tmp_path / "o.tab27")])


def test_the_archive_is_never_fully_extracted(tmp_path, monkeypatch):
    """1.35 GB compressed. extractall would write several GB to disk."""
    called = []
    monkeypatch.setattr(zipfile.ZipFile, "extractall",
                        lambda self, *a, **k: called.append(1))
    _run(tmp_path, [HEADER, _row()])
    assert called == []


def test_it_streams_rather_than_reading_the_member_whole(tmp_path, monkeypatch):
    """`for line in fh` over the member, not fh.read(). The difference is
    several GB of resident memory."""
    real = zipfile.ZipFile.open
    seen = []

    class _NoRead(io.RawIOBase):
        def __init__(self, inner):
            self._inner = inner

        def readable(self):
            return True

        def readinto(self, b):
            return self._inner.readinto(b)

        def read(self, size=-1):
            if size is None or size < 0:
                seen.append("read-all")
                raise AssertionError("member read whole")
            return self._inner.read(size)

    def _open(self, name, mode="r", *a, **k):
        fh = real(self, name, mode, *a, **k)
        if mode != "r":                      # writing the fixture
            return fh
        return io.BufferedReader(_NoRead(fh))

    z = _zip(tmp_path, [HEADER, _row()])     # built BEFORE the patch
    monkeypatch.setattr(zipfile.ZipFile, "open", _open)
    code = fi.main(["--zip", str(z),
                    "--out", str(tmp_path / "intact_viral.tab27")])
    assert code == 0 and seen == []


# --------------------------------------- the driver reads both, labelled
def _releases(tmp: Path):
    """Minimal spine+host so AliasIndex resolves the two interactors."""
    import json
    for name, nodes in (
            ("v0.1-spine", [{"id": "UniProtKB:Q22222", "class": "Protein",
                             "properties": {"is_viral": True,
                                            "taxon_id": "NCBITaxon:64320",
                                            "sequence_hash": "a" * 16,
                                            "reviewed": True}}]),
            ("v0.1-host", [{"id": "UniProtKB:P11111", "class": "Protein",
                            "properties": {"is_viral": False,
                                           "taxon_id": "NCBITaxon:9606",
                                           "sequence_hash": "b" * 16,
                                           "reviewed": True}}])):
        d = tmp / "releases" / name
        d.mkdir(parents=True)
        (d / "nodes.jsonl").write_text(
            "\n".join(json.dumps(n) for n in nodes))
        (d / "edges.jsonl").write_text("")
    return tmp / "releases"


def _drive(tmp: Path, *files: tuple[str, bytes]):
    import io
    from contextlib import redirect_stdout
    raw = tmp / "vhppi"
    raw.mkdir(parents=True)
    for name, body in files:
        (raw / name).write_bytes(body)
    rel = _releases(tmp)
    mod = import_module("ingest_vhppi")
    buf = io.StringIO()
    old = sys.argv
    sys.argv = ["ingest_vhppi.py", "--vhppi-dir", str(raw),
                "--spine", str(rel / "v0.1-spine"),
                "--host", str(rel / "v0.1-host"),
                "--out", str(tmp / "out")]
    try:
        with redirect_stdout(buf):
            mod.main()
    except SystemExit as e:
        return buf.getvalue(), e.code
    finally:
        sys.argv = old
    return buf.getvalue(), 0


def test_an_intact_file_is_read_and_labelled_as_intact(tmp_path):
    """THE WHOLE POINT. The glob was `vhn_*.tab27`, so intact.zip sat in
    the same directory unread while the flavivirus gap was recorded as a
    missing download."""
    import json
    _out, _ = _drive(tmp_path, ("intact_viral.tab27", _row()))
    edges = [json.loads(x) for x in
             (tmp_path / "out" / "edges.jsonl").read_text().splitlines() if x]
    assert edges, "the intact file produced no edge"
    assert all(e["primary_knowledge_source"] == "infores:intact"
               for e in edges), [e["primary_knowledge_source"] for e in edges]


def test_a_vhn_file_keeps_its_own_source(tmp_path):
    """An IntAct row recorded as infores:virhostnet is a false provenance
    claim, and the tier ablation reads these labels."""
    import json
    _out, _ = _drive(tmp_path, ("vhn_64320.tab27", _row()))
    edges = [json.loads(x) for x in
             (tmp_path / "out" / "edges.jsonl").read_text().splitlines() if x]
    assert edges
    assert all(e["primary_knowledge_source"] == "infores:virhostnet"
               for e in edges)


def test_both_families_are_read_in_one_run(tmp_path):
    out, _ = _drive(tmp_path,
                    ("vhn_64320.tab27", _row(a="P11111")),
                    ("intact_viral.tab27", _row(a="P11111")))
    assert "infores:virhostnet" in out and "infores:intact" in out
    assert "2 MITAB files" in out


def test_the_same_fact_from_both_sources_keeps_both_provenances():
    """NOT deduplicated inside the layer, and that is deliberate.

    The first version of this test asserted one edge per pair and failed.
    The code was right: two rows is how both source labels survive to
    assemble_release, which merges on (subject, predicate, object) and
    writes the combined source -- the existing graph already carries
    "infores:string|infores:virhostnet" from exactly this mechanism.
    Deduplicating here would silently drop whichever source came second,
    and the tier ablation reads these labels.
    """
    from kgav.assemble import merge_edges
    a = {"subject": "UniProtKB:P11111", "predicate": "PHYSICALLY_INTERACTS_WITH",
         "object": "UniProtKB:Q22222", "primary_knowledge_source":
         "infores:virhostnet", "first_asserted_date": "2024-01-01",
         "evidence_tier": 1, "publications": ["PMID:12345678"],
         "qualifiers": {}}
    b = dict(a, primary_knowledge_source="infores:intact")
    merged = merge_edges(a, b)
    assert "infores:virhostnet" in merged["primary_knowledge_source"]
    assert "infores:intact" in merged["primary_knowledge_source"]


def test_both_sources_reach_the_layer_for_a_shared_pair(tmp_path):
    """The layer file must contain both, or there is nothing to merge."""
    import json
    _drive(tmp_path,
           ("vhn_64320.tab27", _row()),
           ("intact_viral.tab27", _row()))
    edges = [json.loads(x) for x in
             (tmp_path / "out" / "edges.jsonl").read_text().splitlines() if x]
    assert {e["primary_knowledge_source"] for e in edges} == {
        "infores:virhostnet", "infores:intact"}


def test_an_empty_directory_names_the_filter_script(tmp_path):
    out, code = _drive(tmp_path)
    assert code != 0
    assert "filter_intact.py" in str(out) + str(code)
