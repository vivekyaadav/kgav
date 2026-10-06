import { useMemo, useRef, useState } from "react";
import {
  Activity,
  ArrowDownRight,
  ChevronDown,
  CircleHelp,
  Database,
  Focus,
  GitBranch,
  Info,
  Layers3,
  Maximize2,
  MousePointer2,
  Orbit,
  Play,
  RotateCcw,
  Search,
  SlidersHorizontal,
  Sparkles,
  Target,
  X,
  Zap,
} from "lucide-react";

type NodeType = "Virus" | "Protein" | "Gene" | "Pathway" | "Compound";
type GraphNode = {
  id: string;
  label: string;
  type: NodeType;
  x: number;
  y: number;
  z: number;
  size: number;
  meta: string;
  description: string;
};
type Edge = { source: string; target: string; label: string; strength: number };

const typeColors: Record<NodeType, string> = {
  Virus: "#ff8a72",
  Protein: "#7c9cff",
  Gene: "#d7a7ff",
  Pathway: "#65d9c1",
  Compound: "#f6c86e",
};

const nodes: GraphNode[] = [
  { id: "sars", label: "SARS-CoV-2", type: "Virus", x: 50, y: 48, z: 1, size: 56, meta: "NCBITaxon:2697049", description: "A coronavirus node connected to viral proteins and downstream host biology." },
  { id: "spike", label: "Spike", type: "Protein", x: 31, y: 25, z: .86, size: 33, meta: "UniProt: P0DTC2", description: "Viral surface protein used here to show a cross-entity bridge into host factors." },
  { id: "nsp12", label: "nsp12 / RdRp", type: "Protein", x: 70, y: 27, z: .92, size: 36, meta: "UniProt: P0DTD1", description: "RNA-dependent RNA polymerase. A central viral protein in the repurposing graph." },
  { id: "ace2", label: "ACE2", type: "Protein", x: 20, y: 55, z: .76, size: 30, meta: "UniProt: Q9BYF1", description: "Human host receptor that anchors the illustrative entry path." },
  { id: "tmprss2", label: "TMPRSS2", type: "Protein", x: 21, y: 78, z: .63, size: 25, meta: "UniProt: O15393", description: "Host protease associated with viral entry." },
  { id: "npc1", label: "NPC1", type: "Protein", x: 47, y: 82, z: .74, size: 31, meta: "UniProt: O15118", description: "A host dependency factor surfaced from screen-level evidence." },
  { id: "ly6e", label: "LY6E", type: "Gene", x: 76, y: 73, z: .7, size: 28, meta: "HGNC: 6727", description: "A restriction factor anchor. Its polarity is intentionally visible in the graph." },
  { id: "ifn", label: "Type I IFN", type: "Pathway", x: 87, y: 46, z: .83, size: 31, meta: "Reactome: R-HSA-909733", description: "Pathway context turns individual edges into a biological story." },
  { id: "remdesivir", label: "Remdesivir", type: "Compound", x: 63, y: 59, z: .58, size: 27, meta: "ChEMBL: CHEMBL4065616", description: "A compound node showing how a drug can enter through a target or mechanism path." },
  { id: "pik3c3", label: "PIK3C3", type: "Gene", x: 39, y: 62, z: .52, size: 23, meta: "HGNC: 8974", description: "A dependency-side gene included to demonstrate query expansion." },
];

const edges: Edge[] = [
  { source: "sars", target: "spike", label: "has_protein", strength: 1 },
  { source: "sars", target: "nsp12", label: "has_protein", strength: .95 },
  { source: "spike", target: "ace2", label: "interacts_with", strength: .86 },
  { source: "ace2", target: "tmprss2", label: "coexpressed_with", strength: .66 },
  { source: "nsp12", target: "npc1", label: "host_dependency", strength: .9 },
  { source: "npc1", target: "ly6e", label: "screen_context", strength: .54 },
  { source: "ly6e", target: "ifn", label: "participates_in", strength: .88 },
  { source: "nsp12", target: "remdesivir", label: "targeted_by", strength: .98 },
  { source: "remdesivir", target: "pik3c3", label: "similar_mechanism", strength: .42 },
  { source: "pik3c3", target: "sars", label: "dependency_signal", strength: .6 },
  { source: "ifn", target: "sars", label: "restricts", strength: .48 },
];

const filters: { label: NodeType; count: number }[] = [
  { label: "Virus", count: 7 },
  { label: "Protein", count: 52_763 },
  { label: "Gene", count: 18_204 },
  { label: "Pathway", count: 1_182 },
  { label: "Compound", count: 6_430 },
];

function AppButton({ children, active = false, onClick }: { children: React.ReactNode; active?: boolean; onClick?: () => void }) {
  return <button onClick={onClick} className={`app-button ${active ? "active" : ""}`}>{children}</button>;
}

export default function Home() {
  const [selectedId, setSelectedId] = useState("nsp12");
  const [query, setQuery] = useState("");
  const [activeFilters, setActiveFilters] = useState<NodeType[]>(["Virus", "Protein", "Gene", "Pathway", "Compound"]);
  const [showLabels, setShowLabels] = useState(true);
  const [autoRotate, setAutoRotate] = useState(false);
  const [guided, setGuided] = useState(false);
  const [queryMode, setQueryMode] = useState<"path" | "neighborhood">("path");
  const [zoom, setZoom] = useState(1);
  const [notice, setNotice] = useState("");
  const [rotation, setRotation] = useState({ x: -8, y: 0 });
  const dragRef = useRef({ active: false, x: 0, y: 0, rx: -8, ry: 0 });

  const visibleNodes = useMemo(() => nodes.filter((node) => activeFilters.includes(node.type) && node.label.toLowerCase().includes(query.toLowerCase())), [activeFilters, query]);
  const selected = nodes.find((node) => node.id === selectedId) ?? nodes[0];
  const connected = edges.filter((edge) => edge.source === selectedId || edge.target === selectedId);
  const connectedIds = new Set(connected.flatMap((edge) => [edge.source, edge.target]));
  const visibleEdges = edges.filter((edge) => activeFilters.includes(nodes.find((n) => n.id === edge.source)!.type) && activeFilters.includes(nodes.find((n) => n.id === edge.target)!.type));

  const toggleFilter = (type: NodeType) => setActiveFilters((current) => current.includes(type) ? current.filter((item) => item !== type) : [...current, type]);
  const reset = () => { setQuery(""); setSelectedId("nsp12"); setActiveFilters(filters.map((filter) => filter.label)); setZoom(1); setRotation({ x: -8, y: 0 }); setNotice("View reset"); };
  const handlePointerDown = (event: React.PointerEvent<HTMLDivElement>) => {
    dragRef.current = { active: true, x: event.clientX, y: event.clientY, rx: rotation.x, ry: rotation.y };
    event.currentTarget.setPointerCapture(event.pointerId);
  };
  const handlePointerMove = (event: React.PointerEvent<HTMLDivElement>) => {
    if (!dragRef.current.active) return;
    const nextX = Math.max(-28, Math.min(28, dragRef.current.rx - (event.clientY - dragRef.current.y) * 0.18));
    const nextY = Math.max(-32, Math.min(32, dragRef.current.ry + (event.clientX - dragRef.current.x) * 0.18));
    setRotation({ x: nextX, y: nextY });
  };
  const handlePointerUp = (event: React.PointerEvent<HTMLDivElement>) => {
    dragRef.current.active = false;
    event.currentTarget.releasePointerCapture(event.pointerId);
  };
  const exportGraph = () => {
    const payload = { nodes: visibleNodes, edges: visibleEdges, selected: selected.id, mode: queryMode };
    const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url; link.download = "kgav-graph-view.json"; link.click();
    URL.revokeObjectURL(url); setNotice("Graph view exported");
  };

  return (
    <div className="app-shell">
      <header className="topbar">
        <div className="brand-lockup"><div className="brand-mark"><Orbit size={19} strokeWidth={2.3} /></div><div><div className="brand-name">KGAV <span>/ VISUALISE</span></div><div className="brand-subtitle">Knowledge graph observatory</div></div></div>
        <div className="topbar-actions"><div className="live-pill"><span className="live-dot" /> Live graph <span className="muted">·</span> v1.0</div><button className="icon-button" aria-label="Help" onClick={() => { setGuided(true); setNotice("Guided help opened"); }}><CircleHelp size={18} /></button><button className="avatar" onClick={() => setNotice("KGAV explorer session")}>VA</button></div>
      </header>
      <main className="main-layout">
        <aside className="sidebar">
          <div className="eyebrow"><Sparkles size={13} /> LEARN BY EXPLORING</div>
          <h1>See how<br /><em>knowledge</em><br />connects.</h1>
          <p className="intro">Trace signals from a virus to host biology, pathways, and candidate compounds. Every dot is an entity. Every line is evidence.</p>
          <div className="search-wrap"><Search size={15} /><input value={query} onChange={(event) => setQuery(event.target.value)} placeholder="Find a node..." /><kbd>⌘ K</kbd></div>
          <div className="section-label"><span>Explore layers</span><button className="clear-button" onClick={() => setActiveFilters([])}>Clear</button></div>
          <div className="filter-list">{filters.map((filter) => <button key={filter.label} className={`filter-row ${activeFilters.includes(filter.label) ? "selected" : ""}`} onClick={() => toggleFilter(filter.label)}><span className="filter-check" style={{ background: activeFilters.includes(filter.label) ? typeColors[filter.label] : "transparent", borderColor: typeColors[filter.label] }} /> <span className="filter-name">{filter.label}</span><span className="filter-count">{filter.count.toLocaleString()}</span></button>)}</div>
          <div className="sidebar-divider" />
          <div className="section-label"><span>Graph settings</span><SlidersHorizontal size={14} /></div>
          <label className="setting-row"><span><Layers3 size={15} /> Show labels</span><input type="checkbox" checked={showLabels} onChange={(event) => setShowLabels(event.target.checked)} /><span className="switch" /></label>
          <label className="setting-row"><span><Orbit size={15} /> Auto rotate</span><input type="checkbox" checked={autoRotate} onChange={(event) => setAutoRotate(event.target.checked)} /><span className="switch" /></label>
          <button className={`guided-card ${guided ? "on" : ""}`} onClick={() => setGuided(!guided)}><div className="guided-icon"><Play size={14} fill="currentColor" /></div><div><strong>{guided ? "Guided mode on" : "Take a guided tour"}</strong><span>Learn the graph in 60 seconds</span></div><ArrowDownRight size={17} /></button>
          <div className="sidebar-footer"><span>Data layers synced</span><span className="sync-status"><span className="sync-dot" /> 2m ago</span></div>
        </aside>

        <section className="workspace">
          <div className="workspace-head"><div><div className="breadcrumb"><span>KGAV</span><span>/</span><strong>Graph explorer</strong></div><div className="workspace-title-row"><h2>Antiviral repurposing graph</h2><span className="dataset-chip">DEMO DATASET <span>·</span> 2021</span></div></div><div className="workspace-actions"><AppButton onClick={reset}><RotateCcw size={15} /> Reset view</AppButton><AppButton active onClick={exportGraph}><Database size={15} /> Export <ChevronDown size={14} /></AppButton></div></div>
          <div className="graph-toolbar"><div className="toolbar-left"><button className="toolbar-label toolbar-button" onClick={() => setAutoRotate((value) => !value)}><MousePointer2 size={14} /> {autoRotate ? "Pause motion" : "Resume motion"}</button><span className="toolbar-sep" /><button className="toolbar-label toolbar-button" onClick={() => setZoom(1)}><Focus size={14} /> Fit graph</button></div><div className="toolbar-right"><button className={`mode-button ${queryMode === "path" ? "active" : ""}`} onClick={() => setQueryMode("path")}><GitBranch size={14} /> Path mode</button><button className={`mode-button ${queryMode === "neighborhood" ? "active" : ""}`} onClick={() => setQueryMode("neighborhood")}><Target size={14} /> Neighborhood</button></div></div>
          <div className="graph-card">
            <div className="graph-grid" />
            <div className="graph-badge"><span className="badge-orb"><Activity size={14} /></span><div><strong>{visibleNodes.length} visible nodes</strong><span>{visibleEdges.length} evidence edges</span></div></div>
            <div className={`graph-space ${autoRotate ? "rotating" : ""}`} onPointerDown={handlePointerDown} onPointerMove={handlePointerMove} onPointerUp={handlePointerUp} onPointerCancel={handlePointerUp} style={{ transform: `scale(${zoom}) rotateX(${rotation.x}deg) rotateY(${rotation.y}deg)` }}>
              <svg className="edge-layer" viewBox="0 0 100 100" preserveAspectRatio="none" aria-hidden="true">{visibleEdges.map((edge) => { const a = nodes.find((n) => n.id === edge.source)!; const b = nodes.find((n) => n.id === edge.target)!; const highlighted = connectedIds.has(edge.source) && connectedIds.has(edge.target); return <line key={`${edge.source}-${edge.target}`} x1={a.x} y1={a.y} x2={b.x} y2={b.y} className={highlighted ? "edge highlighted" : "edge"} style={{ stroke: highlighted ? typeColors[selected.type] : "#344563", strokeWidth: highlighted ? 0.65 : 0.35, opacity: highlighted ? .95 : .65 }} />; })}</svg>
              {visibleNodes.map((node) => { const isSelected = node.id === selectedId; const isConnected = connectedIds.has(node.id); return <button key={node.id} className={`graph-node ${isSelected ? "selected" : ""} ${isConnected ? "connected" : ""}`} style={{ left: `${node.x}%`, top: `${node.y}%`, width: node.size, height: node.size, background: typeColors[node.type], zIndex: Math.round(node.z * 10), transform: `translate(-50%, -50%) scale(${node.z})` }} onClick={() => setSelectedId(node.id)} aria-label={`Select ${node.label}`}><span className="node-glow" style={{ background: typeColors[node.type] }} />{isSelected && <span className="node-crosshair" />}{showLabels && <span className="node-label">{node.label}</span>}</button>; })}
            </div>
            <div className="graph-hint"><span>TIP</span> Click any node to inspect its evidence trail</div>
            <div className="zoom-controls"><button onClick={() => setZoom((value) => Math.min(1.3, +(value + .1).toFixed(1)))}>+</button><span>{Math.round(zoom * 100)}%</span><button onClick={() => setZoom((value) => Math.max(.8, +(value - .1).toFixed(1)))}>−</button><button onClick={() => setZoom(1)}><Maximize2 size={14} /></button></div>
          </div>
          <div className="status-strip"><div><span className="status-icon"><Zap size={14} fill="currentColor" /></span><span><strong>Explore a path</strong> to understand how evidence becomes a connection</span></div><button onClick={() => setGuided(true)}>How does this work? <ArrowDownRight size={14} /></button></div>
        </section>

        <aside className="inspector"><div className="inspector-header"><div><div className="eyebrow small">SELECTED NODE</div><h3>Node details</h3></div><button className="icon-button" onClick={() => { setSelectedId("sars"); setNotice("Selection focused on SARS-CoV-2"); }} aria-label="Reset selection"><X size={17} /></button></div><div className="node-profile"><div className="profile-orb" style={{ background: typeColors[selected.type] }}><span>{selected.type === "Protein" ? "P" : selected.type[0]}</span></div><div><div className="node-type" style={{ color: typeColors[selected.type] }}>{selected.type.toUpperCase()}</div><h4>{selected.label}</h4><span className="node-meta">{selected.meta}</span></div></div><p className="node-description">{selected.description}</p><div className="inspector-divider" /><div className="evidence-heading"><span>Evidence trail</span><span className="evidence-count">{connected.length} links</span></div><div className="evidence-list">{connected.map((edge) => { const otherId = edge.source === selected.id ? edge.target : edge.source; const other = nodes.find((node) => node.id === otherId)!; return <button className="evidence-row" key={`${edge.source}-${edge.target}`} onClick={() => setSelectedId(other.id)}><span className="evidence-line" style={{ background: typeColors[other.type] }} /><span className="evidence-copy"><strong>{edge.source === selected.id ? `→ ${edge.label}` : `← ${edge.label}`}</strong><span>{other.label}</span></span><span className="evidence-score">{Math.round(edge.strength * 100)}%</span></button>; })}</div><div className="inspector-divider" /><div className="explain-box"><div className="explain-icon"><Info size={15} /></div><div><strong>How to read this</strong><p>Edges are claims from different source layers. Stronger lines mean more supporting evidence, not certainty.</p></div></div><button className="primary-action" onClick={() => setGuided(true)}><Sparkles size={15} /> Explain this path <ArrowDownRight size={15} /></button></aside>
      </main>
      {notice && <button className="notice-toast" onClick={() => setNotice("")}>{notice}<X size={14} /></button>}
      {guided && <div className="tour-toast"><div className="tour-step">01 <span>/ 03</span></div><div><strong>Start with the orange node</strong><p>That is the virus. Follow its protein edges into host biology, then inspect the evidence trail on the right.</p></div><button onClick={() => setGuided(false)}><X size={15} /></button></div>}
    </div>
  );
}
