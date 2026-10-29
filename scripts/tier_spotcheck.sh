#!/usr/bin/env bash
# Run the agent against ONE VIRUS FROM EACH TIER in section 9 of
# docs/kg_audit_report.txt, and print what it says.
#
# WHY THIS EXISTS. Every hand-check of this agent up to 2026-10-29
# used SARS-CoV-2, which sits in the EVALUABLE tier -- labels, routes
# and reachable positives all present. So the only path ever read by
# eye was the one where everything works. The three tiers below had
# unit tests and had never been read.
#
# The first run found that
#
#     which compounds should be considered for Yellow fever
#
# returned a node profile instead of the refusal written for exactly
# that virus. The refusal text existed and was correct; nothing
# reached it. A unit test cannot find that, because a unit test asks
# the question the author thought to ask.
#
# This is not a test and has no assertions. It prints output for a
# person to read, which is the point: the defect it found was in what
# the agent SAID, and it looked fine in every machine-checkable way.
#
# --no-model: retrieval and guards only, no Ollama needed.
set -u
cd "$(dirname "$0")/.." || exit 1

run () {
  echo
  echo "################################################################"
  echo "## $1  --  $2"
  echo "################################################################"
  python scripts/kgav_ask.py --no-model "$2" 2>&1 || echo "[exit $?]"
}

#          tier                        virus
#          EVALUABLE                   SARS-CoV-2   (the control)
#          ROUTES REACH NO LABEL       DENV-1
#          LABELS, NO FEATURE ROUTE    Yellow fever
#          FEATURES, NOTHING TO SCORE  HCoV-HKU1
for v in "Yellow fever" "DENV-1" "HCoV-HKU1" "SARS-CoV-2"; do
  run "$v" "which compounds should be considered for $v"
  run "$v" "which pathways does $v depend on"
done
