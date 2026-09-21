#!/bin/bash
# Re-run the whole Geant4 analysis (safe on partial data). Usage: analysis/g4/run_all.sh
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
$PY analysis/g4/analyze.py && $PY analysis/g4/figures.py
