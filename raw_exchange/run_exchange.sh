#!/bin/bash
# Restartable driver: one process per job, 1 thread each. Safe to re-run after any interruption (skips finished jobs).
cd "$(dirname "$0")"; source ~/geant4/install/bin/geant4.sh
export G4FORCENUMBEROFTHREADS=1
BIN=~/CHON/recalc_2026-09-19/build/bin/briks; MAC=~/CHON/recalc_2026-09-19/macros/exchange; NPAR=${NPAR:-14}
echo "[$(date -Is)] driver start NPAR=$NPAR binary_sha256=$(sha256sum $BIN | cut -c1-16)" >> driver.log
run_one() {
  n=$1; [ -f "$n/done.flag" ] && exit 0
  rm -rf "$n"; mkdir -p "$n"; cd "$n"
  "$BIN" "$MAC/$n.mac" > run.log 2>&1; rc=$?
  rows=$(( $(cat MAC_results_*.csv 2>/dev/null | wc -l) - 1 ))
  if [ $rc -eq 0 ] && [ $rows -eq 1 ]; then touch done.flag; else echo "[$(date -Is)] FAILED $n rc=$rc rows=$rows" >> ../failed.log; fi
}
export -f run_one; export BIN MAC
xargs -a "$MAC/jobs.txt" -P "$NPAR" -I{} bash -c 'run_one {}'
echo "[$(date -Is)] driver end: done=$(ls -d */ 2>/dev/null | while read d; do [ -f "$d/done.flag" ] && echo x; done | wc -l) of $(wc -l < $MAC/jobs.txt)" >> driver.log
