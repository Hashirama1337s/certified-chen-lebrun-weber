#!/bin/bash
# Finishing production for BAR-MET2/MET3 after amend3 + amend4 (frozen, pinned code only: SHA256SUMS-PROD-amend4.txt).
# Launch at below-normal priority (children inherit). Idempotent and checkpointed: every step skips work already done,
# so the script can simply be started again after an interruption.
cd "$(dirname "$0")"
MC="python"
# same headroom as run_prod2.sh: at most 5 parallel processes x 2 numeric threads
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs/prod3 prodA2 prodB1
S=logs/prod3/status.log
st() { echo "$1 $(date -u +%FT%TZ)" >> $S; }
st "START"

# 0. wait until no amend2 job (prod_B.py, make_cert.py, run_prod2.sh) is alive; two consecutive empty polls. Nothing is killed.
busy() {
  # python jobs prod_B.py / make_cert.py, or the run_prod2.sh driver itself (not shell wrappers that merely mention the names)
  powershell -NoProfile -Command "@(Get-CimInstance Win32_Process | Where-Object { (\$_.Name -like 'python*' -and \$_.CommandLine -match 'prod_B\.py|make_cert\.py') -or (\$_.Name -like 'bash*' -and \$_.CommandLine -match 'run_prod2\.sh' -and \$_.CommandLine -notmatch 'snapshot') }).Count" | tr -d '\r '
}
z=0
while [ $z -lt 2 ]; do
  n=$(busy); if [ "$n" = "0" ]; then z=$((z+1)); else z=0; fi
  [ $z -lt 2 ] && sleep 60
done
st "no amend2 job alive"

# 1. unit tests (stop on failure)
$MC t_unit_safe.py > logs/prod3/unit.log 2>&1; rc=$?; st "unit rc $rc"
[ $rc -ne 0 ] && { st "STOP unit test failed"; exit 1; }
$MC t_amend3.py > logs/prod3/t_amend3.log 2>&1; rc=$?; st "t_amend3 rc $rc"
grep -q "ALL PASS" logs/prod3/t_amend3.log || { st "STOP t_amend3 failed"; exit 1; }

# 2. triage: move (never delete) files that take a new amend3/amend4 path -> prodB1_superseded_amend3/, prodA2_superseded_amend3/
$MC amend3_triage.py B prodB1 --move > logs/prod3/triage_B.log 2>&1; st "triage B: $(tail -2 logs/prod3/triage_B.log | head -1)"
$MC amend3_triage.py A prodA2 --move > logs/prod3/triage_A.log 2>&1; st "triage A: $(tail -2 logs/prod3/triage_A.log | head -1)"

# 3. torus audit of every remaining written file (B first, then A; 5 in parallel); a finished audit log is not redone
for P in B A; do
  D=prodB1; [ $P = A ] && D=prodA2
  for k in 0 1 2 3 4; do
    L=logs/prod3/audit${P}_$k.log
    if grep -q "^audit " $L 2>/dev/null; then continue; fi
    $MC audit_torus.py $P $k $D > $L 2>&1 &
  done
  wait
  for k in 0 1 2 3 4; do
    L=logs/prod3/audit${P}_$k.log
    st "audit $P e$k: $(grep '^audit ' $L || echo 'NO SUMMARY (audit crashed?)')"
    # NaN DROPPED = a crude bound finite in the file but inf now -> move aside, it is recomputed below
    # (logs/prod3/moved.txt makes this run once per file: on a restart a recomputed amend4 file of the same name is NOT moved)
    for f in $(grep "NaN DROPPED" $L | cut -d' ' -f1); do
      grep -qx "$D/$f" logs/prod3/moved.txt 2>/dev/null && continue
      if [ -f $D/$f ]; then mkdir -p ${D}_superseded_amend3; mv $D/$f ${D}_superseded_amend3/; echo "$D/$f" >> logs/prod3/moved.txt; st "  moved $D/$f (NaN dropped)"; fi
    done
    grep "DIFFERENT" $L | sed 's/^/  FLAG (not moved, lead to decide): /' >> $S
  done
done

# 4. part A: recompute missing sub-squares (prod_A skips existing files)
for k in 0 1 2 3 4; do
  n=$(ls prodA2/clw_e${k}_*_*.json 2>/dev/null | wc -l)
  if [ $n -lt 16 ]; then $MC prod_A.py $k prodA2 > logs/prod3/A_$k.log 2>&1 & fi
done
wait
st "A done $(ls prodA2/clw_e*_*_*.json | wc -l) files"

# 5. part B: fill every missing sub-square (prod_B skips existing files)
for k in 0 1 2 3 4; do
  n=$(ls prodB1/clwB_e${k}_*_*.json 2>/dev/null | wc -l)
  if [ $n -lt 16 ]; then $MC prod_B.py $k prodB1 > logs/prod3/B_$k.log 2>&1 & fi
done
wait
st "B done $(ls prodB1/clwB_e*_*_*.json | wc -l) files"

# 6. merge
$MC make_cert.py prodA2 cert.json prodB1 > logs/prod3/make_cert.log 2>&1; st "make_cert rc $? $(grep -o '"certified": [a-z]*' cert.json)"
st "ALL-DONE"
