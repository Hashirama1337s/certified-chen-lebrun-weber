#!/bin/bash
# BAR-MET5 refinement runner (frozen, pinned code only: SHA256SUMS-PROD-refine.txt). Launch at below-normal priority.
# prodA3/ = copy of prodA2/ in which every sub-square with > 0.5% of the total lap2 upper bound is replaced by its 2 x 2 refinement
# (prod_A_refine.py: same degree 28, prec 192, K 4, same a* ball; conservative merge). prodA2/ is never modified. B is unchanged.
# Idempotent and checkpointed: copies never overwrite (cp -n), refined files are skipped, so it can simply be restarted.
cd "$(dirname "$0")"
MC="python"
export NUMBA_NUM_THREADS=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NREF=2
mkdir -p logs/refine prodA3
S=logs/refine/status.log
st() { echo "$1 $(date -u +%FT%TZ)" >> $S; }
st "START (NREF=$NREF)"

# 1. the list is computed from prodA2 (the originals), never from prodA3
LIST=$($MC refine_list.py prodA2 0.005 2> logs/refine/list.log | head -1)
st "list: $LIST"
[ -z "$LIST" ] && { st "STOP empty list"; exit 1; }

# 2. copy prodA2 -> prodA3 without overwriting anything already there
cp -n prodA2/clw_e*_*_*.json prodA3/
st "prodA3 has $(ls prodA3/clw_e*_*_*.json | wc -l) files"

# 3. refine, the list dealt round-robin over at most 5 processes
IFS=',' read -ra L <<< "$LIST"
for p in 0 1 2 3 4; do
  part=""
  for idx in "${!L[@]}"; do [ $((idx % 5)) -eq $p ] && part="${part:+$part,}${L[$idx]}"; done
  [ -n "$part" ] && $MC prod_A_refine.py prodA2 prodA3 "$part" $NREF > logs/refine/refine_$p.log 2>&1 &
done
wait
for p in 0 1 2 3 4; do [ -f logs/refine/refine_$p.log ] && grep -h "merged\|Traceback\|Error" logs/refine/refine_$p.log | sed 's/^/  /' >> $S; done
nref=$($MC -c "import json,glob; print(sum(1 for f in glob.glob('prodA3/clw_e*.json') if json.load(open(f))['A'].get('refined') == $NREF))" | head -1)
st "refined files in prodA3: $nref of ${#L[@]}"
[ "$nref" != "${#L[@]}" ] && { st "STOP refinement incomplete"; exit 1; }

# 4. merge with the unchanged production B
$MC make_cert.py prodA3 cert_refined.json prodB1 > logs/refine/make_cert.log 2>&1
st "make_cert rc $? $(grep -o '"certified": [a-z]*' cert_refined.json) $(grep -o '"lap_rho_L2": "[^"]*"' cert_refined.json)"
st "ALL-DONE"
