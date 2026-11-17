# RUN: recompute the certificate and the assembly

Requirements: Python 3.11 (tested with 3.11.9), python-flint 0.9.0, numpy (tested with 2.4.4), mpmath (tested with 1.3.0).
Production B (`prod_B.py`) also needs numba (tested with 0.66.0). The lab ran every step under an 8 GB memory cap; the
certificate and the assembly need well under 1 GB.

Run from the top folder of this repository (the one with paper.md). Outputs go to `out/`.

## 1. Certificate and assembly from the frozen production data (seconds)

```
mkdir out
cd code/seatE
python make_cert.py ../../data/prodA3 ../../out/cert_refined_a1.json ../../data/prodB1
cd ../seatB2
python assemble.py ../../out/cert_refined_a1.json ../../out/result_refined_a1.json
python assemble.py ../../data/cert_null.json ../../out/result_null.json
cd ../..
cmp out/cert_refined_a1.json data/cert_refined_a1.json
cmp out/result_refined_a1.json data/result_refined_a1.json
```

Expected: `make_cert.py` prints `production B: 80/80 sub-squares ok`; the first assembly ends with
`GATES: 18/18 pass.  VERDICT: ALL GATES HOLD` (full log: `data/logs/assemble_refined_a1.log`); the null assembly fails
gate 10, x < 1 (`x = 1.77478081110`, `GATE x < 1 (K14): FAIL`; it is a sampled dry-run certificate, so the log ends
with `GATES: 9/10 pass.  VERDICT: NOT CLOSED   [inputs SAMPLED, NOT CERTIFIED: a dry run, not a proof]`); both `cmp` lines print nothing (byte-identical files).
`make_cert.py` also reads `code/seatG45/G4_eps_m48.json` (the certified lambda1 lower bound), `code/seatR/sol_clw_B448_M84.json`
(the degree-84 approximant) and `code/seatB2/cert_sampled.json` (only the two settings `w_prime` and `d_target`).

## 2. Containment controls (BAR-CTRL1; about 2.5 minutes)

```
mkdir -p out
cd code/controls
bash run_ctrl1.sh
```
`mkcopy.py` writes a mutated COPY of `seatE/cm.py` to `controls/mut/`; `seatE/cm.py` is hash-checked and never modified.
Expected `TOTAL violations`: 0 (square, degree 28), 0 (square, degree 12), 0 (CLW against the arb evaluator), 350 (inverse
remainder zeroed, degree 28), 810 (tail zeroed, degree 12). Recorded counts: `data/ctrl1/CTRL1.json`.

## 3. Production (hours; not needed to check the certificate)

From `code/seatE/`: `python prod_A.py k prodA2` and `python prod_B.py k prodB1` for k = 0..4 (one process per element), then
`bash run_refine.sh`, which lists, refines (2 x 2) and merges the far-corner sub-squares into `prodA3/` and then runs
`make_cert.py prodA3 cert_refined.json prodB1`. The lab's timings are in paper Section 4. The night runner `run_prod3.sh` is
shipped for the record.

## Hashes

`data/SHA256SUMS-PROD-met6.txt` is the lab's list (paths relative to the lab folder `seatE/`, i.e. `code/seatE/` here).
Files under `code/` match it except the path-line patches in `PATCHES.md` (check with `sha256sum -c MANIFEST.sha256` from the repository root). `data/prod_inputs.sha256` covers `data/prodA3/`
and `data/prodB1/` (run `sha256sum -c prod_inputs.sha256` from `data/`).

## Figures (optional)

The two README figures are drawn from shipped files only: `python figures/make_figures.py .` from the repository root
(Python 3.11, matplotlib 3.11; seconds). It rewrites `figures/*.svg`; two runs give byte-identical SVGs. Captions and
the source of every printed number: `figures/CAPTIONS.md`.
