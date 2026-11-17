# Figure captions: certified-chen-lebrun-weber

All figures are drawn by `figures/make_figures.py <repo root>` from shipped files only (Python 3.11 + matplotlib 3.11;
no other dependency). Two runs give byte-identical SVGs. Printed upper bounds are rounded up.

---

## 1. `figures/checked-patch-by-patch.svg` (hero)

**Point.** The computed metric was checked everywhere, not sampled: the whole space is covered by 80 patches, each with
a certified bound on how far the computed metric misses its equation, and all of them are below 3.81 × 10⁻¹⁷.

**Caption (as it will appear).**
The computed metric, checked patch by patch. The space of the CLW metric is drawn as a pentagon: each inside point
stands for a doughnut surface (2-torus) of points in the 4-dimensional space, so covering the pentagon covers the whole
space. The proof splits it into 80 patches. The colour of each patch is a certified upper bound on the mismatch
|Scal − S|, how far the scalar curvature of the computed metric is from the value its equation demands (log scale,
lighter = larger). Every patch is below 3.81 × 10⁻¹⁷. From such bounds the proof concludes that the true CLW Einstein
metric lies within a factor 1 ± 5.504 × 10⁻⁸ of the rescaled computed one at every point (conditional on the cited
results listed in the paper, one of them a preprint).

**Alt text.** A pentagon (a square with its lower-left corner cut off) divided into 80 curved four-sided patches,
coloured in shades of blue from dark (about 10⁻²²) to light (about 4 × 10⁻¹⁷), with a colour bar from 10⁻²² to 10⁻¹⁶.
Text: every patch below 3.81 × 10⁻¹⁷; result: the true metric lies within a factor 1 ± 5.504 × 10⁻⁸ of the computed one.

**Sources.**

| Plotted / printed | File (field) |
|---|---|
| 80 patches: element k, (s0, t0), side h | `data/prodA3/clw_e{k}_{i}_{j}.json` (`k`, `i`, `j`, `A.s0`, `A.t0`, `A.h`) |
| patch colour: certified upper bound on \|ρ\|, ρ = Scal(g̃) − S | same files (`A.sup.rho`); definition of ρ in `code/seatB2/assemble.py` header ("rho = Scal(g~) - S") |
| every patch below 3.81 × 10⁻¹⁷ (from 3.8046022124832956e-17) | `data/cert_refined_a1.json` (`rho_sup`); the script asserts it equals the maximum of the 80 `A.sup.rho`, rounded up |
| smallest patch bound 4.07e-22 (colour scale 10⁻²² to 10⁻¹⁶) | min / max of `A.sup.rho` over the 80 files |
| 1 ± 5.504 × 10⁻⁸ (from d_h = 5.5036608022109204886e-8) | `data/result_refined_a1.json` (`d_h`); README "At a glance" (Pointwise) |
| pentagon P = [−1, 1]² cut by x₁ + x₂ ≥ −γ, γ = 1.0216003055709606 | `code/seatR/sol_clw_B448_M84.json` (`gamma`); `paper.md` §1 (C3) |
| patch outlines | element map of `code/seatE/common.py` (`ElemMap`, `elem_x`), re-implemented in floating point for drawing only; the script checks that the 80 outlines tile the pentagon (area 3.52137, agreement better than 1e-3) |
| "conditional on cited results, one a preprint" | README "In plain words", "At a glance" (Conditions) and "What we do not claim" (same words) |

The pentagon is the moment polygon of the toric 4-manifold (inside points ↔ 2-tori; edges and corners ↔ circles and
points). Note the colour shows the extremal-equation residual of the computed (extremal Kähler) approximant; the
Einstein metric is its conformal rescaling (README "At a glance").

---

## 2. `figures/the-deciding-test.svg`

**Point.** The whole proof hinges on one number being below 1; the real certificate passes with room to spare, and a
deliberately broken input fails, so the test has teeth.

**Caption.**
The test that decides the proof. The Newton-Kantorovich argument closes only if one number, x, is below 1. With the
certified inputs, x = 0.1773 (pass; 18 of 18 gates hold). In the null control one input (the size of the residual's
Laplacian) is deliberately raised above its budget, using sampled, not certified, inputs as a dry run: x = 1.7748, and
the proof correctly refuses to close.

**Alt text.** Two horizontal bars against a vertical line at 1: a short green bar to 0.1773 labelled PASS (the
certificate), and a long orange bar to 1.7748 labelled FAIL (the null control, a sampled dry run).

**Sources.**

| Plotted / printed | File (field) |
|---|---|
| x = 0.1773 (from 0.17725474960624380611) | `data/result_refined_a1.json` (`x`); `data/logs/assemble_refined_a1.log` line "x = 4 a1 a2^2 a3_Y = 0.177254749606" |
| 18 of 18 gates | `data/result_refined_a1.json` (`gates`, `all_gates`); log line "GATES: 18/18 pass" |
| x = 1.7748 (from 1.77478081110), FAIL | `data/logs/assemble_null.log` (line "x = 4 a1 a2^2 a3_Y" and "GATE x < 1 (K14): FAIL"; the script asserts both) |
| null control: Laplacian of the residual raised above budget; sampled, not certified, a dry run | `data/cert_null.json` (`label`); `data/logs/assemble_null.log` last line; `code/RUN.md` §1 |

## Dropped

- The scalar gap (exactly one eigenvalue below 2Λ, λ₁/Λ ≤ 1.948367728, λ₂/Λ ≥ 2.089704146) would be a strong picture,
  but its numbers come from the companion spectral lane's files (BAR-CLWT-AMEND2.md, REFEREE-CLWT.md), which are not in
  this repository; only the prose in `paper.md` §5 carries them. Not drawn, per the shipped-data rule. If the lane's
  result file is shipped later, a number-line figure (λ₁ window, threshold 2, λ₂ lower bound) is a five-minute add.
- A first-assembly figure (MET4, x = 1.246, not closed) for the same reason: its source log is lab-internal.
