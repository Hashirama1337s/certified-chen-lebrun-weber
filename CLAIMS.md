# CLAIMS: every number in paper.md, with its source file and line

Paths are relative to the lab project folder `PAGE/vectors/v16_clwport/` unless they start with `PAGE/` (the project root) or `lab-data/` (the lab ledger store). `LOG` = `PAGE/LOG.md`.
`LEDGER` = `lab-data/ledger.jsonl` (sealed/settled records). `A6` = `seatB2/logs/assemble_refined_a1.log` (the MET6 run).
Values were copied, not recomputed, except for the hashes in C69 (computed with `sha256sum` while drafting, 2026-11-13).

| id | claim / number in paper.md | source (file:line) |
|---|---|---|
| C1 | p(a), degree 10; a* the unique root in (1, inf); Descartes on p(1+z) | seatG45/NOTES-G45.md:158-164 |
| C2 | a* = 1.95771280518776186068...; FLINT ball radius 1.5e-40 (192 bits) | seatE/NOTES-E.md:252; BAR-MET2.md:26; REFEREE-MET5.md:82 |
| C3 | P = [-1,1]^2 cut by x1 + x2 >= -gamma, gamma = 2/a* | seatG45/NOTES-G45.md:130; SCOPE-EXTREMAL-ROUTE.md:100 |
| C4 | u~ = u_G + f84, f of total degree 84 (exact expansion) | seatG45/NOTES-G45.md:31; seatE/REFEREE-CM.md:22 (item 1c) |
| C5 | H^k sum norms | seatQ2/NOTES-Q2-PROOF.md:8 |
| C6 | Theorem 1: extremal metric for every a in the ball, \|\|psi\|\|_H6 <= 7.568e-10 | REFEREE-MET5.md:82-83; A6:65 (r = 7.56799693326e-10); seatB2/result_refined_a1.json:18 |
| C7 | "complex-gauge sum norm, not a distance to u84"; Scal(g*) = S(x*), c* = 0 | REFEREE-MET5.md:83; seatB2/NOTES-B2.md:66-82 |
| C8 | d_h <= 5.504e-8 (5.50366080221e-8), d_sym 3.19e-8 (3.18593544121e-8), d_cx 3.18593507364e-8 | A6:78-82; REFEREE-MET5.md:87; LEDGER:308 |
| C9 | gauge definitions: complex = identity map, same J; symplectic = same (x, theta); Einstein d_h formula | seatB2/NOTES-B2.md:162-193 |
| C10 | LeBrun Prop. 3 quote | REFEREE-MET5.md:85; SCOPE-EXTREMAL-ROUTE.md:10-11 |
| C11 | Li-Sun Cor. 7.3 (l.4511-4515), Prop. 7.13 (l.5338), preprint | REFEREE-MET5.md:86; REFEREE-EXTREMAL.md:27-30 |
| C12 | S_min = 3.00951670796 > 0 | A6:10; seatE/cert_refined_a1.json:29 |
| C13 | CLW existence known (2008); LeBrun's 2013 "new and rather different proof" | REFEREE-MET5.md:88; SCOPE-EXTREMAL-ROUTE.md:36-43 |
| C14 | "computer-assisted re-proof plus a certified location" | REFEREE-MET5.md:88 |
| C15 | LeBrun's openness step is qualitative; Kantorovich is the quantitative version | SCOPE-EXTREMAL-ROUTE.md:38-44 |
| C16 | LEDGER1 HOLD: Vol(CLW) Lambda^2 = 132.78666547165841..., matches Hall-Murphy 9 digits | LOG:195; LEDGER:252 |
| C17 | Lemma G5 (swap, Euler, det = 2); D(1+z) positive coefficients; D(a*) = 368.81...; family formula from Nifa (6.6)-(6.8) | seatG45/NOTES-G45.md:133-169 |
| C18 | complex gauge g_psi, x_psi, slice X, map G, DG(0), L_t, E | seatB2/NOTES-B2.md:17-29, 14 |
| C19 | Y-norm, w' = 2^-10, \|\|Bo^-1\|\| <= C_reg, kappa_A = 1 + k1 + k2, norm equivalence | seatB2/NOTES-B2.md:37-46, 64 |
| C20 | Pi rho = 0 exactly, so alpha1 = \|\|rho\|\|_H2 | seatB2/NOTES-B2.md:51-52 |
| C21 | Donaldson primary source opened, VERDICT MATCH (factor 2, lattice measure); Abreu Thm 2.8 opened; face-convexity caveat | CITATIONS-ABREU-DONALDSON.md:16-18, 35-39, 57-58, 61-62 |
| C22 | numerical check of Donaldson's identity (strong form; floors to 1e-59 / 1e-70; negative controls 7.02, 14.85 vs 14.04) | seatB2/NOTES-B2.md:102-116 |
| C23 | c* = 0 argument and margin criterion | seatB2/NOTES-B2.md:66-80 |
| C24 | NK used as plain contraction mapping; dP3 Thm 4.1 not load-bearing | seatB2/PROOF-ERRATA.md:82-85 |
| C25 | alpha2 = C_reg/(1 - t), t formula; alpha3_Y = (kappa_A + w') 4 C32 K(2 R_hat) | seatB2/NOTES-B2.md:54-60; A6:45-49, 51-60 |
| C26 | x = 4 a1 a2^2 a3_Y; r = 2 a1 a2/(1 + sqrt(1 - x)) | A6:63-65 |
| C27 | gates 1-7 and their values: mu 0.483347103589, Rm 0, S_min 3.00951670796, V 139.017997713, C_reg 9989.92845835, t 0.00385206247051 | A6:7-11, 27-28, 47-48 |
| C28 | gates 8-13: C9 s = 9.69498e-7; x = 0.177254749606; r <= R_hat; self-map; contraction 0.0929469 | A6:55, 59, 63-68 |
| C29 | C5..C8, C_inf 29.7674049667, C32 30.8202310004, C9 59.5348099334 (dP3 App. B as transcribed, from mu and V only) | A6:30-37 |
| C30 | gates 14-18: b 3.18593507364e-8, c* margin 2.85876742210e-7 (coefficient 377.744262757), d targets 0.006530; 18/18; memcap peak 0.01 GB | A6:71-89 |
| C31 | alpha1 7.19573547909e-14 (pieces 1.5268e-16, 3.2356e-15, 6.8569e-14), budget 4.05953742741e-13, alpha2 10028.5590945, alpha3_Y(R_hat) 6123.30121846, R_hat 1.44325716960e-9, m_R 7.46678054669, kappa_A + w' 4.98350253655 | A6:40-60 |
| C32 | d_target 6.53e-3 is the peer lane's transfer budget | seatT/NOTES-T.md:117; LOG:180; seatB2/NOTES-B2.md:142 |
| C33 | regularity chain lemmas L1', L2', Gj, X0-X2, I, V; F monotone and homogeneous; C_reg := F(1,1,1) | seatX/NOTES-X-PROOFS.md:37-144 |
| C34 | sampled C_reg 1996.217; sampled dry run alpha1 2.385e-20, x 2.116e-9 | seatB2/NOTES-B2.md:246-252 |
| C35 | certified kappa 14.27 vs float 4.40; C_reg 9990 vs float 1996; MET3 ratios <= 1.01 except kappa 3.2x | LEDGER:296, 297; seatE/cert_refined_a1.json:7 |
| C36 | lambda1 = 7.36696546272 (certified lower bound, rounded down) | seatE/cert_refined_a1.json:32; BAR-MET5-amend1.md:8 |
| C37 | Lemma G4; eps' = 1.0663e-6 (1.066268912509878e-06); M30 bound 7.366981173; M84 lower 7.366965462725417 | seatG45/NOTES-G45.md:26, 34-49, 91; seatG45/G4_eps_m48.json:191, 194; seatG45/g45_eps_m48.log:29 |
| C38 | Liu's constant 0.07353 inherited | seatG45/NOTES-G45.md:98; REFEREE-MET5.md:92 |
| C39 | remainder N = n1 + ... + n5; Theorem with K(s), K(0) = m_R + 5/2; inputs I1-I3 | seatQ2/NOTES-Q2-PROOF.md:9-14, 45-72 |
| C40 | Corollary 3 statement; lemma factors 1.00000022015 / 0.999999779854 at d_h | BAR-CLWT-AMEND1.md:10-14 |
| C40a | eta = \|dS\| C_inf r / S_min (1.15886280722e-8); m formula; d_sym = exp(sqrt2 m) - 1 | A6:79-82; seatB2/NOTES-B2.md:185, 189-193 |
| C41 | lambda_1/Lambda <= 1.948367728, lambda_2/Lambda >= 2.089704146, sectors >= 2.268681501 Lambda; descriptive lower end 1.94638098; d_h conservative vs d_sym 3.19e-8 | BAR-CLWT-AMEND2.md:3-6; REFEREE-CLWT.md:11-15 |
| C42 | kappa_min >= 1.9798410293 (raw, nb = 64); \|n\|^2 = 3 absent; 4 x 1.9798410293 = 7.919364 raw = 2.639788 Lambda; 2.06988 retired | REFEREE-SCHED.md:10-14, 30, 34-35, 44-46; LOG:212 |
| C43 | referee's full-cone gradient (-2.18e-87, -2.18e-87, 2.88e-87); null (1.49e-3, 1.49e-3, -1.97e-3) | REFEREE-EXTREMAL.md:22 |
| C44 | 80 sub-squares = 5 elements x 16, side 1/4, chart rules | seatE/NOTES-E.md:74-83; BAR-MET2.md:26 |
| C45 | arb Chebyshev models: degree 28, 192 bits, validated inverses delta < 1/2, exact f expansion, a* ball native | seatE/NOTES-E.md:244-256; BAR-MET2.md:26 |
| C46 | Kato's inequality, three-line proof | seatB2/NOTES-B2.md:265-269 |
| C47 | REFEREE-CM: ACCEPT WITH FIXES, no soundness bug; mutation 350 violations (inverse rem zeroed), 81/81 (tail zeroed, D = 12); fix 1 = add containment controls | seatE/REFEREE-CM.md:8-14, 41-55, 58-66; LOG:214 |
| C48 | float-disc model, (3K + 20) 2.3e-16 slack, inverses, parity pruning; mutation no_errfac 40/80; t_unit_safe exact rationals | seatE/ROUNDING-LEMMA.md:8-53, 83-90, 110-116 |
| C49 | REFEREE-ROUNDING: ACCEPT WITH FIXES; fixes under BAR-MET3-amend2 | seatE/REFEREE-ROUNDING.md:7-13, 77-84; LOG:199-200 |
| C50 | interpolation bound 16 M rho^(1-n)/(rho - 1)^2 (proved) | seatE/ROUNDING-LEMMA.md:55-66 |
| C51 | Rivlin (1.38) opened; first-kind only; even-n note wrong; switch to first-kind nodes | CITATION-LEBESGUE.md:3-12, 24-28; seatE/ROUNDING-LEMMA.md:106-107 |
| C52 | analyticity certificate on the filled polyellipse | seatE/NOTES-E.md:166-168 |
| C53 | timings: A2 START 21:50:16Z, done 23:32:12Z (11-12); B1 68 files 01:35:54Z; amend3 triage (2 A files, 10 B files), audits 0 differ, A done 03:08:36Z, B done 05:01:29Z; controls rc (unit, square, null30, mutation) | seatE/logs/prod2/status.log:1-10; seatE/logs/prod3/status.log:1-21 |
| C54 | MET4 FAILED: x = 1.2464, alpha1 5.06e-13 vs budget 4.06e-13, "~7 orders looser than the true residual"; e3[3,3] mid 1e-41, rad 1e-25 | LEDGER:297; LOG:206 |
| C55 | refine list: 9 sub-squares, shares 64.53% ... 0.51%; criterion > 0.5% | seatE/logs/refine/list.log:1-9; BAR-MET5.md:6-10, 26 |
| C56 | refinement times 458-569 s per parent; START 05:27:06Z, done 05:44:41Z; e3[3,3] lap2 1.525e-25 -> 2.229e-29 | seatE/logs/refine/status.log:1-15 |
| C57 | children tile parents exactly; vol unchanged; merge conservative | REFEREE-MET5.md:40-50; seatE/logs/refine/status.log:4-12 |
| C58 | final certificate: rho_L2 1.5268090057929563e-16, lap_rho_L2 6.856906172320414e-14, rho_H3 2.5105008366062487e-09, mu 0.48334710358870586, detm_inf 0.99892, Gjac_inf 0.24960, vol_sum and V balls | seatE/cert_refined_a1.json:3, 46-48, 53, 56, 62-63 |
| C59 | MET6: lambda1 rounded down, x = 0.177254749606, 18/18; null x = 1.775 fails K14 | LEDGER:308; BAR-MET6.md:4-13; A6:63, 88 |
| C60 | CLWM and CLWN HOLD on M30 (invariant sector; \|n\|^2 <= 2); numbers in the scorecard | LEDGER:218, 226 |
| C61 | methods: exact M/D polynomials (degrees 33, 61) + certified degree-16 surrogate (eps 5e-4); CR separator, C_CR(K) <= 0.1893 h_K (Liu 2024 Table 4.1); strip lemma; Lehmann-Maehly in strong form, no Goerisch dual space | seatM/NOTES-M.md:22-25, 36, 42-47, 97-103; REFEREE-CLWM.md:23-48; LEDGER:218 |
| C62 | CLWT: eps30 <= 1.94e-6; lemma factor; g84 values 1.948367299 / 2.089704606 / 2.268682 | LEDGER:242; BAR-CLWT-AMEND1.md:4-5 |
| C63 | referee verdicts CLWM, CLWN, CLWT, SCHED | REFEREE-CLWM.md:17-48; REFEREE-CLWN.md:12-17; REFEREE-CLWT.md:10-15; REFEREE-SCHED.md:9-15 |
| C64 | scorecard (main table and Appendix A): generated by paper/make_scorecard.py from every settle record in LEDGER whose seal BAR lies under vectors/v16_clwport/ or vectors/v17_geoflow/, plus BAR-LEDGER1. That gives 29 settled seals: 12 HOLD, 16 FAILED, 1 VOID. Verdict, seal id and settle time come from the ledger. The reasons are hand-written (make_scorecard.py REASON, scorecard_reasons.json), each checked against that seal's settle "measured" text (run with --dump). Output: paper/scorecard.md | paper/make_scorecard.py; paper/scorecard.md; LEDGER settle lines 187-308 |
| C65 | amendments sealed before reruns; superseded data unscored | BAR-MET2-amend1.md:41-42; BAR-MET2-amend2.md:15; LEDGER:253, 255, 257, 291, 304; LOG:196-205 |
| C66 | cited-assumptions list | REFEREE-MET5.md:90-97 (plus updates from CITATIONS-ABREU-DONALDSON.md:57-65; CITATION-LEBESGUE.md:1-12) |
| C67 | GEO2 transfer needs 7-jet at P1 within 1e-8, holds iff r <= ~3e-15; not met at r = 7.6e-10 | LOG:159, 201, 208; LEDGER:301 |
| C68 | X = H^6 gives C^3 only; higher jets need an interior bootstrap | SCOPE-EXTREMAL-ROUTE.md:115-122 |
| C69 | hashes: seatE/SHA256SUMS-PROD-met6.txt has 30 entries, and sha256sum -c gives 30/30 OK (rerun 2026-11-13). Its own sha256 is 7b4c5d75c1bc12b5f5fad2d86c1c7ef461945020ef98bfc00b00b82dcb62eaee. make_cert.py is 6105370f0005dd68346952270d771d77e4d3f0a2c62de6db6697e82b946b4633 (pinned under MET6 in LEDGER). result_refined_a1.json is 45c44ca0bdf4b161640fb204d729be886e71c583d901cef157261e9763b61073 (= LEDGER:308) | seatE/SHA256SUMS-PROD-met6.txt; REFEREE-PAPER.md:5 |
| C70 | commands (prod_A/prod_B/make_cert/assemble; run_refine.sh) | seatE/NOTES-E.md:313-319; BAR-MET2-amend2.md:17-18; BAR-MET5-amend1.md:12; seatE/run_refine.sh:16-39; seatB2/NOTES-B2.md:202 |
| C71 | Corollary-3-only cited inputs: X. Liu 2024 book (Table 4.1 C_CR <= 0.1893 h_K; (4.2)-(4.3); Thm 3.2; Cor. 5.1 cited for Lehmann-Maehly); certinert (no pivoting asserted, Higham Thm 9.3); comparison lemma and scheduler (refereed). D. Z. Liu (Zenodo 23598158) is a different author | seatM/NOTES-M.md:103-105; seatN/NOTES-N.md:92-94; REFEREE-CLWM.md:27-29, 46-48; REFEREE-CLWT.md:75; vectors/v13_clw/src/<private>/liu_book2024.txt:1-7; PAGE/sources/liu2026/main.txt:1-3 |
| C72 | no null attacks the refinement merge code (drop one child; sum -> max) | REFEREE-MET5.md:54 |
| C73 | Section 2 citations. dP3 authors, Platt et al. 2611.16908, Wang 2606.05519, Buttsworth-Hodgkinson 2412.01184, 2504.21644, Bunch-Donaldson 0803.0987, Hall-Murphy 1206.5489 Thm 1.1 (2.11 Lambda), Maschler AIP 1093 (2009) 132. CLW: arXiv:0705.0710v2, title and authors. LeBrun: title, J. Reine Angew. Math. issue 678, 2013, DOI 10.1515/crelle.2012.006, arXiv:1009.1270. Page range 69-94 is UNVERIFIED | PRIOR-ART-CERT.md:8, 23-35; vectors/v13_clw/src/0705.0710.txt:1-3; Crossref api.crossref.org/works/10.1515/crelle.2012.006 and arxiv.org/abs/1009.1270 (fetched 2026-11-13) |
| C74 | REFEREE-SCHED item 4: 2.639788 Lambda on M30 -> about 2.63977 Lambda on the true metric (referee arithmetic, not a certificate) | REFEREE-SCHED.md:12, 61-62 |

## Gaps found while drafting
1. The MET6 `make_cert.py` hash is not in any SHA256SUMS file. It is only in the ledger pin. A clean `sha256sum -c SHA256SUMS-PROD-refine.txt` fails on it.
2. The exact working directory and arguments of the MET6 `assemble.py` invocation are not recorded in a log header. The command in paper sec. 9 is reconstructed from BAR-MET5-amend1 and NOTES-B2 sec. 7.
3. The CLW 2008 paper is not opened on disk. The paper cited it by authors and year only at that time. (Since closed: arXiv:0705.0710v2 is cited and was opened.)
4. REFEREE-CM fix 1 (arb point-containment controls in the sealed control list) is not recorded as done in LOG.
5. The lead's brief says "Lehmann-Goerisch". The sources say Lehmann-Maehly in strong form, with no Goerisch dual space (seatM/NOTES-M.md:47; REFEREE-CLWM.md:42). The paper follows the sources.
6. BAR-CLWT-AMEND2 calls the scheduler bound UNREFEREED. REFEREE-SCHED (LOG:212) has since refereed it. The paper cites it as refereed.

## Lead additions (2026-11-13)
- Gap 1 closed: seatE/SHA256SUMS-PROD-met6.txt (refine set with make_cert.py replaced by the MET6 version, plus G4_eps_m48.json, assemble.py, chain.py, cert_refined_a1.json, result_refined_a1.json); `sha256sum -c` gives all OK.
- Gap 2 closed: the exact MET6 commands, run by the lead 2026-11-13T06:00Z (ledger 39be7ec448):
  - from seatE/: `python <memcap.py> 8 python make_cert.py prodA3 cert_refined_a1.json prodB1`
  - from seatB2/: `python <memcap.py> 8 python assemble.py ../seatE/cert_refined_a1.json result_refined_a1.json` (log: logs/assemble_refined_a1.log)
- Gap 4: REFEREE-CM fix 1 (point-containment controls into the sealed list) is NOT yet done. It is listed as an open item; it does not affect MET6 numbers (REFEREE-CM: no soundness bug).
- Gap 4 closed: BAR-CTRL1 (seal 20c4a51c45) settled HOLD 2026-11-13T07:42Z; raw counts in data/ctrl1/CTRL1.json.

## Release note (2026-11-13)
- C73: the Maschler (2009) citation was dropped from the paper (only its abstract could be opened).
- C73: the LeBrun page range 69-94 was verified against the zbMATH record (Zbl 1277.53042) and is now given in the paper.
- Section 2: Bunch-Donaldson (0803.0987), Buttsworth-Hodgkinson (2412.01184), Wang (2504.21644, 2606.05519) were opened
  and the text adjusted to them; the authors of 2611.16908 are Douglas, Gomez-Serrano, Lehmann, Platt and Qi.
- Section 9 was rewritten for this repository (code/, data/); the lab commands are quoted there.
