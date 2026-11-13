# PRIOR-ART-CERT: certified CLW metric (C1), the method (C2), scalar eigenvalue consequence (C3)

Gate run 2026-11-13 by the prior-art seat. Newest dates covered: arXiv listing to 2026-11-12 (2611.17195 and 2611.16908 seen), Zenodo API to 2026-11-12, INSPIRE citers of CLW (recid 750060) to 2026-09-26 (latest citer), Semantic Scholar citers of CLW and of LeBrun 2019.
Downloaded: `v16_clwport/paper/src/2611.16908.pdf` and `.txt` (only new file). Everything else was read from lab copies in `vectors/v13_clw/src/` and `PAGE/sources/liu2026/main.txt`.

## In plain words

Nobody has published a certified CLW metric or a certified count of its scalar eigenvalues. The method is half-taken: in September 2026 Buttsworth, Hadden, Heyes, Platt and Wiseman (2609.32655) certified the dP3 Kahler-Einstein metric with interval arithmetic and a Banach fixed point. They never touch CLW. So the safe wording is "first certified metric for a NON-Kahler, conformally Kahler Einstein metric, through its extremal (Abreu, 4th order) Kahler metric", NOT "first certified non-explicit Einstein metric". On C3, the sharper upper bound is real but small (1.94836773 against Nifa 1.951). The new content is the lower bound for lambda_2 and the "exactly one below 2 Lambda" count. Existing numbers (Liu, Hall-Murphy) agree with ours but are not proofs. Newest neighbour: 2611.16908 (12 Nov 2026, Platt and co-authors) certifies a Calabi-Yau metric on the Fermat quartic. That is not toric and not CLW, but the group is moving fast.

Correction to the brief: LeBrun "Einstein manifolds, self-dual Weyl curvature, and conformally Kahler geometry" is arXiv 1908.01881 (August 2019), not 2013. It is not in INSPIRE (math record).

## Aliases searched

CLW = "Chen-LeBrun-Weber metric" = "Einstein metric on CP2#2(-CP2)" = "conformally Kahler Einstein" = "Hermitian non-Kahler Einstein" = "Einstein-Hermitian" = "toric de Sitter saddle" = "algebraically special toric Einstein" = "Bach-flat extremal Kahler (LeBrun class)" = "two-point blow-up of P2, Einstein".
Method aliases: "computer-assisted proof", "verified numerics", "rigorous numerics", "validated numerics", "interval arithmetic", "a posteriori error bound", "Newton-Kantorovich", "Banach fixed point", "Abreu equation" = "toric extremal Kahler equation" = "symplectic potential, 4th-order".

## Table

| Item | Verdict | Best source | Label | Quoted sentence | Locator |
|---|---|---|---|---|---|
| C1: certified CLW metric (explicit approximant + proven distance to true metric) | **OPEN** | none found | n/a | No arXiv, Zenodo or INSPIRE hit proves a distance bound to CLW. Closest is Hall-Murphy: "there is, as yet, no rigorous theoretical justification for the convergence of our method". | 1407.1272 Sec. 1 (lab copy, lines 75-77) |
| C1 uniqueness hook: "the" CLW metric is well defined | TAKEN (cited, needed) | LeBrun, arXiv 1010.0238 (JDG 90) | OPENED | "or M ~ CP2#2(-CP2) and h is a constant times the metric of [8]." (Theorem A, Einstein-Hermitian metrics) | 1010.0238 p.2 |
| C1 old computer-assisted neighbour | PARTIAL, different object | Maschler, AIP Conf. Proc. 1093 (2009) 132, DOI 10.1063/1.3089199 | ABSTRACT (INSPIRE 816846) | "A computer assisted approach is presented for the cohomological uniqueness of Kaehler metrics conformal to Einstein metrics on the two point blow-up of P2." It concerns the Kahler class, not a certified metric. CLW call it "strong circumstantial evidence". | INSPIRE; CLW 0705.0710 Sec. 9 p.37-38 (OPENED) |
| C1/C2 direct neighbour: certified toric Einstein metric | **TAKEN for dP3 Kahler-Einstein; nothing for CLW** | Buttsworth-Hadden-Heyes-Platt-Wiseman, 2609.32655 (26 Sep 2026) | OPENED | "we use a combination of analytic and computer-assisted techniques to produce an approximate Einstein metric g described explicitly herein, and also prove that the true Einstein metric g_KE is close to g" | abstract; strategy list p.3-4 |
| Does dP3 do C2 or anything for CLW? | No. CLW appears once, as an example of a toric Einstein metric. | same | OPENED | "Another notable example is the Chen-Lebrun-Weber metric on CP2#2(-CP2) which was constructed using the toric assumption in [Chen et al., 2008]." No "extremal", "Bach" or conformal step in the text (grep). Its equation is the complex Monge-Ampere form, solved for a symplectic potential, then a Banach fixed point. It does need a scalar gap: "spectral gap for ... that is bigger than 4". | Sec. 1 p.3 (line 106); Sec. 1 strategy (lines 176-180); line 1752 |
| dP3 own statement of next targets | Calabi-Yau and Hermite-Einstein, not CLW | same | OPENED | "Certifying that one of these numerical candidates is near the genuine Calabi-Yau metric is a central challenge in theoretical physics." | Sec. 1 (line 175) |
| Certified Calabi-Yau (12 Nov 2026), same group | different object | Platt et al., arXiv 2611.16908 | ABSTRACT + text grep (no CLW/extremal) | "We give a rigorous, constructive proof that there exists a Ricci-flat Kahler metric on the Fermat quartic hypersurface, independently of Yau's famous theorem" | abstract |
| Other computer-assisted Einstein existence | TAKEN, other class (ODE / cohomogeneity one) | Buttsworth-Hodgkinson 2412.01184; 2504.21644; Wang 2606.05519 | ABSTRACT | 2504.21644: "The proof involves using rigorous numerics to produce an approximate Einstein metric to high precision ... which is then perturbed to a genuine Einstein metric using fixed-point methods." | abstracts |
| Numerics for extremal toric Kahler (CLW) without proof | TAKEN as numerics only | Bunch-Donaldson 0803.0987; Hall-Murphy 1407.1272 | ABSTRACT / OPENED | "We develop new algorithms for approximating extremal toric Kahler metrics. We focus on an extremal metric on CP2#2(-CP2), which is conformal to an Einstein metric (the Chen-LeBrun-Weber metric)." | 1407.1272 abstract |
| C2: certify non-explicit Einstein metric by certifying a related Kahler metric | **PARTLY TAKEN** | dP3 2609.32655 (direct Kahler case) | OPENED | The reduction "Einstein equation as scalar Kahler equation in symplectic coordinates, then interval arithmetic + Banach fixed point" is dP3's recipe. Not taken: (i) 4th-order Abreu/extremal equation, (ii) LeBrun's conformal step to a non-Kahler Einstein metric, (iii) Newton-Kantorovich on it. No source found for (i)-(iii). | dP3 strategy list; lab search floor below |
| C2 sub-claim "avoids any tensor spectral gap" | Not distinctive vs dP3 | dP3 | OPENED | dP3 also needs no tensor gap, only a scalar gap (finite-element Laplacian, "Steps 1 and 2: FEM and smooth spectral gaps"). Do not present "no tensor gap" as the novelty. | dP3 Sec. 3.4 |
| C3 upper bound lambda_1 | **PARTLY TAKEN (sharpened)** | Nifa, Zenodo 22851442 (2026-09-01) | OPENED (abstract and text) | "the first positive scalar eigenvalue satisfies 4Lambda/3 < lambda_1^sc(g) < 1951Lambda/1000." Ours: lambda_1/Lambda <= 1.94836773. Same units (Ric = Lambda g). 1.951 - 1.94836773 = 0.00263 (0.13%). | Nifa abstract; eq (1.9) line 339 of lab text |
| C3 older rigorous upper bound | TAKEN, weaker | Hall-Murphy 1206.5489 Thm 1.1 | OPENED | "4/3 Lambda < lambda_1 <= 2.11 Lambda" (also quoted as 2.107 in 1407.1272 Sec. 5). Ours is below the 2Lambda threshold; theirs is not. | 1206.5489 p.1; 1407.1272 Sec. 5 |
| C3 lower bound lambda_1 > 4/3 | TAKEN (classical) | Lichnerowicz-Obata | OPENED (HM) | "The lower bound of 4/3 in Theorem 1.1 is just the classical Lichnerowicz-Obata lower bound" | 1206.5489 p.2 |
| C3 numerical lambda_1 | TAKEN as numerics, consistent | Liu, Zenodo 23598158 (2026-05-23) | OPENED | "our converged first eigenvalue, lambda_1 = 1.94738Lambda, lies below the bound 1.951Lambda proved by Nifa". 1.94738 < 1.94836773: consistent. | main.txt lines 160-162 |
| C3 numerical odd-sector Rayleigh quotient | TAKEN as numerics, consistent | Hall-Murphy 1407.1272 Sec. 5 | OPENED | "values for the Rayleigh quotients of ... 2.0940 and ... 1.9481". Even sector 2.0940 is an upper bound on the approximate metric. Ours lambda_2 >= 2.08970414 < 2.0940: consistent. Also "the value 2.0940 is very close to the value 2.0965 found by the authors in [19]". Liu: Rayleigh quotient 1.9480 (his metric), 1.9481 (theirs). | 1407.1272 p.13 (lines 862-870) |
| C3 lower bound for lambda_2 (second scalar eigenvalue) | **OPEN** | none | n/a | Nifa's text (grep for "second", "exactly one", "at most one") has no lambda_2 bound for CLW. Liu and Hall-Murphy give numbers only, no enclosure. | Nifa lab text; Liu main.txt |
| C3 "exactly one scalar eigenvalue in (0, 2 Lambda)" as a theorem | **OPEN** | none | n/a | Liu mentions "the scalar eigenvalues below 2Lambda" only in a numerical check (main.txt line 229). Caveat: our all-sector coverage rests on the lab's own proof, which this gate did not audit. | Liu line 229 |
| Related Nifa result on a different CLW quantity | context | Nifa, Zenodo 22914003 (2026-09-03) | ABSTRACT | "gives, for the Chen-LeBrun-Weber metric at Einstein constant 3, the inequality 207/100 < sup D_k(s^2) < 21/10" (pointwise, on the conformal factor; no certified metric) | Zenodo abstract |
| Instability itself | TAKEN | Biquard-Ozuch, 2310.10109 (JDG 133) | OPENED | "the new case here is the Chen-LeBrun-Weber metric" | p.2 after Thm A |
| Torus-invariant nondegeneracy | TAKEN | Li-Sun 2610.10425v2 Thm 7.4 | OPENED | "All the metrics listed in Corollary 7.3 are torically nondegenerate." Li-Sun cite dP3 as "a computer-assisted approximation of the dP3 metric with rigorous error estimates", no CLW. | Sec. 7; Sec. 1 line 117 |
| LeBrun 2019 citers (1908.01881) | no certified work | Semantic Scholar: 2007.01180; Wu, IMRN 2019 | SNIPPET (titles only) | "Einstein metrics, conformal curvature, and anti-holomorphic involutions"; "Einstein Four-Manifolds With Self-Dual Weyl Curvature of Nonnegative Determinant" | citer list |
| INSPIRE citers of CLW | none certifies CLW | INSPIRE refersto:recid:750060 (20 hits, plausible) | SNIPPET (titles read) | Newest: 2609.32655 (dP3), 2607.21538 (Conformal Kahler rigidity), 2606.24761 (Chen-Teo with Lambda), 2606.10443 (survey). None computes or certifies CLW. | INSPIRE list |

## Verdicts

- C1: OPEN. Wording rule: "first certified metric of a conformally Kahler, non-Kahler Einstein metric via its extremal Kahler metric" is supportable at this search floor. "First certified non-explicit Einstein metric" is TAKEN (dP3, and the ODE-class papers).
- C2: PARTLY TAKEN. dP3 owns the recipe "interval arithmetic, symplectic potential, Banach fixed point" for a Kahler-Einstein metric. The 4th-order Abreu equation plus the conformal step is not found elsewhere. Do not call the lack of a tensor gap the novelty.
- C3: PARTLY TAKEN. The lambda_1 upper bound is a sharpening of Nifa 1951/1000 to 1.94836773 (0.13% better). The lambda_2 lower bound and the exact count are OPEN. All numbers are consistent with Liu and Hall-Murphy numerics (1.94738 < 1.94836773; 2.08970414 < 2.0940).

## Normalisation check

All sources use Ric = Lambda g and the scalar Laplacian eigenvalue in units of Lambda (Nifa abstract, Hall-Murphy Thm 1.1, Liu "Lambda = 3" units divided back out as lambda/Lambda). No rescaling needed. Hall-Murphy 2.0940 is a Rayleigh quotient on a degree-20 approximation, so it bounds the true value only if the approximation were exact. Compare it only as plausibility.

## Search floor (not proof of absence)

- arXiv API: all/abs queries for "Chen-LeBrun-Weber" (6 hits, all known), "conformally Kahler" AND Einstein, "computer-assisted" AND Einstein (7, all read by title; 2609.21797 is a black hole in Einstein-Weyl gravity, unrelated), extremal AND "computer-assisted" (46, none relevant), "Abreu" AND toric AND rigorous (0), "Kantorovich" AND Kahler (1, 2014 optimal transport), "interval arithmetic" AND Einstein (0), "Bach-flat" AND Kahler AND numerical (0); authors Buttsworth, Wiseman, Platt, Ozuch, Hall (titles read). Several exact-phrase combinations returned 0, which is an instrument limit, not evidence.
- Zenodo API: "Chen-LeBrun-Weber" (2: Liu, Nifa), "Chen LeBrun Weber" (1), "conformally Kahler" (hits irrelevant except Nifa), "Abreu equation" (0), "toric Einstein" (0), "Newton-Kantorovich" AND Kahler (0), "extremal Kahler" AND rigorous (0), "Bach-flat" AND Kahler (0), "computer-assisted proof" AND Einstein (5, none relevant), exact constants "1.94836773" (0) and "2.08970414" (0), author Nifa (24+ records; two on CLW: 22851442, 22914003). No Zenodo record certifies a metric.
- INSPIRE: CLW citers (20), date > 2026 sub-query returned 0 (instrument check: the 2026 papers were in the unfiltered list, so the date syntax failed; list was read in full instead).
- Web: three standard and extended searches ("certified extremal Kahler", "validated numerics Abreu equation", "CLW eigenvalue bound"). Nothing new.
- Not covered: paywalled journal pages, MathSciNet, conference talks.
- Reference sweep: dP3, Li-Sun, Nifa, Liu reference lists searched for computer-assisted / rigorous / CLW; leads opened: Maschler (abstract), 2611.16908, Nifa 22914003.

## Most dangerous neighbour

The dP3 group (Buttsworth, Hadden, Heyes, Platt, Wiseman). They hold the pipeline, say certified Calabi-Yau is next, and Platt posted 2611.16908 on 12 Nov 2026. A certified CLW is a natural follow-up: Wiseman co-authored early numerics of this kind, and Hall-Murphy list it as the rigorous gap. Second: Nifa (prolific; already sharpening CLW bounds analytically on Zenodo, so a lambda_1 improvement or a lambda_2 bound by analytic means is plausible). Third: Liu (numerical, could add enclosures).

## Method worth borrowing

- dP3: cancel the facet degeneracy of the inverse Hessian algebraically before interval arithmetic; keep an explicit Banach-fixed-point theorem (their Theorem 4.1) with the three-number test, then compare your Newton-Kantorovich constants against it as a cross-check.
- 2611.16908: its "verified numerics" bounds for eigenvalues of Hermitian matrices (Lemma D.26) and chart-by-chart box covering may port to the polygon.
- Li-Sun Thm 7.4: torus-invariant kernel already ruled out; cite instead of re-proving.

## The gap that is genuinely open

No source proves, with an explicit enclosure, the CLW metric itself (via its extremal Kahler metric and LeBrun's conformal change) or that exactly one scalar Laplace eigenvalue lies below 2 Lambda for it.
