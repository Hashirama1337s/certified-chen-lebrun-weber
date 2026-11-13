# A certified approximation of the Chen-LeBrun-Weber Einstein metric

**Moki & Julio**

Version 1.0. Every number below is listed with its source file and line in `CLAIMS.md` (same folder; those sources are
the lab's working files, most of which are not in this repository). Claim tags like [C12] point to rows of that file.

---

> **In plain words**
>
> The Chen-LeBrun-Weber (CLW) metric is a special shape of a four-dimensional space: an Einstein metric, the kind of
> geometry that solves Einstein's equations in empty space with a cosmological constant (a curved space with no time
> direction). It lives on the manifold
> $`\mathbb{CP}^2\#2\overline{\mathbb{CP}^2}`$. Its existence was proved in 2008, but nobody can write it down by a formula.
>
> We computed a very accurate approximate version of it on a computer. Then we used interval-style arithmetic (every
> number carries a guaranteed error bar) to prove that the true metric sits within a tiny, explicit distance of our
> computed one. At every point, the two metrics differ by a factor of at most $`1 \pm 5.504\times10^{-8}`$ (about 5 parts
> in 100 million). This proof assumes several published results that we cite (listed in Section 7). One of them, the one
> that names the metric as CLW, is a preprint, and one boundary-condition step is not yet fully sourced.
>
> The existence of the metric is not new. What is new is the certified location: an explicit, checkable statement of
> where the metric is. A second result follows from it: the CLW metric has exactly one scalar vibration mode below the
> threshold $`2\Lambda`$, where $`\Lambda`$ is the curvature constant of the space.
>
> We also list every pre-registered test we ran on the way, including the 16 that failed.

---

## 1. Main result

**Setting.** Let $`p(a)`$ be the degree-10 polynomial

```math
p(a)=a^{10}+5a^9+6a^8-12a^7-34a^6-42a^5+20a^3+21a^2-7a-6,
```

and let $`a^*`$ be its unique root in $`(1,\infty)`$ [C1]. The value $`a^* = 1.95771280518776186068\ldots`$ is carried in the
computation as a FLINT-isolated ball of radius $`1.5\times10^{-40}`$ [C2]. Let $`P=P(a^*)`$ be the Delzant pentagon
$`[-1,1]^2\cap\{x_1+x_2\ge-\gamma\}`$ with $`\gamma=2/a^*`$ [C3], and let $`M=\mathbb{CP}^2\#2\overline{\mathbb{CP}^2}`$ be the toric
manifold it defines. Write $`\tilde u=u_G+f_{84}`$ for the computed symplectic potential: the Guillemin potential plus a
degree-84 polynomial correction $`f_{84}`$ [C4]. Write $`\tilde g=g_K(\tilde u)`$ for its toric Kähler metric, $`\tilde x`$
for its moment map, and $`S`$ for Donaldson's affine function of $`P`$. Sobolev norms $`H^k`$ are sum norms
$`\|\psi\|_{H^k}=\sum_{j\le k}\|\nabla^j\psi\|_{L^2(M,\tilde g)}`$ on $`T^2`$-invariant functions [C5].

**Theorem 1 (certified extremal metric; computer-assisted; conditional on the cited inputs of Section 7).** For every
class $`a`$ in the $`a^*`$ ball there is a $`T^2`$-invariant function $`\psi`$ with

```math
\|\psi\|_{H^6(M,\tilde g)}\le 7.568\times10^{-10}
```

such that the Kähler metric $`g_\psi=\tilde g+2(\nabla^2\psi)^+`$ of $`\tilde\omega+dd^c\psi`$ is extremal, with scalar
curvature exactly $`S(x_\psi)`$ [C6, C7, C18]. Here $`dd^c=2i\partial\bar\partial`$, and the radius above is stated for this
normalisation (factor 2).

This is a sum norm of a complex-gauge Kähler potential. It is **not** a distance between symplectic potentials
$`u^*-\tilde u`$ [C7].

**Corollary 2 (location of the CLW metric; conditional on the cited inputs of Section 7).** At the critical class $`a_c=a^*`$ (inside the ball), the metric
$`h=S(x_\psi)^{-2}g_\psi`$ is Einstein. Up to scaling and isometry it is the Chen-LeBrun-Weber metric. Pointwise on all of $`M`$,

```math
(1-d_h)\,\tilde h\;\le\;h\;\le\;(1+d_h)\,\tilde h,\qquad d_h\le 5.504\times10^{-8},
```

where $`\tilde h`$ is the same conformal rescaling of $`\tilde g`$ (complex gauge: same complex structure, identity map). In the
symplectic gauge (same coordinates $`(x,\theta)`$ on the closed polygon) the bound is $`d_{\rm sym}\le 3.19\times10^{-8}`$ [C8, C9].

The three inputs of Corollary 2 are:
- LeBrun, Crelle 678 (2013), Prop. 3: at a critical point $`W`$ of his functional $`\mathcal A`$ on the Kähler cone, an extremal
  Kähler metric with $`s>0`$ is conformal to an Einstein metric via $`h=s^{-2}g`$ [C10];
- criticality of $`a^*`$ on the full cone (Lemma G5, Section 2.6);
- Li-Sun, arXiv:2610.10425v2, Cor. 7.3 and Prop. 7.13: a toric Einstein metric on this manifold is CLW up to scaling and
  isometry. **This is a preprint** [C11]. Independently, LeBrun's Theorem A (arXiv:1010.0238, p. 1) states that a compact
  Hermitian Einstein metric on $`\mathbb{CP}_2\#2\overline{\mathbb{CP}}_2`$ is either Kähler-Einstein or a constant multiple of the
  CLW metric. The manifold carries no Kähler-Einstein metric (standard: by Matsushima's criterion its automorphism group is not reductive), and $`h=S^{-2}g_\psi`$ is Hermitian for $`J`$ and not Kähler (its
  conformal factor is non-constant), so the identification also follows from that theorem.

The positivity $`s>0`$ holds because $`S\ge S_{\min}\ge 3.0095`$ on $`P`$ [C12].

**What is new and what is not.** The existence of the CLW metric is known (Chen, LeBrun and Weber, 2008; a second proof is
LeBrun's 2013 paper [C13]). This paper gives a computer-assisted re-proof of existence, plus a certified location. That
location is an explicit polynomial potential together with explicit $`H^6`$ and pointwise error bars. The extremal metric and
its $`H^6`$ bound hold for every class in the $`a^*`$ ball. The Einstein metric and $`d_h`$ hold at the critical class [C14]. The method is the second contribution: a Newton-Kantorovich argument for the extremal
(Abreu) equation, closed with certified inputs over the whole closed polygon, rods and vertices included.

**Corollary 3 (scalar gap; companion lane, Section 5; conditional on the cited inputs of Section 7).** The scalar Laplacian of the CLW metric has exactly one eigenvalue
in $`(0,2\Lambda)`$, counting all $`T^2`$-sectors. Here
$`\lambda_1/\Lambda\le1.948367728`$ and $`\lambda_2/\Lambda\ge2.089704146`$ [C40-C42].

---

## 2. Prior work and novelty

(Gate: `PRIOR-ART-CERT.md`, arXiv/Zenodo to 2026-11-12, INSPIRE citers of CLW to 2026-09-26. Citations: [C73].)

**Existence.** The Chen-LeBrun-Weber metric exists (X. Chen, C. LeBrun and B. Weber, "On conformally Kähler, Einstein manifolds",
arXiv:0705.0710, published 2008). A second proof is LeBrun's ("Einstein manifolds and extremal Kähler metrics", J. Reine Angew. Math.
678 (2013), DOI 10.1515/crelle.2012.006, arXiv:1009.1270). That proof is a continuity argument whose openness step is a
qualitative inverse-function theorem. We give no new existence theorem. What is new is a certified location of the metric and the
route to it.

**Numerical CLW.** Bunch and Donaldson (arXiv:0803.0987) computed numerical approximations to extremal toric Kähler
metrics, including the pentagon case, whose extremal metric is conformal to CLW. Hall and Murphy computed CLW numerically (arXiv:1407.1272) and note
"there is, as yet, no rigorous theoretical justification for the convergence of our method". No source we found certifies the CLW metric itself.

**Certified Einstein metrics.** Certified non-explicit Einstein metrics already exist:
- the Kähler-Einstein metric on dP3, by Buttsworth, Hadden, Heyes, Platt and Wiseman (arXiv:2609.32655). They use an interval
  a-posteriori error plus a Banach fixed point on a symplectic potential, with only a scalar spectral gap;
- metrics of cohomogeneity one: Buttsworth-Hodgkinson (arXiv:2412.01184) and Wang (arXiv:2504.21644, arXiv:2606.05519);
- very recently, a Calabi-Yau metric on the Fermat quartic, by Douglas, Gómez-Serrano, Lehmann, Platt and Qi (arXiv:2611.16908).

We borrow the form of the dP3 fixed-point theorem and its Appendix-B Sobolev constants (Section 3). Our case
differs in four ways. CLW is Hermitian but not Kähler. We certify its conformally related extremal Kähler metric (the fourth-order
Abreu equation, solved by Newton-Kantorovich with a bordered operator) and pass to the Einstein metric by LeBrun's conformal
criterion at the critical class. To our knowledge, at the search floor of our prior-art gate, this is the first certified non-explicit, compact, conformally Kähler (non-Kähler) Einstein metric (the Page metric is of this type but explicit). We
do not claim it is the first certified non-explicit Einstein metric.

**Spectrum.** Nifa (Zenodo 22851442) proved $`4\Lambda/3<\lambda_1^{sc}<1951\Lambda/1000`$ for CLW. Corollary 3's bound
$`\lambda_1\le1.94836773\Lambda`$ is sharper by 0.13%. The lower bound $`\lambda_2\ge2.0897\Lambda`$ and the statement "exactly one
scalar eigenvalue in $`(0,2\Lambda)`$" were not found in any source. Numerical values agree with ours: D. Z. Liu (Zenodo 23598158; not the X. Liu of the 2024 constants) gives
$`\lambda_1\approx1.94738\Lambda`$, and Hall-Murphy give Rayleigh quotients 1.9481 and 2.0940. These are not proofs.
The earlier rigorous upper bound is Hall-Murphy (arXiv:1206.5489, Thm 1.1): $`4\Lambda/3<\lambda_1\le2.11\Lambda`$, which does not
reach below $`2\Lambda`$.

---

## 3. Method

### 3.1 The equation and the bordered operator
We work in the complex gauge. For $`\psi\in H^6`$ the metric $`g_\psi=\tilde g+2(\nabla^2\psi)^+`$ is the Kähler metric of
$`\tilde\omega+dd^c\psi`$. It lies in the same class, and its moment map is $`x_\psi=\tilde x+d_y\psi`$ [C18]. The gauge directions
are $`A=\mathrm{span}\{1,\tilde x_1,\tilde x_2\}`$, so we use the slice $`X=\{\psi\in H^6:\psi\perp_{L^2(M,\tilde g)}A\}`$ [C18].

The map is

```math
G(\psi,c)=\mathrm{Scal}(g_\psi)-S(x_\psi)+c\cdot e(\tilde x),\qquad (\psi,c)\in X\times\mathbb R^3,\ G\in Y=H^2,
```

with linearisation $`DG(0)(\psi,c)=-L_t\psi+E\psi+c\cdot e`$. Here $`L_t=\Delta^2+2\,\mathrm{Ric}\cdot\nabla^2+\langle ds,d\cdot\rangle`$
and $`E\psi=\langle d\rho,d\psi\rangle`$ [C18]. The residual is $`\rho=\mathrm{Scal}(\tilde g)-S(\tilde x)`$.

- **Bordered inverse.** On $`Y`$ we use the equivalent norm $`\|y\|_Y=\|(1-\Pi)y\|_{H^2}+w'\|\Pi y\|_{L^2}`$, with $`w'=2^{-10}`$.
  Then $`\|B_0^{-1}\|\le C_{\rm reg}`$. The price is the norm equivalence $`\|y\|_Y\le(\kappa_A+w')\|y\|_{H^2}`$, with
  $`\kappa_A=1+k_1+k_2`$ built from Gram matrices of $`e`$ and its derivatives [C19].
- **Donaldson's identity** gives $`\Pi\rho=0`$ exactly, so $`\alpha_1=\|\rho\|_Y=\|\rho\|_{H^2}`$ [C20]. The identity is checked
  against the primary source (arXiv:math/0407486v2, Cor. 1 with Prop. 1 and Def. 1), including the factor 2 and the lattice
  boundary measure in our normalisation [C21]. It is also verified numerically [C22].
- **No Lagrange multiplier at the fixed point.** Donaldson's identity applied to $`g^*`$ forces $`c^*=0`$ once
  $`\|G_P^{-1}\|\,\|\Delta\|_F<1`$ (gate 15) [C23].

### 3.2 The Newton-Kantorovich theorem as used
The fixed point is that of $`T(z)=-DG(0)^{-1}(\rho+N(z))`$ on the closed $`r`$-ball. Our checks are exactly the two
contraction-mapping inequalities

```math
\alpha_2\alpha_1+\alpha_2\alpha_3 r^2\le r,\qquad 2\alpha_2\alpha_3 r<1.
```

So the argument is the plain contraction-mapping theorem and does not rest on the cited form of dP3's Thm 4.1 [C24]. Here
- $`\alpha_1=\|\rho\|_{H^2}`$;
- $`\alpha_2=C_{\rm reg}/(1-t)`$, with $`t=C_{\rm reg}(\kappa_A+w')C_{32}\|\rho\|_{H^3}`$ (moment-map shift) [C25];
- $`\alpha_3`$ is the Lipschitz constant of the remainder (Section 3.4).

The closing quantity is $`x=4\alpha_1\alpha_2^2\alpha_{3,Y}<1`$ (strict), and the radius is
$`r=2\alpha_1\alpha_2/(1+\sqrt{1-x})`$ [C26].

**All 18 gates, with values** (Arb, 256 bits; a gate passes only if Arb proves it) [C27-C30]:

| # | gate | value | result |
|---|---|---|---|
| 1 | $`\mu>0`$ (Ricci lower bound, needed for the Sobolev constants) | $`\mu\ge0.4833471035`$ | PASS |
| 2 | $`R_m\ge0`$, and $`R_m=0`$ only with $`\mu>0`$ | $`R_m=0`$ | PASS |
| 3 | $`\mu\le R_p\le R_{op}`$ not violated | | PASS |
| 4 | $`S_{\min}>0`$ | $`\ge3.0095167079`$ | PASS |
| 5 | $`V\le(2\pi)^2 G_{P,00}`$ | $`V=139.017997713`$ | PASS |
| 6 | $`C_{\rm reg}`$ finite | $`9989.92845835`$ | PASS |
| 7 | $`t<1`$ | $`0.00385206247051`$ | PASS |
| 8 | $`C_9\cdot2r_{\max}<1`$ | $`9.69498\times10^{-7}`$ | PASS |
| 9 | $`C_9\cdot2\hat R<1`$ (hypothesis of the $`\alpha_3`$ theorem) | | PASS |
| 10 | $`x<1`$ strict | $`x=0.177254749606`$ | PASS |
| 11 | $`r\le\hat R`$ ($`\alpha_3`$ evaluated on a larger ball) | $`r=7.568\times10^{-10}`$, $`\hat R=1.443\times10^{-9}`$ | PASS |
| 12 | self-map $`\alpha_2\alpha_1+\alpha_2\alpha_{3,Y}r^2\le r`$ | | PASS |
| 13 | contraction $`2\alpha_2\alpha_{3,Y}r<1`$ | $`0.0929469`$ | PASS |
| 14 | $`b=\sqrt2\,C_\infty r<1`$ ($`g_\psi`$ positive) | $`3.186\times10^{-8}`$ | PASS |
| 15 | $`c^*=0`$: $`\Vert G_P^{-1}\Vert \Vert \Delta\Vert _F<1`$ | $`2.859\times10^{-7}`$ | PASS |
| 16 | $`d_{cx}\le d_{\rm target}=6.53\times10^{-3}`$ | $`3.186\times10^{-8}`$ | PASS |
| 17 | $`d_{\rm sym}\le d_{\rm target}`$ | $`3.186\times10^{-8}`$ | PASS |
| 18 | $`d_h\le d_{\rm target}`$ | $`5.504\times10^{-8}`$ | PASS |

Intermediate values [C31]:
- $`\alpha_1=7.19573547909\times10^{-14}`$, against a budget of $`4.05953742741\times10^{-13}`$;
- $`\alpha_2=10028.5590945`$;
- $`\alpha_{3,Y}(\hat R)=6123.30121846`$;
- $`\hat R=1.44325716960\times10^{-9}`$.

The target $`d_{\rm target}=6.53\times10^{-3}`$ is the companion lane's budget for transferring its spectral bounds (Section 5) [C32].

### 3.3 The regularity chain $`C_{\rm reg}=F(1,1,1)`$
We need $`\|\psi\|_{H^6}\le C_{\rm reg}\|L_t\psi\|_{H^2}`$ on $`X`$. The chain nests two second-order scalar chains, built from:
- (L1') a commutator formula for the rough Laplacian with $`\nabla`$ on tensors of any rank;
- (L2') integrated Bochner;
- (Gj) a telescoped bound for $`\Delta\nabla^j s`$;
- short lemmas X0-X2, I and V, which use $`\lambda_1`$, the contracted Bianchi identity and the trace split of Ric [C33].

Every step is a nondecreasing, positively homogeneous function of earlier bounds. So $`F(w)`$ is monotone and 1-homogeneous,
and $`\|\psi\|_{H^6}\le F(W,W,W)=W\,F(1,1,1)`$. No grid supremum is needed [C33].

Certified value: $`C_{\rm reg}=9989.92845835`$ [C27]. It is about five times the sampled value 1996.2 [C34]. The main cause is the
certified Weitzenböck norm, $`\kappa\le14.27`$ against a float value of $`4.40`$ [C35]. The spectral input is
$`\lambda_1(L_t)\ge7.36696546272`$ on $`A^\perp`$ [C36]. It is obtained from a certified M30 bound $`7.366981173`$ by a
Rayleigh-quotient comparison with $`\epsilon'=1.0663\times10^{-6}`$ (Lemma G4) [C37]. That M30 bound rests on Liu's cited
constant 0.07353 [C38].

### 3.4 The complex-gauge remainder $`\alpha_3`$
With $`B=g^{-1}\,2(\nabla^2\psi)^+`$, $`P=(1+B)^{-1}`$ and $`Z=BP`$, the scalar-curvature remainder is an exact sum
$`N=n_1+\dots+n_5`$ [C39]. Theorem (NOTES-Q2-PROOF sec. 3): for $`s=2r`$ with $`C_9s<1`$ and $`\gamma=1/(1-C_9 s)`$,

```math
\|N(\psi)-N(\psi')\|_{H^2}\le 4C_{32}K(s)\,(\|\psi\|_{H^6}+\|\psi'\|_{H^6})\,\|\psi-\psi'\|_{H^6},
```


```math
K(s)=(m_R+1)\gamma^2+\gamma^3(1+\tfrac{C_{32}s}{4})+\tfrac12\gamma^2(1+C_{32}\gamma s)+\tfrac12C_{32}s\gamma^3(1+\gamma+C_{32}\gamma s),
```

with $`K(0)=m_R+5/2`$ and $`m_R=R_{op}+2Rc_1+Rc_2`$ [C39]. In the $`Y`$ norm, $`\alpha_{3,Y}=(\kappa_A+w')\,4C_{32}K(2\hat R)`$ [C25].

The proof uses only:
- a sup-norm embedding $`\|T\|_\infty\le C_\infty\|T\|_{H^3}`$, $`C_\infty=4C_7C_6C_8C_5`$;
- an $`L^4`$ embedding with constant $`C_5`$;
- Kato's inequality, which has a three-line proof in the notes [C39, C46].

$`C_5..C_8`$ come from dP3 App. B as transcribed, from $`\mu`$ and $`V`$ only [C29]. Values: $`C_\infty=29.7674049667`$,
$`C_{32}=30.8202310004`$, $`C_9=59.5348099334`$, $`m_R=7.46678054669`$ [C29, C31].

### 3.5 The pointwise metric bound
- **Complex gauge:** $`d_{cx}=\sqrt2\,C_\infty r`$, because $`B`$ commutes with $`J`$.
- **Symplectic gauge:** a Moser flow $`V_t=-\mathrm{grad}_{g_t}\psi`$ preserves $`\theta`$ and the moment map. Gronwall then gives
  $`d_{\rm sym}=e^{\sqrt2 m}-1`$, with $`m\le C_\infty r(1+3C_\infty r/(1-b))/(1-b)`$.
- **Einstein rescaling:** $`d_h=(1+d_{cx})/(1-\eta)^2-1`$, with $`\eta=|dS|_{\tilde g}C_\infty r/S_{\min}`$ [C9, C40a].

### 3.6 Criticality on the full cone (Lemma G5)
- $`\mathcal A`$ is invariant under the lattice swap $`\iota(p,q,r)=(q,p,r)`$.
- At $`c(a^*)=(a^*,a^*,a^*+1)`$, three tangent vectors with $`\det=2`$ each have $`d\mathcal A=0`$:
  - the family direction ($`p(a^*)=0`$);
  - the radial direction (Euler, degree-0 homogeneity);
  - the anti-invariant direction (swap symmetry).
- Also $`D(a^*)\ne0`$, since $`D(1+z)`$ has positive coefficients [C17].

A referee checked this independently on the full three-dimensional cone and found
$`\mathrm{grad}\,\mathcal A=(-2.18,-2.18,2.88)\times10^{-87}`$ at truncation level, while the null at $`1.001a^*`$ gives order $`10^{-3}`$ [C43].

### 3.7 The rigorous engine
The constants are certified on 80 sub-squares covering the closed pentagon: 5 vertex elements, each with 16 sub-squares of
side 1/4, in polar, Cartesian-vertex, rod or symplectic charts [C44].

- **Part A, scalar layer (arb Chebyshev models).** Bivariate Chebyshev models with arb coefficients, degree 28, 192 bits. Each
  product is one ball-polynomial multiplication. Truncation tails go into one sup remainder. Inverses are validated a
  posteriori ($`\delta=\|1-qF\|<1/2`$). $`f_{84}`$ is expanded exactly, and the $`a^*`$ ball is carried natively. This layer gives
  $`\rho`$, $`\Delta\rho`$, the Ricci eigenvalue bounds (hence $`\mu`$), $`S`$-derivatives, Gram integrals and admissibility [C45].
  - An independent referee found no soundness bug [C47].
  - A mutation test showed that the original square control cannot detect remainder bugs, while an arb point-containment check
    does (350 violations with the inverse remainder zeroed; 81/81 points violated with the tail zeroed at degree 12) [C47].
  - That check was then sealed as a control (BAR-CTRL1) and run on the pinned production code: 0 violations in 1620 square
    checks (degrees 28 and 12) and 0 in 98 CLW checks against the independent arb evaluator; the two mutations fail as
    required (350 and 810 violations). Raw counts: `data/ctrl1/CTRL1.json`.
- **Part B, tensor constants (float-disc engine).** Complex-disc arithmetic in binary64 with a written rounding model.
  - Each fused sum of $`K`$ products gets slack $`(3K+20)\cdot2.3\times10^{-16}`$.
  - Validated order-0 inverses, exact jet recursions, and parity pruning by the chart's joint reflection [C48].
  - Mutation: dropping the $`(3K+20)u`$ terms breaks containment (40/80), so that guard is load-bearing.
  - The per-operation slack is tested separately by exact-rational unit tests [C48].
  - Refereed: ACCEPT WITH FIXES. The three fixes went in under BAR-MET3-amend2 before production [C49].
- **Interpolation bound** (proved, not cited): for $`F`$ analytic on the polyellipse with $`|F|\le M`$,
  $`|F-IF|\le16M\rho^{1-n}/(\rho-1)^2`$ [C50].
- **Node errors** use the Lebesgue constant $`\Lambda_n\le(2/\pi)\log(n+1)+1`$ for first-kind nodes, Rivlin, *The Chebyshev
  Polynomials*, 2nd ed., eq. (1.38), opened [C51]. A referee caught that second-kind nodes at even $`n`$ fall outside (1.38).
  Production was switched to first-kind nodes [C51].
- **Analyticity** is certified on the whole filled polyellipse by a cover of complex boxes. This bounds $`\mathrm{Re}\,A`$,
  $`\mathrm{Re}\,B`$, $`|l_r|`$ and $`|\det m|`$ below [C52].

---

## 4. The computation

- **Production A** (scalar layer): 80 sub-squares, 5 processes. Code under BAR-MET2-amend2, run 2026-11-12 21:50Z to 23:32Z.
  Two files were recomputed after the amend3 fail-safe, finishing 03:08Z on 11-13 [C53].
- **Production B** (tensor layer): 80 sub-squares. 68 finished in the first night run. The rest finished after fail-safe fixes
  (BAR-MET3-amend3/4): an infinite/NaN crude bound now means "split or switch chart", never a silent value. All 80 were done
  by 05:01Z [C53].
  - The NaN audit found 0 differing files [C53].
  - Certified tensor constants are within $`1.01\times`$ the float node maxima, except $`\kappa`$ ($`3.2\times`$) [C35].
- **First assembly (MET4) did not close:** $`x=1.246`$ [C54]. The cause was enclosure width, not the residual:
  - 64.5% of the $`\int(\Delta\rho)^2`$ bound came from one far-corner sub-square, e3 [3,3];
  - there the bound was set by the model remainder (midpoint $`\sim10^{-41}`$, radius $`\sim10^{-25}`$) [C54, C55].
- **Far-corner refinement (MET5):**
  - The 9 sub-squares holding more than 0.5% each of the $`\int(\Delta\rho)^2`$ bound were re-certified as $`2\times2`$ children
    with the same degree, precision and ball [C55].
  - Each parent took 458-569 s. Wall time was 05:27:06Z to 05:44:41Z [C56].
  - On e3 [3,3] the bound fell from $`1.525\times10^{-25}`$ to $`2.229\times10^{-29}`$ [C56].
  - The children tile their parents exactly, and the volume pieces are unchanged [C57].
- **Final certificate:**
  - $`\|\rho\|_{L^2}\le1.527\times10^{-16}`$, $`\|\Delta\rho\|_{L^2}\le6.857\times10^{-14}`$, $`\|\rho\|_{H^3}\le2.511\times10^{-9}`$;
  - $`\mu\ge0.48334`$;
  - $`\det m\ge0.99891`$ and Jacobian $`\ge0.24960`$ (admissibility) [C58];
  - the volume pieces sum to $`[139.017997713132\pm1.93\times10^{-13}]`$ against
    $`V=[139.017997713132\pm1.49\times10^{-13}]`$ [C58].
- **After the referee's $`\lambda_1`$ fix (MET6):** the input had been rounded up, and is now rounded down to 7.36696546272.
  The assembly gives $`x=0.177254749606`$ and 18/18 gates [C59].

The assembly itself takes seconds (memory peak 0.01 GB) [C30].

---

## 5. Corollary: the scalar gap of the true CLW metric (companion lane)

The following is the companion spectral lane's result. Its records are BAR-CLWT-AMEND1 and -AMEND2, and that lane will
review this section.

**Corollary 3.** The scalar Laplacian of the CLW metric has exactly one eigenvalue in $`(0,2\Lambda)`$, counting all
$`T^2`$-sectors [C40]. With bounds rounded the safe way (referee-recomputed, outward):
- $`\lambda_1/\Lambda\le1.948367728`$ [C41];
- $`\lambda_2/\Lambda\ge2.089704146`$ [C41];
- every non-invariant $`T^2`$-sector has lowest eigenvalue $`\ge2.268681501\,\Lambda`$ [C41];
- sectors with $`|n|^2\ge4`$ are covered by the scheduler bound $`\kappa_{\min}\ge1.9798410293`$ (raw, per $`|n|^2`$, on the M30
  approximant). Since $`|n|^2=3`$ does not occur on this lattice, these sectors sit above $`4\times1.9798410293=7.919364`$ raw,
  which is $`2.639788\Lambda`$ on M30; after the transfer factors it is about $`2.63977\Lambda`$ on the true metric
  (REFEREE-SCHED item 4, referee arithmetic) [C42, C74].

A descriptive lower end is $`\lambda_1/\Lambda\ge1.94638098`$ [C41].

**Chain.**
1. BAR-CLWM and BAR-CLWN (plus their AMEND1s) certify the invariant sector and the $`|n|^2\le2`$ sectors on the M30 approximant [C60].
   The methods:
   - exact polynomial forms $`M=(\prod l)\,\mathrm{adj}\,\mathrm{Hess}\,\tilde u`$ and $`D=(\prod l)\det\mathrm{Hess}\,\tilde u`$,
     plus a certified degree-16 surrogate of the one rational weight;
   - a Crouzeix-Raviart / strip-lemma separator (Liu's constant $`C_{CR}(K)\le0.1893h_K`$, cited);
   - a Lehmann-type lower bound in strong form, a Rayleigh-Ritz argument for $`(L-\rho)^{-1}`$ needing no Goerisch dual
     space [C61].
2. BAR-CLWT carries these to the M84 approximant. It uses the per-sector comparison lemma with factor
   $`[(1-e)^{3/2}/(1+e)^{5/2},\,(1+e)^{3/2}/(1-e)^{5/2}]`$ and the certified $`\epsilon_{30}\le1.94\times10^{-6}`$ [C62].
3. MET5/MET6 (this paper) supply $`e=d_h=5.504\times10^{-8}`$ from M84 to the true metric [C40].
   - The lemma factors at $`d_h`$ are $`1.00000022015`$ (up) and $`0.999999779854`$ (down) [C40].
   - The transfer is conservative: it used the complex-gauge $`d_h`$, while the same-coordinates $`d_{\rm sym}=3.19\times10^{-8}`$
     is smaller, and the lemma accepts either [C41].

**Referees.**
- REFEREE-CLWM: separator and surrogate CORRECT; Lehmann CORRECT WITH FIXES.
- REFEREE-CLWN: new step correct line by line.
- REFEREE-CLWT: chain holds.
- REFEREE-SCHED: scheduler correct [C63].

The figure $`2.06988`$ that circulated earlier is retired, because it cannot be checked [C42].

**Conditions.** The corollary inherits every cited assumption of Section 7, items 1-9, plus its own items 10-12. Its
identification with CLW rests on Li-Sun (preprint).

---

## 6. An honest scorecard

Each test below was sealed before its run. Main table: every sealed test that the theorem chain or Corollary 3 uses. Appendix A: every other settled bar in this line of work, including 14 more FAILED and 1 VOID. Both tables are generated by `make_scorecard.py` from the lab ledger's settle records, in settle order. The reasons are hand-written, one line per row [C64]. Totals: 30 settled seals, 13 HOLD, 16 FAILED, 1 VOID.

**Main table** (9 settled seals: 7 HOLD, 2 FAILED)

| # | bar (seal) | settled (UTC) | verdict | reason |
|---|---|---|---|---|
| 1 | CLWI (2df5647e52) | 2026-11-11 21:57 | ![HOLD](chips/hold.svg) | the approximant family the MET bars build on (calibrated_on/depends_on). Its residual 6.0e-18 at M68 is SAMPLED, not certified. |
| 2 | CLWM (bba6649342) | 2026-11-11 22:56 | ![HOLD](chips/hold.svg) | lambda_1/Lambda in [1.9463965, 1.9483522], lambda_2/Lambda >= 2.0897208 (M30). The referee asked for the lower end to be written rounded down. |
| 3 | CLWN (2459207186) | 2026-11-12 03:50 | ![HOLD](chips/hold.svg) | sectors \|n\|^2 <= 2 (M30): lowest >= 2.2688 Lambda. |
| 4 | CLWT (0da760d766) | 2026-11-12 09:14 | ![HOLD](chips/hold.svg) | transfer M30 -> M84: eps30 <= 1.94e-6; budget left 6.53e-3. |
| 5 | MET2 (92a48677ad) | 2026-11-13 05:02 | ![FAILED](chips/failed.svg) | the thresholds were set from sampled values (rho 1e-18, Delta rho 1e-15, rho_H3 1e-9). The certified values (7.48e-16, 4.86e-13, 2.51e-9) miss them, and with them the fixed point did not close (MET4). The refined values (MET5) do. |
| 6 | MET3 (49a2fe684b) | 2026-11-13 05:02 | ![HOLD](chips/hold.svg) | all 80 sub-squares ok; certified/float ratios <= 1.01, except kappa at 3.2x; the mutation guard fails as required. |
| 7 | MET4 (f58c422d19) | 2026-11-13 05:02 | ![FAILED](chips/failed.svg) | not closed: x = 1.246 > 1. alpha1 = 5.06e-13 is over the 4.06e-13 budget, because one far-corner enclosure was loose. The pre item also failed, since MET2 was FAILED. |
| 8 | MET5 (ea72f3bc07) | 2026-11-13 05:45 | ![HOLD](chips/hold.svg) | x = 0.17725, 18/18 gates. Disclosure: this was a calibrated recomputation (the e3[3,3] test and an x ~ 0.75 estimate came before the seal), and MET4's pre item (MET2 HOLD) was dropped. |
| 9 | MET6 (39be7ec448) | 2026-11-13 06:00 | ![HOLD](chips/hold.svg) | x = 0.177254749606, 18/18, d_h = 5.504e-8, with lambda1 rounded down (referee fix). |

Amendments were sealed before each rerun: MET2-amend1/2; MET3-amend1..4; and MET5-amend1, the $`\lambda_1`$ fix, which was rescored as MET6. Data computed under superseded code was never scored [C65].

---

## 7. Cited assumptions, in full

These are the inputs that this paper cites rather than proves or re-derives, as listed by REFEREE-MET5 (6) [C66]:
1. **dP3 App. B Sobolev constants** ($`C_5..C_8`$, (B.7), (B.8)): arXiv:2609.32655, a preprint. Transcribed from the PDF, where
   the text extraction loses symbols.
2. **Liu's constant 0.07353**, inherited by $`\lambda_1`$ through the M30 certificate and the G4 port.
3. **Rivlin (1.38)**, the Lebesgue constant for first-kind Chebyshev nodes. OPENED.
4. **Donaldson's identity**. REFEREE-MET5 lists it as "secondary sources only". Since then the primary arXiv version has been
   opened and matches the lab's normalisation (factor 2, lattice measure) [C21]. Caveat: that is the arXiv version, not the
   published one.
5. **Abreu/Guillemin boundary class** of the symplectic potential of a smooth toric Kähler metric. Abreu arXiv:math/0004122
   Thm 2.8 has now been opened [C21]. A gap remains: Donaldson's class also asks for strict convexity on each face. That is
   plausible from Abreu's determinant condition, but neither paper states it.
6. **Elliptic bootstrap** (smoothness of $`\psi^*`$) and **Fredholm/onto** for $`L_t`$. These are textbook results, not opened.
7. **Li-Sun** Cor. 7.3 / Prop. 7.13, arXiv:2610.10425v2: a preprint. It is used only to name the Einstein metric as CLW.
8. **LeBrun, Crelle 678 (2013), Prop. 3**. OPENED, at lines 546-549 of the text.
9. **dP3 Thm 4.1** in its Lipschitz form. Read from garbled text, but not load-bearing, since we check the two
   contraction-mapping inequalities directly [C24].

Additional inputs cited only by Corollary 3 (BAR-CLWM/CLWN/CLWT) [C71]:

10. **X. Liu (2024)**, *Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems* (SpringerBriefs).
    The citations are: the Crouzeix-Raviart constant $`C_{CR}(K)\le0.1893\,h_K`$ (Table 4.1); the CR orthogonality ((4.2)-(4.3));
    Thm 3.2 (semi-definite $`a`$). OPENED by REFEREE-CLWM. The Lehmann-Maehly bound is cited there as Liu Cor. 5.1, but the
    seat's actual argument is its own Rayleigh-Ritz proof for $`(L-\rho)^{-1}`$. The referee asks for that proof to be cited as
    such.
11. **The inertia-count certificate** (certinert, inherited from the lab's earlier seat). It assumes that SuperLU does no
    pivoting (asserted in code) and uses Higham, Thm 9.3. Not re-audited by REFEREE-CLWM.
12. **The comparison lemma and the scheduler bound.** The lemma is derived in the lab and refereed (REFEREE-CLWT); the scheduler
    bound is refereed (REFEREE-SCHED). Both are lab results, not citations. They are listed here because they are not part of
    this paper's proof.

---

## 8. Limitations

- **The Killing-tensor count does not transfer yet.** The "exactly 4 quadratic Killing tensors" certificate (GEO2) holds for
  the M84 approximant. Transferring it to the true CLW metric needs a certified 7-jet at the test point to within $`10^{-8}`$.
  The derived (float, unsealed) chain gives that only if $`r\lesssim3\times10^{-15}`$. Our $`r=7.568\times10^{-10}`$ is about
  $`2.5\times10^{5}`$ times too large [C67].
- **$`H^6`$ controls $`C^3`$ only.** Higher jets need an interior bootstrap that is not built [C68].
- **The engine has local, hand-written parts:**
  - the float-disc rounding model, and the Arb/Chebyshev-model code;
  - each was refereed with fixes [C47-C49];
  - the referee's containment controls (REFEREE-CM fix 1) were not yet part of the sealed control list when MET5/MET6 ran [C47].
    They were run afterwards on the same pinned code as BAR-CTRL1, with required-to-fail negatives (Section 3.7).
  - No null yet attacks the refinement merge code; for example, dropping one child or replacing a sum by a max should be
    caught (REFEREE-MET5 (3)) [C72].
- **The certified residual is far above the true residual.** The certified $`\alpha_1=7.2\times10^{-14}`$ is about
  $`3\times10^{6}`$ times the sampled value $`2.4\times10^{-20}`$ (before refinement, MET4's was about 7 orders looser). The
  closing margin ($`x=0.177`$) is real, but enclosure widths, not the approximant, set it [C54, C31, C34].
- **One class versus the ball.** Theorem 1 holds on the whole $`a^*`$ ball, but Corollary 2 only at the exact critical class
  inside it.

---

## 9. Reproduction

The code and frozen data that recompute the final certificate and the assembly are in this repository: `code/` and `data/`.
Commands, expected output and run times are in `code/RUN.md`.

**Hashes.**
- `data/SHA256SUMS-PROD-met6.txt` is the lab's hash list of the frozen code and data (30 entries, paths relative to the lab
  folder `seatE/`). It covers the MET6 `make_cert.py` (`6105370f...4b4633`), `G4_eps_m48.json`, `assemble.py`, `chain.py`,
  `cert_refined_a1.json` and `result_refined_a1.json` (`45c44ca0...61073`, which matches the lab ledger) [C69].
- Four shipped scripts differ from that list in path lines only. `code/PATCHES.md` lists each change with both hashes.
- `data/prod_inputs.sha256` hashes the 160 production files (`data/prodA3/`, `data/prodB1/`).

**Commands** (`python` = Python 3.11 with python-flint 0.9.0, numpy, mpmath; see `code/RUN.md`):
```
# from code/seatE/: production A (scalar) and B (tensor), k = 0..4 (hours; not needed to check the certificate)
python prod_A.py k prodA2                          python prod_B.py k prodB1
# from code/seatE/: far-corner refinement (lists, splits 2x2, merges conservatively)
bash run_refine.sh
# from code/seatE/: certificate with the certified lambda1 (MET6), from the frozen production data
python make_cert.py ../../data/prodA3 ../../out/cert_refined_a1.json ../../data/prodB1
# from code/seatB2/: the assembly
python assemble.py ../../out/cert_refined_a1.json ../../out/result_refined_a1.json
```
The lab ran the last two steps as `make_cert.py prodA3 cert_refined_a1.json prodB1` from `seatE/` and
`assemble.py ../seatE/cert_refined_a1.json result_refined_a1.json` from `seatB2/`, each under an 8 GB memory cap.
A clean-room rerun of both steps from fresh copies gives byte-identical `cert_refined_a1.json` and `result_refined_a1.json`.

The night runners `run_prod2.sh` and `run_prod3.sh` wrapped production A and B with these controls [C53, C70]:
- unit tests;
- the square control;
- the M30 null, which certifies $`\sup|\rho|\ge3.3\times10^{-6}`$;
- the mutation pair;
- the NaN audit.
(`run_prod3.sh` is shipped for the record; it is a Windows/bash script that waits for the previous runner.)

The null certificate `data/cert_null.json` must fail gate 10 ($`x=1.775`$) [C59].

---

## Appendix A: every other settled seal in this line of work

Generated by `make_scorecard.py` (same source and order as Section 6).

| # | bar (seal) | settled (UTC) | verdict | reason |
|---|---|---|---|---|
| 1 | CLWB (c73ed2b2f4) | 2026-11-11 18:16 | ![FAILED](chips/failed.svg) | not closed: the single-domain bound did not carry over to the tensor problem. A new patch lemma (D3) was certified on Page only. |
| 2 | CLWA (07e592fc36) | 2026-11-11 18:59 | ![FAILED](chips/failed.svg) | not closed: no norm pair with certified regularity constants controlled the metric in C^1. Its part A-a (lambda_1 lower bound, later 7.36698) did certify and is used. |
| 3 | FAST (b757d923ad) | 2026-11-11 19:17 | ![HOLD](chips/hold.svg) | two-stage route (separator + Lehmann-Goerisch) certified all four Page classes, in 26.4 s total. |
| 4 | CLWD (96ca1efcbb) | 2026-11-11 19:28 | ![FAILED](chips/failed.svg) | not closed: a C^2 metric bound via the Abreu equation needed regularity constants, curvature sups and residual norms that did not exist yet. |
| 5 | CLWC (2583eeb9db) | 2026-11-11 19:34 | ![FAILED](chips/failed.svg) | nothing certified: the control stayed conditional on float constants, and the CLW case was not built. |
| 6 | CLWE (2ab4b0de4f) | 2026-11-11 20:41 | ![FAILED](chips/failed.svg) | not closed: dP3's fixed-point chain needed an H^3 residual of 1e-21..1e-15, against ~1e-6 for the approximant available then. |
| 7 | CLWG (3806eb5d56) | 2026-11-11 20:45 | ![FAILED](chips/failed.svg) | premise wrong: the rescaling the bar relied on leaves alpha2^2 alpha3 invariant (sealed without checking). Nothing built. |
| 8 | CLWH (0e9ba523a3) | 2026-11-11 21:04 | ![FAILED](chips/failed.svg) | control failed: on S^4 the chain came out up to 822x the exact value at k = 3, against an allowed factor 10. The decisive item passed. |
| 9 | CLWJ (50b789e837) | 2026-11-11 21:22 | ![FAILED](chips/failed.svg) | not closed: the margins were too thin for the inset-FEM comparison, and the mesh it needed exceeded the 8 GB memory cap. |
| 10 | CLWF (afd8a15f2d) | 2026-11-11 21:33 | ![HOLD](chips/hold.svg) | patched two-stage route (Lemma D3) certified on Page (-,+): [0.4660171759, 0.4660172610]. |
| 11 | CLWK (4ee2fab868) | 2026-11-11 22:02 | ![FAILED](chips/failed.svg) | not closed: the bar assumed a route the software did not have (bar-design error). Controls not run. |
| 12 | CLWL (ce982a5a3f) | 2026-11-11 22:22 | ![FAILED](chips/failed.svg) | not closed: no certified lambda_3 separator, so no lower bounds. |
| 13 | CLWM (ac61aa9e3d) | 2026-11-11 22:37 | ![VOID](chips/void.svg) | sealed in error: a name collision with the other lane's BAR-CLWM. No data scored. |
| 14 | MET1 (2e992a73d9) | 2026-11-11 22:58 | ![FAILED](chips/failed.svg) | not closed: ball enclosures of the tensor curvature blew up (~1e10 per unit box radius). This led to the switch to the scalar extremal route. |
| 15 | GEO1 (5648789783) | 2026-11-12 04:29 | ![FAILED](chips/failed.svg) | the control's energy-drift clause missed (2.3e-10 > 1e-10). The decisive item held (93.75% chaotic). |
| 16 | GEO1B (032e1fd606) | 2026-11-12 04:52 | ![HOLD](chips/hold.svg) | Page control at rtol 1e-13: 512/512 regular, energy drift 1.8e-11. Ledger pin warning: result file not found (no result hash recorded). |
| 17 | TT0 (e4ed1591b2) | 2026-11-12 04:59 | ![FAILED](chips/failed.svg) | NO-GO: the projected basis size per class was 11,885-95,117, against a limit of 1500. One control clause also failed as sealed. |
| 18 | GEO2 (fa2546c325) | 2026-11-12 05:17 | ![HOLD](chips/hold.svg) | exactly 4 quadratic Killing tensors on the M84 approximant at two points (Arb rank certificate); approximant only. |
| 19 | TT1P (7f0a75255d) | 2026-11-12 05:35 | ![FAILED](chips/failed.svg) | NO-GO: the best decoupled operator had 12 eigenvalues below 3.5 raw, where the true one has 4. |
| 20 | LEDGER1 (dfbc1c5126) | 2026-11-12 13:09 | ![HOLD](chips/hold.svg) | Vol(CLW) Lambda^2 = 132.786665471658... certified (needs no certified metric); it matches Hall-Murphy's 9 digits. |
| 21 | CTRL1 (20c4a51c45) | 2026-11-13 07:42 | ![HOLD](chips/hold.svg) | point-containment controls: pos-square 0/1620, pos-clw 0/98 (one sub-square, e0 [0.5,0.75]^2), neg-inv 350/810 and neg-tail 810/810 fail as required. The checker script was written after the seal; its hash is in data/ctrl1/CTRL1.json. |

---

## References (opened unless marked)
- X. Chen, C. LeBrun, B. Weber, "On conformally Kähler, Einstein manifolds", arXiv:0705.0710v2 (opened; published 2008, journal
  venue not checked here).
- C. LeBrun, "Einstein manifolds and extremal Kähler metrics", J. Reine Angew. Math. 678 (2013), DOI 10.1515/crelle.2012.006,
  arXiv:1009.1270, pp. 69-94. Prop. 3 opened (text l.546-549). Title, journal, issue, year and DOI were verified against the
  Crossref record and the arXiv abstract page; the page range against the zbMATH record (Zbl 1277.53042).
- Li, Sun, arXiv:2610.10425v2, Cor. 7.3, Prop. 7.13 (opened; preprint).
- Buttsworth, Hadden, Heyes, Platt, Wiseman, arXiv:2609.32655 (dP3), Thm 4.1, App. B (opened, garbled extraction).
- M. R. Douglas, J. Gómez-Serrano, F. Lehmann, D. Platt, Y. Qi, arXiv:2611.16908 (Fermat quartic; opened, abstract p.1).
- T. Buttsworth, L. Hodgkinson, arXiv:2412.01184, abstract p.1 (opened); Q. S. Wang, arXiv:2504.21644, abstract and Sec. 1.1
  p.1-2 (opened); Q. S. Wang, arXiv:2606.05519, abstract p.1 (opened).
- R. S. Bunch, S. K. Donaldson, arXiv:0803.0987, Sec. 1 p.1-2 (opened).
- S. K. Donaldson, arXiv:math/0407486v2 (Collect. Math. 56 (2005) 103-142), Cor. 1, Prop. 1, Def. 1 (opened).
- M. Abreu, arXiv:math/0004122v1, Thm 2.8 (opened).
- Hall-Haslhofer-Siepmann, arXiv:1111.4981, Lemma 7.3 (opened). Hall-Murphy, arXiv:1206.5489 Thm 1.1 and arXiv:1407.1272v2
  (opened).
- T. J. Rivlin, *The Chebyshev Polynomials*, 2nd ed., Wiley 1990, Sec. 1.3, eq. (1.38) (opened).
- X. Liu, *Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems*, SpringerBriefs in Mathematics
  (2024): Table 4.1, Thm 3.2, Cor. 5.1 (opened by REFEREE-CLWM). Source of the constants $`C_{CR}`$ and 0.07353.
- D. Z. Liu, Zenodo 23598158 (opened; numerical $`\lambda_1`$; a different author from X. Liu).
- Nifa, Zenodo 22851442, eqs. (5.2), (6.6)-(6.12) (opened).

---

*Prepared with AI assistance; all results are machine-checked and listed in the scorecard (Section 6).*
