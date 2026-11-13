"""assemble.py -- seat B2: the ONE script the sealed fixed-point bar runs (CLW metric line, extremal route).

Input : a JSON of constants (schema below; NOTES-B2 sec. 7).  Output: every intermediate with its source, the closing
        quantity x = 4 alpha1 alpha2^2 alpha3_Y (strict < 1), the Kantorovich radius r, the c* = 0 margin and the pointwise
        metric bound d in both gauges.  Exit code 0 iff every gate holds.
Arithmetic: Arb balls (python-flint), 256 bits.  Every input is read as a decimal STRING (floats are converted with repr, i.e.
        taken as that exact decimal).  Upper bounds enter as given; a gate "A < B" is passed only if Arb proves it.
        The chain F(1,1,1) is seatX/chain.py imported READ-ONLY, with its module-level sqrt replaced by an Arb-aware one
        (no file of seatX is edited).  F is nondecreasing in every constant and nonincreasing in lambda_1
        (NOTES-X-PROOFS sec. 7, Claim 1), so certified upper bounds (lower for lambda_1) give a certified upper bound for C_reg.
dP3's C12 is NEVER used (REFEREE-EXTREMAL (5), MUST-FIX 10): only C5..C8 (App. B), C_inf = 4 C7 C6 C8 C5 from (B.7),
        C32 = C_inf + 4 C5^2 and C9 = max{2 C_inf, C32, 8 C5^2} ((B.8), re-derived in NOTES-Q2-PROOF Lemma B).

Schema (all sup norms over M for g~ = g_K(u84), for all a in the a* ball; all L2 norms are L2(M, g~) = (2 pi)^2 int_P dx):
  label        free text; printed first.  certified: true/false (false => banner "sampled, not certified").
  mu           LOWER bound, Ric(g~) >= mu g~, must be > 0 (dP3 App. B hypothesis; K4)
  Rm           upper bound of sup max(0, -min eig Ric).  With mu > 0 certified this is 0 (referee fix 1): pass "0".
  Rp, Rop      upper bounds of sup max eig Ric (positive part) and sup |Ric|_op
  kappa        upper bound of the per-slot Weitzenboeck operator norm (NOTES-X sec. 2)
  K            [K_0, K_1, K_2]     sup |nabla^q Riem|        (K_0 is not used once kappa is given; kept for the record)
  Rc           [Rc_1, Rc_2, Rc_3]  sup |nabla^q Ric|, q = 1..3
  R0           [R0_0, R0_1, R0_2]  sup |nabla^q Ric_0|
  S_max, S_min max and min of Donaldson's affine S over P (vertex values, exact on the ball)
  rho_sup      upper bound of sup |rho|, rho = Scal(g~) - S.   sup|s| := S_max + rho_sup (referee fix 2)
  S_d          [S_0, S_1, S_2]     sup |nabla^q ds|, s = Scal(g~)
  gradS_x      Euclidean length of grad_x S (exact from S's coefficients); used with |dS|_g~ <= gradS_x sqrt(U_sup)
  lambda1      LOWER bound of lambda_1(L_t) on affine^perp (K7)
  V            [lo, hi] enclosure of Vol(M, g~) = (2 pi)^2 Area(P)
  kappa_A      upper bound of 1 + k1 + k2, k_j = sqrt(lambda_max(G_P^-1 G_j)), all Gram matrices in dx (K9)
  GPinv_norm   upper bound of ||G_P^-1||_2, G_P = int_P e e^T dx, e = (1, x1, x2)
  GP_diag      upper bounds of [Area(P), int x1^2 dx, int x2^2 dx]
  U_sup        upper bound of sup_P lambda_max((Hess u84)^-1)
  rho_L2, lap_rho_L2, rho_H3   upper bounds of ||rho||_L2(M), ||Delta rho||_L2(M), ||rho||_H3(M) (sum norm)
  w_prime      the weight w = w' of the Y-norm (default 2^-10)
  d_target     the peer lane's bound (default 6.53e-3)
"""
import sys, json, math
from flint import arb, ctx

ctx.prec = 256
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'seatX'))
import chain as CH                                   # read-only import


def _asqrt(x):
    return x.sqrt() if isinstance(x, arb) else math.sqrt(x)


CH.sqrt = _asqrt                                     # Arb-aware sqrt inside chain.step / chain.F (module attribute only)


def A(v):
    if isinstance(v, arb):
        return v
    if isinstance(v, (int,)):
        return arb(v)
    if isinstance(v, float):
        v = repr(v)
    return arb(str(v))


def up(x):
    """rigorous upper bound of an Arb ball, as an exact-point ball"""
    return x.upper()


def amax(*xs):
    m = xs[0]
    for y in xs[1:]:
        m = m.max(y)
    return m


def fmt(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else str(x)


GATES = []


def gate(name, ok, detail=''):
    GATES.append((name, bool(ok)))
    print(f'   GATE {name}: {"PASS" if ok else "FAIL"}  {detail}')


def line(name, val, src):
    print(f'{name:<26s} = {fmt(val):<22s} | {src}')


def sobolev(mu, V):
    """dP3 App. B constants for n = 4 as transcribed in seatA4/sobolev_page.py (PDF p.32-34, checked by REFEREE-EXTREMAL (5)).
    mu is a lower bound (all constants are nonincreasing in mu); V is an enclosure."""
    pi = arb.pi()
    om3, om4 = 2 * pi ** 2, 8 * pi ** 2 / 3
    K42 = (arb(4) / (8 * om4.sqrt())).sqrt()                       # sqrt(4/(n(n-2) om4^(2/n))), n = 4
    C3 = 4 * (om3 / 4).root(4)
    C4 = V.root(4) / arb(2).root(4) * (mu / 3).sqrt() / (arb(2) / 3)   # int_0^{pi/2} cos^3 = 2/3
    t = (om4 / (V * mu ** 2)).root(4)
    Vm14 = 1 / V.root(4)
    C5 = amax(arb(3).sqrt() * K42 * t, Vm14)
    C6 = amax(3 * arb(3).sqrt() * K42 * t, Vm14)
    C7 = amax((arb(11) / 8) ** (arb(11) / 12) * (C3 / C4) * (4 * V).root(6) / om3.root(4), 2 / V.root(12))
    C8 = V.root(12)
    return [up(c) for c in (C5, C6, C7, C8)]


def Kfun(mR, s, C32, C9):
    """NOTES-Q2-PROOF sec. 3 Theorem: alpha3(r) = 4 C32 K(2 r), gamma = 1/(1 - C9 s); needs C9 s < 1."""
    g = 1 / (1 - C9 * s)
    return (mR + 1) * g ** 2 + g ** 3 * (1 + C32 * s / 4) + g ** 2 * (1 + C32 * g * s) / 2 + C32 * s * g ** 3 * (1 + g + C32 * g * s) / 2


def creg(c):
    """C_reg := F(1,1,1) (NOTES-X-PROOFS sec. 7), options split / Rop / kap / interp ON, as in seatX/f111.py."""
    nan = arb('nan')                                  # Rc_0 is never read (Rop replaces it); NaN would poison the result if it were.  K_3 only ever meets the factor m = 0 -> pass 0
    cc = CH.consts(c['lambda1'], c['Rp'], c['Rm'], [nan] + c['Rc'], c['K'] + [arb(0)], c['S_d'], Rop=c['Rop'], R0c=c['R0'], ssup=c['ssup'])
    cc.update(kap=c['kappa'], interp=True)
    tot, a, b = CH.F((arb(1), arb(1), arb(1)), cc, True)
    return up(tot), a, b


def assemble(cert):
    cert_ok = bool(cert.get('certified', False))
    print('=' * 110)
    print('assemble.py (seat B2)   label:', cert.get('label', '(none)'))
    print('INPUT STATUS:', 'CERTIFIED constants' if cert_ok else '*** SAMPLED, NOT CERTIFIED: nothing below is a proof ***')
    print('=' * 110)
    g = lambda k: A(cert[k])
    c = dict(mu=g('mu'), Rm=g('Rm'), Rp=g('Rp'), Rop=g('Rop'), kappa=g('kappa'), K=[A(v) for v in cert['K']],
             Rc=[A(v) for v in cert['Rc']], R0=[A(v) for v in cert['R0']], S_d=[A(v) for v in cert['S_d']], lambda1=g('lambda1'))
    V = A(cert['V'][0]).union(A(cert['V'][1]))
    Smax, Smin, rsup, gSx = g('S_max'), g('S_min'), g('rho_sup'), g('gradS_x')
    kA, GPi, Usup = g('kappa_A'), g('GPinv_norm'), g('U_sup')
    GPd = [A(v) for v in cert['GP_diag']]
    n0, nL, n3 = g('rho_L2'), g('lap_rho_L2'), g('rho_H3')
    wp = A(cert.get('w_prime', '0.0009765625'))
    dtar = A(cert.get('d_target', '6.53e-3'))
    res = {}

    print('\n-- 0. input consistency')
    gate('mu > 0 (dP3 App. B hypothesis, K4)', c['mu'] > 0, f'mu = {fmt(c["mu"])}')
    gate('Rm >= 0 and (Rm = 0 only with mu > 0) (referee fix 1)', (c['Rm'] >= 0) and ((c['Rm'] > 0) or (c['mu'] > 0)), f'Rm = {fmt(c["Rm"])}')
    gate('mu <= Rp <= Rop (not provably violated)', not (c['mu'] > c['Rp']) and not (c['Rp'] > c['Rop']))
    gate('S_min > 0 (K2)', Smin > 0, f'S_min = {fmt(Smin)}')
    gate('V <= (2 pi)^2 GP_diag[0]', V.upper() <= (2 * arb.pi()) ** 2 * GPd[0] * (1 + arb(2) ** -40), f'V = {fmt(V)}')
    c['ssup'] = up(Smax + rsup)
    line('sup|s|', c['ssup'], 'S_max + sup|rho| (referee fix 2)')

    print('\n-- 1. C_reg = F(1,1,1)   [seatX/chain.py imported; NOTES-X-PROOFS sec. 7]')
    Creg, a, b = creg(c)
    for j, v in enumerate(a):
        line(f'  A_{j} (bound of ||nabla^{j} f||)', v, 'chain.F at w = (1,1,1)')
    for j, v in enumerate(b):
        line(f'  B_{j} (bound of ||nabla^{j} Delta f||)', v, 'chain.F at w = (1,1,1)')
    line('C_reg', Creg, 'sum A_j; ||psi||_H6 <= C_reg ||L_t psi||_H2 on X (K8)')
    gate('C_reg finite', Creg.is_finite())

    print('\n-- 2. Sobolev constants   [dP3 App. B as transcribed, from mu and V only; C12 not used]')
    C5, C6, C7, C8 = sobolev(c['mu'], V)
    Cinf = up(4 * C7 * C6 * C8 * C5)
    C32 = up(Cinf + 4 * C5 ** 2)
    C9 = up(amax(2 * Cinf, C32, 8 * C5 ** 2))
    for nme, v, s in (('C5', C5, 'Cor. B.4, H^1 -> L^4'), ('C6', C6, 'Cor. B.4, H^{1,3} -> L^12'), ('C7', C7, 'Cor. B.4, H^{1,12} -> L^inf'),
                      ('C8', C8, 'Prop. B.5, L^4 -> L^3'), ('C_inf', Cinf, '4 C7 C6 C8 C5, (B.7): ||T||_inf <= C_inf ||T||_H3 (Kato)'),
                      ('C32', C32, 'C_inf + 4 C5^2, NOTES-Q2-PROOF Lemma A'), ('C9', C9, 'max{2 C_inf, C32, 8 C5^2}, (B.8) = NOTES-Q2-PROOF Lemma B')):
        line(nme, v, s)

    print('\n-- 3. alpha1   [K10; Pi rho = 0 exactly by Donaldson, so ||rho||_Y = ||rho||_H2]')
    gradr = up((n0 * nL).sqrt())
    hessr = up((nL ** 2 + c['Rm'] * n0 * nL).sqrt())
    alpha1 = up(n0 + gradr + hessr)
    line('||rho||_L2(M)', n0, 'certificate')
    line('||grad rho||_L2', gradr, 'sqrt(||rho|| ||Delta rho||)   (<Delta rho, rho> = ||grad rho||^2)')
    line('||Hess rho||_L2', hessr, 'sqrt(||Delta rho||^2 + Rm ||grad rho||^2)   (integrated Bochner, NOTES-X L2\' at j = 0)')
    line('alpha1 = ||rho||_H2', alpha1, 'sum of the three')

    print('\n-- 4. alpha2   [NOTES-B2 sec. 1c]')
    kw = up(kA + wp)
    t = up(Creg * kw * C32 * n3)
    line('kappa_A + w\'', kw, 'norm equivalence ||y||_Y <= (kappa_A + w\') ||y||_H2')
    line('t', t, 'C_reg (kappa_A + w\') C32 ||rho||_H3   (moment-map shift <d rho, d.>)')
    gate('t < 1 (K11)', t < 1)
    alpha2 = up(Creg / (1 - t))
    line('alpha2', alpha2, 'C_reg / (1 - t)')

    print('\n-- 5. alpha3   [NOTES-Q2-PROOF Theorem: alpha3(r) = 4 C32 K(2r); alpha3_Y = (kappa_A + w\') alpha3(r)]')
    mR = up(c['Rop'] + 2 * c['Rc'][0] + c['Rc'][1])
    line('m_R', mR, 'Rop + 2 Rc_1 + Rc_2')
    a3Y0 = up(kw * 4 * C32 * Kfun(mR, arb(0), C32, C9))
    line('alpha3_Y(0)', a3Y0, '(kappa_A + w\') 4 C32 (m_R + 5/2)')
    rmax = up(1 / (2 * alpha2 * a3Y0))
    line('r_max', rmax, '1 / (2 alpha2 alpha3_Y(0)): no admissible radius exceeds it')
    gate('C9 . 2 r_max < 1', C9 * 2 * rmax < 1, f'C9 s = {fmt(up(C9 * 2 * rmax), 6)}')
    a3Yrmax = up(kw * 4 * C32 * Kfun(mR, 2 * rmax, C32, C9))
    budget = (1 / (4 * alpha2 ** 2 * a3Yrmax)).lower()
    line('alpha3_Y(r_max)', a3Yrmax, 'tail included')
    line('alpha1 budget', budget, '1 / (4 alpha2^2 alpha3_Y(r_max)), lower end')
    Rhat = up(2 * alpha1 * alpha2)
    line('R_hat', Rhat, '2 alpha1 alpha2 (a-priori radius at which alpha3 is evaluated)')
    okR = bool(C9 * 2 * Rhat < 1)
    gate('C9 . 2 R_hat < 1 (Theorem hypothesis)', okR)
    if not okR:
        return finish(res)
    a3Y = up(kw * 4 * C32 * Kfun(mR, 2 * Rhat, C32, C9))
    line('alpha3_Y(R_hat)', a3Y, '(kappa_A + w\') 4 C32 K(2 R_hat)')

    print('\n-- 6. closing   [dP3 Thm 4.1 form; NOTES-B2 K14-K15]')
    x = up(4 * alpha1 * alpha2 ** 2 * a3Y)
    line('x = 4 a1 a2^2 a3_Y', x, 'closing quantity, must be < 1 STRICT')
    okx = bool(x < 1)
    gate('x < 1 (K14)', okx)
    if not okx:
        return finish(res)
    r = up(2 * alpha1 * alpha2 / (1 + (1 - x).sqrt()))
    line('r', r, '(1 - sqrt(1 - x)) / (2 alpha2 alpha3_Y) = 2 a1 a2 / (1 + sqrt(1 - x))')
    gate('r <= R_hat (alpha3 was evaluated on a larger ball)', r <= Rhat)
    gate('alpha2 alpha1 + alpha2 alpha3_Y r^2 <= r (ball mapped into itself)', alpha2 * alpha1 + alpha2 * a3Y * r ** 2 <= r * (1 + arb(2) ** -200))
    gate('2 alpha2 alpha3_Y r < 1 (contraction)', 2 * alpha2 * a3Y * r < 1, f'= {fmt(up(2 * alpha2 * a3Y * r), 6)}')

    print('\n-- 7. ball inside the domain, c* = 0   [NOTES-B2 K16, K17]')
    bcx = up(arb(2).sqrt() * Cinf * r)
    line('b = sqrt2 C_inf r', bcx, '|B|_op <= |B|_F / sqrt2 (J-commuting), ||nabla^2 psi||_inf <= C_inf r')
    gate('b < 1 (g_psi positive, K16)', bcx < 1)
    Ie = [GPd[0], (GPd[0] * GPd[1]).sqrt(), (GPd[0] * GPd[2]).sqrt()]
    coef = up(GPi * (2 * sum((Usup.sqrt() * Cinf * I) ** 2 for I in Ie)).sqrt())
    line('c* coefficient', coef, '||G_P^-1|| ||Delta||_F / r; |Delta_ji| <= sqrt(U_sup) C_inf r int_P|e_j|, int|x_a| <= sqrt(Area G_aa)')
    mc = up(coef * r)
    line('c* margin', mc, '||G_P^-1|| ||Delta||_F, must be < 1')
    gate('c* = 0: margin < 1 (K17)', mc < 1)

    print('\n-- 8. pointwise metric bound   [NOTES-B2 sec. 5, K18]')
    dcx = bcx
    m = up(Cinf * r * (1 + 3 * Cinf * r / (1 - bcx)) / (1 - bcx))
    dsym = up((arb(2).sqrt() * m).expm1())
    eta = up(gSx * Usup.sqrt() * Cinf * r / Smin)
    dh = up((1 + dcx) / (1 - eta) ** 2 - 1)
    line('d_cx', dcx, 'complex gauge: (1 - d) g84 <= g_true <= (1 + d) g84, d = sqrt2 C_inf r')
    line('m', m, 'sup |nabla_t^2 psi|_{g_t} <= C_inf r (1 + 3 C_inf r/(1-b)) / (1-b)')
    line('d_sym', dsym, 'symplectic gauge (same (x, theta)): exp(sqrt2 m) - 1')
    line('eta', eta, '|dS|_g~ C_inf r / S_min, |dS|_g~ <= |grad_x S| sqrt(U_sup) (referee fix 12)')
    line('d_h (complex gauge)', dh, 'Einstein metrics c S^-2 g: (1 + d_cx)/(1 - eta)^2 - 1; symplectic gauge: d_sym')
    gate(f'd_cx <= d_target = {fmt(dtar, 4)}', dcx <= dtar)
    gate(f'd_sym <= d_target = {fmt(dtar, 4)} (K18)', dsym <= dtar)
    gate(f'd_h <= d_target = {fmt(dtar, 4)}', dh <= dtar)
    res.update(C_reg=fmt(Creg, 20), C5=fmt(C5, 20), C6=fmt(C6, 20), C7=fmt(C7, 20), C8=fmt(C8, 20), C_inf=fmt(Cinf, 20), C32=fmt(C32, 20),
               C9=fmt(C9, 20), alpha1=fmt(alpha1, 20), alpha2=fmt(alpha2, 20), m_R=fmt(mR, 20), alpha3_Y_0=fmt(a3Y0, 20),
               alpha3_Y=fmt(a3Y, 20), r_max=fmt(rmax, 20), alpha1_budget=fmt(budget, 20), x=fmt(x, 20), r=fmt(r, 20),
               c_star_margin=fmt(mc, 20), d_cx=fmt(dcx, 20), d_sym=fmt(dsym, 20), d_h=fmt(dh, 20), eta=fmt(eta, 20), t=fmt(t, 20))
    return finish(res, cert_ok)


def finish(res, cert_ok=False):
    allok = all(ok for _, ok in GATES)
    print('\n' + '=' * 110)
    print(f'GATES: {sum(ok for _, ok in GATES)}/{len(GATES)} pass.  VERDICT: {"ALL GATES HOLD" if allok else "NOT CLOSED"}'
          + ('' if cert_ok else '   [inputs SAMPLED, NOT CERTIFIED: a dry run, not a proof]'))
    res.update(gates=[[n, ok] for n, ok in GATES], all_gates=allok, certified_inputs=cert_ok)
    return res


if __name__ == '__main__':
    cert = json.load(open(sys.argv[1]))
    res = assemble(cert)
    if len(sys.argv) > 2:
        json.dump(res, open(sys.argv[2], 'w'), indent=1)
    sys.exit(0 if res['all_gates'] else 1)
