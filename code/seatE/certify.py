"""certify.py -- seat E: certified bounds on one (s,t) sub-square Q = [s0, s0+h] x [t0, t0+h] of vertex element k, for the
Kahler metric g_K:  (a) sup of |nabla^q Riem|^2, |nabla^q Ric|^2, |nabla^q Ric_0|^2, |nabla^q ds|^2, |s|, tr N^8, tr N^16
(kappa), and a certified lower/upper bound of the Ricci eigenvalues;  (b) sup |psi| and sup |Delta psi| (psi = Scal(u) - S).

Method (per analytic function F of (s, t) on Q; local variables u, v in [-1, 1], s = s0 + h (u + 1)/2, ...):
 1. ANALYTICITY on the closed polyellipse E_rho x E_rho: every denominator of the chart formulas is bounded away from 0 on a cover
    of the filled region by complex boxes (conjugation symmetry: half in s): Re A(t), Re B(s) (principal sqrt), |l_r| (r != a, b),
    |det m|, and (polar) |s|, |t|.  The cart chart has det g_K = 16 and the polar chart det g_K = 16 p^2 q^2, so g^-1 adds nothing else.
 2. CRUDE BOUND M_F = max over disc boxes covering the distinguished boundary dE_rho x dE_rho (s-arcs theta in [0, pi], by
    F(conj) = conj F) of the enclosure of |F| (max principle in each variable, given 1).
 3. INTERPOLATION at (n+1)^2 Chebyshev points (2nd kind); |F - p| <= 4 M_F rho^-n (1 + Lambda_n)/(rho - 1) on Q
    (Trefethen ATAP Thm 8.2 in each variable; Lambda_n <= (2/pi) log(n+1) + 1, Thm 15.2).  Node values: float-disc engine (a),
    arb (b); the input points are balls containing the exact Chebyshev nodes.
 4. SUP / INF OF p: Chebyshev coefficients c = A F A^T in disc arithmetic.  (b): |p| <= sum |c|.  (a): p on a grid
    u_g = cos(g pi/(m-1)) plus |p(u) - p(u_g)| <= (sum |c_kl| k^2 + sum |c_kl| l^2) pi/(2(m-1)) (Markov, |T_k'| <= k^2; every
    u lies within pi/(2(m-1)) of a grid point since |d cos| <= |d phi|).
 Ricci: in real dim 4 Kahler, eigenvalues s/4 +- |Ric_0|/2 (multiplicity 2 each).  lambda_min >= c  <=  s/4 - c >= 0 and
    H_c = Q0 - c s/2 + c^2 >= 0 with Q0 = s^2/16 - |Ric_0|^2/4 (= lambda1 lambda2) analytic; interpolation is linear, so
    p_{H_c} = p_Q0 - (c/2) p_s + c^2 and its error <= err(Q0) + (c/2) err(s).  Same with C for lambda_max <= C.
 kappa = max |eig| of N = (g^-1 (x) g^-1) M, M_{(k,a),(i,d)} = R_ikad (symmetric w.r.t. g (x) g):  kappa <= (tr N^{2j})^{1/(2j)}.
"""
import math, time, json
import numpy as np
from flint import ctx, arb
import common
from common import Disc, ElemMap, SAFE, TINY
import fbox, jets, jets2, arbpsi, arbrho

U2 = 2.3e-16
KEYS_A = ['K0', 'K1', 'K2', 'K3', 'Rc0', 'Rc1', 'Rc2', 'Rc3', 'R0c0', 'R0c1', 'R0c2', 'S0', 'S1', 'S2', 's', 'Q0', 'T8', 'T16', 'T32', 'T64']
WANT = ('riem', 'ric', 'abreu', 'kappa')


# ------------------------------------------------------------------ disc helpers
def dmm(Am, Ar, Bm, Br, spec='zij,zjk->zik', K=16):
    mid = np.einsum(spec, Am, Bm)
    big = np.einsum(spec, np.abs(Am) + Ar, np.abs(Bm) + Br)
    sml = np.einsum(spec, np.abs(Am), np.abs(Bm))
    return mid, (big - sml) + (3 * K + 20) * U2 * big + TINY


def kappa_traces(o):
    """tr N^8, tr N^16, tr N^32, tr N^64 as Discs (nb,)."""
    Rm, Rr = o['R0']; gm, gr = o['g0']
    G = jets.inv0(Disc(np.moveaxis(gm, 0, -1).astype(complex), np.moveaxis(gr, 0, -1)))
    Gm, Gr = np.moveaxis(G.m, -1, 0), np.moveaxis(G.r, -1, 0)
    if not np.iscomplexobj(Rm):
        Gm = Gm.real
    # N[z, k, a, i, d] = gi[k, K] gi[a, A] R[i, K, A, d]
    t1m, t1r = dmm(Gm, Gr, Rm, Rr, 'zkK,ziKAd->zkiAd', 4)
    Nm, Nr = dmm(Gm, Gr, t1m, t1r, 'zaA,zkiAd->zkaid', 4)
    nb = Nm.shape[0]
    Nm = Nm.reshape(nb, 16, 16); Nr = Nr.reshape(nb, 16, 16)
    P2 = dmm(Nm, Nr, Nm, Nr); P4 = dmm(*P2, *P2); P8 = dmm(*P4, *P4); P16 = dmm(*P8, *P8); P32 = dmm(*P16, *P16); P64 = dmm(*P32, *P32)
    out = []
    for Pm, Pr in (P8, P16, P32, P64):
        tm = np.einsum('zii->z', Pm); tr = np.einsum('zii->z', Pr) * (1 + SAFE) + 40 * U2 * np.einsum('zii->z', np.abs(Pm))
        out.append(Disc(tm, tr))
    return out


def funcs_a(o):
    F = {}
    for q in range(4):
        F['K%d' % q] = o['riem'][q]; F['Rc%d' % q] = o['Rc'][q]
    for q in range(3):
        F['R0c%d' % q] = o['R0c'][q]; F['S%d' % q] = o['Sd'][q]
    F['s'] = o['s']
    F['Q0'] = o['s'] * o['s'] * (1.0 / 16) - o['R0c'][0] * 0.25
    F['T8'], F['T16'], F['T32'], F['T64'] = kappa_traces(o)
    return F


# ------------------------------------------------------------------ geometry of the sub-square
class Sub:
    def __init__(self, prob, model, k, s0, t0, h, chart):
        self.prob, self.model, self.k, self.s0, self.t0, self.h, self.chart = prob, model, k, s0, t0, h, chart
        self.em = ElemMap(prob, k)

    def st(self, u, v):
        """local -> (s, t) (Disc arrays, u, v Discs)."""
        hs = self.h / 2
        return Disc(self.s0 + hs + hs * u.m, hs * u.r * (1 + SAFE) + 1e-16 * (abs(self.s0) + self.h)), \
            Disc(self.t0 + hs + hs * v.m, hs * v.r * (1 + SAFE) + 1e-16 * (abs(self.t0) + self.h))

    def chartV(self, S, T):
        if self.chart == 'xchart':                     # symplectic chart (interior only): x(s, t) is polynomial
            return list(common.elem_x(self.prob, self.k, S, T))
        p, q = self.em.pq(S, T)
        n = len(S.m)
        z = Disc(np.zeros(n, p.m.dtype))
        if self.chart == 'cart':
            return [p, z, q, z]
        if self.chart == 'polar':
            return [p, q]
        if self.chart == 'rodA':
            return [p, z, q]
        return [p, q, z]

    def run(self, S, T, want=WANT, chunk=512):
        outs = []
        for i in range(0, len(S.m), chunk):
            V = self.chartV(S[i:i + chunk], T[i:i + chunk])
            outs.append(jets2.pipeline(self.model, self.k, self.chart, V, 3, want=want))
        return outs


def ellipse_arcs(rho, n, half):
    """discs covering the ellipse z(th) = (rho e^{i th} + e^{-i th}/rho)/2, th in [0, pi] (half) or [0, 2 pi)."""
    span = np.pi if half else 2 * np.pi
    th = (np.arange(n) + 0.5) * span / n
    z = (rho * np.exp(1j * th) + np.exp(-1j * th) / rho) / 2
    rad = (rho + 1 / rho) / 2 * (span / n) / 2 * (1 + 1e-9)
    return Disc(z, np.full(n, rad))


def filled_cover(rho, d, half):
    a = (rho + 1 / rho) / 2; b = (rho - 1 / rho) / 2
    xs = np.arange(-a + d / 2, a + d / 2, d)
    ys = np.arange(d / 2, b + d / 2, d) if half else np.arange(-b + d / 2, b + d / 2, d)
    Z = (xs[:, None] + 1j * ys[None, :]).ravel()
    return Disc(Z, np.full(len(Z), d / math.sqrt(2) * (1 + 1e-9)))


def prod_boxes(A, B):
    na, nb = len(A.m), len(B.m)
    return Disc(np.repeat(A.m, nb), np.repeat(A.r, nb)), Disc(np.tile(B.m, na), np.tile(B.r, na))


def denominators(sub, rho, d=0.1, chunk=20000):
    """lower bounds of the denominators on a cover of the filled polyellipse (local d -> s units h/2)."""
    U = filled_cover(rho, d, True); V = filled_cover(rho, d, False)
    Uu, Vv = prod_boxes(U, V)
    P = sub.prob.P; k = sub.k
    a_i, b_i = k, (k + 1) % P.d
    Am = [list(P.nu[a_i]), list(P.nu[b_i])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    res = dict(reA=np.inf, reB=np.inf, lr=np.inf, detm=np.inf, abs_s=np.inf, abs_t=np.inf, boxes=len(Uu.m))
    model = sub.model
    for i in range(0, len(Uu.m), chunk):
        S, T = sub.st(Uu[i:i + chunk], Vv[i:i + chunk])
        nb = len(S.m)
        K = lambda j: sub.em.K(j, nb)
        one = Disc(np.ones(nb))
        A = (one - T * T) * K(0) + T * T * K(1)
        B = (one - S * S) * K(2) + S * S * K(3)
        res['reA'] = min(res['reA'], float(np.min(A.m.real - A.r)))
        res['reB'] = min(res['reB'], float(np.min(B.m.real - B.r)))
        res['abs_s'] = min(res['abs_s'], float(np.min(np.abs(S.m) - S.r)))
        res['abs_t'] = min(res['abs_t'], float(np.min(np.abs(T.m) - T.r)))
        la = S * S * A; lb = T * T * B
        LA = la - Disc(np.full(nb, model.lam[a_i][0]), np.full(nb, model.lam[a_i][1]))
        LB = lb - Disc(np.full(nb, model.lam[b_i][0]), np.full(nb, model.lam[b_i][1]))
        x = [LA * float(Ai[i_][0]) + LB * float(Ai[i_][1]) for i_ in range(2)]
        R = [[Disc(np.zeros(nb, complex)) for _ in range(2)] for _ in range(2)]
        for r in range(P.d):
            if r in (a_i, b_i):
                continue
            nvv = [Ai[0][j] * P.nu[r][0] + Ai[1][j] * P.nu[r][1] for j in range(2)]
            lr = x[0] * float(P.nu[r][0]) + x[1] * float(P.nu[r][1]) + Disc(np.full(nb, model.lam[r][0]), np.full(nb, model.lam[r][1]))
            res['lr'] = min(res['lr'], float(np.min(np.abs(lr.m) - lr.r)))
            il = lr.recip()
            for i_ in range(2):
                for j in range(2):
                    if nvv[i_] * nvv[j]:
                        R[i_][j] = R[i_][j] + il * (nvv[i_] * nvv[j] / 2)
        if model.fsrc is not None:
            D = model.fsrc(x[0], x[1])
            H = [[D[(2, 0)], D[(1, 1)]], [D[(1, 1)], D[(0, 2)]]]
            for aa in range(2):
                for bb in range(2):
                    for i_ in range(2):
                        for j in range(2):
                            co = Ai[i_][aa] * Ai[j][bb]
                            if co:
                                R[aa][bb] = R[aa][bb] + H[i_][j] * float(co)
        m11 = one + la * R[0][0] * 2.0; m12 = la * R[0][1] * 2.0; m21 = lb * R[1][0] * 2.0; m22 = one + lb * R[1][1] * 2.0
        dm = m11 * m22 - m12 * m21
        res['detm'] = min(res['detm'], float(np.min(np.abs(dm.m) - dm.r)))
    need = ['reA', 'reB', 'lr', 'detm'] + (['abs_s', 'abs_t'] if sub.chart in ('polar', 'xchart') else []) \
        + (['abs_s'] if sub.chart == 'rodB' else []) + (['abs_t'] if sub.chart == 'rodA' else [])
    res['ok'] = bool(all(res[k_] > 0 for k_ in need))
    res['need'] = need
    return res


def _ub(v):
    """max over the boxes of an upper bound of |F|.  amend3: a NaN anywhere (0 * inf or inf - inf after an infinite radius) means
    the bound is inf; Python's max() would otherwise silently keep the earlier value.  Unchanged value when no NaN occurs."""
    x = float(np.max(v.upper_abs()))
    return math.inf if math.isnan(x) else x


def torus_max(sub, rho, n_arc, chunk=512):
    U = ellipse_arcs(rho, n_arc // 2, True); V = ellipse_arcs(rho, n_arc, False)
    Uu, Vv = prod_boxes(U, V)
    S, T = sub.st(Uu, Vv)
    M = {}
    for o in sub.run(S, T, chunk=chunk):
        F = funcs_a(o)
        F['psi'] = o['psi']; F['lap'] = o['lap']
        for j in range(4):
            F['rho%d' % j] = o['abreu'][j]
        for kk, v in F.items():
            M[kk] = max(M.get(kk, 0.0), _ub(v))
    return M, len(S.m), float(np.max(S.r))


def lebesgue(n):
    return 2 / math.pi * math.log(n + 1) + 1


def interp_err(M, rho, n):
    """2D tensor interpolation error, self-contained (ROUNDING-LEMMA.md sec. 3): |a_kl| <= 4 M rho^-(k+l) (Cauchy on the
    polyellipse), I(T_k T_l) = +-T_k' T_l' or 0 for either node family, so |f - I f| <= 2 sum_{max(k,l) > n} |a_kl|
    <= 16 M rho^(1-n) / (rho - 1)^2."""
    return 16 * M * rho ** (1 - n) / (rho - 1) ** 2 * (1 + 1e-9)


def cheb_nodes(n):
    return np.cos(np.arange(n + 1) * np.pi / n)


def coef_matrix(n):
    j = np.arange(n + 1)
    A = (2.0 / n) * np.cos((np.outer(j, j) % (2 * n)) * np.pi / n)        # reduced: argument in [0, 2 pi), entry error ~1e-15 absolute
    A[:, 0] *= 0.5; A[:, n] *= 0.5
    A[0, :] *= 0.5; A[n, :] *= 0.5
    return A


def cheb1_nodes_arb(n):
    return [arb.cos_pi(arb(2 * j + 1) / (2 * n + 2)) for j in range(n + 1)]


def coef_matrix1(n):
    """first-kind nodes x_j = cos((2j+1) pi/(2n+2)): c_k = (2/(n+1)) sum_j f_j T_k(x_j), c_0 halved."""
    j = np.arange(n + 1)
    A = (2.0 / (n + 1)) * np.cos((np.outer(j, 2 * j + 1) % (4 * n + 4)) * np.pi / (2 * n + 2))   # reduced argument
    A[0, :] *= 0.5
    return A


def cheb_coeffs(Fm, Fr, n, A=None):
    """Fm, Fr: (n+1, n+1) node mids (real) and radii -> coefficient mids and radii (rigorous)."""
    if A is None:
        A = coef_matrix(n)
    cm = A @ Fm @ A.T
    aA = np.abs(A)
    # cr: float error of the computed coefficients of I(Fm) (entries of A to ~1e-15 absolute AFTER argument reduction, sums of <= (n+1)^2 terms);
    # the node radii enter separately and UNIFORMLY: |I(F) - I(Fm)| <= Lambda_n^2 max Fr  (node_err)
    cr = 1e-12 * (aA @ np.abs(Fm) @ aA.T) + 1e-13 * float(np.max(np.abs(Fm)))
    node_err = lebesgue(n) ** 2 * float(np.max(Fr)) * (1 + 1e-9)
    return cm, cr + node_err / (n + 1) ** 2


def sup_inf_grid(cm, cr, m_max=3000, tol=None):
    n = cm.shape[0] - 1
    k = np.arange(n + 1)
    D = (np.abs(cm) * (k[:, None] ** 2 + k[None, :] ** 2)).sum()      # Markov bound for the computed polynomial p~
    tot = (np.abs(cm) + cr).sum()                                    # |p - p~| <= sum cr (uniform on the square)
    m = m_max if tol is None else int(min(m_max, max(64, math.ceil(math.pi * D / (2 * tol)) + 1)))
    g = np.arange(m)
    Tu = np.cos((np.outer(g, k) % (2 * (m - 1))) * np.pi / (m - 1))      # T_k(cos(g pi/(m-1))), reduced argument
    G = Tu @ cm @ Tu.T
    var = D * math.pi / (2 * (m - 1))
    err = cr.sum() + 1e-12 * tot + var
    return float(G.max()) + err, float(G.min()) - err, m, var, G, Tu


def certify(sub, rho, n_arc, n, tol_rel=0.003, do_b=True, b_prec=192, log=print, den=None, M=None, nbox=None, rbox=None, n_rho=None):
    t0 = time.time()
    out = dict(k=sub.k, s0=sub.s0, t0=sub.t0, h=sub.h, chart=sub.chart, rho=rho, n_arc=n_arc, n=n)
    if den is None:
        den = denominators(sub, rho)
    out['denominators'] = den
    log('  analyticity cover: %d boxes, ok=%s %s (%.1fs)' % (den['boxes'], den['ok'], {k_: '%.3g' % den[k_] for k_ in den['need']}, time.time() - t0))
    if not den['ok']:
        out['status'] = 'FAIL analyticity'
        return out
    t1 = time.time()
    if M is None:
        M, nbox, rbox = torus_max(sub, rho, n_arc)
    out['M'] = M; out['torus_boxes'] = nbox; out['torus_box_radius_s'] = rbox
    log('  torus: %d boxes radius %.1e (%.1fs): %s' % (nbox, rbox, time.time() - t1, ' '.join('%s<=%.1e' % kv for kv in M.items())))
    # nodes (a)
    t1 = time.time()
    u = np.cos((2 * np.arange(n + 1) + 1) * np.pi / (2 * n + 2))       # amend2: FIRST-kind nodes for part (a) (Rivlin (1.38))
    UU, VV = np.meshgrid(u, u, indexing='ij')
    eps = 2.0 ** -50
    S, T = sub.st(Disc(UU.ravel(), np.full(UU.size, eps)), Disc(VV.ravel(), np.full(VV.size, eps)))
    Fn = {}
    for o in sub.run(S, T):
        F = funcs_a(o)
        for kk, v in F.items():
            Fn.setdefault(kk, []).append(v)
    Fn = {kk: Disc(np.concatenate([v.m for v in L]), np.concatenate([v.r for v in L])) for kk, L in Fn.items()}
    out['secs_nodes_a'] = time.time() - t1
    log('  nodes (a): %d points (%.1fs)' % (S.m.size, time.time() - t1))
    res = {}
    coef = {}
    A1n = coef_matrix1(n)
    for kk in KEYS_A:
        v = Fn[kk]
        if np.any(np.abs(v.m.imag) > 0):
            raise ValueError('complex node value')
        Fm = v.m.real.reshape(n + 1, n + 1); Fr = v.r.reshape(n + 1, n + 1) + np.abs(v.m.imag).reshape(n + 1, n + 1)
        cm, cr = cheb_coeffs(Fm, Fr, n, A1n)
        coef[kk] = (cm, cr)
        e = interp_err(M[kk], rho, n)
        tol = max(tol_rel * float(np.max(np.abs(Fm))), 1e-300)
        sup, inf, m, var, G, Tu = sup_inf_grid(cm, cr, tol=tol)
        res[kk] = dict(sup=sup + e, inf=inf - e, interp_err=e, grid=m, var=var, node_max=float(np.max(Fm)), node_min=float(np.min(Fm)))
    out['a'] = res
    # derived certified constants (upper bounds of the norms = sqrt of sup of the squares)
    cert = {}
    for kk in ['K0', 'K1', 'K2', 'K3', 'Rc0', 'Rc1', 'Rc2', 'Rc3', 'R0c0', 'R0c1', 'R0c2', 'S0', 'S1', 'S2']:
        cert[kk] = math.sqrt(max(res[kk]['sup'], 0.0)) * (1 + 1e-12)
    cert['s_sup'] = max(abs(res['s']['sup']), abs(res['s']['inf']))
    cert['s_inf'] = res['s']['inf']
    cert['kappa_T8'] = max(res['T8']['sup'], 0.0) ** (1 / 8) * (1 + 1e-12)
    cert['kappa_T16'] = max(res['T16']['sup'], 0.0) ** (1 / 16) * (1 + 1e-12)
    cert['kappa_T32'] = max(res['T32']['sup'], 0.0) ** (1 / 32) * (1 + 1e-12)
    cert['kappa_T64'] = max(res['T64']['sup'], 0.0) ** (1 / 64) * (1 + 1e-12)
    cert['kappa'] = min(cert['kappa_T8'], cert['kappa_T16'], cert['kappa_T32'], cert['kappa_T64'])
    # Ricci eigenvalue certificates: scan c (lower) and C (upper)
    (q0m, q0r), (sm, sr) = coef['Q0'], coef['s']
    eQ, eS = res['Q0']['interp_err'], res['s']['interp_err']

    m = 1500
    nn = q0m.shape[0] - 1
    kk_ = np.arange(nn + 1)
    W2 = kk_[:, None] ** 2 + kk_[None, :] ** 2
    g_ = np.arange(m)
    Tu = np.cos((np.outer(g_, kk_) % (2 * (m - 1))) * np.pi / (m - 1))   # reduced argument
    GQ = Tu @ q0m @ Tu.T; GS = Tu @ sm @ Tu.T
    def Hmin_of(c):
        cm_ = q0m - (c / 2) * sm; cr_ = q0r + (c / 2) * sr
        err = cr_.sum() + 1e-12 * ((np.abs(cm_) + cr_).sum() + c * c) + (np.abs(cm_) * W2).sum() * np.pi / (2 * (m - 1)) + eQ + c / 2 * eS
        return float((GQ - (c / 2) * GS).min()) + c * c - err
    ok_lo = lambda c: Hmin_of(c) >= 0 and res['s']['inf'] / 4 - c >= 0
    ok_hi = lambda c: Hmin_of(c) >= 0 and c - res['s']['sup'] / 4 >= 0
    lo_c, hi_C = None, None
    for c in np.arange(0.0, 3.0, 0.0025):
        if ok_lo(c):
            lo_c = float(c)
        if hi_C is None and ok_hi(c):
            hi_C = float(c)
    if lo_c is not None:                                   # refine (each accepted value is itself a certificate)
        for c in np.arange(lo_c, lo_c + 0.0025, 2.5e-5):
            if ok_lo(c):
                lo_c = float(c)
    if hi_C is not None:
        for c in np.arange(hi_C, hi_C - 0.0025, -2.5e-5):
            if ok_hi(c):
                hi_C = float(c)
    cert['ric_lower'] = lo_c; cert['ric_upper'] = hi_C
    if hi_C is not None and lo_c is not None:
        cert['Ric_op'] = max(abs(hi_C), abs(lo_c))
        cert['m_R'] = cert['Ric_op'] + 2 * cert['Rc1'] + cert['Rc2']
    out['cert_a'] = cert
    log('  (a) certified on this sub-square: ' + ' '.join('%s=%.4g' % (kk, v) for kk, v in cert.items() if v is not None))
    log('      node max (float):            ' + ' '.join('%s=%.4g' % (kk, math.sqrt(max(res[kk]['node_max'], 0))) for kk in ['K0', 'K1', 'K2', 'K3', 'Rc3', 'R0c2', 'S2']))
    # (b)
    if do_b:
        t1 = time.time()
        old = ctx.prec
        ctx.prec = b_prec
        # exact Chebyshev nodes as arb balls (radius ~2^-b_prec): s = s0 + h/2 (1 + cos(j pi/n))
        hs = arb(sub.h) / 2
        ua = [arb.cos_pi(arb(j) / n) for j in range(n + 1)]
        sl = [arb(sub.s0) + hs * (1 + ua[i]) for i in range(n + 1) for j in range(n + 1)]
        tl = [arb(sub.t0) + hs * (1 + ua[j]) for i in range(n + 1) for j in range(n + 1)]
        ps, lp = [], []
        for i in range(0, len(sl), 2000):
            a_, b_ = arbpsi.psi_lap(sub.prob, sub.k, sl[i:i + 2000], tl[i:i + 2000])
            ps += list(a_); lp += list(b_)
        ctx.prec = old
        out['secs_nodes_b'] = time.time() - t1
        rb = {}
        for name, vals in (('psi', ps), ('lap', lp)):
            Fm = np.array([float(v.mid()) for v in vals]).reshape(n + 1, n + 1)
            Fr = np.array([float(v.rad()) + abs(float(v.mid())) * 2.0 ** -52 for v in vals]).reshape(n + 1, n + 1)
            cm, cr = cheb_coeffs(Fm, Fr, n)
            e = interp_err(M[name], rho, n)
            supc = float((np.abs(cm) + cr).sum())
            rb[name] = dict(sup=supc + e, interp_err=e, coef_sum=supc, node_max=float(np.max(np.abs(Fm))))
        out['b'] = rb
        log('  (b) nodes %d arb (%.1fs): sup|psi| <= %.3e (node max %.2e, interp err %.1e);  sup|Delta psi| <= %.3e (node max %.2e, interp err %.1e)' % (
            len(ps), time.time() - t1, rb['psi']['sup'], rb['psi']['node_max'], rb['psi']['interp_err'], rb['lap']['sup'], rb['lap']['node_max'], rb['lap']['interp_err']))
    if do_b and n_rho:
        t1 = time.time()
        old = ctx.prec
        ctx.prec = b_prec
        hs = arb(sub.h) / 2
        ua = cheb1_nodes_arb(n_rho)
        sl = [arb(sub.s0) + hs * (1 + ua[i]) for i in range(n_rho + 1) for j in range(n_rho + 1)]
        tl = [arb(sub.t0) + hs * (1 + ua[j]) for i in range(n_rho + 1) for j in range(n_rho + 1)]
        L = [[], [], [], []]
        for i in range(0, len(sl), 1000):
            N = arbrho.rho_norms(sub.prob, sub.k, sl[i:i + 1000], tl[i:i + 1000], which='psi')
            for j in range(4):
                L[j] += list(N[j])
        ctx.prec = old
        A1 = coef_matrix1(n_rho)
        rr = {}
        for j in range(4):
            Fm = np.array([float(v.mid()) for v in L[j]]).reshape(n_rho + 1, n_rho + 1)
            Fr = np.array([float(v.rad()) + abs(float(v.mid())) * 2.0 ** -52 for v in L[j]]).reshape(n_rho + 1, n_rho + 1)
            cm, cr = cheb_coeffs(Fm, Fr, n_rho, A1)
            e = interp_err(M['rho%d' % j], rho, n_rho)
            supc = float((np.abs(cm) + cr).sum())
            rr['rho%d' % j] = dict(sup_sq=supc + e, sup=math.sqrt(supc + e), interp_err=e, coef_sum=supc, node_max_sq=float(np.max(np.abs(Fm))))
        out['b3'] = rr; out['n_rho'] = n_rho; out['secs_nodes_b3'] = time.time() - t1
        log('  (b3) %d first-kind arb nodes (%.1fs): sup|nabla^j rho| <= %s  (node max %s)' % (
            len(sl), time.time() - t1, ' '.join('%.2e' % rr['rho%d' % j]['sup'] for j in range(4)),
            ' '.join('%.1e' % math.sqrt(rr['rho%d' % j]['node_max_sq']) for j in range(4))))
    out['secs'] = time.time() - t0
    out['status'] = 'ok'
    return out


def choose_n(M, rho, targets, n_max=160):
    """smallest n with interp_err(M_F) <= target_F for every F in targets."""
    for n in range(8, n_max + 1, 2):
        if all(interp_err(M[kk], rho, n) <= tg for kk, tg in targets.items()):
            return n
    return n_max


def presample(sub, k=9):
    """float scale of each function on the real sub-square (k x k points) -- used ONLY to set accuracy targets."""
    u = cheb_nodes(k - 1)
    UU, VV = np.meshgrid(u, u, indexing='ij')
    S, T = sub.st(Disc(UU.ravel()), Disc(VV.ravel()))
    o = sub.run(S, T)[0]
    F = funcs_a(o)
    sc = {kk: float(np.max(np.abs(v.m))) for kk, v in F.items()}
    return sc


def plan(sub, rho, n_arc, tol_rel=0.003, tol_b=1e-14, tol_Q0=1e-4, tol_rho=1e-29, n_max=160, log=print):
    t0 = time.time()
    den = denominators(sub, rho)
    if not den['ok']:
        return dict(den=den, ok=False)
    M, nbox, rbox = torus_max(sub, rho, n_arc)
    sc = presample(sub)
    tg = {kk: tol_rel * sc[kk] for kk in KEYS_A if kk not in ('Q0', 'T16', 'T32', 'T64')}
    tg['s'] = 1e-4
    tg['Q0'] = tol_Q0
    tg['psi'] = tol_b; tg['lap'] = tol_b
    tg = {kk: max(v, 1e-300) for kk, v in tg.items()}
    n = choose_n(M, rho, tg, n_max)
    worst = max(tg, key=lambda kk: interp_err(M[kk], rho, n) / tg[kk])
    tr = {'rho%d' % j: tol_rho for j in range(4)}
    n_rho = choose_n(M, rho, tr, n_max)
    log('  plan: n = %d (binding %s), n_rho = %d (%.1fs)' % (n, worst, n_rho, time.time() - t0))
    return dict(den=den, ok=True, M=M, nbox=nbox, rbox=rbox, n=n, n_rho=n_rho, targets=tg, scale=sc)


def certify_b3(sub, rho, n_arc, tol_rho=1e-24, n_max=200, b_prec=192, log=print):
    """b3-only mode (production A): analyticity cover, torus crude bounds of |nabla^j psi|^2 only, first-kind arb nodes,
    certified sup |nabla^j rho| (j = 0..3) on the sub-square."""
    t0 = time.time()
    den = denominators(sub, rho)
    out = dict(k=sub.k, s0=sub.s0, t0=sub.t0, h=sub.h, chart=sub.chart, rho=rho, n_arc=n_arc, denominators=den)
    if not den['ok']:
        out['status'] = 'FAIL analyticity'
        log('  b3: ANALYTICITY FAIL %s' % den)
        return out
    U = ellipse_arcs(rho, n_arc // 2, True); V = ellipse_arcs(rho, n_arc, False)
    Uu, Vv = prod_boxes(U, V)
    S, T = sub.st(Uu, Vv)
    M = {}
    for o in sub.run(S, T, want=('abreu',)):
        for j in range(4):
            M['rho%d' % j] = max(M.get('rho%d' % j, 0.0), _ub(o['abreu'][j]))
    out['M'] = M; out['torus_boxes'] = len(S.m)
    if not all(np.isfinite(v) for v in M.values()):
        out['status'] = 'FAIL torus (infinite crude bound)'
        log('  b3: torus bound infinite %s' % M)
        return out
    n_rho = choose_n(M, rho, {'rho%d' % j: tol_rho for j in range(4)}, n_max)
    old = ctx.prec
    ctx.prec = b_prec
    hs = arb(sub.h) / 2
    ua = cheb1_nodes_arb(n_rho)
    sl = [arb(sub.s0) + hs * (1 + ua[i]) for i in range(n_rho + 1) for j in range(n_rho + 1)]
    tl = [arb(sub.t0) + hs * (1 + ua[j]) for i in range(n_rho + 1) for j in range(n_rho + 1)]
    L = [[], [], [], []]
    for i in range(0, len(sl), 1000):
        N = arbrho.rho_norms(sub.prob, sub.k, sl[i:i + 1000], tl[i:i + 1000], which='psi')
        for j in range(4):
            L[j] += list(N[j])
    ctx.prec = old
    A1 = coef_matrix1(n_rho)
    rr = {}
    for j in range(4):
        Fm = np.array([float(v.mid()) for v in L[j]]).reshape(n_rho + 1, n_rho + 1)
        Fr = np.array([float(v.rad()) + abs(float(v.mid())) * 2.0 ** -52 for v in L[j]]).reshape(n_rho + 1, n_rho + 1)
        cm_, cr = cheb_coeffs(Fm, Fr, n_rho, A1)
        e = interp_err(M['rho%d' % j], rho, n_rho)
        supc = float((np.abs(cm_) + cr).sum())
        rr['rho%d' % j] = dict(sup=math.sqrt(supc + e) * (1 + 1e-12), interp_err=e, node_max=math.sqrt(float(np.max(np.abs(Fm)))))
    out['b3'] = rr; out['n_rho'] = n_rho; out['secs'] = time.time() - t0; out['status'] = 'ok'
    log('  b3: n_rho %d, %d torus boxes (%.0fs): sup|nabla^j rho| <= %s' % (n_rho, len(S.m), out['secs'], ' '.join('%.2e' % rr['rho%d' % j]['sup'] for j in range(4))))
    return out
