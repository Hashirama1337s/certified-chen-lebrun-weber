"""cmtoric.py -- seat E: certified SCALAR quantities of the toric Kahler metric g_K (u = u_G + f) on (s,t) sub-squares of the
vertex elements, by Chebyshev models with arb coefficients (cm.py), in the vertex y-chart y = (l_a, l_b).
Everything is polynomial algebra in models except the validated inverses 1/l_r (r != a, b) and 1/det m.

Problem data are arb BALLS: a* = the root of Liu's polynomial in (1, 2) (isolated by FLINT), gamma = 2/a*, polygon, Donaldson's
affine S (closed 3x3 system with exact Gauss integrals of polynomials, in arb).  f = seat R's coefficients (exact dyadic rationals).

Formulas (y-chart, U = (Hess_y u)^-1 = m^-1 Dg, Dg = diag(2y1, 2y2), m = I + Dg R, R = sum_{r != a,b} (1/2) n_r n_r^T / l_r
+ Ai^T Hess_x f Ai; all point-checked against the tensor engine in t_ricU.py, t_scalarU.py):
  s = Scal = -U^{ij}_{,ij}                       (Abreu);  rho = s - S
  A^j_i = -(1/2) U^{jl}_{,il}                    Ricci endomorphism (eigenvalues = Ricci eigenvalues, each twice); tr A = s/2
  |Ric|^2 = 2 tr A^2;  |Ric_0|^2 = |Ric|^2 - s^2/4;  lambda_min >= mu  <=>  s/4 - mu >= 0 and det(A - mu) >= 0
  |Riem|^2 = U^{ij}_{,kl} U^{kl}_{,ij}
  Delta phi = -(U^{ij} phi_j)_{,i}  (positive Laplacian, invariant phi);  |d phi|^2 = U^{ij} phi_i phi_j
  |dS|^2 = c^T U c  (c = dS/dy, constant);  |nabla dS|^2 = (1/2) (U^{ij} q_j)_{,i} + c_i A^i_j U^{jk} c_k,  q = |dS|^2
  |nabla^2 dS|^2 <= |Riem|^2 |dS|^2   (nabla^2 K = Riem . K for the Killing field K = J grad S)
  admissibility (K1): m11 > 0, det m > 0 (then Hess_y u = Dg^-1 m is positive definite in the interior, and the Guillemin
  boundary behaviour is that of u_G since m is smooth and nondegenerate up to the boundary)
  dV = dx dtheta: int_M F dV = (2 pi)^2 int F |dy/d(s,t)| ds dt,  |dy/d(s,t)| = 4 s t (A B - s^2 t^2 dA dB) >= 0.
"""
import math, time, json
from math import factorial
import numpy as np
from flint import arb, arb_mat, fmpz_poly, ctx
import common
import cm
from cm import CM, Jet2
import fbox

LIU_POLY = [1, 5, 6, -12, -34, -42, 0, 20, 21, -7, -6]      # high -> low (seat R clwr.py / Nifa-Liu)
PAGE_POLY = [1, 4, -6, 12, -3]                                # nu^4 + 4 nu^3 - 6 nu^2 + 12 nu - 3


def root_ball(coeffs_hi_lo, lo, hi):
    p = fmpz_poly(list(reversed(coeffs_hi_lo)))
    rts = [z.real for z, m in p.complex_roots() if z.imag.contains(0) and z.real > lo and z.real < hi]
    assert len(rts) == 1
    return rts[0]


class APoly:
    """Delzant polygon with arb data (same conventions as seat R's clwr.Poly)."""

    def __init__(self, nus, lams):
        self.nu = [tuple(int(v) for v in n) for n in nus]
        self.lam = [arb(l) for l in lams]
        self.d = d = len(self.nu)
        V = []
        for k in range(d):
            a, b = self.nu[k], self.nu[(k + 1) % d]
            det = a[0] * b[1] - a[1] * b[0]
            assert abs(det) == 1
            r1, r2 = -self.lam[k], -self.lam[(k + 1) % d]
            V.append(((r1 * b[1] - a[1] * r2) / det, (a[0] * r2 - b[0] * r1) / det))
        self.V = V
        self.O = (sum((v[0] for v in V), arb(0)) / d, sum((v[1] for v in V), arb(0)) / d)
        self.Mid = [((V[r - 1][0] + V[r][0]) / 2, (V[r - 1][1] + V[r][1]) / 2) for r in range(d)]

    def ell(self, r, x1, x2):
        return self.nu[r][0] * x1 + self.nu[r][1] * x2 + self.lam[r]

    def elem(self, k):
        a, b = k, (k + 1) % self.d
        return dict(k=k, a=a, b=b, c=self.V[k], Mb=self.Mid[b], Ma=self.Mid[a], O=self.O)


def gauss01(n):
    xs, ws = [], []
    for k in range(n):
        x, w = arb.legendre_p_root(n, k, weight=True)
        xs.append((x + 1) / 2); ws.append(w / 2)
    return xs, ws


def S_affine(P, nq=6):
    """Donaldson's S = a + b x1 + c x2: int_P S l dx = 2 oint l dsigma for l in {1, x1, x2} (Gauss exact for these polynomials)."""
    xs, ws = gauss01(nq)
    fs = [lambda x1, x2: arb(1), lambda x1, x2: x1, lambda x1, x2: x2]

    def pint(F):
        tot = arb(0)
        for k in range(P.d):
            E = P.elem(k)
            c, Mb, O, Ma = E['c'], E['Mb'], E['O'], E['Ma']
            for sg, wa in zip(xs, ws):
                for ta, wb in zip(xs, ws):
                    w00, w10, w11, w01 = (1 - sg) * (1 - ta), sg * (1 - ta), sg * ta, (1 - sg) * ta
                    x1 = w00 * c[0] + w10 * Mb[0] + w11 * O[0] + w01 * Ma[0]
                    x2 = w00 * c[1] + w10 * Mb[1] + w11 * O[1] + w01 * Ma[1]
                    dxs = [(1 - ta) * (Mb[i] - c[i]) + ta * (O[i] - Ma[i]) for i in range(2)]
                    dxt = [(1 - sg) * (Ma[i] - c[i]) + sg * (O[i] - Mb[i]) for i in range(2)]
                    J = dxs[0] * dxt[1] - dxs[1] * dxt[0]
                    assert J > 0 or J < 0
                    tot += wa * wb * (J if J > 0 else -J) * F(x1, x2)
        return tot

    def bint(F):
        tot = arb(0)
        for r in range(P.d):
            A, B = P.V[r - 1], P.V[r]
            n = P.nu[r]
            L = ((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2).sqrt() / arb(n[0] ** 2 + n[1] ** 2).sqrt()
            for x, w in zip(xs, ws):
                tot += w * L * F(A[0] + x * (B[0] - A[0]), A[1] + x * (B[1] - A[1]))
        return tot
    Am = arb_mat(3, 3, [pint(lambda a, b, i=i, j=j: fs[i](a, b) * fs[j](a, b)) for i in range(3) for j in range(3)])
    rhs = arb_mat(3, 1, [2 * bint(fs[i]) for i in range(3)])
    sol = Am.solve(rhs)
    area = pint(fs[0])
    return [sol[0, 0], sol[1, 0], sol[2, 0]], area


class AProb:
    """case: 'clw' (pentagon at the a* ball), 'page' (seat C's trapezoid at the nu ball), 'square', 'triangle'."""

    def __init__(self, case, solfile=None, prec=256, Cfloat=None, B_float=448):
        ctx.prec = prec
        self.case = case
        if case == 'clw':
            self.astar = root_ball(LIU_POLY, 1, 2)
            self.P = APoly([(1, 0), (1, 1), (0, 1), (-1, 0), (0, -1)], [arb(1), 2 / self.astar, arb(1), arb(1), arb(1)])
        elif case == 'page':
            nu = root_ball(PAGE_POLY, 0, 1)
            gam = 2 - 3 * ((1 - nu) / (1 + nu)) ** 2
            nus = [(1, 0), (1, 1), (0, 1), (-1, -1)]
            lam0 = [arb(1), gam, arb(1), arb(1)]
            self.P = APoly(nus, [(lam0[r] + arb(nus[r][0] + nus[r][1]) / 2) / (arb(3) / 2) for r in range(4)])
        elif case == 'square':
            self.P = APoly([(1, 0), (0, 1), (-1, 0), (0, -1)], [1, 1, 1, 1])
        else:
            self.P = APoly([(1, 0), (0, 1), (-1, -1)], [1, 1, 1])
        self.lam = self.P.lam
        self.Sabc, self.area = S_affine(self.P)
        self.c = arb(1)                     # conformal constant: unused for g_K (cpow = 0)
        self.Cint = None
        if solfile is not None:
            import clwr
            C, d = clwr.load(solfile)
            self.Cint, self.B, self.M = C, int(d['B']), C.shape[0] - 1
            self.sol_meta = {k: d[k] for k in ('gamma', 'Sabc') if k in d}
        elif Cfloat is not None:
            # float coefficient matrix (exact binary64 values) rounded to 2^-B_float: a fixed polynomial (documented approximation of the input)
            self.B = B_float
            self.Cint = np.array([[int(round(float(v) * 2.0 ** 60)) << (B_float - 60) for v in row] for row in Cfloat], dtype=object)
            self.M = self.Cint.shape[0] - 1

    def Cmat(self, cplx=False):
        M = self.M
        return arb_mat(M + 1, M + 1, [int(v) for v in self.Cint.ravel()]) * (arb(2) ** (-self.B))


def _opcoef(Ai, aa, bb, mu):
    """integer coefficients {(g1, g2): c} of the x-derivative operator d_y^mu (Ai d_x)_aa (Ai d_x)_bb (x = Ai (y - lam))."""
    def pmul(p, q):
        out = {}
        for (a1, a2), c1 in p.items():
            for (b1, b2), c2 in q.items():
                k = (a1 + b1, a2 + b2); out[k] = out.get(k, 0) + c1 * c2
        return out
    op = {(0, 0): 1}
    lin = lambda m_: {(1, 0): Ai[0][m_], (0, 1): Ai[1][m_]}
    for _ in range(mu[0]):
        op = pmul(op, lin(0))
    for _ in range(mu[1]):
        op = pmul(op, lin(1))
    op = pmul(pmul(op, lin(aa)), lin(bb))
    return {k: v for k, v in op.items() if v}


def build(prob, k, s0, t0, h, K=4, m_taylor=50, R_f=0.6, log=None):
    """U jets (order K in dy, CM coefficients) on the sub-square [s0, s0+h] x [t0, t0+h] of element k, plus chart data."""
    P = prob.P
    E = P.elem(k)
    a, b = E['a'], E['b']
    Am = [list(P.nu[a]), list(P.nu[b])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    al1, al2 = P.ell(a, *E['Mb']), P.ell(a, *P.O)
    be1, be2 = P.ell(b, *E['Ma']), P.ell(b, *P.O)
    u, v = CM.var(0), CM.var(1)
    hs = arb(h) / 2
    s = u.scale(hs) + (arb(s0) + hs); t = v.scale(hs) + (arb(t0) + hs)
    s2, t2 = s * s, t * t
    A = t2.scale(al2 - al1) + al1; Bm = s2.scale(be2 - be1) + be1
    y = [s2 * A, t2 * Bm]
    Gjac = A * Bm - s2 * t2 * ((al2 - al1) * (be2 - be1))
    jac = (s * t * Gjac).scale(4)          # dy/d(s,t) = 4 s t Gjac; s, t >= 0 on [0,1]^2 and Gjac > 0 (certified) => jac >= 0
    x = [(y[0] - P.lam[a]).scale(Ai[i][0]) + (y[1] - P.lam[b]).scale(Ai[i][1]) for i in range(2)]
    R = [[None, None], [None, None]]

    def jadd(J, o):
        return o if J is None else J + o
    for r in range(P.d):
        if r in (a, b):
            continue
        nv = [Ai[0][i] * P.nu[r][0] + Ai[1][i] * P.nu[r][1] for i in range(2)]
        l0 = x[0].scale(P.nu[r][0]) + x[1].scale(P.nu[r][1]) + P.lam[r]
        inv = l0.inv()
        pw = [None, inv]
        for _ in range(K):
            pw.append(pw[-1] * inv)
        c = {}
        for m1 in range(K + 1):
            for m2 in range(K + 1 - m1):
                n_ = m1 + m2
                co = (-1) ** n_ * math.comb(n_, m1) * nv[0] ** m1 * nv[1] ** m2
                if co:
                    c[(m1, m2)] = pw[n_ + 1].scale(co)
        il = Jet2(c, K)
        for i in range(2):
            for j in range(2):
                if nv[i] * nv[j]:
                    R[i][j] = jadd(R[i][j], il.scale(arb(nv[i] * nv[j]) / 2))
    if prob.Cint is not None:
        c0 = (float(x[0].p.coeffs()[cm.D + cm.D * cm.W].mid()), float(x[1].p.coeffs()[cm.D + cm.D * cm.W].mid()))
        dx = [x[0] - arb(c0[0]), x[1] - arb(c0[1])]
        rp = 1.02 * max(dx[0].sup(), dx[1].sup()) + 1e-6
        # f is a polynomial of total degree M: its Taylor expansion at c0 to order M is EXACT (no Cauchy remainder).
        m = prob.M
        pa = fbox.FPatch(prob, c0, rp, m=m, R=max(R_f, 2 * rp), kd=K + 2)
        p1 = [CM.const(1)]; p2 = [CM.const(1)]
        for _ in range(m):
            p1.append(p1[-1] * dx[0]); p2.append(p2[-1] * dx[1])
        s1, s2 = arb(dx[0].sup()), arb(dx[1].sup())
        n1 = [s1 ** i for i in range(m + 1)]; n2 = [s2 ** i for i in range(m + 1)]
        SKIP = arb(10) ** (-40)
        F = {}
        for g1 in range(K + 3):
            for g2 in range(K + 3 - g1):
                if g1 + g2 < 2:
                    continue
                mm = m - g1 - g2
                tot = None
                skipped = arb(0)
                for b1 in range(mm + 1):
                    inner = None
                    for b2 in range(mm + 1 - b1):
                        co = pa.a[g1 + b1, g2 + b2] * (factorial(g1 + b1) // factorial(b1)) * (factorial(g2 + b2) // factorial(b2))
                        bnd = co.abs_upper() * n1[b1] * n2[b2]          # |term| <= |co| ||dx1||^b1 ||dx2||^b2 on the square
                        if bnd < SKIP:
                            skipped += bnd
                            continue
                        tt = p2[b2].scale(co)
                        inner = tt if inner is None else inner + tt
                    if inner is None:
                        continue
                    tt = p1[b1] * inner
                    tot = tt if tot is None else tot + tt
                tot.rem = (tot.rem + skipped).abs_upper()
                F[(g1, g2)] = tot                                   # exact expansion: nothing to add
        for aa in range(2):
            for bb in range(aa, 2):
                c = {}
                for m1 in range(K + 1):
                    for m2 in range(K + 1 - m1):
                        op = _opcoef(Ai, aa, bb, (m1, m2))
                        tot = None
                        for g, co in op.items():
                            tt = F[g].scale(arb(co) / (factorial(m1) * factorial(m2)))
                            tot = tt if tot is None else tot + tt
                        if tot is not None:
                            c[(m1, m2)] = tot
                Jf = Jet2(c, K)
                R[aa][bb] = jadd(R[aa][bb], Jf)
                if aa != bb:
                    R[bb][aa] = jadd(R[bb][aa], Jf)
    zero = Jet2({}, K)
    R = [[zero if R[i][j] is None else R[i][j] for j in range(2)] for i in range(2)]

    def ymul(i, J, fac=2):
        """(fac * Y_i) * J with Y_i = y_i + dy_i."""
        c = {}
        for mth, vv in J.c.items():
            tt = (y[i] * vv).scale(fac)
            c[mth] = (c[mth] + tt) if mth in c else tt
            m2_ = (mth[0] + 1, mth[1]) if i == 0 else (mth[0], mth[1] + 1)
            if m2_[0] + m2_[1] <= J.K:
                t2_ = vv.scale(fac)
                c[m2_] = (c[m2_] + t2_) if m2_ in c else t2_
        return Jet2(c, J.K)
    m11 = ymul(0, R[0][0]) + 1; m12 = ymul(0, R[0][1]); m21 = ymul(1, R[1][0]); m22 = ymul(1, R[1][1]) + 1
    detm = m11 * m22 - m12 * m21
    idet = detm.recip()
    U00 = ymul(0, m22 * idet); U01 = -ymul(1, m12 * idet); U11 = ymul(1, m11 * idet)
    U = [[U00, U01], [U01, U11]]
    cS = [prob.Sabc[1] * Ai[0][i] + prob.Sabc[2] * Ai[1][i] for i in range(2)]
    S0 = x[0].scale(prob.Sabc[1]) + x[1].scale(prob.Sabc[2]) + prob.Sabc[0]
    return dict(U=U, m11=m11.val(), m22=m22.val(), detm=detm.val(), jac=jac, Gjac=Gjac, y=y, x=x, cS=cS, S0=S0, K=K, Ai=Ai, h=h)


def quantities(B_):
    """CM scalars from the U jets."""
    U, K, cS = B_['U'], B_['K'], B_['cS']
    d2 = {}
    for i in range(2):
        for j in range(2):
            for kk in range(2):
                for l in range(2):
                    d2[(i, j, kk, l)] = U[i][j].d(kk).d(l)                # jets of U^{ij}_{,kl} (order K-2)
    v = lambda J: J.val()
    Q = {}
    # scalar curvature jet (order K-2) and rho jet
    sJ = -(d2[(0, 0, 0, 0)] + d2[(0, 1, 0, 1)].scale(2) + d2[(1, 1, 1, 1)])
    Sj = Jet2({(0, 0): B_['S0'], (1, 0): CM.const(cS[0]), (0, 1): CM.const(cS[1])}, sJ.K)
    rJ = sJ - Sj
    Q['s'] = v(sJ); Q['rho'] = v(rJ)
    A = [[-(v(d2[(j, 0, i, 0)]) + v(d2[(j, 1, i, 1)])).scale(arb(1) / 2) for i in range(2)] for j in range(2)]      # A[j][i] = A^j_i
    Q['A'] = A
    Q['trA2'] = A[0][0] * A[0][0] + (A[0][1] * A[1][0]).scale(2) + A[1][1] * A[1][1]
    Q['detA'] = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    Q['Ric2'] = Q['trA2'].scale(2)
    Q['Ric0_2'] = Q['Ric2'] - (Q['s'] * Q['s']).scale(arb(1) / 4)
    T = None
    for i in range(2):
        for j in range(2):
            for kk in range(2):
                for l in range(2):
                    tt = v(d2[(i, j, kk, l)]) * v(d2[(kk, l, i, j)])
                    T = tt if T is None else T + tt
    Q['Riem2'] = T
    # rho derivatives
    Uv = [[v(U[i][j]) for j in range(2)] for i in range(2)]
    dU = lambda i, j, kk: v(U[i][j].d(kk))

    def lap_pos(J):
        """positive Laplacian -(U^ij phi_j)_,i of an invariant function given by a jet (order >= 2)."""
        tot = None
        for i in range(2):
            for j in range(2):
                tt = Uv[i][j] * v(J.d(i).d(j)) + dU(i, j, i) * v(J.d(j))
                tot = tt if tot is None else tot + tt
        return -tot

    def grad2(J):
        tot = None
        for i in range(2):
            for j in range(2):
                tt = Uv[i][j] * (v(J.d(i)) * v(J.d(j)))
                tot = tt if tot is None else tot + tt
        return tot
    Q['lap_rho'] = lap_pos(rJ)
    Q['drho2'] = grad2(rJ)
    # S (affine) quantities
    qJ = None
    for i in range(2):
        for j in range(2):
            tt = U[i][j].scale(cS[i] * cS[j])
            qJ = tt if qJ is None else qJ + tt
    Q['dS2'] = v(qJ)
    ric = None
    for i in range(2):
        for j in range(2):
            for kk in range(2):
                tt = (A[i][j] * Uv[j][kk]).scale(cS[i] * cS[kk])
                ric = tt if ric is None else ric + tt
    Q['ndS2'] = (-lap_pos(qJ)).scale(arb(1) / 2) + ric
    Q['n2dS2_bound'] = Q['Riem2'] * Q['dS2']
    # tr of U in x coordinates (U_x = Ai U_y Ai^T >= 0, so lambda_max(U_x) <= tr U_x)
    Ai = B_['Ai']
    trU = None
    for a_ in range(2):
        for b_ in range(2):
            co = Ai[0][a_] * Ai[0][b_] + Ai[1][a_] * Ai[1][b_]
            if co:
                tt = Uv[a_][b_].scale(co)
                trU = tt if trU is None else trU + tt
    Q['trUx'] = trU
    Q['S'] = B_['S0']
    # x-chart Gram integrands (x = Ai (y - lam)): (U_x)_ab = Ai_ak Ai_bl U^kl ; (A U)_x sym: G2 integrand
    # |nabla dx_a|^2 = (1/2) L(U^aa) + Ric(grad x_a, grad x_a) (Bochner for the Killing field J grad x_a); the divergence term
    # integrates to 0 over P (U nu = 0 on the facets), so int_P <nabla dx_a, nabla dx_b> dx = int_P (A U)_x,(ab) dx (symmetrised).
    AU = [[None, None], [None, None]]
    for i in range(2):
        for j in range(2):
            tt = None
            for l in range(2):
                u_ = A[i][l] * Uv[l][j]
                tt = u_ if tt is None else tt + u_
            AU[i][j] = tt
    Ux = {}; AUx = {}
    for a_, b_ in ((0, 0), (0, 1), (1, 1)):
        acc1 = None; acc2 = None
        for kk in range(2):
            for l in range(2):
                co = Ai[a_][kk] * Ai[b_][l]
                if not co:
                    continue
                u1 = Uv[kk][l].scale(co)
                u2 = (AU[kk][l] + AU[l][kk]).scale(arb(co) / 2)
                acc1 = u1 if acc1 is None else acc1 + u1
                acc2 = u2 if acc2 is None else acc2 + u2
        Ux[(a_, b_)] = acc1 if acc1 is not None else CM.const(0)
        AUx[(a_, b_)] = acc2 if acc2 is not None else CM.const(0)
    Q['Ux'] = Ux; Q['AUx'] = AUx
    return Q


def certify_sub(prob, k, s0, t0, h, K=4, log=print, mus=None):
    t_0 = time.time()
    B_ = build(prob, k, s0, t0, h, K=K)
    t_1 = time.time()
    Q = quantities(B_)
    out = dict(k=k, s0=s0, t0=t0, h=h, K=K, degree=cm.D, prec=ctx.prec)
    bd = {}
    for name in ('m11', 'm22', 'detm'):
        bd[name] = B_[name].bounds()
    for name in ('s', 'Riem2', 'Ric2', 'Ric0_2', 'dS2', 'ndS2', 'n2dS2_bound', 'detA', 'trUx', 'S'):
        bd[name] = Q[name].bounds()
    out['bounds'] = {kk: [float(vv[0]), float(vv[1])] for kk, vv in bd.items()}
    out['rem'] = {name: float(Q[name].rem.abs_upper()) for name in ('s', 'rho', 'lap_rho', 'Riem2', 'Ric2', 'ndS2', 'drho2')}
    out['rem']['detm'] = float(B_['detm'].rem.abs_upper())
    # Ricci eigenvalues: lambda_min >= mu <= inf(s/4 - mu) >= 0 and inf(detA - mu trA + mu^2) >= 0 (tr A = s/2);
    # lambda_max <= C <= inf(C - s/4) >= 0 and the same quadratic >= 0.  Bisection; every accepted value is a certificate.
    half_s = Q['s'].scale(arb(1) / 2)

    def okq(mu):
        Hm = Q['detA'] - half_s.scale(arb(mu)) + arb(mu) ** 2
        return Hm.bounds()[0] >= 0
    best_lo, best_hi = None, None
    lo, hi = 0.0, bd['s'][0] / 4
    if hi > 0 and okq(0.0):
        best_lo = 0.0
        for _ in range(24):
            mid = (lo + hi) / 2
            if okq(mid):
                best_lo = lo = mid
            else:
                hi = mid
    lo, hi = bd['s'][1] / 4, max(bd['s'][1], 1.0) * 4
    if okq(hi):
        best_hi = hi
        for _ in range(24):
            mid = (lo + hi) / 2
            if mid >= bd['s'][1] / 4 and okq(mid):
                best_hi = hi = mid
            else:
                lo = mid
    out['ric_lower'], out['ric_upper'] = best_lo, best_hi
    # integrals over the sub-square piece of M: (2 pi)^2 (h/2)^2 int F jac du dv
    fac = (2 * arb.pi()) ** 2 * (arb(h) / 2) ** 2
    ints = {}
    out['Gjac_inf'] = B_['Gjac'].bounds()[0]
    if not out['Gjac_inf'] > 0:
        raise ValueError('Jacobian factor not certified positive')
    terms = [('vol', CM.const(1)), ('rho2', Q['rho'] * Q['rho']), ('lap2', Q['lap_rho'] * Q['lap_rho']), ('drho2', Q['drho2'])]
    terms += [('Ux_%d%d' % ab, Q['Ux'][ab]) for ab in ((0, 0), (0, 1), (1, 1))]
    terms += [('AUx_%d%d' % ab, Q['AUx'][ab]) for ab in ((0, 0), (0, 1), (1, 1))]
    fac_dx = (arb(h) / 2) ** 2                       # dx = dy = jac (h/2)^2 du dv (no (2 pi)^2)
    for name, F in terms:
        I, r4 = (F * B_['jac']).integral()
        f_ = fac if name in ('vol', 'rho2', 'lap2', 'drho2') else fac_dx
        ball = f_ * I
        rad = float((f_ * r4).abs_upper()) + float(ball.rad())
        mf = float(ball.mid())
        ints[name] = dict(mid=repr(mf), rad=(rad + abs(mf) * 2.0 ** -52) * (1 + 1e-9) + 1e-300,       # exact value in [mid - rad, mid + rad]
                          upper=float((ball.abs_upper() + f_ * r4).abs_upper()))
    out['integrals'] = ints
    out['rho_sup_lower'] = Q['rho'].sup_lower()
    out['sup'] = dict(rho=Q['rho'].sup(), lap_rho=Q['lap_rho'].sup(), drho=math.sqrt(max(Q['drho2'].bounds()[1], 0.0)))
    cert = dict(K0=math.sqrt(max(bd['Riem2'][1], 0)), Rc0=math.sqrt(max(bd['Ric2'][1], 0)), R0_0=math.sqrt(max(bd['Ric0_2'][1], 0)),
                s_inf=bd['s'][0], s_sup=bd['s'][1], dS=math.sqrt(max(bd['dS2'][1], 0)), ndS=math.sqrt(max(bd['ndS2'][1], 0)),
                n2dS=math.sqrt(max(bd['n2dS2_bound'][1], 0)), detm_inf=bd['detm'][0], m11_inf=bd['m11'][0], m22_inf=bd['m22'][0],
                ric_lower=best_lo, ric_upper=best_hi, U_sup=bd['trUx'][1], S_min=bd['S'][0], S_max=bd['S'][1])
    out['cert'] = cert
    out['secs_build'] = t_1 - t_0; out['secs'] = time.time() - t_0
    log('  e%d [%.3f,%.3f]x[%.3f,%.3f] (%.0fs build, %.0fs total): ' % (k, s0, s0 + h, t0, t0 + h, t_1 - t_0, out['secs'])
        + ' '.join('%s=%.5g' % (kk, vv) for kk, vv in cert.items() if vv is not None))
    log('      sup|rho| <= %.2e  sup|Delta rho| <= %.2e  sup|d rho| <= %.2e ;  int rho^2 <= %.2e  int (Delta rho)^2 <= %.2e  vol piece %.6f'
        % (out['sup']['rho'], out['sup']['lap_rho'], out['sup']['drho'], ints['rho2']['upper'], ints['lap2']['upper'], float(ints['vol']['upper'])))
    log('      sup|rho| >= %.3e (certified lower)  Gjac >= %.4f' % (out['rho_sup_lower'], out['Gjac_inf']))
    return out


def gram_P(P, nq=6):
    """G_P = int_P e e^T dx, e = (1, x1, x2), exact Gauss (degree <= 5 per variable on the bilinear elements), arb 3x3."""
    xs, ws = gauss01(nq)
    G = arb_mat(3, 3)
    for k in range(P.d):
        E = P.elem(k)
        c, Mb, O, Ma = E['c'], E['Mb'], E['O'], E['Ma']
        for sg, wa in zip(xs, ws):
            for ta, wb in zip(xs, ws):
                w00, w10, w11, w01 = (1 - sg) * (1 - ta), sg * (1 - ta), sg * ta, (1 - sg) * ta
                x1 = w00 * c[0] + w10 * Mb[0] + w11 * O[0] + w01 * Ma[0]
                x2 = w00 * c[1] + w10 * Mb[1] + w11 * O[1] + w01 * Ma[1]
                dxs = [(1 - ta) * (Mb[i] - c[i]) + ta * (O[i] - Ma[i]) for i in range(2)]
                dxt = [(1 - sg) * (Ma[i] - c[i]) + sg * (O[i] - Mb[i]) for i in range(2)]
                J = dxs[0] * dxt[1] - dxs[1] * dxt[0]
                J = J if J > 0 else -J
                e = [arb(1), x1, x2]
                for i in range(3):
                    for j in range(3):
                        G[i, j] += wa * wb * J * e[i] * e[j]
    return G
