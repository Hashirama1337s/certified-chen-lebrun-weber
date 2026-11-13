"""arbpsi.py -- seat E: high-precision (arb) node values of the Abreu residual psi = Scal(u) - S and of its positive Laplacian
Delta psi = -d_i(U^ij d_j psi), at real points of a vertex element (s, t), vectorised over points (object arrays).

Vertex chart y = (l_a, l_b) (affine, unimodular: x = Ai (y - (lam_a, lam_b))).  Hess_y u = diag(1/(2 l_a), 1/(2 l_b)) + R,
R = sum_{r != a,b} (1/2) n_r n_r^T / l_r + Ai^T Hess_x f Ai  (n_r = Ai^T nu_r), so
U = (Hess_y u)^-1 = m^-1 Dg,  Dg = diag(2 l_a, 2 l_b), m = I + Dg R   (analytic through the rods l_a = 0, l_b = 0 and the vertex).
Scal = -d_i d_j U^ij (Abreu; affine invariant), Delta psi = -d_i(U^ij d_j psi) (measure dx = dy).  Truncated bivariate Taylor jets in
dy at each point (seat R's met1_core.Jet, arb balls): U to order 4, psi to order 2.  f derivatives to order 6 at the exact point
(met1_core.fderivs, point recurrence).  Outputs are arb balls (rigorous for the exact point (s, t) given as a binary64 float).
Same method family as seat X3 (x3.py), written independently here."""
import numpy as np
from math import factorial
from flint import arb, ctx
import common
import met1_core as mc
from met1_core import Jet, obj, monos


def _var(v0, which, K):
    n = len(v0)
    ms, ix, _ = monos(K)
    c = [None] * len(ms)
    c[0] = v0
    c[ix[(1, 0)] if which == 0 else ix[(0, 1)]] = obj([arb(1)] * n)
    return Jet(c, K)


def psi_lap(prob, k, s, t, K=4, rad=0.0):
    """s, t: lists of floats.  Returns (psi, lap) as object arrays of arb."""
    P = prob.P
    E = P.elem(k)
    a, b = E['a'], E['b']
    Am = [list(P.nu[a]), list(P.nu[b])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    if len(s) and isinstance(s[0], arb):
        sa, ta = obj(list(s)), obj(list(t))                     # exact arb balls (e.g. Chebyshev nodes)
    else:
        sa = obj([arb(float(v), rad) for v in s]); ta = obj([arb(float(v), rad) for v in t])
    al1, al2 = mc.A(P.ell(a, *E['Mb'])), mc.A(P.ell(a, *P.O))
    be1, be2 = mc.A(P.ell(b, *E['Ma'])), mc.A(P.ell(b, *P.O))
    la = sa * sa * ((1 - ta * ta) * al1 + ta * ta * al2)
    lb = ta * ta * ((1 - sa * sa) * be1 + sa * sa * be2)
    lam = prob.lam
    x1 = (la - lam[a]) * Ai[0][0] + (lb - lam[b]) * Ai[0][1]
    x2 = (la - lam[a]) * Ai[1][0] + (lb - lam[b]) * Ai[1][1]
    Y1, Y2 = _var(la, 0, K), _var(lb, 1, K)
    n = len(s)
    Z = obj([arb(0)] * n)
    R = [[Jet.const(Z, K), Jet.const(Z, K)], [Jet.const(Z, K), Jet.const(Z, K)]]
    for r in range(P.d):
        if r in (a, b):
            continue
        nv = [Ai[0][i] * P.nu[r][0] + Ai[1][i] * P.nu[r][1] for i in range(2)]
        l0 = x1 * P.nu[r][0] + x2 * P.nu[r][1] + lam[r]
        lr = Jet.const(l0, K)
        ms, ix, _ = monos(K)
        lr.c[ix[(1, 0)]] = obj([arb(nv[0])] * n); lr.c[ix[(0, 1)]] = obj([arb(nv[1])] * n)
        il = lr.recip()
        for i in range(2):
            for j in range(2):
                if nv[i] * nv[j]:
                    R[i][j] = R[i][j] + il * (arb(nv[i] * nv[j]) / 2)
    if prob is not None and getattr(prob, 'Cint', None) is not None:
        D = mc.fderivs(prob, x1, x2, K + 2, False)
        # dx = Ai dy (linear): monomial coefficients of dx_i
        ms, ix, _ = monos(K)
        dx = []
        for i in range(2):
            c = [None] * len(ms)
            c[ix[(1, 0)]] = obj([arb(Ai[i][0])] * n); c[ix[(0, 1)]] = obj([arb(Ai[i][1])] * n)
            dx.append(Jet(c, K))
        one = Jet.const(obj([arb(1)] * n), K)
        p1 = [one]; p2 = [one]
        for _ in range(K):
            p1.append(p1[-1] * dx[0]); p2.append(p2[-1] * dx[1])
        Hx = [[None, None], [None, None]]
        for i in range(2):
            for j in range(i, 2):
                al = [0, 0]; al[i] += 1; al[j] += 1
                acc = Jet.const(Z, K)
                for a1 in range(K + 1):
                    for a2 in range(K + 1 - a1):
                        acc = acc + (p1[a1] * p2[a2]) * (D[(al[0] + a1, al[1] + a2)] * (arb(1) / (factorial(a1) * factorial(a2))))
                Hx[i][j] = acc
        Hx[1][0] = Hx[0][1]
        for aa in range(2):
            for bb in range(2):
                for i in range(2):
                    for j in range(2):
                        co = Ai[i][aa] * Ai[j][bb]
                        if co:
                            R[aa][bb] = R[aa][bb] + Hx[i][j] * arb(co)
    D1, D2 = Y1 * 2, Y2 * 2
    m11 = D1 * R[0][0] + 1; m12 = D1 * R[0][1]; m21 = D2 * R[1][0]; m22 = D2 * R[1][1] + 1
    idet = (m11 * m22 - m12 * m21).recip()
    U = [[m22 * D1 * idet, -(m12 * D2 * idet)], [-(m21 * D1 * idet), m11 * D2 * idet]]
    Sc = -(U[0][0].d(0).d(0) + U[0][1].d(0).d(1) + U[1][0].d(1).d(0) + U[1][1].d(1).d(1))
    Sg = [prob.Sabc[1] * Ai[0][i] + prob.Sabc[2] * Ai[1][i] for i in range(2)]
    S0 = x1 * prob.Sabc[1] + x2 * prob.Sabc[2] + prob.Sabc[0]
    Ss = Jet.const(S0, Sc.K)
    ms2, ix2, _ = monos(Sc.K)
    Ss.c[ix2[(1, 0)]] = obj([Sg[0]] * n); Ss.c[ix2[(0, 1)]] = obj([Sg[1]] * n)
    psi = Sc - Ss                                    # order 2
    lap = None
    for i in range(2):
        for j in range(2):
            dj = psi.d(j)
            tt = U[i][j].c[0] * dj.d(i).c[0] + U[i][j].d(i).c[0] * dj.c[0]
            lap = tt if lap is None else lap + tt
    return psi.c[0], -lap


class CtrlProb:
    """control polygons (f = 0) in the same interface: P (clwr.Poly), lam (arb), Sabc (arb), Cint None."""

    def __init__(self, kind, prec=320):
        import clwr, mpmath as mp
        ctx.prec = prec
        mp.mp.prec = prec + 100
        if kind == 'square':
            self.P = clwr.Poly([(1, 0), (0, 1), (-1, 0), (0, -1)], [1, 1, 1, 1], 'square')
        else:
            self.P = clwr.Poly([(1, 0), (0, 1), (-1, -1)], [1, 1, 1], 'triangle')
        self.lam = [mc.A(v) for v in self.P.lam]
        self.Sabc = [arb(4), arb(0), arb(0)]
        self.c = arb(1)
        self.Cint = None
        self.prec = prec
