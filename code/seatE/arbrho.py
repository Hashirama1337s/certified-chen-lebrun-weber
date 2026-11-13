"""arbrho.py -- seat E: arb (high-precision) covariant norms |nabla^j phi|^2_{g_K}, j = 0..3, of T^2-invariant functions phi on the
toric Kahler manifold, in the vertex y-chart (y1, y2, th1, th2), y = (l_a, l_b):
  g_K = H dy dy + U dth dth,   H = Hess_y u = diag(1/(2 y1), 1/(2 y2)) + R,   U = H^-1 = m^-1 Dg  (m = I + Dg R, Dg = diag(2y1, 2y2)).
phi = psi = Scal(u) - S (Abreu, Scal = -d_i d_j U^ij) or phi = Scal (validation: |nabla^j Scal| = |nabla^{j-1} ds|).
The chart is singular at the rods (H ~ 1/y), so it is used only at INTERIOR points (first-kind Chebyshev nodes); the high
precision absorbs the 1/y cancellations.  Christoffel symbols and the covariant derivative loop follow seat R's met1_core
(generic coordinate formulas, theta derivatives zero, entries with an odd number of theta indices exactly 0 by the joint
reflection isometry).  Jets: met1_core.Jet (arb balls), order K = 5 (U to order 5 -> phi to order 3)."""
import numpy as np
from math import factorial
from flint import arb, ctx
import common
import met1_core as mc
from met1_core import Jet, obj, monos, jadd, jmul, nth, inv_blocks
import arbpsi


def _lin(v0, c1, c2, K, n):
    ms, ix, _ = monos(K)
    J = Jet.const(v0, K)
    J.c = list(J.c)
    J.c[ix[(1, 0)]] = obj([arb(c1)] * n); J.c[ix[(0, 1)]] = obj([arb(c2)] * n)
    return J


def build(prob, k, sa, ta, K=5):
    """metric jets g (4x4 list), psi jet, Scal jet at the arb points (sa, ta) of element k."""
    P = prob.P
    E = P.elem(k)
    a, b = E['a'], E['b']
    Am = [list(P.nu[a]), list(P.nu[b])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    n = len(sa)
    al1, al2 = common.to_arb(P.ell(a, *E['Mb'])), common.to_arb(P.ell(a, *P.O))
    be1, be2 = common.to_arb(P.ell(b, *E['Ma'])), common.to_arb(P.ell(b, *P.O))
    la = sa * sa * ((1 - ta * ta) * al1 + ta * ta * al2)
    lb = ta * ta * ((1 - sa * sa) * be1 + sa * sa * be2)
    lam = prob.lam
    x1 = (la - lam[a]) * Ai[0][0] + (lb - lam[b]) * Ai[0][1]
    x2 = (la - lam[a]) * Ai[1][0] + (lb - lam[b]) * Ai[1][1]
    Y1, Y2 = _lin(la, 1, 0, K, n), _lin(lb, 0, 1, K, n)
    Z = obj([arb(0)] * n)
    R = [[Jet.const(Z, K) for _ in range(2)] for _ in range(2)]
    for r in range(P.d):
        if r in (a, b):
            continue
        nv = [Ai[0][i] * P.nu[r][0] + Ai[1][i] * P.nu[r][1] for i in range(2)]
        l0 = x1 * P.nu[r][0] + x2 * P.nu[r][1] + lam[r]
        il = _lin(l0, nv[0], nv[1], K, n).recip()
        for i in range(2):
            for j in range(2):
                if nv[i] * nv[j]:
                    R[i][j] = R[i][j] + il * (arb(nv[i] * nv[j]) / 2)
    if getattr(prob, 'Cint', None) is not None:
        D = mc.fderivs(prob, x1, x2, K + 2, False)
        dx = [_lin(Z, Ai[i][0], Ai[i][1], K, n) for i in range(2)]
        for t in dx:
            t.c[0] = None
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
    H = [[Y1.recip() * arb(0.5) + R[0][0], R[0][1]], [R[1][0], Y2.recip() * arb(0.5) + R[1][1]]]
    Sc = -(U[0][0].d(0).d(0) + U[0][1].d(0).d(1) + U[1][0].d(1).d(0) + U[1][1].d(1).d(1))
    Sg = [prob.Sabc[1] * Ai[0][i] + prob.Sabc[2] * Ai[1][i] for i in range(2)]
    S0 = x1 * prob.Sabc[1] + x2 * prob.Sabc[2] + prob.Sabc[0]
    ms2, ix2, _ = monos(Sc.K)
    Ss = Jet.const(S0, Sc.K); Ss.c = list(Ss.c)
    Ss.c[ix2[(1, 0)]] = obj([Sg[0]] * n); Ss.c[ix2[(0, 1)]] = obj([Sg[1]] * n)
    g = [[H[0][0], H[0][1], None, None], [H[1][0], H[1][1], None, None],
         [None, None, U[0][0], U[0][1]], [None, None, U[1][0], U[1][1]]]
    return g, Sc - Ss, Sc


def christoffel(g):
    K = g[0][0].K
    gi = inv_blocks(g)
    dg = {}
    for m in (0, 1):
        for a in range(4):
            for b in range(4):
                if g[a][b] is not None:
                    dg[(m, a, b)] = g[a][b].d(m)
    DG = lambda m, a, b: dg.get((m, a, b)) if m < 2 else None
    Gam = {}
    for l in range(4):
        for m in range(4):
            for nn in range(m, 4):
                if nth((l, m, nn)) % 2:
                    continue
                acc = None
                for kk in range(4):
                    if gi[l][kk] is None:
                        continue
                    t = jadd(jadd(DG(m, kk, nn), DG(nn, kk, m)), None if DG(kk, m, nn) is None else -DG(kk, m, nn))
                    if t is not None:
                        acc = jadd(acc, gi[l][kk].trunc(K - 1) * t)
                if acc is not None:
                    acc = acc * arb(0.5)
                    Gam[(l, m, nn)] = acc; Gam[(l, nn, m)] = acc
    return Gam, gi


def nabla_norms(phi, Gam, gi, jmax=3):
    """|nabla^j phi|^2 (j = 0..jmax) at the points (order-0 arb arrays)."""
    G = lambda l, m, nn: Gam.get((l, m, nn))
    gi0 = [[None if gi[a][b] is None else gi[a][b].val() for b in range(4)] for a in range(4)]
    T = {(): phi}
    out = [mc.sqnorm(T, gi0, 0)]
    rank = 0
    Kt = phi.K
    for j in range(jmax):
        Kt -= 1
        T2 = {}
        idxs = [()]
        for _ in range(rank + 1):
            idxs = [i + (c,) for i in idxs for c in range(4)]
        for idx in idxs:
            if nth(idx) % 2:
                continue
            m, rest = idx[0], idx[1:]
            acc = None
            if m < 2 and rest in T:
                acc = T[rest].d(m)
            for i in range(rank):
                for e in range(4):
                    gg = G(e, m, rest[i])
                    if gg is None:
                        continue
                    src = rest[:i] + (e,) + rest[i + 1:]
                    if src not in T:
                        continue
                    acc = jadd(acc, -(gg.trunc(Kt) * T[src].trunc(Kt)))
            if acc is not None:
                T2[idx] = acc.trunc(Kt)
        T = T2
        rank += 1
        out.append(mc.sqnorm(T, gi0, rank))
    return out


def rho_norms(prob, k, sa, ta, which='psi'):
    g, psi, Sc = build(prob, k, obj(list(sa)), obj(list(ta)))
    Gam, gi = christoffel(g)
    phi = psi if which == 'psi' else Sc
    return nabla_norms(phi, Gam, gi, 3)
