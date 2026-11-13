"""met1_core.py -- seat R, BAR-MET1.  Ball-arithmetic (python-flint arb / acb) evaluation of the Einstein residual
E = Ric(g) - Lambda g and its covariant derivatives nabla^j E (j <= 3) of g = c S^-2 g_K, vectorised over points.

Coordinates per vertex element k: (p, q, theta_a, theta_b), p = s sqrt(A(t)), q = t sqrt(B(s)) (seat W / NOTES-R sec. 0).
Order-K bivariate jets in (p, q) (K = 5 for j <= 3).  All arithmetic is arb (real points) or acb (complex boxes); inputs may be
balls, so outputs enclose the exact values for every point of the input box.  Components with an odd number of theta indices are
set to exact zero (the reflection theta -> -theta is an isometry: no dp dtheta terms in g).
"""
import numpy as np
import mpmath as mp
from math import factorial
from flint import arb, acb, arb_mat, acb_mat, ctx

LAM = 3


def obj(a):
    return np.array(a, dtype=object)


# ------------------------------------------------------------------ bivariate jets
_MONO = {}


def monos(K):
    if K not in _MONO:
        ms = [(i, n - i) for n in range(K + 1) for i in range(n, -1, -1)]
        ix = {m: k for k, m in enumerate(ms)}
        tab = []
        for a, (i1, j1) in enumerate(ms):
            for b, (i2, j2) in enumerate(ms):
                if i1 + j1 + i2 + j2 <= K:
                    tab.append((a, b, ix[(i1 + i2, j1 + j2)]))
        _MONO[K] = (ms, ix, tab)
    return _MONO[K]


class Jet:
    """scalar jet: c[k] coefficient of monomial ms[k] (dp^i dq^j), each an object array over points (or None = 0)."""
    __slots__ = ('c', 'K')

    def __init__(self, c, K):
        self.c, self.K = c, K

    @staticmethod
    def const(v, K):
        ms, _, _ = monos(K)
        c = [None] * len(ms); c[0] = v
        return Jet(c, K)

    def trunc(self, K):
        if K == self.K:
            return self
        n = len(monos(K)[0])
        return Jet(self.c[:n], K)

    def __add__(self, o):
        if o is None:
            return self
        if not isinstance(o, Jet):
            c = list(self.c); c[0] = o if c[0] is None else c[0] + o
            return Jet(c, self.K)
        K = min(self.K, o.K)
        a, b = self.trunc(K).c, o.trunc(K).c
        return Jet([x if y is None else (y if x is None else x + y) for x, y in zip(a, b)], K)

    __radd__ = __add__

    def __neg__(self):
        return Jet([None if x is None else -x for x in self.c], self.K)

    def __sub__(self, o):
        return self + (-o if isinstance(o, Jet) else -o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet([None if x is None else x * o for x in self.c], self.K)
        K = min(self.K, o.K)
        _, _, tab = monos(K)
        out = [None] * len(monos(K)[0])
        a, b = self.c, o.c
        for i, j, k in tab:
            if a[i] is None or b[j] is None:
                continue
            v = a[i] * b[j]
            out[k] = v if out[k] is None else out[k] + v
        return Jet(out, K)

    __rmul__ = __mul__

    def recip(self):
        # exact recursion: y_0 = 1/a_0, y_m = -(1/a_0) sum_{m' != 0} a_{m'} y_{m - m'} (monomials in graded order)
        ms, ix, _ = monos(self.K)
        a = self.c
        y = [None] * len(ms)
        inv0 = 1 / a[0]
        y[0] = inv0
        for n, (i, j) in enumerate(ms):
            if n == 0:
                continue
            acc = None
            for i1 in range(i + 1):
                for j1 in range(j + 1):
                    if i1 + j1 == 0:
                        continue
                    am = a[ix[(i1, j1)]]
                    yy = y[ix[(i - i1, j - j1)]]
                    if am is None or yy is None:
                        continue
                    v = am * yy
                    acc = v if acc is None else acc + v
            y[n] = None if acc is None else -(acc * inv0)
        return Jet(y, self.K)

    def d(self, k):
        ms, ix, _ = monos(self.K)
        K = self.K - 1
        ms2, ix2, _ = monos(K)
        out = [None] * len(ms2)
        for n, (i, j) in enumerate(ms2):
            src = ix[(i + 1, j)] if k == 0 else ix[(i, j + 1)]
            v = self.c[src]
            if v is not None:
                out[n] = v * ((i + 1) if k == 0 else (j + 1))
        return Jet(out, K)

    def val(self):
        return self.c[0]


def jadd(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return a + b


def jmul(a, b):
    if a is None or b is None:
        return None
    return a * b


def nth(idx):
    return sum(1 for i in idx if i >= 2)


# ------------------------------------------------------------------ problem data (exact, arb)
class Problem:
    """polygon (mp exact data converted to arb), S, c, coefficients C (exact rationals -> arb balls)."""

    def __init__(self, P, Sabc_arb, c_arb, Cint, B, prec):
        self.P = P
        self.prec = prec
        ctx.prec = prec
        self.Sabc = Sabc_arb
        self.c = c_arb
        self.lam = [A(v) for v in P.lam]
        self.M = Cint.shape[0] - 1
        self.Cint, self.B = Cint, B

    def Cmat(self, cplx):
        M = self.M
        sc = arb(2) ** (-self.B)
        if cplx:
            return acb_mat(M + 1, M + 1, [acb(int(v)) * sc for v in self.Cint.ravel()])
        return arb_mat(M + 1, M + 1, [int(v) for v in self.Cint.ravel()]) * sc


def A(v):
    return arb(mp.nstr(v, int(ctx.prec * 0.3011) + 20, strip_zeros=False))


def cheb_tables(x, M, kmax):
    n = len(x)
    one = x[0] * 0 + 1
    zero = x[0] * 0
    T = [[None] * (M + 1) for _ in range(kmax + 1)]
    for k in range(kmax + 1):
        T[k][0] = obj([one if k == 0 else zero] * n)
        T[k][1] = x.copy() if k == 0 else obj([one if k == 1 else zero] * n)
    for i in range(1, M):
        for k in range(kmax + 1):
            v = 2 * x * T[k][i] - T[k][i - 1]
            if k >= 1:
                v = v + 2 * k * T[k - 1][i]
            T[k][i + 1] = v
    return T


def cheb_tables_closed(x, M, kmax):
    """T_i^{(k)}(x) by closed forms (no three-term recurrence, whose ball radii grow like 2.4^i):
    T_n = chebyshev_t(n), T_n^{(k)} = n 2^(k-1) (k-1)! C^{(k)}_{n-k} (Gegenbauer), k >= 1."""
    n = len(x)
    zero = x[0] * 0
    T = [[None] * (M + 1) for _ in range(kmax + 1)]
    for i in range(M + 1):
        T[0][i] = obj([v.chebyshev_t(i) for v in x])
        for k in range(1, kmax + 1):
            if i < k:
                T[k][i] = obj([zero] * n)
            else:
                f = i * 2 ** (k - 1) * factorial(k - 1)
                T[k][i] = obj([v.gegenbauer_c(i - k, k) * f for v in x])
    return T


CLOSED = {'on': False}


def fderivs(prob, x1, x2, kmax, cplx):
    M = prob.M
    n = len(x1)
    tab = cheb_tables_closed if CLOSED['on'] else cheb_tables
    T1 = tab(x1, M, kmax)
    T2 = tab(x2, M, kmax)
    Cm = prob.Cmat(cplx)
    MAT = acb_mat if cplx else arb_mat
    D = {}
    for a1 in range(kmax + 1):
        A1 = MAT(n, M + 1, [T1[a1][i][p] for p in range(n) for i in range(M + 1)])
        Pl = np.array((A1 * Cm).tolist(), dtype=object)
        for a2 in range(kmax + 1 - a1):
            if a1 + a2 < 2:
                continue
            T2m = np.array([T2[a2][j] for j in range(M + 1)], dtype=object).T
            D[(a1, a2)] = (Pl * T2m).sum(axis=1)
    return D


def metric_jet(prob, k, p0, q0, K, cplx):
    """g[mu][nu] as Jets of order K around (p0, q0) (object arrays), coordinates (p, q, theta_a, theta_b)."""
    P = prob.P
    E = P.elem(k)
    a, b = E['a'], E['b']
    n = len(p0)
    Am = [list(P.nu[a]), list(P.nu[b])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    lab = [prob.lam[a], prob.lam[b]]
    ms, ix, _ = monos(K)
    one = p0[0] * 0 + 1
    onearr = obj([one] * n)

    def var(v0, which):
        c = [None] * len(ms); c[0] = v0; c[ix[(1, 0)] if which == 0 else ix[(0, 1)]] = onearr
        return Jet(c, K)
    pj, qj = var(p0, 0), var(q0, 1)
    p2, q2 = pj * pj, qj * qj
    X = [(p2 - lab[0]) * Ai[i][0] + (q2 - lab[1]) * Ai[i][1] for i in range(2)]
    x0 = [X[0].val(), X[1].val()]
    R = [[None, None], [None, None]]
    for r in range(P.d):
        if r in (a, b):
            continue
        nv = [Ai[0][i] * P.nu[r][0] + Ai[1][i] * P.nu[r][1] for i in range(2)]
        lr = X[0] * P.nu[r][0] + X[1] * P.nu[r][1] + prob.lam[r]
        il = lr.recip()
        for i in range(2):
            for j in range(2):
                if nv[i] * nv[j] != 0:
                    R[i][j] = jadd(R[i][j], il * (arb(nv[i] * nv[j]) / 2))
    D = fderivs(prob, x0[0], x0[1], K + 2, cplx)
    dx = [X[0] - x0[0], X[1] - x0[1]]
    for t in dx:
        t.c[0] = None
    pw = {}
    for a1 in range(K + 1):
        for a2 in range(K + 1 - a1):
            t = Jet.const(onearr, K)
            for _ in range(a1):
                t = t * dx[0]
            for _ in range(a2):
                t = t * dx[1]
            pw[(a1, a2)] = t * (arb(1) / (factorial(a1) * factorial(a2)))

    def H(i, j):
        al = [0, 0]; al[i] += 1; al[j] += 1
        acc = None
        for (a1, a2), t in pw.items():
            acc = jadd(acc, t * D[(al[0] + a1, al[1] + a2)])
        return acc
    Hf = [[H(0, 0), H(0, 1)], [None, H(1, 1)]]
    Hf[1][0] = Hf[0][1]
    for aa in range(2):
        for bb in range(2):
            for i in range(2):
                for j in range(2):
                    co = Ai[i][aa] * Ai[j][bb]
                    if co != 0:
                        R[aa][bb] = jadd(R[aa][bb], Hf[i][j] * co)
    gpp = 2 + 4 * p2 * R[0][0]
    gpq = 4 * pj * qj * R[0][1]
    gqq = 2 + 4 * q2 * R[1][1]
    m11 = 1 + 2 * p2 * R[0][0]; m12 = 2 * p2 * R[0][1]; m21 = 2 * q2 * R[0][1]; m22 = 1 + 2 * q2 * R[1][1]
    idet = (m11 * m22 - m12 * m21).recip()
    Gaa = m22 * 2 * p2 * idet
    Gab = -(m12 * 2 * q2) * idet
    Gbb = m11 * 2 * q2 * idet
    S = X[0] * prob.Sabc[1] + X[1] * prob.Sabc[2] + prob.Sabc[0]
    Si = S.recip()
    cf = Si * Si * prob.c
    g = [[cf * gpp, cf * gpq, None, None], [cf * gpq, cf * gqq, None, None],
         [None, None, cf * Gaa, cf * Gab], [None, None, cf * Gab, cf * Gbb]]
    return g, x0


def inv_blocks(g):
    """g^-1 as jets (block diagonal 2x2 + 2x2)."""
    gi = [[None] * 4 for _ in range(4)]
    for (u, v) in ((0, 1), (2, 3)):
        a, b, d = g[u][u], g[u][v], g[v][v]
        idt = (a * d - b * b).recip()
        gi[u][u] = d * idt; gi[v][v] = a * idt; gi[u][v] = gi[v][u] = -(b * idt)
    return gi


def residual_tensors(g, jmax=3):
    """returns list [E, nabla E, ..., nabla^jmax E] as dicts index-tuple -> value arrays (order-0), and g0, gi0 values."""
    K = g[0][0].K
    gi = inv_blocks(g)
    dg = {}
    for m in (0, 1):
        for a in range(4):
            for b in range(4):
                if g[a][b] is not None:
                    dg[(m, a, b)] = g[a][b].d(m)

    def DG(m, a, b):
        return dg.get((m, a, b)) if m < 2 else None
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
    G = lambda l, m, nn: Gam.get((l, m, nn))
    Kr = K - 2
    Ric = {}
    for s in range(4):
        for nn in range(s, 4):
            if nth((s, nn)) % 2:
                continue
            acc = None
            for r in (0, 1):
                if G(r, nn, s) is not None:
                    acc = jadd(acc, G(r, nn, s).d(r))
            if nn < 2:
                for r in range(4):
                    if G(r, r, s) is not None:
                        acc = jadd(acc, -G(r, r, s).d(nn))
            for r in range(4):
                for l in range(4):
                    acc = jadd(acc, jmul(G(r, r, l), G(l, nn, s)))
                    t = jmul(G(r, nn, l), G(l, r, s))
                    if t is not None:
                        acc = jadd(acc, -t)
            Ric[(s, nn)] = acc.trunc(Kr)
            Ric[(nn, s)] = Ric[(s, nn)]
    E = {}
    for (s, nn), v in Ric.items():
        gv = g[s][nn]
        E[(s, nn)] = v - (gv.trunc(Kr) * LAM if gv is not None else 0)
    out = [E]
    T = E
    rank = 2
    for j in range(jmax):
        Kt = Kr - j - 1
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
        out.append(T)
    gi0 = [[None if gi[a][b] is None else gi[a][b].val() for b in range(4)] for a in range(4)]
    return out, gi0


def sqnorm(T, gi0, rank):
    """|T|^2_g = g^{a1 b1} ... T_a T_b at order 0 (values)."""
    V = {k: v.val() for k, v in T.items() if v.val() is not None}
    cur = dict(V)
    for slot in range(rank):
        nxt = {}
        for idx, v in cur.items():
            a = idx[slot]
            for b in range(4):
                gg = gi0[a][b]
                if gg is None:
                    continue
                nidx = idx[:slot] + (b,) + idx[slot + 1:]
                nxt[nidx] = gg * v if nidx not in nxt else nxt[nidx] + gg * v
        cur = nxt
    tot = None
    for idx, v in cur.items():
        if idx in V:
            t = v * V[idx]
            tot = t if tot is None else tot + t
    return tot


def element_pq(prob, k, s, t):
    """(p, q) = (s sqrt(A(t)), t sqrt(B(s))) for object arrays s, t (arb or acb)."""
    P = prob.P
    E = P.elem(k)
    a, b = E['a'], E['b']
    al1, al2 = A(P.ell(a, *E['Mb'])), A(P.ell(a, *P.O))
    be1, be2 = A(P.ell(b, *E['Ma'])), A(P.ell(b, *P.O))
    Aa = (1 - t * t) * al1 + t * t * al2
    Bb = (1 - s * s) * be1 + s * s * be2
    p = s * obj([v.sqrt() for v in Aa])
    q = t * obj([v.sqrt() for v in Bb])
    return p, q, Aa, Bb


def norms_at(prob, k, s, t, cplx=False, jmax=3):
    p, q, _, _ = element_pq(prob, k, s, t)
    g, x0 = metric_jet(prob, k, p, q, jmax + 2, cplx)
    Ts, gi0 = residual_tensors(g, jmax)
    return [sqnorm(T, gi0, j + 2) for j, T in enumerate(Ts)], x0


def weight(prob, k, s, t, g, p, q, Sval):
    """volume density per ds dt: dV = c^2 S^-4 dx dtheta, dx = d(l_a) d(l_b) = 4 p q dp dq, and
    p q d(p,q)/d(s,t) = s t (A B - s^2 t^2 dA dB): w = (2 pi)^2 4 c^2 S^-4 s t (A B - s^2 t^2 dA dB)  (polynomial x S^-4)."""
    P = prob.P
    E = P.elem(k); a, b = E['a'], E['b']
    al1, al2 = A(P.ell(a, *E['Mb'])), A(P.ell(a, *P.O))
    be1, be2 = A(P.ell(b, *E['Ma'])), A(P.ell(b, *P.O))
    Aa = (1 - t * t) * al1 + t * t * al2
    Bb = (1 - s * s) * be1 + s * s * be2
    pqJ = s * t * (Aa * Bb - s * s * t * t * (al2 - al1) * (be2 - be1))
    Si = 1 / Sval
    return (2 * arb.pi()) ** 2 * 4 * prob.c ** 2 * Si ** 4 * pqJ


def G_at(prob, k, s, t, cplx=False, jmax=3):
    """G_j = |nabla^j E|^2_g * w  (j = 0..jmax) and the sqrt argument h (must avoid the branch cut)."""
    p, q, _, _ = element_pq(prob, k, s, t)
    g, x0 = metric_jet(prob, k, p, q, jmax + 2, cplx)
    Ts, gi0 = residual_tensors(g, jmax)
    F = [sqnorm(T, gi0, j + 2) for j, T in enumerate(Ts)]
    Sval = x0[0] * prob.Sabc[1] + x0[1] * prob.Sabc[2] + prob.Sabc[0]
    w = weight(prob, k, s, t, g, p, q, Sval)
    return [f * w for f in F], F, w
