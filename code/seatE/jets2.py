"""jets2.py -- seat E: fast certified curvature pipeline on ball tensor jets (box-first layout, numba kernels, parity pruning).

KJ: arrays m (nb, Nmono, *tshape) (float64 at real points/boxes, complex128 on complex boxes) and r (radii, same shape).
Products: C[o] = sum_t s_t A[ia_t] B[ib_t] over statically generated triples (monomial pairs x index assignments of one or
several signed einsum specs), one fused kernel per product computing mid, big = sum (|a| + ra)(|b| + rb), small = sum |a||b|,
and radius = (big - small) + (3 K_o + 20) u big + tiny  (K_o = number of terms of entry o; s_t in {+-1, +-1/2}, exact).
(big - small) is the exact disc-product bound summed; the second term dominates binary64 rounding of the recursive sums, of
the complex/real products and of the final difference (each relative error <= K_o u, complex product <= sqrt5 u).
Parity pruning: every object of the curvature pipeline (g, g^-1, dg, Gamma, Riem, nabla^q Riem, Ric, E, Scal, psi) is built
naturally from the metric and so is even under the chart's JOINT reflection isometry J (polar: (ta, tb) -> (-ta, -tb);
cart: (Y, W) -> (-Y, -W); rodA: (Y, tb) -> (-Y, -tb); rodB: (W, ta) -> (-W, -ta)).  [A single-angle flip is NOT an isometry:
the g_{ta tb} term changes sign.]  J is LINEAR in the chart coordinates and fixes the evaluation points (Y = W = 0 exactly).
Hence a jet entry (monomial e, indices I) with odd  #{i : I_i in flipped coordinates} + sum_{v flipped} e[v]  is exactly
zero; such entries are never formed (exact 0, radius 0).  Metric construction (which uses J-odd scalars such as Y) runs
unpruned (sym = NONE) and the pruned entries of g are set to exact 0 (NOTES-E sec. 3).
Chart formulas: see jets.py docstring and NOTES-E sec. 2 (identical formulas; jets.py is the dense reference engine).
"""
import itertools
from math import factorial
import numpy as np
from numba import njit, prange
import common
from common import Disc, SAFE, TINY
import jets
from jets import mtab, CHARTS

U2 = 2.3e-16
LET = 'abcdefghijklmnopqrstuvwxy'
SYM = {'cart': ((1, 3), (1, 3)), 'polar': ((2, 3), ()), 'xchart': ((2, 3), ()), 'rodA': ((1, 3), (1,)), 'rodB': ((2, 3), (2,))}   # (flipped indices, vars)
NONE = ((), ())


@njit(cache=True)
def _sub_probe(x):
    return x / 2.0


# amend2: the rounding model needs round-to-nearest and gradual underflow (no flush-to-zero), in numpy AND in numba code
assert 1.0 + 2.0 ** -53 == 1.0 and 1.0 + 2.0 ** -52 != 1.0 and (1.0 + 3 * 2.0 ** -53) == 1.0 + 2.0 ** -51, 'not round-to-nearest'
assert np.float64(2.0 ** -1022) / 2.0 != 0.0, 'flush-to-zero in numpy'
assert _sub_probe(2.0 ** -1022) != 0.0, 'flush-to-zero in numba'


@njit(parallel=True, cache=True)
def _kr(Am, Ar, Bm, Br, ia, ib, io, sg, nout, ef):
    nb = Am.shape[0]
    Om = np.zeros((nb, nout)); R = np.zeros((nb, nout))
    for z in prange(nb):
        Ob = np.zeros(nout); Os = np.zeros(nout)
        for t in range(ia.shape[0]):
            a = Am[z, ia[t]]; b = Bm[z, ib[t]]
            aa = abs(a); ba = abs(b)
            o = io[t]; s = sg[t]
            Om[z, o] += s * (a * b)
            w = abs(s)
            Ob[o] += w * ((aa + Ar[z, ia[t]]) * (ba + Br[z, ib[t]]))
            Os[o] += w * (aa * ba)
        for o in range(nout):
            if ef[o] > 0:
                R[z, o] = (Ob[o] - Os[o]) + ef[o] * Ob[o] + 1e-290
    return Om, R


@njit(parallel=True, cache=True)
def _kc(Am, Ar, Bm, Br, ia, ib, io, sg, nout, ef):
    nb = Am.shape[0]
    Om = np.zeros((nb, nout), np.complex128); R = np.zeros((nb, nout))
    for z in prange(nb):
        Ob = np.zeros(nout); Os = np.zeros(nout)
        for t in range(ia.shape[0]):
            a = Am[z, ia[t]]; b = Bm[z, ib[t]]
            aa = abs(a); ba = abs(b)
            o = io[t]; s = sg[t]
            Om[z, o] += s * (a * b)
            w = abs(s)
            Ob[o] += w * ((aa + Ar[z, ia[t]]) * (ba + Br[z, ib[t]]))
            Os[o] += w * (aa * ba)
        for o in range(nout):
            if ef[o] > 0:
                R[z, o] = (Ob[o] - Os[o]) + ef[o] * Ob[o] + 1e-290
    return Om, R


def _parity_mono(nv, N, sym):
    mons = mtab(nv, N)['mons']
    return np.array([sum(e[v] for v in sym[1]) for e in mons], int)


_TR = {}


def triples(nv, N, Na, Nb, specs, sym, block=None):
    """specs: tuple of (spec, sign); all specs share the output letters.  block = (d, d1min) or None."""
    key = (nv, N, Na, Nb, specs, sym, block)
    if key in _TR:
        return _TR[key]
    Ta, Tb, To = mtab(nv, Na), mtab(nv, Nb), mtab(nv, N)
    pma, pmb = _parity_mono(nv, Na, sym), _parity_mono(nv, Nb, sym)
    pairs = []
    for i, e1 in enumerate(Ta['mons']):
        d1 = sum(e1)
        if d1 > N:
            break
        for j, e2 in enumerate(Tb['mons']):
            d2 = sum(e2)
            if d1 + d2 > N:
                break
            if block is not None and (d1 + d2 != block[0] or d1 < block[1]):
                continue
            pairs.append((i, j, To['ix'][tuple(x + y for x, y in zip(e1, e2))]))
    pairs = np.array(pairs, int).reshape(-1, 3)
    fl = np.array(sym[0], int)
    IA, IB, IO, SG = [], [], [], []
    for spec, sign in specs:
        sa, sbs = spec.split('->')[0].split(',')
        so = spec.split('->')[1]
        U = sorted(set(sa + sbs))
        vals = np.array(list(itertools.product(range(4), repeat=len(U)))).reshape(-1, len(U)) if U else np.zeros((1, 0), int)
        col = {c: i for i, c in enumerate(U)}

        def flat(s):
            if not s:
                return np.zeros(len(vals), int)
            return sum(vals[:, col[c]] * 4 ** (len(s) - 1 - i) for i, c in enumerate(s))
        offa, offb, offo = flat(sa), flat(sbs), flat(so)
        pa = np.isin(vals[:, [col[c] for c in sa]], fl).sum(1) if sa else np.zeros(len(vals), int)
        pb = np.isin(vals[:, [col[c] for c in sbs]], fl).sum(1) if sbs else np.zeros(len(vals), int)
        okA = (pma[pairs[:, 0]][:, None] + pa[None]) % 2 == 0
        okB = (pmb[pairs[:, 1]][:, None] + pb[None]) % 2 == 0
        P_, S_ = np.nonzero(okA & okB)
        sza, szb, szo = 4 ** len(sa), 4 ** len(sbs), 4 ** len(so)
        IA.append(pairs[P_, 0] * sza + offa[S_]); IB.append(pairs[P_, 1] * szb + offb[S_])
        if block is None:
            IO.append(pairs[P_, 2] * szo + offo[S_])
        else:
            IO.append((pairs[P_, 2] - To['starts'][block[0]]) * szo + offo[S_])
        SG.append(np.full(len(P_), float(sign)))
    so = specs[0][0].split('->')[1]
    szo = 4 ** len(so)
    nout = (To['n'] if block is None else To['starts'][block[0] + 1] - To['starts'][block[0]]) * szo
    ia, ib, io, sg = (np.concatenate(v) for v in (IA, IB, IO, SG))
    order = np.argsort(io, kind='stable')                     # output-ordered triples (cache locality)
    ia, ib, io, sg = ia[order].astype(np.int64), ib[order].astype(np.int64), io[order].astype(np.int64), sg[order]
    K = np.bincount(io, minlength=nout)
    ef = np.where(K > 0, (3 * K + 20) * U2, 0.0)
    _TR[key] = (ia, ib, io, sg, nout, ef, len(so))
    return _TR[key]


def _run(a, b, tr):
    ia, ib, io, sg, nout, ef, _ = tr
    nb = a.nb
    Am, Ar = a.m.reshape(nb, -1), a.r.reshape(nb, -1)
    Bm, Br = b.m.reshape(nb, -1), b.r.reshape(nb, -1)
    if np.iscomplexobj(Am) or np.iscomplexobj(Bm):
        Om, R = _kc(np.ascontiguousarray(Am, np.complex128), np.ascontiguousarray(Ar), np.ascontiguousarray(Bm, np.complex128),
                    np.ascontiguousarray(Br), ia, ib, io, sg, nout, ef)
    else:
        Om, R = _kr(np.ascontiguousarray(Am), np.ascontiguousarray(Ar), np.ascontiguousarray(Bm), np.ascontiguousarray(Br),
                    ia, ib, io, sg, nout, ef)
    bad = ~np.isfinite(Om) | ~np.isfinite(R)
    if bad.any():
        R = np.where(bad, np.inf, R); Om = np.where(bad, 0, Om)
    return Om, R


class KJ:
    __slots__ = ('m', 'r', 'nv', 'N', 'sym')

    def __init__(self, m, r, nv, N, sym):
        self.m, self.r, self.nv, self.N, self.sym = m, r, nv, N, sym

    @property
    def rank(self):
        return self.m.ndim - 2

    @property
    def nb(self):
        return self.m.shape[0]

    @staticmethod
    def const(d, nv, N, sym=NONE):
        n = mtab(nv, N)['n']
        m = np.zeros((len(d.m), n), d.m.dtype); r = np.zeros((len(d.m), n))
        m[:, 0] = d.m; r[:, 0] = d.r
        return KJ(m, r, nv, N, sym)

    @staticmethod
    def var(d, k, nv, N):
        J = KJ.const(d, nv, N)
        e = [0] * nv; e[k] = 1
        J.m[:, mtab(nv, N)['ix'][tuple(e)]] = 1
        return J

    def trunc(self, N):
        if N >= self.N:
            return self
        n = mtab(self.nv, N)['n']
        return KJ(self.m[:, :n], self.r[:, :n], self.nv, N, self.sym)

    def __add__(self, o):
        if not isinstance(o, KJ):
            m = self.m.copy(); r = self.r.copy()
            m[:, 0] = m[:, 0] + o
            r[:, 0] = r[:, 0] * (1 + SAFE) + SAFE * np.abs(m[:, 0]) + TINY
            return KJ(m, r, self.nv, self.N, self.sym)
        N = min(self.N, o.N)
        a, b = self.trunc(N), o.trunc(N)
        m = a.m + b.m
        r = (a.r + b.r) * (1 + SAFE) + SAFE * np.abs(m) + TINY
        return KJ(m, r, self.nv, N, self.sym)

    __radd__ = __add__

    def __neg__(self):
        return KJ(-self.m, self.r, self.nv, self.N, self.sym)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, c):
        """real constant c (the intended constant may differ from the float c by <= u|c|: covered by SAFE |m c|)."""
        m = self.m * c
        return KJ(m, self.r * abs(c) * (1 + SAFE) + SAFE * np.abs(m) + TINY, self.nv, self.N, self.sym)

    def bscale(self, d):
        """multiply by a per-box scalar Disc (nb,)."""
        dm = d.m if np.any(d.m.imag != 0) or np.iscomplexobj(self.m) else d.m.real
        sh = (-1,) + (1,) * (self.m.ndim - 1)
        dm = dm.reshape(sh); dr = d.r.reshape(sh)
        am = np.abs(self.m)
        m = self.m * dm
        r = (am * dr + self.r * np.abs(dm) + self.r * dr) * (1 + SAFE) + SAFE * np.abs(m) + TINY
        return KJ(m, r, self.nv, self.N, self.sym)

    def __mul__(self, o):
        if isinstance(o, KJ):
            assert self.rank == 0 or o.rank == 0
            la = LET[:self.rank]; lb = LET[:o.rank]
            return kmul(self, o, '%s,%s->%s' % (la, lb, la or lb))
        if isinstance(o, Disc):
            return self.bscale(o)
        return self.scale(o)

    __rmul__ = __mul__

    def recip(self):
        """scalar: y0 = 1/a0, y_d = -y0 sum_{d1 >= 1} a_{d1} y_{d - d1} (exact recursion)."""
        st = mtab(self.nv, self.N)['starts']
        y0 = Disc(self.m[:, 0], self.r[:, 0]).recip()
        y0m = y0.m if np.iscomplexobj(self.m) else y0.m.real
        Y = KJ.const(Disc(y0m, y0.r) if True else y0, self.nv, self.N, self.sym)
        Y.m = Y.m.astype(self.m.dtype)
        for d in range(1, self.N + 1):
            bm, br = kblock(self, Y, ',->', self.N, d, 1)
            blk = KJ(bm, br, self.nv, 0, self.sym).bscale(y0)
            Y.m[:, st[d]:st[d + 1]] = -blk.m
            Y.r[:, st[d]:st[d + 1]] = blk.r
        return Y

    def perm(self, spec):
        a, b = spec.split('->')
        ax = [2 + a.index(c) for c in b]
        return KJ(np.transpose(self.m, [0, 1] + ax), np.transpose(self.r, [0, 1] + ax), self.nv, self.N, self.sym)

    def val(self):
        return self.trunc(0)

    def d(self, k):
        T = mtab(self.nv, self.N)
        N = self.N - 1
        n = mtab(self.nv, N)['n']
        if k >= self.nv:
            z = np.zeros((self.nb, n) + self.m.shape[2:], self.m.dtype)
            return KJ(z, np.zeros(z.shape), self.nv, N, self.sym)
        src, fac = T['D'][k]
        f = fac.reshape((1, -1) + (1,) * self.rank)
        m_ = self.m[:, src] * f
        return KJ(m_, self.r[:, src] * f * (1 + SAFE) + SAFE * np.abs(m_) + TINY, self.nv, N, self.sym)      # amend2: rounding guard


def kmul(a, b, spec, extra=()):
    """a * b by spec (plus extra (spec, sign) terms on the same operands, summed in one kernel)."""
    N = min(a.N, b.N)
    specs = ((spec, 1.0),) + tuple(extra)
    tr = triples(a.nv, N, a.N, b.N, specs, a.sym)
    Om, R = _run(a, b, tr)
    sh = (a.nb, mtab(a.nv, N)['n']) + (4,) * tr[6]
    return KJ(Om.reshape(sh), R.reshape(sh), a.nv, N, a.sym)


def kmulx(a, b, specs):
    N = min(a.N, b.N)
    tr = triples(a.nv, N, a.N, b.N, tuple(specs), a.sym)
    Om, R = _run(a, b, tr)
    sh = (a.nb, mtab(a.nv, N)['n']) + (4,) * tr[6]
    return KJ(Om.reshape(sh), R.reshape(sh), a.nv, N, a.sym)


def kblock(a, b, spec, N, d, d1min):
    tr = triples(a.nv, N, a.N, b.N, ((spec, 1.0),), a.sym, (d, d1min))
    Om, R = _run(a, b, tr)
    st = mtab(a.nv, N)['starts']
    sh = (a.nb, st[d + 1] - st[d]) + (4,) * tr[6]
    return Om.reshape(sh), R.reshape(sh)


def _realify(d):
    return Disc(d.m.real, d.r) if np.all(d.m.imag == 0) else d


# ------------------------------------------------------------------ metric in charts (KJ, unpruned)
def metric_x(model, V, N):
    """symplectic chart (x1, x2, th1, th2), interior only: g_K = H dx dx + H^-1 dth dth, H = (1/2) sum_r nu_r nu_r^T / l_r + Hess f."""
    V = [_realify(v) for v in V]
    cplx = any(np.iscomplexobj(v.m) for v in V)
    if cplx:
        V = [Disc(v.m.astype(complex), v.r) for v in V]
    nb = len(V[0].m); nv = 2

    def K(v):
        d = Disc(np.full(nb, v[0]), np.full(nb, v[1]))
        return d if cplx else Disc(d.m.real, d.r)
    x = [KJ.var(V[i], i, nv, N) for i in range(2)]
    zero = KJ.const(Disc(np.zeros(nb, V[0].m.dtype)), nv, N)
    H = [[zero, zero], [zero, zero]]
    for r in range(model.d):
        nu = model.nu[r]
        lr = x[0] * nu[0] + x[1] * nu[1] + KJ.const(K(model.lam[r]), nv, N)
        il = lr.recip()
        for i in range(2):
            for j in range(2):
                if nu[i] * nu[j]:
                    H[i][j] = H[i][j] + il * (nu[i] * nu[j] / 2)
    if model.fsrc is not None:
        D = model.fsrc(Disc(V[0].m.astype(complex), V[0].r), Disc(V[1].m.astype(complex), V[1].r))
        D = {al: (Disc(dd.m.real, dd.r + np.abs(dd.m.imag)) if not cplx else dd) for al, dd in D.items()}
        dx = [x[0] - KJ.const(V[0], nv, N), x[1] - KJ.const(V[1], nv, N)]
        for t in dx:
            t.m[:, 0] = 0; t.r[:, 0] = 0
        one = KJ.const(Disc(np.ones(nb, V[0].m.dtype)), nv, N)
        p1 = [one]; p2 = [one]
        for _ in range(N):
            p1.append(p1[-1] * dx[0]); p2.append(p2[-1] * dx[1])
        for i in range(2):
            for j in range(2):
                al = [0, 0]; al[i] += 1; al[j] += 1
                acc = zero
                for a1 in range(N + 1):
                    for a2 in range(N + 1 - a1):
                        acc = acc + ((p1[a1] * p2[a2]) * (1.0 / (factorial(a1) * factorial(a2)))).bscale(D[(al[0] + a1, al[1] + a2)])
                H[i][j] = H[i][j] + acc
    idt = (H[0][0] * H[1][1] - H[0][1] * H[1][0]).recip()
    G = [[H[1][1] * idt, -(H[0][1] * idt)], [-(H[1][0] * idt), H[0][0] * idt]]
    n = mtab(nv, N)['n']
    m = np.zeros((nb, n, 4, 4), V[0].m.dtype); r = np.zeros((nb, n, 4, 4))
    for i in range(2):
        for j in range(2):
            m[:, :, i, j] = H[i][j].m; r[:, :, i, j] = H[i][j].r
            m[:, :, 2 + i, 2 + j] = G[i][j].m; r[:, :, 2 + i, 2 + j] = G[i][j].r
    g = KJ(m, r, nv, N, SYM['xchart'])
    S = x[0] * K(model.Sabc[1]) + x[1] * K(model.Sabc[2]) + KJ.const(K(model.Sabc[0]), nv, N)
    S.sym = SYM['xchart']
    assert model.cpow == 0
    return g, [V[0], V[1]], S


def metric(model, k, chart, V, N):
    if chart == 'xchart':
        return metric_x(model, V, N)
    nv = CHARTS[chart]
    V = [_realify(v) for v in V]
    cplx = any(np.iscomplexobj(v.m) for v in V)
    nb = len(V[0].m)

    def K(v):
        d = Disc(np.full(nb, v[0]), np.full(nb, v[1]))
        return d if cplx else Disc(d.m.real, d.r)
    if cplx:
        V = [Disc(v.m.astype(complex), v.r) for v in V]
    a_i, b_i = k, (k + 1) % model.d
    Am = [list(model.nu[a_i]), list(model.nu[b_i])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    v = [KJ.var(V[i], i, nv, N) for i in range(nv)]
    if chart == 'polar':
        p, q = v; a = p * p; b = q * q
    elif chart == 'cart':
        X, Y, Z, W = v; a = X * X + Y * Y; b = Z * Z + W * W
    elif chart == 'rodA':
        X, Y, q = v; a = X * X + Y * Y; b = q * q
    else:
        p, Z, W = v; a = p * p; b = Z * Z + W * W
    A_ = a - KJ.const(K(model.lam[a_i]), nv, N); B_ = b - KJ.const(K(model.lam[b_i]), nv, N)
    Xx = [A_ * Ai[i][0] + B_ * Ai[i][1] for i in range(2)]
    x0 = [Disc(Xx[0].m[:, 0], Xx[0].r[:, 0]), Disc(Xx[1].m[:, 0], Xx[1].r[:, 0])]
    zero = KJ.const(Disc(np.zeros(nb, V[0].m.dtype)), nv, N)
    R = [[zero, zero], [zero, zero]]
    for r in range(model.d):
        if r in (a_i, b_i):
            continue
        nvv = [Ai[0][i] * model.nu[r][0] + Ai[1][i] * model.nu[r][1] for i in range(2)]
        lr = Xx[0] * model.nu[r][0] + Xx[1] * model.nu[r][1] + KJ.const(K(model.lam[r]), nv, N)
        il = lr.recip()
        for i in range(2):
            for j in range(2):
                if nvv[i] * nvv[j]:
                    R[i][j] = R[i][j] + il * (nvv[i] * nvv[j] / 2)
    if model.fsrc is not None:
        D = model.fsrc(Disc(x0[0].m.astype(complex), x0[0].r), Disc(x0[1].m.astype(complex), x0[1].r))
        D = {al: (Disc(dd.m.real, dd.r + np.abs(dd.m.imag)) if not cplx else dd) for al, dd in D.items()}
        dx = [Xx[0] - KJ.const(x0[0], nv, N), Xx[1] - KJ.const(x0[1], nv, N)]
        for t in dx:
            t.m[:, 0] = 0; t.r[:, 0] = 0
        one = KJ.const(Disc(np.ones(nb, V[0].m.dtype)), nv, N)
        p1 = [one]; p2 = [one]
        for _ in range(N):
            p1.append(p1[-1] * dx[0]); p2.append(p2[-1] * dx[1])
        pw = {(a1, a2): (p1[a1] * p2[a2]) * (1.0 / (factorial(a1) * factorial(a2))) for a1 in range(N + 1) for a2 in range(N + 1 - a1)}

        def H(i, j):
            al = [0, 0]; al[i] += 1; al[j] += 1
            acc = zero
            for (a1, a2), t in pw.items():
                acc = acc + t.bscale(D[(al[0] + a1, al[1] + a2)])
            return acc
        Hf = [[H(0, 0), H(0, 1)], [None, H(1, 1)]]
        Hf[1][0] = Hf[0][1]
        for aa in range(2):
            for bb in range(2):
                for i in range(2):
                    for j in range(2):
                        co = Ai[i][aa] * Ai[j][bb]
                        if co:
                            R[aa][bb] = R[aa][bb] + Hf[i][j] * co
    R00, R01, R11 = R[0][0], R[0][1], R[1][1]
    m11 = 1 + (a * R00) * 2; m12 = (a * R01) * 2; m21 = (b * R01) * 2; m22 = 1 + (b * R11) * 2
    idet = (m11 * m22 - m12 * m21).recip()
    one = KJ.const(Disc(np.ones(nb, V[0].m.dtype)), nv, N)
    if chart == 'polar':
        gpp = 2 + (a * R00) * 4; gpq = (p * q * R01) * 4; gqq = 2 + (b * R11) * 4
        Gaa = (a * m22 * idet) * 2; Gab = -((a * b * R01 * idet) * 4); Gbb = (b * m11 * idet) * 2
        comps = [[gpp, gpq, None, None], [gpq, gqq, None, None], [None, None, Gaa, Gab], [None, None, Gab, Gbb]]
    else:
        Pa = ((R00 * m22) * -2 + (b * R01 * R01) * 4) * idet
        Pb = ((R11 * m11) * -2 + (a * R01 * R01) * 4) * idet
        Qd = (R01 * idet) * -4
        if chart == 'cart':
            al = [X, Y, None, None]; alp = [-Y, X, None, None]; be = [None, None, Z, W]; bep = [None, None, -W, Z]
            terms = [(R00 * 4, al, al), (R11 * 4, be, be), (R01 * 4, al, be), (R01 * 4, be, al),
                     (Pa * 2, alp, alp), (Pb * 2, bep, bep), (Qd, alp, bep), (Qd, bep, alp)]
            diag = [2, 2, 2, 2]
        elif chart == 'rodA':
            al = [X, Y, None, None]; alp = [-Y, X, None, None]; dq = [None, None, one, None]; dtb = [None, None, None, one]
            gqq = 2 + (b * R11) * 4; Gbb = (b * m11 * idet) * 2; c1 = (R01 * q) * 4; c2 = (b * R01 * idet) * -4
            terms = [(R00 * 4, al, al), (Pa * 2, alp, alp), (c1, al, dq), (c1, dq, al), (gqq, dq, dq), (Gbb, dtb, dtb),
                     (c2, alp, dtb), (c2, dtb, alp)]
            diag = [2, 2, 0, 0]
        else:
            be = [None, Z, W, None]; bep = [None, -W, Z, None]; dp = [one, None, None, None]; dta = [None, None, None, one]
            gpp = 2 + (a * R00) * 4; Gaa = (a * m22 * idet) * 2; c1 = (R01 * p) * 4; c2 = (a * R01 * idet) * -4
            terms = [(R11 * 4, be, be), (Pb * 2, bep, bep), (c1, dp, be), (c1, be, dp), (gpp, dp, dp), (Gaa, dta, dta),
                     (c2, bep, dta), (c2, dta, bep)]
            diag = [0, 2, 2, 0]
        comps = [[None] * 4 for _ in range(4)]
        for i in range(4):
            if diag[i]:
                comps[i][i] = one * float(diag[i])
        for co, u, w in terms:
            for i in range(4):
                if u[i] is None:
                    continue
                for j in range(4):
                    if w[j] is None:
                        continue
                    t = co * u[i] * w[j]
                    comps[i][j] = t if comps[i][j] is None else comps[i][j] + t
    S = Xx[0] * K(model.Sabc[1]) + Xx[1] * K(model.Sabc[2]) + KJ.const(K(model.Sabc[0]), nv, N)
    if model.cpow == -2:
        Si = S.recip()
        cf = (Si * Si) * K(model.c)
    else:
        cf = one
    n = mtab(nv, N)['n']
    m = np.zeros((nb, n, 4, 4), V[0].m.dtype); r = np.zeros((nb, n, 4, 4))
    for i in range(4):
        for j in range(4):
            if comps[i][j] is not None:
                t = cf * comps[i][j]
                m[:, :, i, j] = t.m; r[:, :, i, j] = t.r
    g = KJ(m, r, nv, N, SYM[chart])
    # exact zeros of J-odd entries (symmetry argument)
    pm = _parity_mono(nv, N, SYM[chart])
    fl = np.array(SYM[chart][0])
    pt = (np.isin(np.arange(4), fl)[:, None].astype(int) + np.isin(np.arange(4), fl)[None, :].astype(int))
    odd = ((pm[:, None, None] + pt[None]) % 2 == 1)
    g.m[:, odd] = 0; g.r[:, odd] = 0
    S.sym = SYM[chart]
    return g, x0, S


# ------------------------------------------------------------------ curvature pipeline
def kinv(G):
    T = mtab(G.nv, G.N); st = T['starts']
    g0 = Disc(np.moveaxis(G.m[:, 0], 0, -1).astype(complex), np.moveaxis(G.r[:, 0], 0, -1))
    Y0 = jets.inv0(g0)
    Y0m = np.moveaxis(Y0.m, -1, 0); Y0r = np.moveaxis(Y0.r, -1, 0)
    if not np.iscomplexobj(G.m):
        Y0m = Y0m.real
    m = np.zeros(G.m.shape, Y0m.dtype); r = np.zeros(G.m.shape)
    m[:, 0] = Y0m; r[:, 0] = Y0r
    Y = KJ(m, r, G.nv, G.N, G.sym)
    Y0j = KJ(Y0m[:, None], Y0r[:, None], G.nv, 0, G.sym)
    for d in range(1, G.N + 1):
        am, ar = kblock(G, Y, 'ij,jk->ik', G.N, d, 1)
        s = 'zij,zPjk->zPik'
        mm = np.einsum(s, Y0m, am)
        bb = np.einsum(s, np.abs(Y0m) + Y0r, np.abs(am) + ar)
        ss = np.einsum(s, np.abs(Y0m), np.abs(am))
        Y.m[:, st[d]:st[d + 1]] = -mm
        Y.r[:, st[d]:st[d + 1]] = (bb - ss) + (3 * 4 + 20) * U2 * bb + TINY
    return Y


def dcoord(T):
    parts = [T.d(k) for k in range(4)]
    m = np.stack([p.m for p in parts], axis=2); r = np.stack([p.r for p in parts], axis=2)
    return KJ(m, r, T.nv, T.N - 1, T.sym)


def christoffel(g, gi):
    dg = dcoord(g)                                   # [m, k, n] = d_m g_kn
    Cc = dg.perm('mkn->kmn') + dg.perm('nkm->kmn') - dg
    return kmul(gi.trunc(Cc.N), Cc, 'lk,kmn->lmn').scale(0.5)


def riemann(g, Gam):
    dG = dcoord(Gam)                                 # [c, l, m, n]
    Gt = Gam.trunc(dG.N)
    Rup = (dG.perm('mrns->rsmn') - dG.perm('nrms->rsmn')) + kmul(Gt, Gt, 'rml,lns->rsmn', extra=(('rnl,lms->rsmn', -1.0),))
    return kmul(g.trunc(Rup.N), Rup, 'rl,lsmn->rsmn')


def nabla(T, Gam):
    rank = T.rank
    out = dcoord(T)
    if rank == 0:
        return out
    idx = LET[:rank]
    specs = []
    for i in range(rank):
        src = idx[:i] + 'z' + idx[i + 1:]
        specs.append(('zm%s,%s->m%s' % (idx[i], src, idx), 1.0))
    corr = kmulx(Gam.trunc(T.N - 1), T.trunc(T.N - 1), specs)
    return out - corr


def sqnorm(T, gi):
    rank = T.rank
    T0 = T.val(); G0 = gi.val()
    v = T0
    idx = list(LET[:rank])
    for s in range(rank):
        o = idx.copy(); o[s] = 'z'
        v = kmul(G0, v, '%sz,%s->%s' % (idx[s], ''.join(idx), ''.join(o)))
    s = ''.join(idx)
    out = kmul(v, T0, '%s,%s->' % (s, s))
    return Disc(out.m[:, 0], out.r[:, 0])


CH = dict(CHARTS); CH['xchart'] = 2


def pipeline(model, k, chart, V, qmax=3, want=('riem', 'E')):
    """'riem': |nabla^q Riem|^2_g (q <= qmax); 'E': |nabla^j E|^2_g, E = Ric - 3 g; 'abreu' (model with cpow = 0):
    psi = Scal(g_K) - S(x) and |nabla^j psi|^2_{g_K} (j <= qmax).  All Disc arrays over the boxes (bilinear contractions)."""
    N = qmax + 2
    g, x0, S = metric(model, k, chart, V, N)
    gi = kinv(g)
    Gam = christoffel(g, gi)
    Rdn = riemann(g, Gam)
    out = {'x0': x0}
    if 'riem' in want:
        T = Rdn; L = []
        for q in range(qmax + 1):
            L.append(sqnorm(T, gi))
            if q < qmax:
                T = nabla(T, Gam)
        out['riem'] = L
    if 'kappa' in want:
        out['R0'] = (Rdn.m[:, 0], Rdn.r[:, 0]); out['g0'] = (g.m[:, 0], g.r[:, 0])
    if 'E' in want or 'abreu' in want or 'ric' in want:
        Ric = kmul(gi.trunc(Rdn.N), Rdn, 'rm,rsmn->sn')
    if 'ric' in want:
        # Rc_q = |nabla^q Ric|^2 (q <= qmax), R0_q = |nabla^q Ric_0|^2 (q <= qmax-1), S_q = |nabla^q ds|^2 (q <= qmax-1),
        # s value; Kahler: Ric eigenvalues are s/4 +- |Ric_0|/2 (each of multiplicity 2)
        sc = kmul(gi.trunc(Ric.N), Ric, 'sn,sn->')
        out['s'] = Disc(sc.m[:, 0], sc.r[:, 0])
        T = Ric; L = []
        for q in range(qmax + 1):
            L.append(sqnorm(T, gi))
            if q < qmax:
                T = nabla(T, Gam)
        out['Rc'] = L
        R0 = Ric - kmul(sc, g.trunc(sc.N), ',ab->ab').scale(0.25)
        T = R0; L = []
        for q in range(qmax):
            L.append(sqnorm(T, gi))
            if q < qmax - 1:
                T = nabla(T, Gam)
        out['R0c'] = L
        T = nabla(sc, Gam); L = []
        for q in range(qmax):
            L.append(sqnorm(T, gi))
            if q < qmax - 1:
                T = nabla(T, Gam)
        out['Sd'] = L
    if 'E' in want:
        T = Ric - g.trunc(Ric.N).scale(3.0); L = []
        for j in range(qmax + 1):
            L.append(sqnorm(T, gi))
            if j < qmax:
                T = nabla(T, Gam)
        out['E'] = L
    if 'abreu' in want:
        scal = kmul(gi.trunc(Ric.N), Ric, 'sn,sn->')
        psi = scal - S.trunc(scal.N)
        out['psi'] = Disc(psi.m[:, 0], psi.r[:, 0])
        T = psi; L = []
        for j in range(qmax + 1):
            L.append(sqnorm(T, gi))
            if j == 2:                                    # positive Laplacian  Delta psi = -g^ab nabla_a nabla_b psi
                lp = kmul(gi.trunc(T.N), T, 'ab,ab->').scale(-1.0)
                out['lap'] = Disc(lp.m[:, 0], lp.r[:, 0])
            if j < qmax:
                T = nabla(T, Gam)
        out['abreu'] = L
    return out


# ------------------------------------------------------------------ output jets (for centred / Taylor-model forms)
def sqnorm_jet(T, gi, L):
    """|T|^2_g as a jet of order L (bilinear)."""
    rank = T.rank
    T0 = T.trunc(L); G0 = gi.trunc(L)
    v = T0
    idx = list(LET[:rank])
    for s in range(rank):
        o = idx.copy(); o[s] = 'z'
        v = kmul(G0, v, '%sz,%s->%s' % (idx[s], ''.join(idx), ''.join(o)))
    s = ''.join(idx)
    return kmul(v, T0, '%s,%s->' % (s, s))


def pipeline_jet(model, k, chart, V, qmax, L, want=('riem',)):
    """as pipeline, but every output is a scalar KJ of order L in the chart variables (Taylor coefficients at the input
    points/boxes).  'riem': |nabla^q Riem|^2; 'abreu': psi and Delta psi (positive Laplacian)."""
    N = qmax + 2 + L
    g, x0, S = metric(model, k, chart, V, N)
    gi = kinv(g)
    Gam = christoffel(g, gi)
    Rdn = riemann(g, Gam)
    out = {}
    if 'riem' in want:
        T = Rdn; Ls = []
        for q in range(qmax + 1):
            Ls.append(sqnorm_jet(T, gi, L))
            if q < qmax:
                T = nabla(T, Gam)
        out['riem'] = Ls
    if 'abreu' in want:
        Ric = kmul(gi.trunc(Rdn.N), Rdn, 'rm,rsmn->sn')
        scal = kmul(gi.trunc(Ric.N), Ric, 'sn,sn->')
        psi = scal - S.trunc(scal.N)
        out['psi'] = psi.trunc(L)
        H2 = nabla(nabla(psi, Gam), Gam)
        out['lap'] = kmul(gi.trunc(H2.N), H2, 'ab,ab->').scale(-1.0).trunc(L)
    return out


def centred_upper(Jc, Jb, L, r):
    """upper bound of |F| over the box |d_i| <= r (componentwise, active vars): Jc = jet of F at the centre (order L-1 used),
    Jb = jet of F over the box (order-L coefficients used).  Returns array (nb,)."""
    T = mtab(Jc.nv, L)
    st = T['starts']
    tot = np.zeros(Jc.nb)
    for d in range(L):
        blk = slice(st[d], st[d + 1])
        tot = tot + (np.abs(Jc.m[:, blk]) + Jc.r[:, blk]).sum(1) * r ** d
    blk = slice(st[L], st[L + 1])
    tot = tot + (np.abs(Jb.m[:, blk]) + Jb.r[:, blk]).sum(1) * r ** L
    return tot * (1 + 1e-12)
