"""jets.py -- seat E: vectorised ball (complex-disc) tensor jets and the curvature pipeline in four chart types.

TJet: truncated Taylor polynomial of order N in nv chart variables, tensor-valued; coefficients stored as two arrays
  m (complex mids) and r (radii), shape (Nmono, *tshape, Nbox).  Inputs may be boxes (discs), so every coefficient encloses the
  corresponding Taylor coefficient at EVERY point of the box (analytic functions; standard ball-jet argument).
Rounding: see common.py.  Products are formed blockwise by degree with numpy einsum; the radius of a product block is
  sum (|am| + ar)(|bm| + br) - sum |am||bm|   (exact identity for the disc product)  plus  errfac * sum (|am| + ar)(|bm| + br),
  errfac = (3 K + 20) * 2.3e-16 with K >= the number of float products summed into any one output entry (includes the
  monomial pairs), which dominates the binary64 rounding of mids, big and small sums and their difference.

Charts (coordinates: active variables first, the rest cyclic Killing angles):
  'polar' (p, q, ta, tb)         p = sqrt l_a, q = sqrt l_b              (seat W / seat R vertex chart)
  'cart'  (X, Y, Z, W)            X + iY = p e^{i ta}, Z + iW = q e^{i tb}  (vertex chart, smooth at rods and vertex)
  'rodA'  (X, Y, q, tb)           Cartesian at rod a only
  'rodB'  (p, Z, W, ta)           Cartesian at rod b only
Metric g_K (seat W's blocks) and its Cartesian form (derived in NOTES-E.md sec. 2):
  R = sum_{r != a,b} (1/2) n_r n_r^T / l_r + Ai^T Hess f Ai,  a = l_a, b = l_b,  m11 = 1 + 2aR00, m12 = 2aR01, m21 = 2bR01,
  m22 = 1 + 2bR11, det = m11 m22 - m12 m21,  Pa = (-2R00 m22 + 4bR01^2)/det, Pb = (-2R11 m11 + 4aR01^2)/det, Q = -4R01/det:
  g_K = 2 dX^2 + 2 dY^2 + 2 dZ^2 + 2 dW^2 + 4R00 al al + 4R11 be be + 4R01 (al be + be al) + 2Pa al' al' + 2Pb be' be'
        + Q (al' be' + be' al'),     al = X dX + Y dY, al' = X dY - Y dX, be = Z dZ + W dW, be' = Z dW - W dZ.
  g = c S^cpow g_K (cpow = -2: the CLW Einstein candidate; cpow = 0: the Kahler metric itself, for Abreu's Scal).
"""
import itertools
from math import factorial
import numpy as np
import common
from common import Disc, SAFE

U2 = 2.3e-16
LET = 'abcdefghijklmnopqrstuvwxy'


# ------------------------------------------------------------------ monomials
_MT = {}


def mtab(nv, N):
    key = (nv, N)
    if key in _MT:
        return _MT[key]
    mons, starts = [], []
    for d in range(N + 1):
        starts.append(len(mons))
        blk = [e for e in itertools.product(range(d + 1), repeat=nv) if sum(e) == d]
        blk.sort(reverse=True)
        mons += blk
    starts.append(len(mons))
    ix = {e: i for i, e in enumerate(mons)}
    # scatter matrices for block products (d1, d2)
    S = {}
    for d1 in range(N + 1):
        for d2 in range(N + 1 - d1):
            B1 = mons[starts[d1]:starts[d1 + 1]]; B2 = mons[starts[d2]:starts[d2 + 1]]
            o0 = starts[d1 + d2]; no = starts[d1 + d2 + 1] - o0
            tg = np.array([ix[tuple(x + y for x, y in zip(e1, e2))] - o0 for e1 in B1 for e2 in B2])
            order = np.argsort(tg, kind='stable')
            st_ = np.flatnonzero(np.r_[True, np.diff(tg[order]) != 0])
            assert len(st_) == no
            ident = len(tg) == no and np.all(order == np.arange(no))
            S[(d1, d2)] = (order, st_, no, ident)
    # derivative maps
    D = []
    for k in range(nv):
        src, fac = [], []
        for e in mons:
            if sum(e) == N:
                break
            e2 = list(e); e2[k] += 1
            src.append(ix[tuple(e2)]); fac.append(e[k] + 1)
        D.append((np.array(src, dtype=int), np.array(fac, dtype=float)))
    _MT[key] = dict(mons=mons, starts=starts, ix=ix, S=S, D=D, n=len(mons))
    return _MT[key]


def _spec3(spec):
    """'ij,jk->ik' -> (lhs a, lhs b, out, K = product of contracted dims placeholder letters)"""
    lhs, out = spec.split('->')
    sa, sb = lhs.split(',')
    return sa, sb, out


def disc_einsum(spec, am, ar, bm, br, dims=4):
    """einsum over tensor letters with a trailing box axis (letter 'Z' appended). returns mid, rad (unfinished: big, sml)."""
    sa, sb, so = _spec3(spec)
    s = sa + 'Z,' + sb + 'Z->' + so + 'Z'
    mid = np.einsum(s, am, bm)
    big = np.einsum(s, np.abs(am) + ar, np.abs(bm) + br)
    sml = np.einsum(s, np.abs(am), np.abs(bm))
    K = dims ** len(set(sa + sb) - set(so))
    return mid, big, sml, K


class TJet:
    __slots__ = ('m', 'r', 'nv', 'N', 'ts')

    def __init__(self, m, r, nv, N):
        self.m, self.r, self.nv, self.N = m, r, nv, N
        self.ts = m.shape[1:-1]

    @property
    def nb(self):
        return self.m.shape[-1]

    @staticmethod
    def zeros(nv, N, ts, nb):
        n = mtab(nv, N)['n']
        return TJet(np.zeros((n,) + tuple(ts) + (nb,), complex), np.zeros((n,) + tuple(ts) + (nb,)), nv, N)

    @staticmethod
    def const(d, nv, N, ts=()):
        """d: Disc with shape ts + (nb,)"""
        J = TJet.zeros(nv, N, ts, d.m.shape[-1])
        J.m[0] = d.m; J.r[0] = d.r
        return J

    @staticmethod
    def var(d, k, nv, N):
        """coordinate jet v_k = d + dv_k (d: Disc (nb,))."""
        J = TJet.const(d, nv, N)
        e = [0] * nv; e[k] = 1
        J.m[mtab(nv, N)['ix'][tuple(e)]] = 1
        return J

    def trunc(self, N):
        if N >= self.N:
            return self
        n = mtab(self.nv, N)['n']
        return TJet(self.m[:n], self.r[:n], self.nv, N)

    def val(self):
        return Disc(self.m[0], self.r[0])

    def __add__(self, o):
        if isinstance(o, TJet):
            N = min(self.N, o.N)
            a, b = self.trunc(N), o.trunc(N)
            m = a.m + b.m
            r = (a.r + b.r) * (1 + SAFE) + SAFE * np.abs(m) + common.TINY
            return TJet(m, r, self.nv, N)
        if isinstance(o, Disc):
            return self + TJet.const(o, self.nv, self.N, self.ts)
        m = self.m.copy(); r = self.r.copy()
        m[0] = m[0] + o
        r[0] = r[0] * (1 + SAFE) + SAFE * np.abs(m[0]) + common.TINY
        return TJet(m, r, self.nv, self.N)

    __radd__ = __add__

    def __neg__(self):
        return TJet(-self.m, self.r, self.nv, self.N)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, d):
        """multiply by a box scalar Disc (nb,) or a python number."""
        if not isinstance(d, Disc):
            d = Disc(np.full(self.nb, complex(d)))
        dm, dr = d.m, d.r
        am = np.abs(self.m)
        m = self.m * dm
        r = (am * dr + self.r * np.abs(dm) + self.r * dr) * (1 + SAFE) + SAFE * am * np.abs(dm) + common.TINY
        return TJet(m, r, self.nv, self.N)

    def mul(self, o, spec=None):
        if not isinstance(o, TJet):
            return self.scale(o)
        if spec is None:
            la = LET[:len(self.ts)]; lb = LET[:len(o.ts)]
            if len(self.ts) and len(o.ts) and self.ts != o.ts:
                raise ValueError('need spec')
            spec = ('%s,%s->%s' % (la, lb, la or lb)) if (not la or not lb or la == lb) else None
        N = min(self.N, o.N)
        T = mtab(self.nv, N)
        st = T['starts']
        a, b = self.trunc(N), o.trunc(N)
        outm = outb = outs = None
        for d in range(N + 1):
            accm, accb, accs = _bsum(a, b, d, 0, spec, T)
            if outm is None:
                outm = np.zeros((T['n'],) + accm.shape[1:], complex)
                outb = np.zeros(outm.shape); outs = np.zeros(outm.shape)
            outm[st[d]:st[d + 1]] = accm; outb[st[d]:st[d + 1]] = accb; outs[st[d]:st[d + 1]] = accs
        return TJet(*_finish(outm, outb, outs, spec, T), self.nv, N)

    def __mul__(self, o):
        if isinstance(o, TJet) and (not self.ts or not o.ts):
            la = LET[:len(self.ts)]; lb = LET[:len(o.ts)]
            return self.mul(o, '%s,%s->%s' % (la, lb, la or lb))
        return self.mul(o)

    __rmul__ = __mul__

    def recip(self):
        """scalar jet reciprocal: y0 = 1/a0, y_d = -y0 sum_{d1 >= 1} a_{d1} y_{d - d1} (exact recursion, ball arithmetic)."""
        assert not self.ts
        T = mtab(self.nv, self.N); st = T['starts']
        y0 = self.val().recip()
        Y = TJet.const(y0, self.nv, self.N)
        for d in range(1, self.N + 1):
            m_, b_, s_ = _bsum(self, Y, d, 1, ',->', T)
            m_, r_ = _finish(m_, b_, s_, ',->', T)
            blk = TJet(m_, r_, self.nv, self.N).scale(y0)        # (block-shaped: scale is elementwise)
            Y.m[st[d]:st[d + 1]] = -blk.m
            Y.r[st[d]:st[d + 1]] = blk.r
        return Y

    def d(self, k):
        T = mtab(self.nv, self.N)
        N = self.N - 1
        n = mtab(self.nv, N)['n']
        if k >= self.nv:
            return TJet.zeros(self.nv, N, self.ts, self.nb)
        src, fac = T['D'][k]
        sh = (-1,) + (1,) * (self.m.ndim - 1)
        f = fac.reshape(sh)
        return TJet(self.m[src] * f, self.r[src] * f * (1 + SAFE), self.nv, N)


def _bsum(a, b, d, d1min, spec, T):
    """degree-d block of a*b summed over d1 = d1min..d (d2 = d - d1): returns mid, big, sml (unfinished)."""
    st = T['starts']
    sa, sb, so = _spec3(spec)
    s = 'P' + sa + 'Z,Q' + sb + 'Z->PQ' + so + 'Z'
    accm = accb = accs = None
    for d1 in range(d1min, d + 1):
        d2 = d - d1
        A = slice(st[d1], st[d1 + 1]); B = slice(st[d2], st[d2 + 1])
        am_, bm_ = a.m[A], b.m[B]
        aa, ba = np.abs(am_), np.abs(bm_)
        em = np.einsum(s, am_, bm_)
        eb = np.einsum(s, aa + a.r[A], ba + b.r[B])
        es = np.einsum(s, aa, ba)
        order, sts, no, ident = T['S'][(d1, d2)]
        sh = em.shape
        if ident:
            f = lambda x: x.reshape((no,) + sh[2:])
        else:
            f = lambda x: np.add.reduceat(x.reshape((sh[0] * sh[1],) + sh[2:])[order], sts, axis=0)
        em, eb, es = f(em), f(eb), f(es)
        accm = em if accm is None else accm + em
        accb = eb if accb is None else accb + eb
        accs = es if accs is None else accs + es
    return accm, accb, accs


def _finish(m, big, sml, spec, T):
    sa, sb, so = _spec3(spec)
    K = 4 ** len(set(sa + sb) - set(so)) * T['n']
    errfac = (3 * K + 20) * U2
    r = (big - sml) + errfac * big + common.TINY
    bad = ~np.isfinite(m) | ~np.isfinite(r)
    if bad.any():
        r = np.where(bad, np.inf, r); m = np.where(bad, 0, m)
    return m, r


def _block(J, d):
    """jet keeping only the degree-d homogeneous part."""
    T = mtab(J.nv, J.N); st = T['starts']
    m = np.zeros_like(J.m); r = np.zeros_like(J.r)
    m[st[d]:st[d + 1]] = J.m[st[d]:st[d + 1]]; r[st[d]:st[d + 1]] = J.r[st[d]:st[d + 1]]
    return TJet(m, r, J.nv, J.N)


def stack(comps, nv, N, nb):
    """nested list (shape ts) of scalar TJets / None -> tensor TJet."""
    arr = np.array(comps, dtype=object)
    ts = arr.shape
    J = TJet.zeros(nv, N, ts, nb)
    for idx in np.ndindex(*ts):
        e = arr[idx]
        if e is None:
            continue
        e = e.trunc(N)
        J.m[(slice(None),) + idx] = e.m
        J.r[(slice(None),) + idx] = e.r
    return J


def comp(J, *idx):
    return TJet(J.m[(slice(None),) + idx], J.r[(slice(None),) + idx], J.nv, J.N)


# ------------------------------------------------------------------ verified matrix inverse (order 0) and jet inverse
def inv0(G):
    """G: Disc mats with shape (n, n, nb).  Returns Disc enclosing G^-1 for every matrix in the disc set.
    X = float inverse; E = I - G X (disc); e = ||E||_inf < 1 => G^-1 = X + X E (I - E)^-1, |entry| <= ||X||_inf e / (1 - e)."""
    n = G.m.shape[0]
    Gm = np.moveaxis(G.m, -1, 0)
    bad = None
    try:
        X = np.linalg.inv(Gm)                                 # (nb, n, n)
    except np.linalg.LinAlgError:
        # amend3 fail-safe: a singular (or non-finite) float mid in some box.  Invert box by box; every box whose mid is
        # non-finite or singular gets mid 0 and radius inf (the enclosure is then the whole plane: trivially valid).
        bad = ~np.isfinite(Gm).all(axis=(1, 2))
        X = np.zeros_like(Gm)
        for z in np.nonzero(~bad)[0]:
            try:
                X[z] = np.linalg.inv(Gm[z])
            except np.linalg.LinAlgError:
                bad[z] = True
    if not np.isfinite(X).all():                              # amend3: NaN/inf float inverse (no exception) -> same rule
        bad = (bad if bad is not None else np.zeros(X.shape[0], bool)) | ~np.isfinite(X).all(axis=(1, 2))
    if bad is not None:
        X = np.where(bad[:, None, None], 0, X)
    Xd = np.moveaxis(X, 0, -1)
    em = np.einsum('ijZ,jkZ->ikZ', G.m, Xd)
    eb = np.einsum('ijZ,jkZ->ikZ', np.abs(G.m) + G.r, np.abs(Xd))
    es = np.einsum('ijZ,jkZ->ikZ', np.abs(G.m), np.abs(Xd))
    Er = (eb - es) + (3 * n + 20) * U2 * eb
    Em = np.eye(n)[:, :, None] - em
    e = (np.abs(Em) + Er).sum(axis=1).max(axis=0) * (1 + SAFE)
    xn = np.abs(Xd).sum(axis=1).max(axis=0) * (1 + SAFE)
    with np.errstate(divide='ignore', invalid='ignore'):
        b = np.where(e < 1, xn * e / (1 - e) * (1 + SAFE), np.inf)
    if bad is not None:
        b = np.where(bad, np.inf, b)                          # amend3
    return Disc(Xd, np.broadcast_to(b, Xd.shape).copy())


def inv_jet(G):
    """matrix-valued jet inverse: Y0 = inv0(G0), Y_d = -Y0 sum_{d1>=1} G_{d1} Y_{d-d1} (ball arithmetic, exact recursion)."""
    T = mtab(G.nv, G.N); st = T['starts']
    Y0 = inv0(G.val())
    Y = TJet.const(Y0, G.nv, G.N, G.ts)
    s = 'ijZ,PjkZ->PikZ'
    for d in range(1, G.N + 1):
        m_, b_, s_ = _bsum(G, Y, d, 1, 'ij,jk->ik', T)
        am, ar = _finish(m_, b_, s_, 'ij,jk->ik', T)
        mm = np.einsum(s, Y0.m, am)
        bb = np.einsum(s, np.abs(Y0.m) + Y0.r, np.abs(am) + ar)
        ss = np.einsum(s, np.abs(Y0.m), np.abs(am))
        Y.m[st[d]:st[d + 1]] = -mm
        Y.r[st[d]:st[d + 1]] = (bb - ss) + (3 * 4 + 20) * U2 * bb + common.TINY
    return Y


# ------------------------------------------------------------------ model and metric in charts
class Model:
    """polygon data as Discs; fsrc(x1, x2) -> dict alpha -> Disc of D^alpha f over the boxes (or None: f = 0)."""

    def __init__(self, nus, lam, Sabc, c, cpow=-2, fsrc=None):
        self.nu = [tuple(int(v) for v in n) for n in nus]
        self.d = len(self.nu)
        self.lam = lam          # list of (complex mid, rad)
        self.Sabc = Sabc
        self.c = c
        self.cpow = cpow
        self.fsrc = fsrc

    def K(self, v, nb):
        return Disc(np.full(nb, v[0]), np.full(nb, v[1]))


def model_from_prob(prob, fsrc, cpow=-2):
    from flint import arb
    lam = [common.disc_from_acb(v) for v in prob.lam]
    S = [common.disc_from_acb(v) for v in prob.Sabc]
    c = common.disc_from_acb(prob.c)
    return Model(prob.P.nu, lam, S, c, cpow, fsrc)


CHARTS = {'polar': 2, 'cart': 4, 'rodA': 3, 'rodB': 3}


def metric(model, k, chart, V, N):
    """V: list of Disc (nb,) values of the active coordinates at the boxes.  Returns g (4,4) TJet, x0 (Disc pair), S jet."""
    nv = CHARTS[chart]
    nb = V[0].m.shape[0]
    a_i, b_i = k, (k + 1) % model.d
    Am = [list(model.nu[a_i]), list(model.nu[b_i])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    v = [TJet.var(V[i], i, nv, N) for i in range(nv)]
    if chart == 'polar':
        p, q = v; a = p * p; b = q * q
    elif chart == 'cart':
        X, Y, Z, W = v; a = X * X + Y * Y; b = Z * Z + W * W
    elif chart == 'rodA':
        X, Y, q = v; a = X * X + Y * Y; b = q * q
    else:
        p, Z, W = v; a = p * p; b = Z * Z + W * W
    la, lb = model.K(model.lam[a_i], nb), model.K(model.lam[b_i], nb)
    A_ = a - TJet.const(la, nv, N); B_ = b - TJet.const(lb, nv, N)
    Xx = [A_.scale(Ai[i][0]) + B_.scale(Ai[i][1]) for i in range(2)]
    x0 = [Xx[0].val(), Xx[1].val()]
    R = [[None, None], [None, None]]
    for r in range(model.d):
        if r in (a_i, b_i):
            continue
        nvv = [Ai[0][i] * model.nu[r][0] + Ai[1][i] * model.nu[r][1] for i in range(2)]
        lr = Xx[0].scale(model.nu[r][0]) + Xx[1].scale(model.nu[r][1]) + TJet.const(model.K(model.lam[r], nb), nv, N)
        il = lr.recip()
        for i in range(2):
            for j in range(2):
                if nvv[i] * nvv[j]:
                    t = il.scale(nvv[i] * nvv[j] / 2)
                    R[i][j] = t if R[i][j] is None else R[i][j] + t
    if model.fsrc is not None:
        D = model.fsrc(x0[0], x0[1])
        dx = [Xx[0] - TJet.const(x0[0], nv, N), Xx[1] - TJet.const(x0[1], nv, N)]
        for t in dx:
            t.m[0] = 0; t.r[0] = 0
        one = TJet.const(Disc(np.ones(nb)), nv, N)
        pw = {}
        p1 = [one]
        for _ in range(N):
            p1.append(p1[-1] * dx[0])
        p2 = [one]
        for _ in range(N):
            p2.append(p2[-1] * dx[1])
        for a1 in range(N + 1):
            for a2 in range(N + 1 - a1):
                pw[(a1, a2)] = (p1[a1] * p2[a2]).scale(1.0 / (factorial(a1) * factorial(a2)))

        def H(i, j):
            al = [0, 0]; al[i] += 1; al[j] += 1
            acc = None
            for (a1, a2), t in pw.items():
                u = t.scale(D[(al[0] + a1, al[1] + a2)])
                acc = u if acc is None else acc + u
            return acc
        Hf = [[H(0, 0), H(0, 1)], [None, H(1, 1)]]
        Hf[1][0] = Hf[0][1]
        for aa in range(2):
            for bb in range(2):
                for i in range(2):
                    for j in range(2):
                        co = Ai[i][aa] * Ai[j][bb]
                        if co:
                            t = Hf[i][j].scale(co)
                            R[aa][bb] = t if R[aa][bb] is None else R[aa][bb] + t
    zero = TJet.zeros(nv, N, (), nb)
    R = [[zero if R[i][j] is None else R[i][j] for j in range(2)] for i in range(2)]
    R00, R01, R11 = R[0][0], R[0][1], R[1][1]
    m11 = 1 + (a * R00).scale(2); m12 = (a * R01).scale(2); m21 = (b * R01).scale(2); m22 = 1 + (b * R11).scale(2)
    idet = (m11 * m22 - m12 * m21).recip()
    if chart == 'polar':
        gpp = 2 + (a * R00).scale(4); gpq = (p * q * R01).scale(4); gqq = 2 + (b * R11).scale(4)
        Gaa = (a * m22 * idet).scale(2); Gab = -(a * b * R01 * idet).scale(4); Gbb = (b * m11 * idet).scale(2)
        Z0 = None
        comps = [[gpp, gpq, Z0, Z0], [gpq, gqq, Z0, Z0], [Z0, Z0, Gaa, Gab], [Z0, Z0, Gab, Gbb]]
    else:
        Pa = ((R00 * m22).scale(-2) + (b * R01 * R01).scale(4)) * idet
        Pb = ((R11 * m11).scale(-2) + (a * R01 * R01).scale(4)) * idet
        Qd = (R01 * idet).scale(-4)
        if chart == 'cart':
            al = [X, Y, None, None]; alp = [-Y, X, None, None]
            be = [None, None, Z, W]; bep = [None, None, -W, Z]
            terms = [(R00.scale(4), al, al), (R11.scale(4), be, be), (R01.scale(4), al, be), (R01.scale(4), be, al),
                     (Pa.scale(2), alp, alp), (Pb.scale(2), bep, bep), (Qd, alp, bep), (Qd, bep, alp)]
            diag = [2, 2, 2, 2]
        elif chart == 'rodA':
            # coords (X, Y, q, tb)
            al = [X, Y, None, None]; alp = [-Y, X, None, None]
            dq = [None, None, one_(nv, N, nb), None]; dtb = [None, None, None, one_(nv, N, nb)]
            gqq = 2 + (b * R11).scale(4); Gbb = (b * m11 * idet).scale(2)
            terms = [(R00.scale(4), al, al), (Pa.scale(2), alp, alp), ((R01 * q).scale(4), al, dq), ((R01 * q).scale(4), dq, al),
                     (gqq, dq, dq), (Gbb, dtb, dtb), ((b * R01 * idet).scale(-4), alp, dtb), ((b * R01 * idet).scale(-4), dtb, alp)]
            diag = [2, 2, 0, 0]
        else:
            # coords (p, Z, W, ta)
            be = [None, Z, W, None]; bep = [None, -W, Z, None]
            dp = [one_(nv, N, nb), None, None, None]; dta = [None, None, None, one_(nv, N, nb)]
            gpp = 2 + (a * R00).scale(4); Gaa = (a * m22 * idet).scale(2)
            terms = [(R11.scale(4), be, be), (Pb.scale(2), bep, bep), ((R01 * p).scale(4), dp, be), ((R01 * p).scale(4), be, dp),
                     (gpp, dp, dp), (Gaa, dta, dta), ((a * R01 * idet).scale(-4), bep, dta), ((a * R01 * idet).scale(-4), dta, bep)]
            diag = [0, 2, 2, 0]
        comps = [[None] * 4 for _ in range(4)]
        for i in range(4):
            if diag[i]:
                comps[i][i] = TJet.const(Disc(np.full(nb, float(diag[i]))), nv, N)
        for co, u, w in terms:
            for i in range(4):
                if u[i] is None:
                    continue
                for j in range(4):
                    if w[j] is None:
                        continue
                    t = co * u[i] * w[j]
                    comps[i][j] = t if comps[i][j] is None else comps[i][j] + t
    S = Xx[0].scale(model.K(model.Sabc[1], nb)) + Xx[1].scale(model.K(model.Sabc[2], nb)) + TJet.const(model.K(model.Sabc[0], nb), nv, N)
    cK = model.K(model.c, nb)
    if model.cpow == -2:
        Si = S.recip()
        cf = (Si * Si).scale(cK)
    else:                                                    # cpow = 0: the Kahler metric g_K itself (no c)
        cf = TJet.const(Disc(np.ones(nb)), nv, N)
    comps = [[None if e is None else cf * e for e in row] for row in comps]
    g = stack(comps, nv, N, nb)
    return g, x0, S


def one_(nv, N, nb):
    return TJet.const(Disc(np.ones(nb)), nv, N)


# ------------------------------------------------------------------ curvature pipeline
def dcoord(T):
    """d_m T with the derivative index first (m = 0..3; zero for cyclic coordinates)."""
    parts = [T.d(k) for k in range(4)]
    m = np.stack([p.m for p in parts], axis=1)
    r = np.stack([p.r for p in parts], axis=1)
    return TJet(m, r, T.nv, T.N - 1)


def christoffel(g, gi):
    dg = dcoord(g)                                           # [m, a, b]
    Cm = np.einsum('xmkn...->xkmn...', dg.m) + np.einsum('xnkm...->xkmn...', dg.m) - dg.m
    Cr = (np.einsum('xmkn...->xkmn...', dg.r) + np.einsum('xnkm...->xkmn...', dg.r) + dg.r) * (1 + 2 * SAFE)
    Cc = TJet(Cm, Cr + SAFE * np.abs(Cm), g.nv, dg.N)
    return gi.mul(Cc, 'lk,kmn->lmn').scale(0.5)


def tperm(T, spec):
    return TJet(np.einsum('x' + spec.replace('->', '...->x') + '...', T.m), np.einsum('x' + spec.replace('->', '...->x') + '...', T.r), T.nv, T.N)


def riemann(g, gi, Gam):
    dG = dcoord(Gam)                                         # [c, l, m, n]
    A = tperm(dG, 'mrns->rsmn'); B = tperm(dG, 'nrms->rsmn')
    Rup = (A - B) + Gam.mul(Gam, 'rml,lns->rsmn') - Gam.mul(Gam, 'rnl,lms->rsmn')
    return g.mul(Rup, 'rl,lsmn->rsmn')                       # R_{rsmn}


def nabla(T, Gam):
    """covariant derivative, new index first."""
    rank = len(T.ts)
    dT = dcoord(T)
    out = dT
    idx = LET[:rank]
    for i in range(rank):
        src = idx[:i] + 'z' + idx[i + 1:]
        spec = 'zm%s,%s->m%s' % (idx[i], src, idx)
        out = out - Gam.mul(T, spec)
    return out


def sqnorm0(T, gi):
    """|T|^2_g at order 0 (Disc over boxes), bilinear contraction with g^-1 (analytic continuation)."""
    rank = len(T.ts)
    Tm, Tr = T.m[0], T.r[0]
    G = gi.val()
    vm, vr = Tm, Tr
    for s in range(rank):
        idx = list(LET[:rank])
        o = idx.copy(); o[s] = 'z'
        spec = '%szZ,%sZ->%sZ' % (idx[s], ''.join(idx), ''.join(o))
        mid = np.einsum(spec, G.m, vm)
        big = np.einsum(spec, np.abs(G.m) + G.r, np.abs(vm) + vr)
        sml = np.einsum(spec, np.abs(G.m), np.abs(vm))
        vr = (big - sml) + (3 * 4 + 20) * U2 * big
        vm = mid
    s = ''.join(LET[:rank])
    spec = '%sZ,%sZ->Z' % (s, s)
    mid = np.einsum(spec, vm, Tm)
    big = np.einsum(spec, np.abs(vm) + vr, np.abs(Tm) + Tr)
    sml = np.einsum(spec, np.abs(vm), np.abs(Tm))
    K = 4 ** rank
    rad = (big - sml) + (3 * K + 20) * U2 * big
    return Disc(*common._clean(mid, rad))


def curvature(model, k, chart, V, qmax=3, want=('riem', 'E', 'scal')):
    """returns dict: 'riem': [|nabla^q Riem|^2], 'E': [|nabla^j E|^2] (E = Ric - 3 g), 'scal': scalar curvature value,
    'g0': metric value, x0."""
    N = qmax + 2
    g, x0, S = metric(model, k, chart, V, N)
    gi = inv_jet(g)
    Gam = christoffel(g, gi)
    Rdn = riemann(g, gi, Gam)
    out = {'x0': x0}
    Ric = gi.trunc(Rdn.N).mul(Rdn, 'rm,rsmn->sn')
    Ric0 = Ric.val(); gi0 = gi.val()
    sm = np.einsum('snZ,snZ->Z', Ric0.m, gi0.m)
    sb = np.einsum('snZ,snZ->Z', np.abs(Ric0.m) + Ric0.r, np.abs(gi0.m) + gi0.r)
    ss = np.einsum('snZ,snZ->Z', np.abs(Ric0.m), np.abs(gi0.m))
    out['scal'] = Disc(sm, (sb - ss) + 80 * U2 * sb)
    out['Ric'] = Ric
    out['g'] = g; out['gi'] = gi; out['Gam'] = Gam
    if 'riem' in want:
        T = Rdn; L = []
        for q in range(qmax + 1):
            L.append(sqnorm0(T, gi))
            if q < qmax:
                T = nabla(T, Gam.trunc(T.N - 1))
        out['riem'] = L
    if 'E' in want:
        E = Ric - g.trunc(Ric.N).scale(3)
        T = E; L = []
        for j in range(qmax + 1):
            L.append(sqnorm0(T, gi))
            if j < qmax:
                T = nabla(T, Gam.trunc(T.N - 1))
        out['E'] = L
    return out
