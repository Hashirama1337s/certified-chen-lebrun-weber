"""cm.py -- seat E: bivariate Chebyshev models with arb coefficients (rigorous by construction; no float rounding model).

A model on the square (u, v) in [-1, 1]^2 is  F = P + E,  P a polynomial of degree <= D in each variable, |E| <= rem on the square.
Storage: with u = (z + 1/z)/2, v = (w + 1/w)/2, T_k(u) = (z^k + z^-k)/2, so P is a symmetric Laurent polynomial
P = sum_{|i|,|j| <= D} p_ij z^i w^j (p_ij = p_|i||j|), packed as ONE arb_poly: coefficient of X^{(i + D) + (j + D) W}, W = 4D + 1.
Products are single arb_poly multiplications (Kronecker packing: i-sums lie in [0, 4D] < W, no overlap), which are ball
arithmetic, so every coefficient is an arb ball containing the exact one.  On |z| = |w| = 1: |P| <= ||P|| := sum |p_ij|.
Truncation of a product to degree D moves the dropped coefficients' absolute sum into rem.  Product rule for remainders:
(P + E)(P' + E') = PP' + (terms bounded by ||P|| rem' + ||P'|| rem + rem rem').
Inverse (validated a posteriori): q = float interpolant of 1/P; e = 1 - q F as a model, delta = ||e|| + rem(e) < 1  =>
F does not vanish on the square and 1/F = q/(1 - e) = q + q e/(1 - e), |q e/(1 - e)| <= ||q|| delta/(1 - delta).
Integral over the square: exact on P (int T_k = 2/(1 - k^2), k even) plus 4 rem.
"""
import math
import numpy as np
from flint import arb, arb_poly, ctx

D = 28
W = 4 * D + 1
NEWTON = 2
ZERO = None


def set_degree(d):
    global D, W
    D = d; W = 4 * d + 1


def _up(x):
    """arb upper bound of |x| as an arb (exact point)."""
    return x.abs_upper()


class CM:
    __slots__ = ('p', 'rem', 'cst', 'deg', '_norm')

    def __init__(self, p, rem, cst=None, deg=None):
        self.p, self.rem, self.cst, self.deg = p, rem, cst, (D if deg is None else deg)
        self._norm = None

    # ---- constructors
    @staticmethod
    def const(c):
        c = arb(c)
        L = [arb(0)] * (D + D * W) + [c]
        return CM(arb_poly(L), arb(0), cst=c, deg=0)

    @staticmethod
    def var(which):
        L = [arb(0)] * ((2 * D + 1) + 2 * D * W)
        h = arb(1) / 2
        if which == 0:
            L[(D - 1) + D * W] = h; L[(D + 1) + D * W] = h
        else:
            L[D + (D - 1) * W] = h; L[D + (D + 1) * W] = h
        return CM(arb_poly(L), arb(0), deg=1)

    @staticmethod
    def from_cheb_float(Q):
        """Q[i, j] float Chebyshev coefficients (exact binary64 values)."""
        n1, n2 = Q.shape
        L = [arb(0)] * ((2 * D + 1) + 2 * D * W)
        deg = 0
        for i in range(n1):
            for j in range(n2):
                q = float(Q[i, j])
                if q == 0.0:
                    continue
                c = arb(q) / ((2 if i else 1) * (2 if j else 1))
                for si in ((i, -i) if i else (0,)):
                    for sj in ((j, -j) if j else (0,)):
                        L[(si + D) + (sj + D) * W] = c
                deg = max(deg, i, j)
        return CM(arb_poly(L), arb(0), deg=deg)

    # ---- norms and bounds
    def norm(self):
        if self._norm is None:
            s = arb(0)
            for x in self.p.coeffs():
                s += x.abs_upper()
            self._norm = s.abs_upper()
        return self._norm

    def sup(self):
        """float upper bound of |F| on the square."""
        return float((self.norm() + self.rem).abs_upper())

    def cheb(self):
        """Chebyshev coefficient matrix as arb (list of lists, (D+1) x (D+1))."""
        L = self.p.coeffs()
        n = len(L)
        out = [[arb(0)] * (D + 1) for _ in range(D + 1)]
        for j in range(D + 1):
            for i in range(D + 1):
                k = (i + D) + (j + D) * W
                if k < n:
                    out[i][j] = L[k] * ((2 if i else 1) * (2 if j else 1))
        return out

    def cheb_float(self):
        """(mid matrix, err) : float mids and a float bound of sum |c - mid| (ball radii + float rounding)."""
        C = self.cheb()
        M = np.zeros((D + 1, D + 1)); err = 0.0
        for i in range(D + 1):
            for j in range(D + 1):
                c = C[i][j]
                m = float(c.mid())
                M[i, j] = m
                err += float(c.rad()) + abs(m) * 2.0 ** -52
        return M, err * (1 + 1e-9)

    # ---- arithmetic
    def __add__(self, o):
        if not isinstance(o, CM):
            o = CM.const(o)
        c = (self.cst + o.cst) if (self.cst is not None and o.cst is not None) else None
        return CM(self.p + o.p, (self.rem + o.rem).abs_upper(), cst=c, deg=max(self.deg, o.deg))

    __radd__ = __add__

    def __neg__(self):
        return CM(-self.p, self.rem, cst=None if self.cst is None else -self.cst, deg=self.deg)

    def __sub__(self, o):
        return self + (-o if isinstance(o, CM) else CM.const(-arb(o)))

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, c):
        c = arb(c)
        return CM(self.p * c, (self.rem * c.abs_upper()).abs_upper(), cst=None if self.cst is None else self.cst * c, deg=self.deg)

    def __mul__(self, o):
        if not isinstance(o, CM):
            return self.scale(o)
        if self.cst is not None:
            return o.scale(self.cst)
        if o.cst is not None:
            return self.scale(o.cst)
        prod = self.p * o.p
        rem = self.norm() * o.rem + o.norm() * self.rem + self.rem * o.rem
        dsum = self.deg + o.deg
        if dsum <= D:
            p = prod.right_shift(D + D * W)
            return CM(p, rem.abs_upper(), deg=dsum)
        L = prod.coeffs()
        n = len(L)
        need = (4 * D + 1) * W
        if n < need:
            L = L + [arb(0)] * (need - n)
        keep = []
        tail = arb(0)
        pad = [arb(0)] * (W - (2 * D + 1))
        lo, hi = 2 * D - dsum, 2 * D + dsum          # possibly nonzero index range in each variable
        for J in range(max(lo, 0), min(hi, 4 * D) + 1):
            row = L[J * W:J * W + 4 * D + 1]
            if D <= J <= 3 * D:
                for x in row[max(lo, 0):D]:
                    tail += x.abs_upper()
                for x in row[3 * D + 1:min(hi, 4 * D) + 1]:
                    tail += x.abs_upper()
            else:
                for x in row[max(lo, 0):min(hi, 4 * D) + 1]:
                    tail += x.abs_upper()
        for J in range(D, 3 * D + 1):
            keep += L[J * W + D:J * W + 3 * D + 1]
            keep += pad
        return CM(arb_poly(keep), (rem + tail).abs_upper(), deg=D)

    __rmul__ = __mul__

    def inv(self):
        """validated reciprocal (see module docstring)."""
        if self.cst is not None:
            return CM.const(1 / self.cst)
        C, _ = self.cheb_float()
        k = np.arange(D + 1)
        T = np.cos(np.outer(k, k) * np.pi / D)          # T[g, k] = T_k(u_g), u_g = cos(g pi/D)
        V = T @ C @ T.T
        A = (2.0 / D) * np.cos(np.outer(k, k) * np.pi / D)
        A[:, 0] *= 0.5; A[:, D] *= 0.5; A[0, :] *= 0.5; A[D, :] *= 0.5
        Q = A @ (1.0 / V) @ A.T
        q = CM.from_cheb_float(Q)
        # Newton refinement of the CANDIDATE only (remainders discarded: q is not trusted, the validation below is the proof)
        P_ = CM(self.p, arb(0), deg=self.deg)
        for _ in range(NEWTON):
            q = q * (CM.const(2) - P_ * q)
            q = CM(q.p, arb(0), deg=q.deg)
        e = CM.const(1) - q * self
        delta = (e.norm() + e.rem).abs_upper()
        if not (delta < arb(1) / 2):
            raise ValueError('CM.inv: validation failed, delta = %s' % delta)
        return CM(q.p, (q.norm() * delta / (1 - delta)).abs_upper(), deg=q.deg)

    def integral(self):
        """arb enclosure-free pair: (integral of the polynomial part as arb ball, 4 rem upper)."""
        C = self.cheb()
        tot = arb(0)
        for i in range(0, D + 1, 2):
            for j in range(0, D + 1, 2):
                tot += C[i][j] * (arb(2) / (1 - i * i)) * (arb(2) / (1 - j * j))
        return tot, (4 * self.rem).abs_upper()

    def sup_lower(self, m=257):
        """rigorous LOWER bound of sup |F| on the square: max over grid points of |P(u_g, v_g)| minus (coefficient conversion
        error + float evaluation slack + rem).  Each grid point is an actual point of the square, so no variation term."""
        C, err = self.cheb_float()
        k = np.arange(D + 1)
        g = np.arange(m)
        Tu = np.cos((np.outer(g, k) % (2 * (m - 1))) * np.pi / (m - 1))     # reduced argument (amend2)
        G = Tu @ C @ Tu.T
        slack = (err + 1e-12 * np.abs(C).sum() + float(self.rem.abs_upper())) * (1 + 1e-9)
        return float(np.abs(G).max()) - slack

    def bounds(self, m_max=1500, tol_rel=1e-3):
        """rigorous (inf, sup) of F on the square: float grid evaluation of P at u_g = cos(g pi/(m-1)) + Markov variation
        (|T_k'| <= k^2) + coefficient conversion error + float evaluation slack (1e-12 ||P||) + rem."""
        C, err = self.cheb_float()
        k = np.arange(D + 1)
        Dk = (np.abs(C) * (k[:, None] ** 2 + k[None, :] ** 2)).sum()
        tot = np.abs(C).sum()
        tol = max(tol_rel * tot, 1e-300)
        m = int(min(m_max, max(64, math.ceil(math.pi * Dk / (2 * tol)) + 1)))
        g = np.arange(m)
        Tu = np.cos((np.outer(g, k) % (2 * (m - 1))) * np.pi / (m - 1))     # reduced argument (amend2)
        G = Tu @ C @ Tu.T
        slack = Dk * math.pi / (2 * (m - 1)) + err + 1e-12 * tot + float(self.rem.abs_upper())
        slack *= (1 + 1e-9)
        return float(G.min()) - slack, float(G.max()) + slack


class Jet2:
    """truncated Taylor jet in (dy1, dy2) of order K with CM coefficients: c[(i, j)] = coefficient of dy1^i dy2^j."""
    __slots__ = ('c', 'K')

    def __init__(self, c, K):
        self.c, self.K = c, K

    @staticmethod
    def const(cm_, K):
        return Jet2({(0, 0): cm_}, K)

    def get(self, m):
        return self.c.get(m)

    def __add__(self, o):
        if not isinstance(o, Jet2):
            o = Jet2.const(o if isinstance(o, CM) else CM.const(o), self.K)
        K = min(self.K, o.K)
        c = {m: v for m, v in self.c.items() if m[0] + m[1] <= K}
        for m, v in o.c.items():
            if m[0] + m[1] <= K:
                c[m] = (c[m] + v) if m in c else v
        return Jet2(c, K)

    __radd__ = __add__

    def __neg__(self):
        return Jet2({m: -v for m, v in self.c.items()}, self.K)

    def __sub__(self, o):
        return self + (-o if isinstance(o, Jet2) else -(o if isinstance(o, CM) else CM.const(o)))

    def scale(self, s):
        return Jet2({m: v * s for m, v in self.c.items()}, self.K)

    def __mul__(self, o):
        if not isinstance(o, Jet2):
            return self.scale(o)
        K = min(self.K, o.K)
        c = {}
        for (i1, j1), a in self.c.items():
            for (i2, j2), b in o.c.items():
                if i1 + j1 + i2 + j2 > K:
                    continue
                t = a * b
                m = (i1 + i2, j1 + j2)
                c[m] = (c[m] + t) if m in c else t
        return Jet2(c, K)

    __rmul__ = __mul__

    def d(self, k):
        c = {}
        for (i, j), v in self.c.items():
            if (i if k == 0 else j) == 0:
                continue
            c[(i - 1, j) if k == 0 else (i, j - 1)] = v.scale(i if k == 0 else j)
        return Jet2(c, self.K - 1)

    def val(self):
        return self.c.get((0, 0), CM.const(0))

    def recip(self):
        """1/J: y0 = inv(c00) (validated), y_m = -y0 sum_{m' != 0} c_{m'} y_{m - m'}."""
        y0 = self.c[(0, 0)].inv()
        y = {(0, 0): y0}
        for n in range(1, self.K + 1):
            for i in range(n, -1, -1):
                j = n - i
                acc = None
                for (i1, j1), a in self.c.items():
                    if (i1, j1) == (0, 0) or i1 > i or j1 > j:
                        continue
                    yy = y.get((i - i1, j - j1))
                    if yy is None:
                        continue
                    t = a * yy
                    acc = t if acc is None else acc + t
                if acc is not None:
                    y[(i, j)] = -(acc * y0)
        return Jet2(y, self.K)
