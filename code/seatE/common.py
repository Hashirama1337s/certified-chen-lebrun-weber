"""common.py -- seat E: problem data (read-only import of seat R's validated setup) and a rigorous float64 complex-disc
arithmetic used for fast vectorised enclosures.

Disc arithmetic (Disc): a value is a set {z : |z - mid| <= rad}, mid complex128 array, rad float64 array (>= 0).
Rounding model: IEEE binary64, round-to-nearest, unit roundoff u = 2^-53.  Every operation returns mid = fl(op(mids)) and a
radius that bounds (exact op over the input discs) - mid, computed as  (analytic bound) * (1 + SAFE) + SAFE * |magnitude|,
where SAFE = 2^-40 dominates the accumulated relative rounding error of the few float operations used to form the bound and
the mid (each <= ~10 u = 1.1e-15 << 9.1e-13).  Overflow/NaN -> rad = inf (a failed enclosure, never a silent one).
"""
import sys, os
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEATR = os.path.join(HERE, '..', 'seatR')
SEATA4 = os.path.join(HERE, '..', 'seatA4')
for p in (SEATR,):
    if p not in sys.path:
        sys.path.insert(0, p)

SOL84 = os.path.join(SEATR, 'sol_clw_B448_M84.json')
SAFE = 2.0 ** -40
TINY = 1e-290


def problem(prec=640, sol=SOL84, case='clw'):
    import met1_setup
    return met1_setup.problem(case, sol, prec)


def _clean(mid, rad):
    bad = ~np.isfinite(mid) | ~np.isfinite(rad)
    if np.any(bad):
        rad = np.where(bad, np.inf, rad)
        mid = np.where(bad, 0.0, mid)
    return mid, rad


class Disc:
    __slots__ = ('m', 'r')
    __array_priority__ = 1000

    def __init__(self, m, r=None):
        self.m = np.asarray(m, dtype=np.complex128)
        self.r = np.zeros(self.m.shape) if r is None else np.asarray(r, dtype=np.float64)

    @property
    def shape(self):
        return self.m.shape

    def mag(self):
        return np.abs(self.m) + self.r

    def __add__(self, o):
        if not isinstance(o, Disc):
            o = Disc(o)
        m = self.m + o.m
        r = (self.r + o.r) * (1 + SAFE) + SAFE * np.abs(m) + TINY
        return Disc(*_clean(m, r))
    __radd__ = __add__

    def __neg__(self):
        return Disc(-self.m, self.r)

    def __sub__(self, o):
        return self + (-o if isinstance(o, Disc) else Disc(-np.asarray(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Disc):
            o = Disc(o)
        am, bm = np.abs(self.m), np.abs(o.m)
        m = self.m * o.m
        r = (am * o.r + self.r * bm + self.r * o.r) * (1 + SAFE) + SAFE * am * bm + TINY
        return Disc(*_clean(m, r))
    __rmul__ = __mul__

    def recip(self):
        a = np.abs(self.m)
        lo = a * (1 - SAFE) - self.r * (1 + SAFE)          # rigorous lower bound of |z| on the disc
        ok = lo > 0
        with np.errstate(divide='ignore', invalid='ignore'):
            m = np.where(ok, 1 / self.m, 0)
            r = np.where(ok, self.r / (a * lo) * (1 + SAFE) + SAFE / lo, np.inf)
        return Disc(*_clean(m, r))

    def __truediv__(self, o):
        if not isinstance(o, Disc):
            o = Disc(o)
        return self * o.recip()

    def __rtruediv__(self, o):
        return Disc(o) * self.recip()

    def __getitem__(self, k):
        return Disc(self.m[k], self.r[k])

    def upper_abs(self):
        return np.abs(self.m) * (1 + SAFE) + self.r

    def contains(self, z):
        """float test (diagnostic only): |z - mid| <= rad."""
        return np.abs(np.asarray(z) - self.m) <= self.r * (1 + 1e-9) + 1e-300


def disc_from_acb(v):
    """acb/arb ball -> Disc enclosing it (outward)."""
    from flint import acb, arb
    if isinstance(v, arb):
        v = acb(v)
    re, im = v.real, v.imag
    mr, mi = float(re.mid()), float(im.mid())
    rad = float(re.rad()) + float(im.rad()) + abs(float(re.mid()) - mr) + abs(float(im.mid()) - mi)
    # float(arb.mid()) rounds to nearest: |exact mid - float| <= u |float|; add it
    rad = rad * (1 + SAFE) + (abs(mr) + abs(mi)) * 2.0 ** -52 + TINY
    return complex(mr, mi), rad


def disc_sqrt(d):
    """principal sqrt over a disc lying in the open right half-plane (Re mid - rad > 0):
    |sqrt z - sqrt c| = |z - c| / |sqrt z + sqrt c| <= r / Re sqrt c  (both roots have Re > 0)."""
    m = np.asarray(d.m, complex)
    if np.any(m.real * (1 - SAFE) - d.r * (1 + SAFE) <= 0):
        raise ValueError('disc_sqrt: disc not in the right half-plane')
    s = np.sqrt(m)
    rad = d.r / (s.real * (1 - SAFE)) * (1 + SAFE) + SAFE * np.abs(s) + TINY
    return Disc(s, rad)


class ElemMap:
    """seat R / seat W vertex element k: (s, t) in [0,1]^2 -> chart (p, q) = (s sqrt(A(t)), t sqrt(B(s))),
    A(t) = (1 - t^2) al1 + t^2 al2, B(s) = (1 - s^2) be1 + s^2 be2 (al1 = l_a(Mb), al2 = l_a(O), be1 = l_b(Ma), be2 = l_b(O))."""

    def __init__(self, prob, k):
        P = prob.P
        E = P.elem(k)
        a, b = E['a'], E['b']
        self.c = [disc_from_acb(to_arb(v)) for v in (P.ell(a, *E['Mb']), P.ell(a, *P.O), P.ell(b, *E['Ma']), P.ell(b, *P.O))]

    def K(self, i, nb):
        return Disc(np.full(nb, self.c[i][0]), np.full(nb, self.c[i][1]))

    def pq(self, s, t):
        nb = len(s.m)
        one = Disc(np.ones(nb))
        A = (one - t * t) * self.K(0, nb) + t * t * self.K(1, nb)
        B = (one - s * s) * self.K(2, nb) + s * s * self.K(3, nb)
        return s * disc_sqrt(A), t * disc_sqrt(B)


def elem_x(prob, k, s, t):
    """element (s,t) -> x discs, polynomial (no sqrt): l_a = s^2 A(t), l_b = t^2 B(s), x = Ai (l_a - lam_a, l_b - lam_b)."""
    em = ElemMap(prob, k)
    P = prob.P
    a, b = k, (k + 1) % P.d
    Am = [list(P.nu[a]), list(P.nu[b])]
    det = Am[0][0] * Am[1][1] - Am[0][1] * Am[1][0]
    Ai = [[Am[1][1] // det, -Am[0][1] // det], [-Am[1][0] // det, Am[0][0] // det]]
    nb = len(s.m)
    one = Disc(np.ones(nb))
    A = (one - t * t) * em.K(0, nb) + t * t * em.K(1, nb)
    B = (one - s * s) * em.K(2, nb) + s * s * em.K(3, nb)
    la = disc_from_acb(prob.lam[a]); lb = disc_from_acb(prob.lam[b])
    LA = s * s * A - Disc(np.full(nb, la[0]), np.full(nb, la[1]))
    LB = t * t * B - Disc(np.full(nb, lb[0]), np.full(nb, lb[1]))
    return LA * float(Ai[0][0]) + LB * float(Ai[0][1]), LA * float(Ai[1][0]) + LB * float(Ai[1][1])


def to_arb(v):
    """arb data pass through (a* ball route); mp values (legacy seat R setup) are converted at the current precision."""
    from flint import arb
    if isinstance(v, arb):
        return v
    import met1_core as mc
    return mc.A(v)
