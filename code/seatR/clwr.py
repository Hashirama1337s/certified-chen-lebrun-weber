"""clwr.py -- seat R, Boss 5 (BAR-CLWI).  High-precision Donaldson/Abreu solver on toric polygons.

Fixed-point integer arithmetic (B fractional bits) for everything touching the solution; FLINT fmpz_mat for products;
Arb (approximate LU at prec >= 256) for the chord-Newton matrix.  Conventions follow seat W (NOTES-W.md sec. 2-3);
see NOTES-R.md sec. 0 for the design.
"""
import time, json
import numpy as np
import mpmath as mp
from flint import fmpz_mat, arb_mat, arb, ctx

LIU_POLY = [1, 5, 6, -12, -34, -42, 0, 20, 21, -7, -6]


def astar(dps=160):
    mp.mp.dps = dps
    f = lambda a: mp.polyval(LIU_POLY, a)
    r = mp.findroot(f, mp.mpf('1.9577128'))
    eps = mp.mpf(10) ** (-(dps - 20))
    assert f(1) * f(2) < 0 and 1 < r < 2 and f(r - eps) * f(r + eps) < 0
    return r


def page_nu(dps=160):
    mp.mp.dps = dps
    p = lambda v: v**4 + 4 * v**3 - 6 * v**2 + 12 * v - 3
    r = mp.findroot(p, mp.mpf('0.2817'))
    eps = mp.mpf(10) ** (-(dps - 20))
    assert p(r - eps) * p(r + eps) < 0
    return r


class Poly:
    """Delzant polygon l_r = nu_r . x + lam_r, cyclic; vertex k = edge k cap edge k+1 (exact mp data)."""

    def __init__(self, nus, lams, name=''):
        self.nu = [tuple(int(v) for v in n) for n in nus]
        self.lam = [mp.mpf(l) for l in lams]
        self.d = d = len(self.nu)
        self.name = name
        V = []
        for k in range(d):
            a, b = self.nu[k], self.nu[(k + 1) % d]
            det = a[0] * b[1] - a[1] * b[0]
            assert abs(det) == 1, 'not Delzant'
            r1, r2 = -self.lam[k], -self.lam[(k + 1) % d]
            V.append(((r1 * b[1] - a[1] * r2) / det, (a[0] * r2 - b[0] * r1) / det))
        self.V = V
        self.O = (sum(v[0] for v in V) / d, sum(v[1] for v in V) / d)
        self.Mid = [((V[r - 1][0] + V[r][0]) / 2, (V[r - 1][1] + V[r][1]) / 2) for r in range(d)]

    def ell(self, r, x1, x2):
        return self.nu[r][0] * x1 + self.nu[r][1] * x2 + self.lam[r]

    def elem(self, k):
        d = self.d
        a, b = k, (k + 1) % d
        return dict(k=k, a=a, b=b, c=self.V[k], Mb=self.Mid[b], Ma=self.Mid[a], O=self.O)

    def swap_orbits(self):
        """Element orbits under (x1,x2) -> (x2,x1); returns [(rep, multiplicity)]."""
        d = self.d
        ep = []
        for r in range(d):
            sn = (self.nu[r][1], self.nu[r][0])
            cand = [s for s in range(d) if self.nu[s] == sn and abs(self.lam[s] - self.lam[r]) < mp.mpf(10) ** -30]
            assert len(cand) == 1, 'polygon not swap-symmetric'
            ep.append(cand[0])
        vp = []
        for k in range(d):
            e = {ep[k], ep[(k + 1) % d]}
            j = [j for j in range(d) if {j, (j + 1) % d} == e]
            vp.append(j[0])
        seen, out = set(), []
        for k in range(d):
            if k in seen:
                continue
            orb = {k, vp[k]}
            seen |= orb
            out.append((min(orb), len(orb)))
        return out


def pentagon(gamma):
    return Poly([(1, 0), (1, 1), (0, 1), (-1, 0), (0, -1)], [1, gamma, 1, 1, 1], 'pentagon')


def page_trapezoid(gamma):
    """Seat C's coordinates: lam' = (lam0 + nus.sh)/sc, sh = (1/2,1/2), sc = 3/2."""
    nus = [(1, 0), (1, 1), (0, 1), (-1, -1)]
    lam0 = [mp.mpf(1), gamma, mp.mpf(1), mp.mpf(1)]
    lams = [(lam0[r] + (nus[r][0] + nus[r][1]) / mp.mpf(2)) / mp.mpf('1.5') for r in range(4)]
    return Poly(nus, lams, 'trapezoid')


def elem_point(E, sg, ta):
    c, Mb, O, Ma = E['c'], E['Mb'], E['O'], E['Ma']
    w00, w10, w11, w01 = (1 - sg) * (1 - ta), sg * (1 - ta), sg * ta, (1 - sg) * ta
    x1 = w00 * c[0] + w10 * Mb[0] + w11 * O[0] + w01 * Ma[0]
    x2 = w00 * c[1] + w10 * Mb[1] + w11 * O[1] + w01 * Ma[1]
    dxs = [(1 - ta) * (Mb[i] - c[i]) + ta * (O[i] - Ma[i]) for i in range(2)]
    dxt = [(1 - sg) * (Ma[i] - c[i]) + sg * (O[i] - Mb[i]) for i in range(2)]
    J = abs(dxs[0] * dxt[1] - dxs[1] * dxt[0])
    return x1, x2, J


_GL = {}


def gauss01(n):
    """Gauss-Legendre on [0,1] at current mp precision (Newton-polished from numpy)."""
    key = (n, mp.mp.prec)
    if key in _GL:
        return _GL[key]
    z0, _ = np.polynomial.legendre.leggauss(n)
    xs, ws = [], []
    for z in z0:
        x = mp.mpf(z)
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(2) ** (-mp.mp.prec + 8):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append((x + 1) / 2)
        ws.append(1 / ((1 - x * x) * dp * dp))
    _GL[key] = (xs, ws)
    return xs, ws


def polygon_integral(P, F, nq, elems=None):
    """int_P F(x1,x2) dx by element Gauss (mp)."""
    xs, ws = gauss01(nq)
    tot = mp.mpf(0)
    for k in (range(P.d) if elems is None else elems):
        E = P.elem(k)
        for i in range(nq):
            for j in range(nq):
                x1, x2, J = elem_point(E, xs[i], xs[j])
                tot += ws[i] * ws[j] * J * F(x1, x2)
    return tot


def boundary_integral(P, F, nq):
    xs, ws = gauss01(nq)
    tot = mp.mpf(0)
    for r in range(P.d):
        A, B = P.V[r - 1], P.V[r]
        n = P.nu[r]
        L = mp.sqrt((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2) / mp.sqrt(n[0] ** 2 + n[1] ** 2)
        for i in range(nq):
            tot += ws[i] * L * F(A[0] + xs[i] * (B[0] - A[0]), A[1] + xs[i] * (B[1] - A[1]))
    return tot


def S_affine(P):
    """Donaldson's S = a + b x1 + c x2 : int S l = 2 oint l dsigma, l in {1, x1, x2} (exact quadrature)."""
    fs = [lambda x1, x2: 1, lambda x1, x2: x1, lambda x1, x2: x2]
    A = mp.matrix(3, 3)
    rhs = mp.matrix(3, 1)
    for i in range(3):
        for j in range(3):
            A[i, j] = polygon_integral(P, lambda x1, x2: fs[i](x1, x2) * fs[j](x1, x2), 6)
        rhs[i] = 2 * boundary_integral(P, fs[i], 6)
    s = mp.lu_solve(A, rhs)
    return [s[0], s[1], s[2]]


def c_polygon(P, Sabc, LAM=3, nq=80):
    S = lambda x1, x2: Sabc[0] + Sabc[1] * x1 + Sabc[2] * x2
    num = 2 * boundary_integral(P, lambda a, b: S(a, b) ** -2, nq)
    den = 4 * LAM * polygon_integral(P, lambda a, b: S(a, b) ** -4, nq)
    return num / den


# ------------------------------------------------------------------ fixed point helpers
def fix(v, B):
    return int(mp.nint(v * mp.mpf(2) ** B))


def obj(a):
    return np.array(a, dtype=object)


def cheb_tables(X, M, B, kmax):
    """X: object array of fixed-point ints (scale 2^B). Returns list T[k] of object arrays (M+1, n), scale 2^B."""
    n = len(X)
    ONE = 1 << B
    T = [np.empty((M + 1, n), dtype=object) for _ in range(kmax + 1)]
    for k in range(kmax + 1):
        T[k][0] = obj([ONE if k == 0 else 0] * n)
        if M >= 1:
            T[k][1] = obj([X[i] if k == 0 else (ONE if k == 1 else 0) for i in range(n)])
    X2 = 2 * X
    for i in range(1, M):
        for k in range(kmax + 1):
            v = ((X2 * T[k][i]) >> B) - T[k][i - 1]
            if k >= 1:
                v = v + 2 * k * T[k - 1][i]
            T[k][i + 1] = v
    return T


def to_fm(A):
    """object 2D array -> fmpz_mat"""
    A = np.asarray(A)
    return fmpz_mat(A.shape[0], A.shape[1], [int(v) for v in A.ravel()])


def from_fm(F):
    return np.array([[int(v) for v in row] for row in F.tolist()], dtype=object)


def sym_index(M):
    return [(i, j) for i in range(M + 1) for j in range(i + 1) if 2 <= i + j <= M]


class Solver:
    def __init__(self, P, M, B, nq=None, BH=None, log=print):
        self.P, self.M, self.B = P, M, B
        self.BH = BH or B
        self.nq = nq or (M + 24)
        self.log = log
        mp.mp.prec = max(B, self.BH) + 120
        t0 = time.time()
        self.Sabc = S_affine(P)
        xs, ws = gauss01(self.nq)
        X1, X2, W, ELL, SS = [], [], [], [[] for _ in range(P.d)], []
        for (k, mult) in P.swap_orbits():
            E = P.elem(k)
            for i in range(self.nq):
                for j in range(self.nq):
                    x1, x2, J = elem_point(E, xs[i], xs[j])
                    X1.append(fix(x1, B)); X2.append(fix(x2, B)); W.append(fix(ws[i] * ws[j] * J * mult, B))
                    for r in range(P.d):
                        ELL[r].append(fix(P.ell(r, x1, x2), B))
                    SS.append(fix(self.Sabc[0] + self.Sabc[1] * x1 + self.Sabc[2] * x2, B))
        self.X1, self.X2, self.W, self.SS = obj(X1), obj(X2), obj(W), obj(SS)
        self.ELL = [obj(e) for e in ELL]
        self.npts = len(W)
        one2 = 1 << (2 * B)
        # Hess u_G = 1/2 sum nu nu^T / l
        h11 = h12 = h22 = 0
        for r in range(P.d):
            il = one2 // self.ELL[r]       # scale B
            n = P.nu[r]
            h11 = h11 + n[0] * n[0] * il; h12 = h12 + n[0] * n[1] * il; h22 = h22 + n[1] * n[1] * il
        self.HG = (h11 >> 1, h12 >> 1, h22 >> 1)
        self.T1 = cheb_tables(self.X1, M, B, 2)
        self.T2 = cheb_tables(self.X2, M, B, 2)
        self.T1f = [to_fm(self.T1[k].T) for k in range(3)]        # npts x (M+1)
        self.T1fT = [to_fm(self.T1[k]) for k in range(3)]         # (M+1) x npts
        # boundary term bnd[i][j] = oint T_i(x1) T_j(x2) dsigma (exact Gauss)
        nqb = M // 2 + 6
        bx, bw = gauss01(nqb)
        bnd = [[mp.mpf(0)] * (M + 1) for _ in range(M + 1)]
        for r in range(P.d):
            A_, B_ = P.V[r - 1], P.V[r]
            n = P.nu[r]
            L = mp.sqrt((B_[0] - A_[0]) ** 2 + (B_[1] - A_[1]) ** 2) / mp.sqrt(n[0] ** 2 + n[1] ** 2)
            for g in range(nqb):
                y1 = A_[0] + bx[g] * (B_[0] - A_[0]); y2 = A_[1] + bx[g] * (B_[1] - A_[1])
                t1 = [mp.mpf(1), y1]; t2 = [mp.mpf(1), y2]
                for i in range(1, M):
                    t1.append(2 * y1 * t1[i] - t1[i - 1]); t2.append(2 * y2 * t2[i] - t2[i - 1])
                for i in range(M + 1):
                    for j in range(M + 1 - i):
                        bnd[i][j] += bw[g] * L * t1[i] * t2[j]
        self.bnd = [[fix(bnd[i][j], B) if i + j <= M else 0 for j in range(M + 1)] for i in range(M + 1)]
        self.idx = sym_index(M)
        self.nb = len(self.idx)
        self.log(f'[setup] M={M} B={B} nq={self.nq} npts={self.npts} nb={self.nb} secs={time.time() - t0:.1f}')

    # ---------------- coefficients: C full symmetric matrix (object ints, scale 2^B)
    def C_from_sym(self, c):
        M = self.M
        C = np.zeros((M + 1, M + 1), dtype=object)
        C[:] = 0
        for v, (i, j) in zip(c, self.idx):
            C[i, j] = v; C[j, i] = v
        return C

    def sym_from_full(self, G):
        return obj([G[i, j] + G[j, i] if i != j else G[i, i] for (i, j) in self.idx])

    def hess_f(self, C):
        B = self.B
        Cf = to_fm(C)
        out = []
        for (a, b) in ((2, 0), (1, 1), (0, 2)):
            P_ = from_fm(self.T1f[a] * Cf)           # npts x (M+1), scale 2B
            v = (P_ * self.T2[b].T).sum(axis=1)      # scale 3B
            out.append(v >> (2 * B))
        return out

    def U_of(self, C):
        B = self.B
        F20, F11, F02 = self.hess_f(C)
        h11, h12, h22 = self.HG[0] + F20, self.HG[1] + F11, self.HG[2] + F02
        det = (h11 * h22 - h12 * h12) >> B
        assert all(d > 0 for d in det), 'Hess u not positive'
        return (h22 << B) // det, -((h12 << B) // det), (h11 << B) // det

    def grad(self, C):
        B = self.B
        U11, U12, U22 = self.U_of(C)
        W = self.W
        c20 = (W * U11) >> B; c11 = (2 * W * U12) >> B; c02 = (W * U22) >> B; c00 = (W * self.SS) >> B
        Y0 = c02[:, None] * self.T2[2].T + c00[:, None] * self.T2[0].T      # scale 2B
        Y1 = c11[:, None] * self.T2[1].T
        Y2 = c20[:, None] * self.T2[0].T
        G = self.T1fT[0] * to_fm(Y0) + self.T1fT[1] * to_fm(Y1) + self.T1fT[2] * to_fm(Y2)   # scale 3B
        G = from_fm(G)
        Gm = np.empty((self.M + 1, self.M + 1), dtype=object)
        for i in range(self.M + 1):
            for j in range(self.M + 1):
                Gm[i, j] = -(G[i, j] >> (2 * B)) + 2 * self.bnd[i][j]
        return self.sym_from_full(Gm), (U11, U12, U22)

    def hessian(self, U, chunk=6000):
        """H_ab = sum_pts w tr(U Phi_a U Phi_b), symmetric basis; returns arb_mat (scaled to true values)."""
        B, BH = self.B, self.BH
        sh = B - BH
        U11, U12, U22 = [u >> sh if sh > 0 else u for u in U]
        W = self.W >> sh if sh > 0 else self.W
        ii = np.array([i for i, j in self.idx]); jj = np.array([j for i, j in self.idx])
        diag = ii == jj
        H = None
        t0 = time.time()
        for s0 in range(0, self.npts, chunk):
            sl = slice(s0, min(self.npts, s0 + chunk))
            T1 = [(self.T1[k][:, sl] >> sh) if sh > 0 else self.T1[k][:, sl] for k in range(3)]
            T2 = [(self.T2[k][:, sl] >> sh) if sh > 0 else self.T2[k][:, sl] for k in range(3)]
            Phi = []
            for (a, b) in ((2, 0), (1, 1), (0, 2)):
                A = T1[a][ii] * T2[b][jj] + T1[a][jj] * T2[b][ii]       # nb x n, scale 2BH
                A[diag] = T1[a][ii[diag]] * T2[b][jj[diag]]
                Phi.append((A >> BH).T)                                    # n x nb, scale BH
            u11, u12, u22 = U11[sl], U12[sl], U22[sl]
            w = W[sl]
            Mq = [[u11 * u11, 2 * u11 * u12, u12 * u12],
                  [2 * u11 * u12, 2 * (u11 * u22 + u12 * u12), 2 * u12 * u22],
                  [u12 * u12, 2 * u12 * u22, u22 * u22]]               # scale 2BH
            for m in range(3):
                Y = None
                for n_ in range(3):
                    cf = (w * Mq[m][n_]) >> (2 * BH)                       # scale BH
                    t = cf[:, None] * Phi[n_]                                # scale 2BH
                    Y = t if Y is None else Y + t
                Y = Y >> BH
                prod = to_fm(Phi[m].T) * to_fm(Y)                            # scale 2BH
                H = prod if H is None else H + prod
            self.log(f'  [hess] chunk {s0}/{self.npts} secs={time.time() - t0:.1f}')
        Ha = arb_mat(H) * arb(2) ** (-2 * BH)
        return Ha

    def chord_solve(self, C0, maxit=12, tol_bits=None, prec=None, Hinv=None, chunk=4000):
        """chord Newton from C0 (full symmetric object matrix, scale 2^B). Returns C, info (info['Hinv'] = chord inverse)."""
        B = self.B
        ctx.prec = prec or max(256, self.BH + 64)
        t0 = time.time()
        g, U = self.grad(C0)
        cond = None
        if Hinv is None:
            H = self.hessian(U, chunk=chunk)
            self.log(f'  [hess] assembled secs={time.time() - t0:.1f}')
            I = arb_mat(self.nb, self.nb, [1 if i == j else 0 for i in range(self.nb) for j in range(self.nb)])
            Hinv = H.solve(I, algorithm='approx')
            # condition estimate (max-abs norms)
            hn = max(abs(float(H[i, j].mid())) for i in range(self.nb) for j in range(self.nb))
            hin = max(abs(float(Hinv[i, j].mid())) for i in range(self.nb) for j in range(self.nb))
            cond = hn * hin * self.nb
            del H
            self.log(f'[chord] H built+inverted secs={time.time() - t0:.1f}  cond~{cond:.3e} (log2 {np.log2(cond):.1f}, B={B}, BH={self.BH})')
        else:
            self.log(f'[chord] reusing chord inverse (B={B})')
        C = C0.copy()
        hist = []
        scale = arb(2) ** (-B)
        for it in range(maxit):
            gv = arb_mat(self.nb, 1, [int(v) for v in g]) * scale
            d = Hinv * gv
            dec2 = sum((d[i, 0] * gv[i, 0]).mid() for i in range(self.nb))
            dint = []
            for i in range(self.nb):
                m, e = d[i, 0].mid().man_exp()
                m, e = int(m), int(e) + B
                dint.append(-(m << e) if e >= 0 else -(m >> (-e)))
            dec = float(abs(dec2)) ** 0.5
            hist.append(dict(it=it, decrement=dec, gmax=max(abs(int(v)) for v in g) / 2.0 ** B, secs=time.time() - t0))
            self.log(f'  chord {it}: decrement={dec:.3e} gmax={hist[-1]["gmax"]:.3e} secs={time.time() - t0:.1f}')
            if dec < 2.0 ** (-(tol_bits or min(B - 20, 140))):
                break
            if it >= 2 and dec > 0.5 * hist[-2]['decrement']:
                self.log('  chord stagnated (precision floor)')
                C = C + self.C_from_sym(dint)
                g, U = self.grad(C)
                break
            C = C + self.C_from_sym(dint)
            g, U = self.grad(C)
        self.Hinv_last = Hinv
        return C, dict(cond=cond, hist=hist)


def prolong(C, Mto):
    M = C.shape[0] - 1
    out = np.zeros((Mto + 1, Mto + 1), dtype=object)
    out[:] = 0
    m = min(M, Mto)
    for i in range(m + 1):
        for j in range(m + 1 - i):
            out[i, j] = C[i, j]
    return out


def rescale(C, Bfrom, Bto):
    if Bto >= Bfrom:
        return np.vectorize(lambda v: int(v) << (Bto - Bfrom), otypes=[object])(C)
    return np.vectorize(lambda v: int(v) >> (Bfrom - Bto), otypes=[object])(C)


def from_float_npz(path, B):
    d = np.load(path)
    a, M = d['a'], int(d['M'])
    idx = [(i, j) for i in range(M + 1) for j in range(M + 1) if 2 <= i + j <= M]
    C = np.zeros((M + 1, M + 1))
    for v, (i, j) in zip(a, idx):
        C[i, j] = v
    C = (C + C.T) / 2
    return np.vectorize(lambda v: int(round(v * 2.0 ** 60)) << (B - 60) if B >= 60 else int(round(v * 2.0 ** B)), otypes=[object])(C)


def save(path, C, B, meta):
    M = C.shape[0] - 1
    d = dict(meta, B=B, M=M, C=[[str(int(C[i, j])) for j in range(M + 1)] for i in range(M + 1)])
    json.dump(d, open(path, 'w'))


def load(path):
    d = json.load(open(path))
    C = np.array([[int(v) for v in row] for row in d['C']], dtype=object)
    return C, d
