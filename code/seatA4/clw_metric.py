"""clw_metric.py -- seat W, Boss 2 (CLW).  One self-contained module for the BACKGROUND.

Toric extremal Kahler metric on a Delzant polygon by Donaldson's functional / Abreu's equation,
symplectic potential u = u_G + f, u_G = 1/2 sum_r l_r log l_r (Guillemin), f = sum a_ij T_i(x1) T_j(x2)
(Chebyshev, 2 <= i+j <= M, affine part excluded: it is the kernel).  Einstein metric g = c S^-2 g_K.
All conventions are documented in NOTES-W.md.  Plain numpy / mpmath; no symbolic algebra.

Geometry is evaluated per VERTEX ELEMENT k (vertex between edges a=k and b=k+1) in the coordinates
(p, q) = (sqrt l_a, sqrt l_b), angles (theta_a, theta_b) of the vertex lattice basis (nu_a, nu_b):
phi = theta_a nu_a + theta_b nu_b.  In these coordinates every metric coefficient is smooth up to both rods.
Element parameter (s, t) in [0,1]^2:  x = bilinear corners (c, M_b, O, M_a) at (sigma, tau) = (s^2, t^2),
so l_a = s^2 A(t), l_b = t^2 B(s) exactly (A, B > 0 affine in t^2, s^2).
"""
import numpy as np
import mpmath as mp
from numpy.polynomial import chebyshev as Ch

# ------------------------------------------------------------------ exact class
LIU_POLY = [1, 5, 6, -12, -34, -42, 0, 20, 21, -7, -6]   # a^10 + 5a^9 + ... - 6 (Liu SM S4)


def astar(dps=40):
    """Root in (1,2) of Liu's degree-10 polynomial (LeBrun's Bach-flat class). gamma = 2/a*."""
    mp.mp.dps = dps
    f = lambda a: mp.polyval(LIU_POLY, a)
    r = mp.findroot(f, mp.mpf('1.9577128'))
    assert f(1) * f(2) < 0 and 1 < r < 2
    return r


def gamma_bach(dps=40):
    return 2 / astar(dps)


# ------------------------------------------------------------------ polygons
class Polygon:
    """Edges l_r(x) = nu_r . x + lam_r >= 0 in cyclic (counter-clockwise) order; vertex k = edge k cap edge k+1."""

    def __init__(self, nus, lams, name=''):
        self.nu = np.array(nus, float)
        self.lam = np.array(lams, float)
        self.d = len(lams)
        self.name = name
        d = self.d
        V = []
        for k in range(d):
            A = np.array([self.nu[k], self.nu[(k + 1) % d]])
            assert abs(round(np.linalg.det(A)) - np.linalg.det(A)) < 1e-12 and abs(np.linalg.det(A)) == 1, 'not Delzant'
            V.append(np.linalg.solve(A, -np.array([self.lam[k], self.lam[(k + 1) % d]])))
        self.V = np.array(V)
        self.O = self.V.mean(axis=0)              # centre used for the element interfaces
        # edge r runs from vertex r-1 to vertex r
        self.Mid = np.array([(self.V[r - 1] + self.V[r]) / 2 for r in range(d)])

    def ell(self, x):   # x: (..., 2) -> (..., d)
        return x @ self.nu.T + self.lam

    def elem(self, k):
        """Data of vertex element k: edges a=k, b=k+1."""
        d = self.d
        a, b = k, (k + 1) % d
        A = np.array([self.nu[a], self.nu[b]])          # y = A x + lab
        lab = np.array([self.lam[a], self.lam[b]])
        Ainv = np.linalg.inv(A)
        c, Mb, Ma, O = self.V[k], self.Mid[b], self.Mid[a], self.O
        la = lambda x: self.nu[a] @ x + self.lam[a]
        lb = lambda x: self.nu[b] @ x + self.lam[b]
        return dict(k=k, a=a, b=b, A=A, Ainv=Ainv, lab=lab, c=c, Mb=Mb, Ma=Ma, O=O,
                    al1=la(Mb), al2=la(O), be1=lb(Ma), be2=lb(O),
                    others=[r for r in range(d) if r not in (a, b)])


def pentagon(gamma):
    """CLW pentagon (Liu SM S4): l = 1+x1, gamma+x1+x2, 1+x2, 1-x1, 1-x2."""
    return Polygon([(1, 0), (1, 1), (0, 1), (-1, 0), (0, -1)], [1, gamma, 1, 1, 1], 'pentagon')


def square():
    return Polygon([(1, 0), (0, 1), (-1, 0), (0, -1)], [1, 1, 1, 1], 'square')


def triangle():
    return Polygon([(1, 0), (0, 1), (-1, -1)], [1, 1, 1], 'triangle')


# ------------------------------------------------------------------ element maps
def elem_x(E, s, t):
    """x(s,t) for element E, plus (p,q) and Jacobian d(p,q)/d(s,t)."""
    sg, ta = s * s, t * t
    x = ((1 - sg) * (1 - ta))[..., None] * E['c'] + (sg * (1 - ta))[..., None] * E['Mb'] \
        + (sg * ta)[..., None] * E['O'] + ((1 - sg) * ta)[..., None] * E['Ma']
    Aa = (1 - ta) * E['al1'] + ta * E['al2']
    Bb = (1 - sg) * E['be1'] + sg * E['be2']
    p, q = s * np.sqrt(Aa), t * np.sqrt(Bb)
    ps, pt = np.sqrt(Aa), s * t * (E['al2'] - E['al1']) / np.sqrt(Aa)
    qt, qs = np.sqrt(Bb), t * s * (E['be2'] - E['be1']) / np.sqrt(Bb)
    J = np.array([[ps, pt], [qs, qt]])               # d(p,q)/d(s,t)
    return x, p, q, J


def elem_dxdst(E, s, t):
    """|det dx/d(s,t)| (for integrals over the polygon in x)."""
    sg, ta = s * s, t * t
    dxds = 2 * s[..., None] * (-(1 - ta)[..., None] * E['c'] + (1 - ta)[..., None] * E['Mb'] + ta[..., None] * E['O'] - ta[..., None] * E['Ma'])
    dxdt = 2 * t[..., None] * (-(1 - sg)[..., None] * E['c'] - sg[..., None] * E['Mb'] + sg[..., None] * E['O'] + (1 - sg)[..., None] * E['Ma'])
    return np.abs(dxds[..., 0] * dxdt[..., 1] - dxds[..., 1] * dxdt[..., 0])


# ------------------------------------------------------------------ Chebyshev tools
def basis_index(M):
    return [(i, j) for i in range(M + 1) for j in range(M + 1) if 2 <= i + j <= M]


def cheb_derivs(x, deg, kmax=4):
    """T[k, i, :] = T_i^{(k)}(x), k = 0..kmax."""
    out = np.zeros((kmax + 1, deg + 1) + np.shape(x))
    for i in range(deg + 1):
        e = np.zeros(i + 1); e[i] = 1
        for k in range(kmax + 1):
            out[k, i] = Ch.chebval(x, Ch.chebder(e, k)) if k <= i else 0.0
    return out


def f_derivs(coef, x, kmax=4):
    """All partial derivatives d^(a1,a2) f, a1+a2 <= kmax, of f = sum coef[i,j] T_i(x1) T_j(x2) at x (...,2)."""
    deg = coef.shape[0] - 1
    T1 = cheb_derivs(x[..., 0], deg, kmax)
    T2 = cheb_derivs(x[..., 1], deg, kmax)
    D = {}
    for a1 in range(kmax + 1):
        for a2 in range(kmax + 1 - a1):
            D[(a1, a2)] = np.einsum('ij,i...,j...->...', coef, T1[a1], T2[a2])
    return D


def coef_matrix(avec, M):
    C = np.zeros((M + 1, M + 1))
    for c, (i, j) in zip(avec, basis_index(M)):
        C[i, j] = c
    return C


# ------------------------------------------------------------------ jets in (p,q): [v, d_p, d_q, d_pp, d_pq, d_qq]
class Jet:
    __slots__ = ('c',)

    def __init__(self, c):
        self.c = c

    @staticmethod
    def const(v, shape):
        c = np.zeros((6,) + shape); c[0] = v
        return Jet(c)

    def __add__(self, o):
        return Jet(self.c + (o.c if isinstance(o, Jet) else np.concatenate([[o * np.ones_like(self.c[0])], np.zeros_like(self.c[1:])])))
    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.c)

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.c * o)
        a, b = self.c, o.c
        return Jet(np.array([a[0] * b[0], a[1] * b[0] + a[0] * b[1], a[2] * b[0] + a[0] * b[2],
                             a[3] * b[0] + 2 * a[1] * b[1] + a[0] * b[3],
                             a[4] * b[0] + a[1] * b[2] + a[2] * b[1] + a[0] * b[4],
                             a[5] * b[0] + 2 * a[2] * b[2] + a[0] * b[5]]))
    __rmul__ = __mul__

    def fn(self, f0, f1, f2):
        a = self.c
        return Jet(np.array([f0, f1 * a[1], f1 * a[2], f2 * a[1] ** 2 + f1 * a[3], f2 * a[1] * a[2] + f1 * a[4], f2 * a[2] ** 2 + f1 * a[5]]))

    def recip(self):
        v = self.c[0]
        return self.fn(1 / v, -1 / v ** 2, 2 / v ** 3)

    def __truediv__(self, o):
        return self * (o.recip() if isinstance(o, Jet) else 1.0 / o)

    def sqrt(self):
        v = self.c[0]; r = np.sqrt(v)
        return self.fn(r, 0.5 / r, -0.25 / (r * v))

    def pow(self, e):
        v = self.c[0]
        return self.fn(v ** e, e * v ** (e - 1), e * (e - 1) * v ** (e - 2))


def compose(F, gF, hF, X):
    """Jet of F(x(p,q)): F value, gF[i] gradient, hF[i][j] Hessian (arrays), X = [Jet x1, Jet x2]."""
    c = np.zeros((6,) + np.shape(F))
    c[0] = F
    xc = [X[0].c, X[1].c]
    for i in range(2):
        c[1] += gF[i] * xc[i][1]; c[2] += gF[i] * xc[i][2]
        c[3] += gF[i] * xc[i][3]; c[4] += gF[i] * xc[i][4]; c[5] += gF[i] * xc[i][5]
        for j in range(2):
            c[3] += hF[i][j] * xc[i][1] * xc[j][1]
            c[4] += hF[i][j] * xc[i][1] * xc[j][2]
            c[5] += hF[i][j] * xc[i][2] * xc[j][2]
    return Jet(c)


# ------------------------------------------------------------------ the model
class Model:
    """Toric Kahler potential u = u_G + f on polygon P; S affine (Donaldson); conformal power w: g = c S^w g_K."""

    def __init__(self, P, coef=None, M=0, Sabc=None, cpow=-2, cscale=1.0):
        self.P, self.M = P, M
        self.coef = coef if coef is not None else np.zeros((max(M, 1) + 1, max(M, 1) + 1))
        self.Sabc = Sabc if Sabc is not None else S_affine(P)
        self.cpow, self.cscale = cpow, cscale

    # ---- R_y = Hess_y u - D  (smooth part, vertex coordinates) as value only, at x
    def Ry_value(self, E, x, with_f=True):
        P = self.P
        Ai = E['Ainv']
        ell = P.ell(x)
        R = np.zeros((2, 2) + x.shape[:-1])
        for r in E['others']:
            n = Ai.T @ P.nu[r]
            R += 0.5 * np.outer(n, n)[..., None] / ell[..., r]
        if with_f:
            D = f_derivs(self.coef, x, 2)
            Hf = np.array([[D[(2, 0)], D[(1, 1)]], [D[(1, 1)], D[(0, 2)]]])
            R += np.einsum('ai,ij...,jb->ab...', Ai.T, Hf, Ai)
        return R

    def Ux(self, E, x, with_f=True):
        """u^{ij} (inverse Hessian in x coords) at points x of element E, via the vertex form (accurate at rods)."""
        y = x @ E['A'].T + E['lab']
        la, lb = y[..., 0], y[..., 1]
        R = self.Ry_value(E, x, with_f)
        m11 = 1 + 2 * la * R[0, 0]; m12 = 2 * la * R[0, 1]; m21 = 2 * lb * R[0, 1]; m22 = 1 + 2 * lb * R[1, 1]
        det = m11 * m22 - m12 * m21
        G = np.array([[m22 * 2 * la, -m12 * 2 * lb], [-m21 * 2 * la, m11 * 2 * lb]]) / det   # (I+E R)^-1 E
        G = 0.5 * (G + G.transpose(1, 0, *range(2, G.ndim)))
        Ai = E['Ainv']
        return np.einsum('ia,ab...,jb->ij...', Ai, G, Ai)

    def S_of(self, x):
        a, b, c = self.Sabc
        return a + b * x[..., 0] + c * x[..., 1]

    # ---- the 4D metric as jets in (p,q) at element points
    def metric_jets(self, E, p, q):
        """Return list-of-lists of Jets g[mu][nu], coordinates (p, q, theta_a, theta_b)."""
        P = self.P
        shp = np.shape(p)
        pj = Jet(np.array([p, np.ones(shp), np.zeros(shp), np.zeros(shp), np.zeros(shp), np.zeros(shp)]))
        qj = Jet(np.array([q, np.zeros(shp), np.ones(shp), np.zeros(shp), np.zeros(shp), np.zeros(shp)]))
        p2, q2 = pj * pj, qj * qj
        Ai, lab = E['Ainv'], E['lab']
        X = [ (p2 - lab[0]) * Ai[i, 0] + (q2 - lab[1]) * Ai[i, 1] for i in range(2)]
        xv = np.stack([X[0].c[0], X[1].c[0]], -1)
        # R_y jets
        R = [[Jet.const(0.0, shp) for _ in range(2)] for _ in range(2)]
        for r in E['others']:
            n = Ai.T @ P.nu[r]
            lr = X[0] * P.nu[r][0] + X[1] * P.nu[r][1] + P.lam[r]
            il = lr.recip()
            for i in range(2):
                for j in range(2):
                    R[i][j] = R[i][j] + il * (0.5 * n[i] * n[j])
        if np.any(self.coef):
            D = f_derivs(self.coef, xv, 4)
            def H(i, j):       # jet of f_{ij}(x(p,q))
                al = [0, 0]; al[i] += 1; al[j] += 1
                def d(*e):
                    b = list(al)
                    for k in e: b[k] += 1
                    return D[tuple(b)]
                return compose(d(), [d(0), d(1)], [[d(0, 0), d(0, 1)], [d(1, 0), d(1, 1)]], X)
            Hf = [[H(0, 0), H(0, 1)], [H(0, 1), H(1, 1)]]
            for a in range(2):
                for b in range(2):
                    acc = Jet.const(0.0, shp)
                    for i in range(2):
                        for j in range(2):
                            co = Ai[i, a] * Ai[j, b]
                            if co != 0:
                                acc = acc + Hf[i][j] * co
                    R[a][b] = R[a][b] + acc
        # Kahler metric blocks
        gpp = 2 + 4 * p2 * R[0][0]
        gpq = 4 * pj * qj * R[0][1]
        gqq = 2 + 4 * q2 * R[1][1]
        m11 = 1 + 2 * p2 * R[0][0]; m12 = 2 * p2 * R[0][1]; m21 = 2 * q2 * R[0][1]; m22 = 1 + 2 * q2 * R[1][1]
        idet = (m11 * m22 - m12 * m21).recip()
        Gaa = m22 * 2 * p2 * idet
        Gab = -(m12 * 2 * q2) * idet
        Gbb = m11 * 2 * q2 * idet
        # conformal factor c S^w
        a0, b0, c0 = self.Sabc
        S = X[0] * b0 + X[1] * c0 + a0
        cf = S.pow(self.cpow) * self.cscale
        Z = Jet.const(0.0, shp)
        g = [[cf * gpp, cf * gpq, Z, Z], [cf * gpq, cf * gqq, Z, Z], [Z, Z, cf * Gaa, cf * Gab], [Z, Z, cf * Gab, cf * Gbb]]
        return g, S, xv


def S_affine(P, nq=40):
    """Donaldson's affine S: int S l dx = 2 oint l dsigma for affine l (dsigma = Lebesgue / |nu_r|)."""
    pts, w = polygon_quad(P, nq)
    L = np.stack([np.ones(len(w)), pts[:, 0], pts[:, 1]], 1)
    A = (L * w[:, None]).T @ L
    bpts, bw = boundary_quad(P, nq)
    Lb = np.stack([np.ones(len(bw)), bpts[:, 0], bpts[:, 1]], 1)
    rhs = 2 * (Lb * bw[:, None]).sum(0)
    return np.linalg.solve(A, rhs)


def gauss01(n):
    z, w = np.polynomial.legendre.leggauss(n)
    return (z + 1) / 2, w / 2


def polygon_quad(P, nq, return_elem=False):
    s1, w1 = gauss01(nq)
    S_, T_ = np.meshgrid(s1, s1, indexing='ij'); W_ = np.outer(w1, w1)
    pts, ws, ks = [], [], []
    for k in range(P.d):
        E = P.elem(k)
        x, p, q, J = elem_x(E, S_.ravel(), T_.ravel())
        pts.append(x); ws.append(W_.ravel() * elem_dxdst(E, S_.ravel(), T_.ravel())); ks.append(np.full(len(x), k))
    if return_elem:
        return np.concatenate(pts), np.concatenate(ws), np.concatenate(ks)
    return np.concatenate(pts), np.concatenate(ws)


def boundary_quad(P, nq):
    s1, w1 = gauss01(nq)
    pts, ws = [], []
    for r in range(P.d):
        A, B = P.V[r - 1], P.V[r]
        L = np.linalg.norm(B - A) / np.linalg.norm(P.nu[r])
        pts.append(A[None] + s1[:, None] * (B - A)[None]); ws.append(w1 * L)
    return np.concatenate(pts), np.concatenate(ws)


# ------------------------------------------------------------------ Donaldson / Abreu solver
class AbreuSolver:
    def __init__(self, P, M, nq=None):
        self.P, self.M = P, M
        self.idx = basis_index(M)
        nq = nq or (2 * M + 16)     # x is quadratic in (s,t): degree doubles; M+14 under-integrates (stall at 4e-4)
        self.nq = nq
        self.pts, self.w, self.kel = polygon_quad(P, nq, return_elem=True)
        self.Sabc = S_affine(P)
        S = self.Sabc[0] + self.Sabc[1] * self.pts[:, 0] + self.Sabc[2] * self.pts[:, 1]
        T1 = cheb_derivs(self.pts[:, 0], M, 2); T2 = cheb_derivs(self.pts[:, 1], M, 2)
        nb = len(self.idx)
        self.Phi = np.zeros((3, nb, len(self.w)))     # (f11, f12, f22) of each basis function
        self.phi0 = np.zeros((nb, len(self.w)))
        for n, (i, j) in enumerate(self.idx):
            self.Phi[0, n] = T1[2, i] * T2[0, j]; self.Phi[1, n] = T1[1, i] * T2[1, j]; self.Phi[2, n] = T1[0, i] * T2[2, j]
            self.phi0[n] = T1[0, i] * T2[0, j]
        bp, bw = boundary_quad(P, nq)
        B1 = cheb_derivs(bp[:, 0], M, 0)[0]; B2 = cheb_derivs(bp[:, 1], M, 0)[0]
        bvec = np.array([(B1[i] * B2[j] * bw).sum() for (i, j) in self.idx])
        self.lin = 2 * bvec - self.phi0 @ (S * self.w)      # linear part of F
        self.UG = self._U(np.zeros(nb), with_f=False)

    def _U(self, a, with_f=True):
        mod = Model(self.P, coef_matrix(a, self.M), self.M, self.Sabc)
        U = np.zeros((2, 2, len(self.w)))
        for k in range(self.P.d):
            m = self.kel == k
            U[:, :, m] = mod.Ux(self.P.elem(k), self.pts[m], with_f)
        return U

    def F_grad_hess(self, a):
        U = self._U(a)
        Hf = np.einsum('knq,n->kq', self.Phi, a)
        # F~(a) = -int log det(I + U_G Hf) + lin.a
        H2 = np.array([[Hf[0], Hf[1]], [Hf[1], Hf[2]]])
        UG = self.UG
        Mx = np.eye(2)[:, :, None] + np.einsum('ijq,jkq->ikq', UG, H2)
        det = Mx[0, 0] * Mx[1, 1] - Mx[0, 1] * Mx[1, 0]
        F = (-np.log(det) * self.w).sum() + self.lin @ a if np.all(det > 0) else np.inf
        trUPhi = U[0, 0] * self.Phi[0] + 2 * U[0, 1] * self.Phi[1] + U[1, 1] * self.Phi[2]
        grad = -trUPhi @ self.w + self.lin
        # Hessian: int tr(U Phi_a U Phi_b) = sum_kl Phi_k M_kl Phi_l
        U11, U12, U22 = U[0, 0], U[0, 1], U[1, 1]
        Mq = np.array([[U11 * U11, 2 * U11 * U12, U12 * U12],
                       [2 * U11 * U12, 2 * (U11 * U22 + U12 * U12), 2 * U12 * U22],
                       [U12 * U12, 2 * U12 * U22, U22 * U22]])
        Y = np.einsum('klq,lnq->knq', Mq, self.Phi) * self.w
        Hs = np.einsum('kmq,knq->mn', self.Phi, Y)
        return F, grad, Hs, det.min() if np.all(np.isfinite(det)) else -1

    def solve(self, a0=None, tol=1e-13, maxit=60, verbose=True):
        a = np.zeros(len(self.idx)) if a0 is None else a0.copy()
        hist = []
        for it in range(maxit):
            F, g, H, dmin = self.F_grad_hess(a)
            lam_, Q_ = np.linalg.eigh(H)
            kp = lam_ > 1e-13 * lam_.max()
            gn = np.abs(Q_[:, kp] @ (Q_[:, kp].T @ g)).max()      # gradient on the resolved subspace
            hist.append((it, F, gn, np.abs(g).max()))
            if verbose:
                print(f'  newton {it:2d}  F={F:.15f}  |grad|={gn:.3e}', flush=True)
            if gn < tol:
                break
            # truncated-spectrum Newton step: the Chebyshev basis of the square is nearly dependent on the pentagon
            # at high M (cond(H) ~ 1e16 at M = 26 made plain Newton fail); drop directions below 1e-13 of the top.
            lam, Q = np.linalg.eigh(H)
            keep = lam > 1e-13 * lam.max()
            step = -Q[:, keep] @ ((Q[:, keep].T @ g) / lam[keep])
            alpha = 1.0
            while alpha > 1e-6:
                Fn = self.F_grad_hess_F(a + alpha * step)
                if Fn <= F + 1e-4 * alpha * (g @ step) or (gn < 1e-8 and np.isfinite(Fn)):
                    break
                alpha *= 0.5
            a = a + alpha * step
        return a, hist

    def F_grad_hess_F(self, a):
        Hf = np.einsum('knq,n->kq', self.Phi, a)
        H2 = np.array([[Hf[0], Hf[1]], [Hf[1], Hf[2]]])
        Mx = np.eye(2)[:, :, None] + np.einsum('ijq,jkq->ikq', self.UG, H2)
        det = Mx[0, 0] * Mx[1, 1] - Mx[0, 1] * Mx[1, 0]
        if not np.all(det > 0):
            return np.inf
        return (-np.log(det) * self.w).sum() + self.lin @ a


def prolong(a, Mfrom, Mto):
    C = coef_matrix(a, Mfrom)
    out = []
    for (i, j) in basis_index(Mto):
        out.append(C[i, j] if i <= Mfrom and j <= Mfrom else 0.0)
    return np.array(out)


# ------------------------------------------------------------------ 4D geometry at element points
def geometry(mod, k, s, t):
    """Everything the spectral code needs at element-k points (s,t) (1D arrays of equal length).
    Coordinates (p, q, theta_a, theta_b).  Returns dict with
      x, p, q, J (d(p,q)/d(s,t)), sqrtg, e[A,mu] (orthonormal frame), Gam[C,A,D] = g(nabla_{e_C} e_A, e_D),
      RF[A,B,C,D] frame Riemann (MTW sign: sphere positive R_1212), Ric[A,B], scal, S (affine), E (element data)."""
    E = mod.P.elem(k)
    x, p, q, J = elem_x(E, s, t)
    gJ, Sj, xv = mod.metric_jets(E, p, q)
    n = len(p)
    G = np.zeros((6, 4, 4, n))
    for m in range(4):
        for nn in range(4):
            G[:, m, nn] = gJ[m][nn].c
    g = G[0]; dg = np.zeros((4, 4, 4, n)); dg[0] = G[1]; dg[1] = G[2]          # dg[alpha, mu, nu]
    ddg = np.zeros((4, 4, 4, 4, n)); ddg[0, 0] = G[3]; ddg[0, 1] = G[4]; ddg[1, 0] = G[4]; ddg[1, 1] = G[5]
    gm = np.moveaxis(g, -1, 0)
    ginv = np.moveaxis(np.linalg.inv(gm), 0, -1)
    # Christoffel Gamma^l_{mn} and derivatives
    Cc = np.einsum('mkn...->kmn...', dg) + np.einsum('nkm...->kmn...', dg) - dg      # C[k,m,n] = d_m g_kn + d_n g_km - d_k g_mn
    Gam = 0.5 * np.einsum('lk...,kmn...->lmn...', ginv, Cc)
    dginv = -np.einsum('ia...,cab...,bj...->cij...', ginv, dg, ginv)                     # d_c g^{ij}
    dC = np.einsum('cmkn...->ckmn...', ddg) + np.einsum('cnkm...->ckmn...', ddg) - np.einsum('ckmn...->ckmn...', ddg)
    # dC[c,k,m,n] = d_c d_m g_kn + d_c d_n g_km - d_c d_k g_mn
    dGam = 0.5 * (np.einsum('clk...,kmn...->clmn...', dginv, Cc) + np.einsum('lk...,ckmn...->clmn...', ginv, dC))
    # Riemann R^r_{s m n} = d_m Gam^r_{n s} - d_n Gam^r_{m s} + Gam^r_{m l} Gam^l_{n s} - Gam^r_{n l} Gam^l_{m s}
    Rup = (np.einsum('mrns...->rsmn...', dGam) - np.einsum('nrms...->rsmn...', dGam)
           + np.einsum('rml...,lns...->rsmn...', Gam, Gam) - np.einsum('rnl...,lms...->rsmn...', Gam, Gam))
    Rdn = np.einsum('rl...,lsmn...->rsmn...', g, Rup)
    # orthonormal frame (jets): base e2 = d_q/|d_q|, e1 = (d_p - g_pq/g_qq d_q)/norm; fibre e3 = d_a/|d_a|, e4 = (d_b - G_ab/G_aa d_a)/norm
    gpp, gpq, gqq = gJ[0][0], gJ[0][1], gJ[1][1]
    Gaa, Gab, Gbb = gJ[2][2], gJ[2][3], gJ[3][3]
    iq = gqq.sqrt().recip()
    rq = gpq / gqq
    n1 = (gpp - gpq * rq).sqrt().recip()
    ia = Gaa.sqrt().recip()
    ra = Gab / Gaa
    n4 = (Gbb - Gab * ra).sqrt().recip()
    eJ = np.zeros((6, 4, 4, n))     # eJ[:, A, mu]
    eJ[:, 0, 0] = n1.c; eJ[:, 0, 1] = (-(rq * n1)).c
    eJ[:, 1, 1] = iq.c
    eJ[:, 2, 2] = ia.c
    eJ[:, 3, 3] = n4.c; eJ[:, 3, 2] = (-(ra * n4)).c
    e = eJ[0]; de = np.zeros((4, 4, 4, n)); de[0] = eJ[1]; de[1] = eJ[2]      # de[alpha, A, mu]
    # nabla_{e_C} e_A = e_C^m (d_m e_A^k + Gam^k_{m l} e_A^l)
    nab = np.einsum('Cm...,mAk...->CAk...', e, de) + np.einsum('Cm...,kml...,Al...->CAk...', e, Gam, e)
    FG = np.einsum('CAk...,kn...,Dn...->CAD...', nab, g, e)       # Gam[C,A,D] = g(nabla_C e_A, e_D)
    RF = np.einsum('rsmn...,Ar...,Bs...,Cm...,Dn...->ABCD...', Rdn, e, e, e, e)
    Ric = np.einsum('CACB...->AB...', RF)
    scal = np.einsum('AA...->...', Ric)
    sqrtg = np.sqrt(np.linalg.det(gm))
    return dict(x=x, p=p, q=q, J=J, sqrtg=sqrtg, e=e, Gam=FG, RF=RF, Ric=Ric, scal=scal, S=Sj.c[0], E=E)
