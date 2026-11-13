"""fbox.py -- seat E stage 1: vectorised enclosures of D^alpha f (f = sum C_ij T_i(x1) T_j(x2), seat R's M84 solution)
over many complex boxes per call.

Patch = (real centre c0, polydisc radius rp).  In arb (prec bits):
  a_beta = Taylor coefficients of f at c0 for |beta| <= m (Chebyshev Taylor tables by t_{n+1,k} = 2c t_{n,k} + 2 t_{n,k-1}
  - t_{n-1,k} at the exact centre, then a = T1^T C T2 as arb matrices).
Remainder (Cauchy): f is entire; on the polydisc |z_i - c0_i| <= R,  |f| <= Mf = sum |C_ij| rho1^i rho2^j  with
  rho_i = A_i + sqrt(A_i^2 - 1),  A_i = (|c0_i - 1| + |c0_i + 1|)/2 + R  (the polydisc lies inside E_rho1 x E_rho2,
  |T_n| <= rho^n there), so |a_beta| <= Mf R^-|beta|, and for |d_i| <= rp < R
  |D^alpha (f - P_m)(c0 + d)| <= eps_alpha = Mf sum_{beta >= alpha, |beta| > m} beta!/(beta - alpha)! R^-|beta| rp^|beta - alpha|
  = Mf [ prod_i alpha_i! R^-alpha_i (1 - rp/R)^-(alpha_i + 1)  -  finite part |beta| <= m ]   (evaluated in arb).
Evaluation over boxes (float64 discs, common.Disc rounding model): D^alpha P_m(c0 + d) = sum_g A_alpha[g] d1^g1 d2^g2 with
  mid = complex matmul of float mids, rad = sum (|A| + rad A) U1^g1 U2^g2 - sum |A| |d1m|^g1 |d2m|^g2 + SAFE2 * magnitude,
  U_i = |d_i mid| + rad d_i  (|d^k - dm^k| <= U^k - |dm|^k), SAFE2 = 2^-36 >> (2m + 10) u for m <= 120.
Result: Disc enclosing D^alpha f over every point of the box, for 0 <= |alpha| <= kd.
"""
import math, time
from math import factorial
import numpy as np
from flint import arb, acb, arb_mat, ctx
import common
from common import Disc, SAFE

SAFE2 = 2.0 ** -36


class FPatch:
    def __init__(self, prob, c0, rp, m=56, R=0.45, kd=7):
        self.c0 = (float(c0[0]), float(c0[1]))
        self.rp, self.m, self.R, self.kd = rp, m, R, kd
        M = prob.M
        cs = [arb(v) for v in self.c0]          # exact binary64 centres
        T = []
        for c in cs:
            t = [[arb(0)] * (m + 1) for _ in range(M + 1)]
            t[0][0] = arb(1)
            t[1][0] = c; t[1][1] = arb(1)
            for n in range(1, M):
                for k in range(m + 1):
                    v = 2 * c * t[n][k] - t[n - 1][k]
                    if k:
                        v += 2 * t[n][k - 1]
                    t[n + 1][k] = v
            T.append(arb_mat(t))                 # (M+1) x (m+1)
        Cm = prob.Cmat(False)
        a = T[0].transpose() * Cm * T[1]         # (m+1) x (m+1), a[k1,k2]
        self.a = a
        # Mf
        Cabs = np.array([[abs(int(v)) for v in row] for row in prob.Cint], dtype=object)
        rho = []
        for c in cs:
            A = ((c - 1).abs_upper() + (c + 1).abs_upper()) / 2 + arb(R)
            rho.append(A + (A * A - 1).sqrt())
        pw = [[r ** i for i in range(M + 1)] for r in rho]
        Mf = arb(0)
        for i in range(M + 1):
            for j in range(M + 1):
                if Cabs[i][j]:
                    Mf += arb(int(Cabs[i][j])) * pw[0][i] * pw[1][j]
        Mf = Mf * arb(2) ** (-prob.B)
        self.Mf = float(Mf.upper())
        # derivative tables and remainders
        self.alphas = [(a1, a2) for n in range(kd + 1) for a1 in range(n, -1, -1) for a2 in [n - a1]]
        self.A = {}
        self.eps = {}
        Ra, rpa = arb(R), arb(rp)
        for (a1, a2) in self.alphas:
            mm = m - a1 - a2
            Am = np.zeros((mm + 1, mm + 1), dtype=np.complex128)
            Ar = np.zeros((mm + 1, mm + 1))
            for g1 in range(mm + 1):
                for g2 in range(mm + 1 - g1):
                    v = a[a1 + g1, a2 + g2] * (factorial(a1 + g1) // factorial(g1)) * (factorial(a2 + g2) // factorial(g2))
                    Am[g1, g2], Ar[g1, g2] = common.disc_from_acb(v)
            self.A[(a1, a2)] = (Am, Ar, np.abs(Am))
            full = arb(1)
            for ai in (a1, a2):
                full *= arb(factorial(ai)) * Ra ** (-ai) * (1 - rpa / Ra) ** (-(ai + 1))
            fin = arb(0)
            for b1 in range(a1, m + 1):
                for b2 in range(a2, m + 1 - b1):
                    fin += arb(factorial(b1) // factorial(b1 - a1) * (factorial(b2) // factorial(b2 - a2))) * Ra ** (-(b1 + b2)) * rpa ** (b1 + b2 - a1 - a2)
            self.eps[(a1, a2)] = float((Mf * (full - fin)).upper())

    def eval(self, x1, x2, alphas=None):
        """x1, x2: Disc arrays (complex boxes).  Returns dict alpha -> Disc over the boxes."""
        m1 = x1.m - self.c0[0]; m2 = x2.m - self.c0[1]
        d1 = Disc(m1, x1.r + 2.0 ** -52 * np.abs(m1))            # amend2: rounding of the shift (u |result|)
        d2 = Disc(m2, x2.r + 2.0 ** -52 * np.abs(m2))
        U1 = np.abs(d1.m) + d1.r; U2 = np.abs(d2.m) + d2.r
        if np.any(np.maximum(U1, U2) > self.rp * (1 - 1e-12)):
            raise ValueError('box outside patch polydisc')
        m = self.m
        k = np.arange(m + 1)
        P1 = d1.m[:, None] ** k; P2 = d2.m[:, None] ** k                 # mids of powers
        Q1 = np.abs(d1.m)[:, None] ** k; Q2 = np.abs(d2.m)[:, None] ** k
        V1 = U1[:, None] ** k; V2 = U2[:, None] ** k
        out = {}
        for al in (alphas or self.alphas):
            Am, Ar, Aa = self.A[al]
            mm = Am.shape[0]
            mid = ((P1[:, :mm] @ Am) * P2[:, :mm]).sum(1)
            big = ((V1[:, :mm] @ (Aa + Ar)) * V2[:, :mm]).sum(1)
            sml = ((Q1[:, :mm] @ Aa) * Q2[:, :mm]).sum(1)
            rad = (big - sml) + SAFE2 * big + self.eps[al] + common.TINY
            out[al] = Disc(*common._clean(mid, rad))
        return out


def patch_grid(prob, h=0.1, rp=0.075, **kw):
    """real centres on a grid of spacing h over [-1.1, 1.1]^2 (covers the pentagon with margin)."""
    cs = np.arange(-1.1, 1.1 + 1e-9, h)
    return {(i, j): FPatch(prob, (cs[i], cs[j]), rp, **kw) for i in range(len(cs)) for j in range(len(cs))}, cs


class PatchBank:
    """lazily built patches on a real grid (spacing h); a box goes to the nearest centre; polydisc radius rp must contain
    it (|x - c0| + rad <= rp componentwise), else ValueError (no silent fallback)."""

    def __init__(self, prob, h=0.1, rp=0.12, m=60, R=0.5, kd=7):
        self.prob, self.h, self.rp, self.m, self.R, self.kd = prob, h, rp, m, R, kd
        self.P = {}

    def patch(self, key):
        if key not in self.P:
            self.P[key] = FPatch(self.prob, (key[0] * self.h, key[1] * self.h), self.rp, m=self.m, R=self.R, kd=self.kd)
        return self.P[key]

    def __call__(self, x1, x2):
        k1 = np.rint(x1.m.real / self.h).astype(int); k2 = np.rint(x2.m.real / self.h).astype(int)
        nb = len(k1)
        out = None
        for key in set(zip(k1.tolist(), k2.tolist())):
            sel = np.nonzero((k1 == key[0]) & (k2 == key[1]))[0]
            D = self.patch(key).eval(x1[sel], x2[sel])
            if out is None:
                out = {al: Disc(np.zeros(nb, complex), np.zeros(nb)) for al in D}
            for al, v in D.items():
                out[al].m[sel] = v.m; out[al].r[sel] = v.r
        return out
