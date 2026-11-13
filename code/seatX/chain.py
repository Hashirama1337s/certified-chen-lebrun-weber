"""chain.py -- seatX: the H^2 -> H^6 regularity chain for L = Delta^2 + 2 Ric^ab nabla_a nabla_b + <ds, d.> (= seat A's toric
operator d_i d_j(u^ik f_kl u^lj), ddstar.py) on T^2-invariant f perp affine, as a float function of the pointwise constants
(NOTES-X.md sec. 2).  a_j = ||nabla^j f||, b_j = ||nabla^j v||, v = Delta f, w_j = ||nabla^j L f||, all L2(g_K).
Constants: lam1 (gap of L on affine^perp), Rp = sup max(0, max eig Ric), Rm = sup max(0, -min eig Ric), Rc[q] = sup|nabla^q Ric|
(q<=3), K[q] = sup|nabla^q Riem| (q<=3), S[q] = sup|nabla^q ds| (q<=2).  Kd[q] := 2 Rc[q+1] bounds |nabla^q div Riem| (Bianchi).
  X0  a0 <= w0/lam1
  X1  b0 <= Rp a0 + sqrt(Rp^2 a0^2 + w0 a0)            (<Lf,f> = ||Delta f||^2 - 2 int Ric(grad f, grad f), a1^2 <= b0 a0)
  X2  a1 <= sqrt(b0 a0)
  Gj  scalar 2nd-order step (source sigma_i = ||nabla^i Delta s||):
      D_j = sigma_j + sum_{m<j} sum_{q<=j-1-m} C(j-1-m,q) [(Rc[q] + 2 m K[q]) a_{j-q} + m Kd[q] a_{j-1-q}]
      a_{j+2} <= sqrt(D_j^2 + (Rm + 2 j K[0] + j Kd[0]/2) a_{j+1}^2 + (j Kd[0]/2) a_j^2)
  V   Delta v = L f - B f, B f = 2 Ric.nabla^2 f + <ds, df>:  ||nabla^i B f|| <= E_i = sum_q C(i,q)(2 Rc[q] a_{i-q+2} + S[q] a_{i-q+1})
      b1 <= sqrt((w0 + E0) b0); b_{j+2} by Gj with sigma_j = w_j + E_j (j <= 2)
  f   a_{j+2} by Gj with sigma_j = b_j (j <= 4).
C_reg = sup over w >= 0, w0+w1+w2 = 1 of sum_{j<=6} a_j (grid; the map is homogeneous of degree 1)."""
from math import comb, sqrt
import itertools


def step(j, sig, a, c):
    Rc, K, Rm = c['Rc'], c['K'], c['Rm']
    Kd = [2 * Rc[q + 1] for q in range(3)]
    D = sig[j]
    for m in range(j):
        for q in range(j - m):
            D += comb(j - 1 - m, q) * (((c.get('Rop', Rc[0]) if q == 0 else Rc[q]) + 2 * m * (c.get('kap', K[0]) if q == 0 else K[q])) * a[j - q] + (m * Kd[q] * a[j - 1 - q] if m else 0.0))
    if c.get('interp') and j >= 1:            # a_{j+1}^2 = <Delta nabla^j s, nabla^j s> <= D_j a_j  (DERIVED)
        a[j + 1] = min(a[j + 1], sqrt(D * a[j]))
    return sqrt(D * D + (Rm + 2 * j * c.get('kap', K[0]) + j * Kd[0] / 2) * a[j + 1] ** 2 + (j * Kd[0] / 2) * a[j] ** 2)


def F(w, c, detail=False):
    lam1, Rp, Rc, S = c['lam1'], c['Rp'], c['Rc'], c['S']
    a = [w[0] / lam1]
    b = [Rp * a[0] + sqrt(Rp ** 2 * a[0] ** 2 + w[0] * a[0])]
    a.append(sqrt(b[0] * a[0]))
    if c.get('split'):        # 2 Ric.nabla^2 f = 2 Ric_0.nabla^2 f - (s/2) v ;  sq[q] = sup|nabla^q s| (sq[0] = sup|s|)
        R0, sq = c['R0c'], [c['ssup']] + list(S)
        E = lambda i: sum(comb(i, q) * (2 * R0[q] * a[i - q + 2] + 0.5 * sq[q] * b[i - q] + S[q] * a[i - q + 1]) for q in range(i + 1))
    else:
        E = lambda i: sum(comb(i, q) * (2 * Rc[q] * a[i - q + 2] + S[q] * a[i - q + 1]) for q in range(i + 1))
    a.append(step(0, b, a, c))                       # a2
    sv = [w[0] + E(0)]
    b.append(sqrt(sv[0] * b[0]))                      # b1
    b.append(step(0, sv, b, c))                       # b2
    a.append(step(1, b, a, c))                       # a3
    a.append(step(2, b, a, c))                       # a4
    sv.append(w[1] + E(1)); b.append(step(1, sv, b, c))   # b3
    sv.append(w[2] + E(2)); b.append(step(2, sv, b, c))   # b4
    a.append(step(3, b, a, c))                       # a5
    a.append(step(4, b, a, c))                       # a6
    return (sum(a), a, b) if detail else sum(a)


def creg(c, n=40):
    best, arg = 0.0, None
    for i in range(n + 1):
        for k in range(n + 1 - i):
            w = (i / n, k / n, (n - i - k) / n)
            v = F(w, c)
            if v > best:
                best, arg = v, w
    return best, arg


def consts(lam1, Rp, Rm, Rc, K, S, Rop=None, R0c=None, ssup=None):
    d = dict(lam1=lam1, Rp=Rp, Rm=Rm, Rc=list(Rc), K=list(K), S=list(S))
    if Rop is not None:
        d['Rop'] = Rop
    if R0c is not None:
        d.update(split=True, R0c=list(R0c), ssup=ssup)
    return d


if __name__ == '__main__':
    import json
    out = {}
    ctl = {'square': consts(8.0, 1.0, 0.0, [2, 0, 0, 0], [8 ** 0.5, 0, 0, 0], [0, 0, 0]),
           'triangle': consts(160 / 9, 1.0, 0.0, [2, 0, 0, 0], [4 / 3 ** 0.5, 0, 0, 0], [0, 0, 0]),
           'square_split': consts(8.0, 1.0, 0.0, [2, 0, 0, 0], [8 ** 0.5, 0, 0, 0], [0, 0, 0], Rop=1.0, R0c=[0, 0, 0], ssup=4.0) | {'kap': 1.0},
           'triangle_split': consts(160 / 9, 1.0, 0.0, [2, 0, 0, 0], [4 / 3 ** 0.5, 0, 0, 0], [0, 0, 0], Rop=1.0, R0c=[0, 0, 0], ssup=4.0) | {'kap': 1.0}}
    exact = json.load(open('results_controls.json'))
    for nm in list(ctl):
        ctl[nm + '_interp'] = dict(ctl[nm], interp=True)
    for name, c in ctl.items():
        C, w = creg(c)
        tot, a, b = F(w, c, True)
        Rex = max(m['R'] for m in exact[name.split('_')[0] + '_modes'])
        m0 = exact[name.split('_')[0] + '_modes'][0]
        print('  chain a', [round(v, 4) for v in a], ' b', [round(v, 4) for v in b])
        print('  exact a (lowest mode, w0=1)', [round((r ** 0.5) / m0['Lam'], 4) for r in m0['a2_over_a0']])
        out[name] = dict(C_reg=C, argmax_w=w, a=a, b=b, exact_lower=Rex, factor=C / Rex)
        print(name, 'C_reg', round(C, 3), 'at w', w, 'exact >=', round(Rex, 4), 'factor <=', round(C / Rex, 2), flush=True)
    json.dump(out, open('results_chain_ctrl.json', 'w'), indent=1)
