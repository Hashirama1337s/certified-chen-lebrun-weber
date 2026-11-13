"""t_amend3.py -- unit checks of the amend3 fail-safes (jets.inv0, certify._ub).
1. a batch with no singular / non-finite box: amend3 inv0 is bit-identical to the amend2 body (reproduced verbatim below).
2. a batch with a singular box, a NaN box and a good box: no exception; bad boxes get mid 0, radius inf; the good box encloses
   the true inverse.  3. _ub: NaN -> inf; finite unchanged."""
import numpy as np, math
import jets, certify
from common import Disc, SAFE
U2 = 2.3e-16


def inv0_amend2(G):
    n = G.m.shape[0]
    Gm = np.moveaxis(G.m, -1, 0)
    X = np.linalg.inv(Gm)
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
    return Disc(Xd, np.broadcast_to(b, Xd.shape).copy())


rng = np.random.default_rng(1)
A = rng.normal(size=(4, 4, 300)) + 1j * rng.normal(size=(4, 4, 300)) + 4 * np.eye(4)[:, :, None]
G = Disc(A, np.full(A.shape, 1e-9))
o, n_ = inv0_amend2(G), jets.inv0(G)
ok1 = np.array_equal(o.m, n_.m) and np.array_equal(o.r, n_.r)
B = A[:, :, :3].copy()
B[:, :, 0] = 0; B[0, 0, 0] = 1                      # singular
B[1, 2, 1] = np.nan                                 # NaN
G2 = Disc(B, np.full(B.shape, 1e-9))
with np.errstate(all='ignore'):
    r = jets.inv0(G2)
true = np.linalg.inv(B[:, :, 2])
ok2 = bool(np.all(np.isinf(r.r[:, :, 0])) and np.all(np.isinf(r.r[:, :, 1])) and np.all(r.m[:, :, :2] == 0)
           and np.all(np.abs(r.m[:, :, 2] - true) <= r.r[:, :, 2]) and np.all(np.isfinite(r.r[:, :, 2])))
ok3 = certify._ub(Disc(np.array([1.0, np.nan]), np.zeros(2))) == math.inf and certify._ub(Disc(np.array([1.0, -3.0]), np.zeros(2))) == 3.0 * (1 + SAFE)
print('1 bit-identical on a regular batch: %s\n2 singular/NaN boxes -> inf, good box encloses: %s\n3 _ub NaN -> inf: %s' % (ok1, ok2, ok3))
print('ALL PASS' if ok1 and ok2 and ok3 else 'FAIL')
