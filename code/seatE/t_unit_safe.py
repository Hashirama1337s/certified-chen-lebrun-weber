"""exact-rational unit test of the per-operation guards (referee recommendation): the disc must contain the EXACT result.
python t_unit_safe.py [mutate]   (mutate: SAFE = TINY = 0, must FAIL)"""
import sys
from fractions import Fraction as Fr
import numpy as np
import common
if len(sys.argv) > 1:
    common.SAFE = 0.0; common.TINY = 0.0
from common import Disc
ok = True
cases = []
a, b = 1.0, 3 * 2.0 ** -54
cases.append(('Disc(1) + Disc(3*2^-54)', Disc(a) + Disc(b), Fr(a) + Fr(b)))
c = 1.0 + 2.0 ** -27
cases.append(('Disc(1+2^-27)^2', Disc(c) * Disc(c), Fr(c) * Fr(c)))
cases.append(('1/Disc(3)', Disc(3.0).recip(), Fr(1, 3)))
for name, D, exact in cases:
    mid = Fr(float(D.m.real.item())); rad = Fr(float(D.r.item()))
    good = abs(exact - mid) <= rad
    ok &= good
    print('%-26s exact-mid = %.3e  rad = %.3e  contained: %s' % (name, float(abs(exact - mid)), float(rad), good))
print('ALL CONTAINED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
