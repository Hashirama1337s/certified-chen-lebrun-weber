"""exact problem data in arb for BAR-MET1: a* / nu balls, polygon, Donaldson S and c by closed-form edge integrals."""
import mpmath as mp
from flint import arb, ctx
import clwr, met1_core as mc


def problem(case, solfile, prec=320):
    ctx.prec = prec
    mp.mp.prec = prec + 200
    if case == 'page':
        nu = clwr.page_nu(int(prec * 0.3) + 80)
        gam = 2 - 3 * ((1 - nu) / (1 + nu)) ** 2
        P = clwr.page_trapezoid(gam)
    else:
        a = clwr.astar(int(prec * 0.3) + 80)
        if case == 'null':
            a = a * mp.mpf('1.001')
        P = clwr.pentagon(2 / a)
    Sabc = clwr.S_affine(P)
    c = clwr.c_polygon(P, Sabc, nq=120)
    C, d = clwr.load(solfile)
    pr = mc.Problem(P, [mc.A(v) for v in Sabc], mc.A(c), C, d['B'], prec)
    return pr, P, Sabc, c
