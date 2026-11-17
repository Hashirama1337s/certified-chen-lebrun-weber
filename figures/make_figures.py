"""Draw the README figures of certified-chen-lebrun-weber from the data the repository ships.

Usage:  python figures/make_figures.py <repo root> [--png]
Needs Python 3.11 and matplotlib 3.11 (nothing else). Writes figures/*.svg next to this script
(or into ./figures when the script does not sit in a folder named figures). --png also writes PNG previews.

Every plotted number is read from shipped files:
  data/prodA3/clw_e{k}_{i}_{j}.json   the 80 sub-squares of the certificate: element k, offsets (s0, t0), side h, and
                                      the certified upper bound A.sup.rho of |rho| = |Scal(g~) - S| on the sub-square
  data/cert_refined_a1.json           rho_sup (the maximum over the 80 sub-squares, as used by the assembly)
  data/result_refined_a1.json         x (the closing quantity, must be < 1) and d_h
  data/logs/assemble_null.log         x of the null control (deliberately broken, sampled inputs; must fail)
  code/seatR/sol_clw_B448_M84.json    gamma = 2/a*, which fixes the pentagon P = [-1,1]^2 cut by x1 + x2 >= -gamma
The sub-square outlines are drawn with the element map of code/seatE/common.py (ElemMap / elem_x):
  l_a = s^2 A(t), l_b = t^2 B(s), A(t) = (1-t^2) l_a(M_b) + t^2 l_a(O), B(s) = (1-s^2) l_b(M_a) + s^2 l_b(O),
re-implemented here in floating point for drawing only (no plotted value depends on it).
"""
import glob
import json
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import FancyBboxPatch, Polygon

# ---------------------------------------------------------------- brand style (figures/STYLE.md)
BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
W_IN = 12.5  # 900 pt viewBox width

plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "axes.facecolor": BG, "figure.facecolor": "none", "savefig.facecolor": "none",
    "axes.edgecolor": AXC, "axes.labelcolor": TXT, "xtick.color": AXC, "ytick.color": AXC,
    "text.color": TXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.fontsize": 14, "lines.linewidth": 2.0, "axes.axisbelow": True,
})
PNG = "--png" in sys.argv
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
# sequential, one hue (blues from #3987e5): dark = small mismatch, light = larger mismatch
BLUES = LinearSegmentedColormap.from_list("blues", ["#15263d", "#1f4f8a", "#3987e5", "#8fbcf2", "#e3eefc"])


def new_fig(h_in):
    fig = plt.figure(figsize=(W_IN, h_in))
    fig.patches.append(FancyBboxPatch((0.004, 0.006), 0.992, 0.988, boxstyle="round,pad=0,rounding_size=0.015",
                                      transform=fig.transFigure, facecolor=BG, edgecolor=BORDER, linewidth=1.5,
                                      zorder=-10))
    return fig


def title(fig, main, sub):
    h = fig.get_figheight()
    fig.text(0.035, 1 - 0.42 / h, main, fontsize=21, fontweight="bold", color=TXT, va="top")
    fig.text(0.035, 1 - 0.88 / h, sub, fontsize=15, color=AXC, va="top")


def save(fig, out, name):
    path = os.path.join(out, name)
    fig.savefig(path, format="svg", metadata={"Date": None, "Creator": None})
    with open(path, encoding="utf-8") as fh:  # mathtext pieces name one font only; give them the full stack
        svg = fh.read().replace("font-family: 'Segoe UI';", "font-family: 'Segoe UI', 'Helvetica', 'Arial', sans-serif;")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(svg)
    if PNG:
        prev = os.path.join(out, os.pardir, "_preview")
        os.makedirs(prev, exist_ok=True)
        fig.savefig(os.path.join(prev, name.replace(".svg", ".png")), format="png", dpi=144)
    plt.close(fig)
    print("wrote figures/" + name)


def load(root, rel):
    with open(os.path.join(root, rel), encoding="utf-8") as fh:
        return json.load(fh)


def sci(x, nd=1):
    m, e = ("%.*e" % (nd, float(x))).split("e")
    return "%s × 10%s" % (m, str(int(e)).translate(SUP))


def sci_up(x, nd=1):
    """like sci, but the mantissa is rounded up (for upper bounds)."""
    e = math.floor(math.log10(x))
    m = math.ceil(x / 10 ** e * 10 ** nd - 1e-9) / 10 ** nd
    if m >= 10:
        m, e = m / 10, e + 1
    return "%.*f × 10%s" % (nd, m, str(e).translate(SUP))


# ---------------------------------------------------------------- the pentagon and its 5 x 16 sub-squares
class Pentagon:
    """l_r = nu_r . x + lam_r >= 0; vertex k = edge k cap edge k+1 (code/seatR/clwr.py, pentagon())."""

    def __init__(self, gamma):
        self.nu = [(1, 0), (1, 1), (0, 1), (-1, 0), (0, -1)]
        self.lam = [1.0, gamma, 1.0, 1.0, 1.0]
        d = 5
        V = []
        for k in range(d):
            a, b = self.nu[k], self.nu[(k + 1) % d]
            det = a[0] * b[1] - a[1] * b[0]
            r1, r2 = -self.lam[k], -self.lam[(k + 1) % d]
            V.append(((r1 * b[1] - a[1] * r2) / det, (a[0] * r2 - b[0] * r1) / det))
        self.V = V
        self.O = (sum(v[0] for v in V) / d, sum(v[1] for v in V) / d)
        self.Mid = [((V[r - 1][0] + V[r][0]) / 2, (V[r - 1][1] + V[r][1]) / 2) for r in range(d)]

    def ell(self, r, x):
        return self.nu[r][0] * x[0] + self.nu[r][1] * x[1] + self.lam[r]

    def elem_x(self, k, s, t):
        a, b = k, (k + 1) % 5
        Ma, Mb = self.Mid[a], self.Mid[b]
        al1, al2 = self.ell(a, Mb), self.ell(a, self.O)
        be1, be2 = self.ell(b, Ma), self.ell(b, self.O)
        la = s * s * ((1 - t * t) * al1 + t * t * al2)
        lb = t * t * ((1 - s * s) * be1 + s * s * be2)
        n1, n2 = self.nu[a], self.nu[b]
        det = n1[0] * n2[1] - n1[1] * n2[0]
        r1, r2 = la - self.lam[a], lb - self.lam[b]
        return ((r1 * n2[1] - n1[1] * r2) / det, (n1[0] * r2 - n2[0] * r1) / det)

    def area(self):
        V = self.V
        return abs(sum(V[i][0] * V[(i + 1) % 5][1] - V[(i + 1) % 5][0] * V[i][1] for i in range(5))) / 2


def patch_outline(P, k, s0, t0, h, n=16):
    pts = []
    for i in range(n):
        pts.append(P.elem_x(k, s0 + h * i / n, t0))
    for i in range(n):
        pts.append(P.elem_x(k, s0 + h, t0 + h * i / n))
    for i in range(n):
        pts.append(P.elem_x(k, s0 + h - h * i / n, t0 + h))
    for i in range(n):
        pts.append(P.elem_x(k, s0, t0 + h - h * i / n))
    return pts


def shoelace(pts):
    return abs(sum(pts[i][0] * pts[i - 1][1] - pts[i - 1][0] * pts[i][1] for i in range(len(pts)))) / 2


# ---------------------------------------------------------------- figure 1: the mismatch map (hero)
def fig_map(root, out):
    gamma = float(load(root, "code/seatR/sol_clw_B448_M84.json")["gamma"])
    P = Pentagon(gamma)
    files = sorted(glob.glob(os.path.join(root, "data", "prodA3", "clw_e*_*_*.json")))
    subs = []
    for f in files:
        d = json.load(open(f, encoding="utf-8"))
        A = d["A"]
        subs.append((d["k"], d["i"], d["j"], A["s0"], A["t0"], A["h"], A["sup"]["rho"]))
    keys = {(k, i, j) for k, i, j, *_ in subs}
    assert keys == {(k, i, j) for k in range(5) for i in range(4) for j in range(4)}
    rho_sup = float(load(root, "data/cert_refined_a1.json")["rho_sup"])
    rmax = max(s[6] for s in subs)
    assert rmax <= rho_sup <= rmax * (1 + 1e-11)  # the certificate's rho_sup is this maximum, rounded up
    res = load(root, "data/result_refined_a1.json")
    d_h = float(res["d_h"])

    fig = new_fig(8.6)
    title(fig, "The computed metric, checked patch by patch",
          "How far its curvature misses the value the equation demands, certified on each of 80 patches")
    ax = fig.add_axes([0.06, 0.08, 0.56, 0.76])
    lo = 10 ** math.floor(math.log10(min(s[6] for s in subs)))
    hi = 10 ** math.ceil(math.log10(rmax))
    norm = LogNorm(vmin=lo, vmax=hi)
    tot = 0.0
    for k, i, j, s0, t0, h, r in subs:
        pts = patch_outline(P, k, s0, t0, h)
        tot += shoelace(pts)
        ax.add_patch(Polygon(pts, closed=True, facecolor=BLUES(norm(r)), edgecolor=BG, lw=0.8))
    assert abs(tot - P.area()) < 1e-3 * P.area()  # the 80 patches tile the pentagon
    ax.add_patch(Polygon(P.V, closed=True, fill=False, edgecolor=AXC, lw=1.6))
    ax.set_aspect("equal")
    ax.set_xlim(-1.12, 1.12)
    ax.set_ylim(-1.12, 1.30)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("x₁")
    ax.set_ylabel("x₂", rotation=0, labelpad=12)
    i_max = max(range(len(subs)), key=lambda q: subs[q][6])
    k, i, j, s0, t0, h, r = subs[i_max]
    cx, cy = P.elem_x(k, s0 + h / 2, t0 + h / 2)
    ax.annotate("largest patch bound: %s" % sci_up(r, 2), (cx, cy), xytext=(-0.98, 1.20), color=TXT,
                fontsize=13, va="center", arrowprops=dict(arrowstyle="-", color=TXT, lw=1.2))

    cax = fig.add_axes([0.66, 0.44, 0.025, 0.30])
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=BLUES)
    cb = fig.colorbar(sm, cax=cax)
    cb.outline.set_edgecolor(BORDER)
    cb.ax.tick_params(colors=AXC, labelsize=13)
    fig.text(0.66, 0.76, "certified bound on the mismatch\nin each patch (log scale)", fontsize=14, color=TXT,
             va="bottom")
    fig.text(0.75, 0.59, "every patch:\nbelow %s" % sci_up(rho_sup, 2), fontsize=15, color=TXT, va="center")
    fig.text(0.66, 0.38, "Each inside point of the pentagon stands for\na doughnut surface (2-torus) of points in the\n"
             "4-dimensional space, so the 80 patches\ncover all of it.", fontsize=13, color=AXC, va="top")
    fig.text(0.66, 0.20, "Result: the true CLW Einstein metric lies\nwithin a factor 1 ± %s of the\n"
             "rescaled computed one, at every point\n(conditional on cited results, one a\npreprint)." % sci_up(d_h, 3), fontsize=14,
             color=ACC, va="top")
    save(fig, out, "checked-patch-by-patch.svg")
    return dict(n=len(subs), rho_min=min(s[6] for s in subs), rho_max=rmax, rho_sup=rho_sup, gamma=gamma,
                area_sum=tot, area=P.area(), largest=(k, i, j), d_h=d_h,
                printed=(sci_up(rho_sup, 2), sci_up(d_h, 3)))


# ---------------------------------------------------------------- figure 2: the closing test, and a control that fails
def fig_closing(root, out):
    x = float(load(root, "data/result_refined_a1.json")["x"])
    log = open(os.path.join(root, "data/logs/assemble_null.log"), encoding="utf-8").read()
    m = re.search(r"^x = 4 a1 a2\^2 a3_Y\s*=\s*([0-9.eE+-]+)", log, re.M)
    x_null = float(m.group(1))
    assert "GATE x < 1 (K14): FAIL" in log and x < 1 < x_null
    fig = new_fig(5.6)
    title(fig, "The test that decides the proof",
          "The proof closes only if one number, x, is below 1. A deliberately broken input must fail it.")
    ax = fig.add_axes([0.30, 0.17, 0.66, 0.52])
    rows = [("the certificate\n(certified inputs)", x, ACC), ("null control\n(broken input, sampled\ndry run; must fail)", x_null, CAT[1])]
    for yi, (lab, v, col) in enumerate(rows):
        y = 1 - yi
        ax.barh(y, v, height=0.5, color=col, alpha=0.30, edgecolor=col, lw=1)
        ax.text(v + 0.03, y, "x = %.4f   %s" % (v, "PASS" if v < 1 else "FAIL"), color=col, fontsize=15,
                va="center", fontweight="bold")
        ax.text(-0.03, y, lab, color=TXT, fontsize=14, va="center", ha="right",
                transform=matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData))
    ax.axvline(1, color=TXT, lw=2)
    ax.text(1.02, 1.55, "must be below 1", color=TXT, fontsize=14, va="center")
    ax.set_xlim(0, 2.4)
    ax.set_ylim(-0.8, 1.8)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("closing quantity x")
    save(fig, out, "the-deciding-test.svg")
    return dict(x=x, x_null=x_null)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if len(args) != 1:
        sys.exit(__doc__)
    root = args[0]
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")
    os.makedirs(out, exist_ok=True)
    info = dict(map=fig_map(root, out), closing=fig_closing(root, out))
    if "--info" in sys.argv:
        print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
