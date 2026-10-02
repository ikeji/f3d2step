import math
import unittest

from f3d2step import nurbs
from f3d2step.model import BSpline, Helix, Cursor, helix_curve
from f3d2step import sab
from f3d2step.step import fmt


def evaluate(bs, t):
    """De Boor evaluation of a (possibly rational) clamped B-spline curve."""
    U = nurbs._expand(bs.uknots, bs.umults)
    p = bs.udeg
    n = len(bs.poles)
    k = p
    while k < n - 1 and not (U[k] <= t < U[k + 1]):
        k += 1
    w = bs.weights or [1.0] * n
    d = [nurbs._homog(bs.poles[j + k - p], w[j + k - p]) for j in range(p + 1)]
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            a = (t - U[j + k - p]) / (U[j + 1 + k - r] - U[j + k - p])
            d[j] = tuple((1 - a) * x + a * y for x, y in zip(d[j - 1], d[j]))
    return nurbs._dehomog(d[p])[0]


def closed_curve():
    """Closed cubic through a rounded square, clamped with coincident ends."""
    bs = BSpline()
    bs.udeg = 3
    pts = [(1, 0, 0), (1, 1, 0), (0, 1, 0.5), (-1, 1, 0), (-1, 0, 0), (-1, -1, 0), (0, -1, -0.5), (1, -1, 0), (1, 0, 0)]
    bs.poles = [tuple(float(x) for x in p) for p in pts]
    bs.uknots = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    bs.umults = [4, 1, 1, 1, 1, 1, 4]
    return bs


class ReseamTest(unittest.TestCase):
    def check(self, rational):
        ref = closed_curve()
        bs = closed_curve()
        if rational:
            for c in (ref, bs):
                c.rational = True
                c.weights = [1.0, 2.0, 0.5, 1.0, 3.0, 1.0, 0.7, 1.5, 1.0]
        nurbs.reseam_curve(bs, 4.3, 8.9)
        self.assertAlmostEqual(bs.uknots[0], 4.3)
        self.assertAlmostEqual(bs.uknots[-1], 10.3)
        self.assertEqual(sum(bs.umults), len(bs.poles) + bs.udeg + 1)
        for i in range(60):
            t = 4.3 + 5.999 * i / 59.0
            src = t if t < 6.0 else t - 6.0
            for a, b in zip(evaluate(bs, t), evaluate(ref, src)):
                self.assertAlmostEqual(a, b, places=9)

    def test_polynomial(self):
        self.check(False)

    def test_rational(self):
        self.check(True)

    def test_untouched_when_range_fits(self):
        bs = closed_curve()
        nurbs.reseam_curve(bs, 1.0, 5.0)
        self.assertEqual(bs.uknots, closed_curve().uknots)

    def test_open_curve_untouched(self):
        bs = closed_curve()
        bs.poles[-1] = (2.0, 0.0, 0.0)
        nurbs.reseam_curve(bs, 4.3, 8.9)
        self.assertEqual(bs.uknots[0], 0.0)


class HelixTest(unittest.TestCase):
    def test_approximation_error(self):
        D = sab.TAG_DOUBLE
        toks = [(sab.TAG_TRUE, None), (D, 1.0), (sab.TAG_TRUE, None), (D, 30.0)]
        for v in ((0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 0.5)):
            toks.append((sab.TAG_VECTOR, tuple(float(x) for x in v)))
        toks += [(D, 0.0), (sab.TAG_VECTOR, (0.0, 0.0, 1.0))]
        bs = helix_curve(toks, 0, len(toks))
        hx = Helix(Cursor(toks, 0))
        worst = 0.0
        for i in range(500):
            t = 1.0 + 28.999 * i / 499.0
            worst = max(worst, math.dist(evaluate(bs, t), hx.pos(t)))
        self.assertLess(worst, 1e-4)


class FormatTest(unittest.TestCase):
    def test_reals(self):
        self.assertEqual(fmt(1), "1.")
        self.assertEqual(fmt(0.5), "0.5")
        self.assertEqual(fmt(1e-20), "0.")
        self.assertEqual(fmt(1.5e-5), "1.5E-05")
        self.assertEqual(fmt(1e22), "1.E+22")


if __name__ == "__main__":
    unittest.main()
