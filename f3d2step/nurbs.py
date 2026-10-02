"""B-spline helpers: moving the seam of closed (periodic) curves and surfaces.

ASM stores periodic splines as clamped splines whose ends coincide, and lets
edges and faces run across that seam.  STEP readers treat such splines as
plain open ones, so the seam is moved to where the edge/face starts.
"""

TOL = 1e-9


def _expand(knots, mults):
    out = []
    for k, m in zip(knots, mults):
        out.extend([k] * m)
    return out


def _compress(expanded):
    knots, mults = [], []
    for k in expanded:
        if knots and k == knots[-1]:
            mults[-1] += 1
        else:
            knots.append(k)
            mults.append(1)
    return knots, mults


def _insert(U, rows, deg, t):
    """Boehm knot insertion of t into every pole row (homogeneous poles)."""
    k = 0
    for i in range(len(U) - 1):
        if U[i] <= t < U[i + 1]:
            k = i
            break
    alphas = [(t - U[i]) / (U[i + deg] - U[i]) for i in range(k - deg + 1, k + 1)]
    out = []
    for P in rows:
        Q = P[: k - deg + 1]
        for a, i in zip(alphas, range(k - deg + 1, k + 1)):
            Q.append(tuple(a * x + (1.0 - a) * y for x, y in zip(P[i], P[i - 1])))
        Q.extend(P[k:])
        out.append(Q)
    return U[: k + 1] + [t] + U[k + 1 :], out


def _rotate_rows(knots, mults, deg, rows, t):
    """Restart closed clamped pole rows at parameter t.

    rows: list of homogeneous pole lists sharing one knot vector.
    Returns (knots, mults, rows) for the domain [t, t + period].
    """
    U = _expand(knots, mults)
    k0, k1 = U[0], U[-1]
    period = k1 - k0
    # snap to an existing knot to avoid slivers
    for k in knots:
        if abs(k - t) < TOL * max(1.0, abs(period)):
            t = k
    have = U.count(t)
    for _ in range(deg - have):
        U, rows = _insert(U, rows, deg, t)
    r = max(i for i, k in enumerate(U) if k == t)
    i = r - deg  # pole index of the curve point at t
    inner_hi = [k for k in U if t < k < k1]
    inner_lo = [k + period for k in U if k0 < k < t]
    newU = [t] * (deg + 1) + inner_hi + [k1] * deg + inner_lo + [t + period] * (deg + 1)
    new_rows = [P[i:] + P[1 : i + 1] for P in rows]
    nk, nm = _compress(newU)
    return nk, nm, new_rows


def _wrap_start(lo, hi, k0, k1):
    """Parameter at which a closed spline over [k0, k1] must restart so that
    [lo, hi] does not cross its seam; None if no change is needed."""
    period = k1 - k0
    if lo is None or hi is None or period <= 0:
        return None
    eps = TOL * max(1.0, abs(period))
    if hi - lo > period + eps:
        return None
    if lo >= k0 - eps and hi <= k1 + eps:
        return None
    t = lo - period * ((lo - k0) // period)
    if t - k0 < eps or k1 - t < eps or t + (hi - lo) <= k1 + eps:
        return None
    return t


def _close(a, b, scale):
    return all(abs(x - y) <= 1e-7 * scale for x, y in zip(a, b))


def _scale(points):
    lo = [min(p[j] for p in points) for j in range(3)]
    hi = [max(p[j] for p in points) for j in range(3)]
    return max(1e-12, max(hi[j] - lo[j] for j in range(3)))


def _homog(p, w):
    return (p[0] * w, p[1] * w, p[2] * w, w)


def _dehomog(q):
    w = q[3]
    return (q[0] / w, q[1] / w, q[2] / w), w


def reseam_curve(bs, lo, hi):
    """Move the seam of a closed curve so that [lo, hi] is contiguous."""
    if len(bs.poles) < 2 or not _close(bs.poles[0], bs.poles[-1], _scale(bs.poles)):
        return
    t = _wrap_start(lo, hi, bs.uknots[0], bs.uknots[-1])
    if t is None:
        return
    w = bs.weights or [1.0] * len(bs.poles)
    row = [_homog(p, x) for p, x in zip(bs.poles, w)]
    bs.uknots, bs.umults, rows = _rotate_rows(bs.uknots, bs.umults, bs.udeg, [row], t)
    pw = [_dehomog(q) for q in rows[0]]
    bs.poles = [p for p, _ in pw]
    if bs.rational:
        bs.weights = [x for _, x in pw]


def reseam_surface(bs, urange, vrange):
    """Move the seams of a closed surface so the given box is contiguous."""
    nu, nv = len(bs.poles), len(bs.poles[0])
    flat = [p for row in bs.poles for p in row]
    scale = _scale(flat)

    def weight(i, j):
        return bs.weights[i][j] if bs.rational else 1.0

    if nu > 1 and all(_close(bs.poles[0][j], bs.poles[-1][j], scale) for j in range(nv)):
        t = _wrap_start(urange[0], urange[1], bs.uknots[0], bs.uknots[-1])
        if t is not None:
            rows = [[_homog(bs.poles[i][j], weight(i, j)) for i in range(nu)] for j in range(nv)]
            bs.uknots, bs.umults, rows = _rotate_rows(bs.uknots, bs.umults, bs.udeg, rows, t)
            nu = len(rows[0])
            pw = [[_dehomog(rows[j][i]) for j in range(nv)] for i in range(nu)]
            bs.poles = [[p for p, _ in r] for r in pw]
            if bs.rational:
                bs.weights = [[x for _, x in r] for r in pw]
    if nv > 1 and all(_close(bs.poles[i][0], bs.poles[i][-1], scale) for i in range(nu)):
        t = _wrap_start(vrange[0], vrange[1], bs.vknots[0], bs.vknots[-1])
        if t is not None:
            rows = [[_homog(bs.poles[i][j], weight(i, j)) for j in range(nv)] for i in range(nu)]
            bs.vknots, bs.vmults, rows = _rotate_rows(bs.vknots, bs.vmults, bs.vdeg, rows, t)
            pw = [[_dehomog(q) for q in r] for r in rows]
            bs.poles = [[p for p, _ in r] for r in pw]
            if bs.rational:
                bs.weights = [[x for _, x in r] for r in pw]
