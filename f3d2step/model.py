"""Entity model on top of the raw SAB records: topology traversal and geometry."""

import math

from . import nurbs, sab

NULL = -1


class Unsupported(Exception):
    pass


class Cursor(object):
    """Sequential typed access to a record's tokens."""

    def __init__(self, tokens, i=0, end=None):
        self.t = tokens
        self.i = i
        self.end = len(tokens) if end is None else end
        self.pending = []  # components of a partially consumed position/vector

    def peek(self):
        return self.t[self.i] if self.i < self.end else (None, None)

    def next(self):
        if self.i >= self.end:
            raise Unsupported("unexpected end of record")
        tok = self.t[self.i]
        self.i += 1
        return tok

    def f(self):
        if self.pending:
            return self.pending.pop(0)
        tag, val = self.next()
        if tag in (sab.TAG_DOUBLE, sab.TAG_FLOAT, sab.TAG_LONG, sab.TAG_SHORT):
            return float(val)
        if tag in (sab.TAG_POSITION, sab.TAG_VECTOR, sab.TAG_VECTOR2D):
            self.pending = list(val[1:])
            return val[0]
        raise Unsupported("expected a number, got tag 0x%02X" % tag)

    def vec(self):
        return (self.f(), self.f(), self.f())

    def int(self):
        tag, val = self.next()
        if tag not in (sab.TAG_LONG, sab.TAG_SHORT, sab.TAG_ENUM, sab.TAG_CHAR):
            raise Unsupported("expected an integer, got tag 0x%02X" % tag)
        return val

    def bool(self):
        tag, val = self.next()
        if tag == sab.TAG_TRUE:
            return True
        if tag == sab.TAG_FALSE:
            return False
        raise Unsupported("expected a boolean, got tag 0x%02X" % tag)

    def ref(self):
        tag, val = self.next()
        if tag != sab.TAG_REF:
            raise Unsupported("expected a reference, got tag 0x%02X" % tag)
        return val.index

    def limit(self):
        """Interval end: a 'finite' flag optionally followed by the value."""
        if self.bool():
            return self.f()
        return None


class BSpline(object):
    """Clamped B-spline curve (v* unused) or surface."""

    def __init__(self):
        self.rational = False
        self.udeg = self.vdeg = 0
        self.uknots = self.vknots = ()
        self.umults = self.vmults = ()
        self.poles = []  # curve: [pt]; surface: [u][v] -> pt
        self.weights = None


_STR_TAGS = (sab.TAG_IDENT, sab.TAG_STR8, sab.TAG_STR16, sab.TAG_STR32, sab.TAG_SUBIDENT)
_CLOSURE_WORDS = ("open", "closed", "periodic")


def _knots(c, n, deg):
    knots, mults = [], []
    for _ in range(n):
        knots.append(c.f())
        mults.append(c.int())
    if n < 2:
        raise Unsupported("degenerate knot vector")
    # ACIS omits one knot at each end; STEP wants fully clamped ends.
    if mults[0] == deg and mults[-1] == deg:
        mults[0] = mults[-1] = deg + 1
    elif not (mults[0] == deg + 1 and mults[-1] == deg + 1):
        raise Unsupported("unclamped (periodic-form) knot vector")
    return knots, mults, sum(mults) - deg - 1


def _enum(c):
    tag, val = c.next()
    if tag in (sab.TAG_ENUM, sab.TAG_LONG) or tag in _STR_TAGS or tag in (sab.TAG_TRUE, sab.TAG_FALSE):
        return val
    raise Unsupported("expected an enum, got tag 0x%02X" % tag)


def read_bs_curve(tokens, i, end, dim=3):
    """Parse 'nubs|nurbs deg closure nknots ...' starting at index i."""
    c = Cursor(tokens, i, end)
    s = BSpline()
    s.rational = c.next()[1] == "nurbs"
    s.udeg = c.int()
    _enum(c)
    n = c.int()
    s.uknots, s.umults, npoles = _knots(c, n, s.udeg)
    if s.rational:
        s.weights = []
    for _ in range(npoles):
        s.poles.append(tuple(c.f() for _ in range(dim)))
        if s.rational:
            s.weights.append(c.f())
    return s, c.i


def read_bs_surface(tokens, i, end):
    """Parse 'nubs|nurbs degu degv [both|u|v] closures singularities nu nv ...'."""
    c = Cursor(tokens, i, end)
    s = BSpline()
    s.rational = c.next()[1] == "nurbs"
    s.udeg = c.int()
    s.vdeg = c.int()
    tag, val = c.peek()
    if tag in _STR_TAGS and val in ("both", "u", "v"):
        c.next()
    for _ in range(4):
        _enum(c)
    nu = c.int()
    nv = c.int()
    s.uknots, s.umults, npu = _knots(c, nu, s.udeg)
    s.vknots, s.vmults, npv = _knots(c, nv, s.vdeg)
    s.poles = [[None] * npv for _ in range(npu)]
    if s.rational:
        s.weights = [[1.0] * npv for _ in range(npu)]
    for v in range(npv):
        for u in range(npu):
            s.poles[u][v] = c.vec()
            if s.rational:
                s.weights[u][v] = c.f()
    return s, c.i


HELIX_SPAN = math.pi / 12.0  # 15 degrees per cubic span


def _helix_nodes(lo, hi):
    """Span boundaries on a global grid, so that helical curves and the
    surfaces they bound get the same approximation."""
    nodes = [lo]
    k = math.floor(lo / HELIX_SPAN) + 1
    while k * HELIX_SPAN < hi:
        u = k * HELIX_SPAN
        if u - nodes[-1] > 1e-6 and hi - u > 1e-6:
            nodes.append(u)
        k += 1
    nodes.append(hi)
    return nodes


def _hermite(nodes, pos, tan):
    """Poles of the C1 piecewise cubic through pos/tan at the nodes."""
    poles = [pos(nodes[0])]
    for i in range(len(nodes) - 1):
        u0, u1 = nodes[i], nodes[i + 1]
        h = (u1 - u0) / 3.0
        p0, t0, p1, t1 = pos(u0), tan(u0), pos(u1), tan(u1)
        poles.append(tuple(p0[j] + t0[j] * h for j in range(3)))
        poles.append(tuple(p1[j] - t1[j] * h for j in range(3)))
        poles.append(p1)
    return poles


class Helix(object):
    """Analytic helix as stored in 'helix_int_cur' data."""

    def __init__(self, c):
        self.lo, self.hi = c.limit(), c.limit()
        self.center, self.major, self.minor, self.pitch = c.vec(), c.vec(), c.vec(), c.vec()
        taper = c.f()
        self.axis = c.vec()
        if abs(taper) > 1e-12:
            raise Unsupported("tapered helix")

    def frame(self, u, x, y, z):
        """x*major + y*minor rotated by angle u about the axis, plus z*axis."""
        cu, su = math.cos(u), math.sin(u)
        a, b = x * cu - y * su, x * su + y * cu
        M, N, A = self.major, self.minor, self.axis
        return tuple(a * M[j] + b * N[j] + z * A[j] for j in range(3))

    def pos(self, u):
        r = self.frame(u, 1.0, 0.0, 0.0)
        k = u / (2.0 * math.pi)
        return tuple(self.center[j] + r[j] + self.pitch[j] * k for j in range(3))

    def tan(self, u):
        t = self.frame(u, 0.0, 1.0, 0.0)
        k = 1.0 / (2.0 * math.pi)
        return tuple(t[j] + self.pitch[j] * k for j in range(3))


def _skip_version(c):
    if c.peek()[0] == sab.TAG_LONG:
        c.int()


def helix_curve(tokens, start, end):
    """Cubic B-spline approximation of a 'helix_int_cur' (which is stored
    analytically, without a NURBS approximation)."""
    c = Cursor(tokens, start, end)
    _skip_version(c)
    hx = Helix(c)
    if hx.lo is None or hx.hi is None or not hx.hi > hx.lo:
        raise Unsupported("unbounded helix")
    nodes = _helix_nodes(hx.lo, hx.hi)
    s = BSpline()
    s.udeg = 3
    s.uknots = nodes
    s.umults = [4] + [3] * (len(nodes) - 2) + [4]
    s.poles = _hermite(nodes, hx.pos, hx.tan)
    return s


def helix_line_surface(tokens, start, end, urange, vrange):
    """B-spline for 'helix_spl_line': a straight line swept along a helix,
    S(u, v) = helix(v) + u * line_direction(v)."""
    c = Cursor(tokens, start, end)
    _skip_version(c)
    u0, u1, v0, v1 = c.limit(), c.limit(), c.limit(), c.limit()
    hx = Helix(c)
    for _ in range(4):  # two null surfaces, two null pcurves
        tag, val = c.next()
        if tag not in _STR_TAGS or not val.startswith("null"):
            raise Unsupported("helix surface with support data")
    d = c.vec()
    if urange[0] is not None and urange[1] is not None:
        u0, u1 = urange
    if vrange[0] is not None and vrange[1] is not None:
        v0, v1 = vrange
    if None in (u0, u1, v0, v1) or not (u1 > u0 and v1 > v0):
        raise Unsupported("unbounded helix surface")
    nodes = _helix_nodes(v0, v1)
    s = BSpline()
    s.udeg, s.vdeg = 1, 3
    s.uknots, s.umults = [u0, u1], [2, 2]
    s.vknots = nodes
    s.vmults = [4] + [3] * (len(nodes) - 2) + [4]
    for u in (u0, u1):
        def pos(v):
            p, r = hx.pos(v), hx.frame(v, d[0], d[1], d[2])
            return tuple(p[j] + u * r[j] for j in range(3))

        def tan(v):
            t, r = hx.tan(v), hx.frame(v, -d[1], d[0], 0.0)
            return tuple(t[j] + u * r[j] for j in range(3))

        s.poles.append(_hermite(nodes, pos, tan))
    return s


class Transform(object):
    def __init__(self, rows=None, trans=(0.0, 0.0, 0.0), scale=1.0):
        self.rows = rows or ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
        self.trans = trans
        self.scale = scale
        self.mirror = False

    def point(self, p):
        r, t, s = self.rows, self.trans, self.scale
        return tuple((p[0] * r[0][j] + p[1] * r[1][j] + p[2] * r[2][j]) * s + t[j] for j in range(3))

    def direction(self, d):
        r = self.rows
        return tuple(d[0] * r[0][j] + d[1] * r[1][j] + d[2] * r[2][j] for j in range(3))


def helix_circ_surface(tokens, start, end, vrange):
    """B-spline for 'helix_spl_circ' (coil): a circle lying in the plane
    through the helix axis, swept along the helix.  u is the angle on the
    circle (exact, rational quadratic), v the helix angle."""
    c = Cursor(tokens, start, end)
    _skip_version(c)
    c.limit(), c.limit()
    v0, v1 = c.limit(), c.limit()
    c.f()
    hx = Helix(c)
    for _ in range(4):  # two null surfaces, two null pcurves
        tag, val = c.next()
        if tag not in _STR_TAGS or not val.startswith("null"):
            raise Unsupported("helix surface with support data")
    radius = c.f()
    if vrange[0] is not None and vrange[1] is not None:
        v0, v1 = vrange
    if v0 is None or v1 is None or not v1 > v0:
        raise Unsupported("unbounded helix surface")
    rmaj = math.sqrt(sum(x * x for x in hx.major))
    if rmaj < 1e-300:
        raise Unsupported("degenerate helix")
    nodes = _helix_nodes(v0, v1)
    s = BSpline()
    s.rational = True
    s.udeg, s.vdeg = 2, 3
    hp = math.pi / 2.0
    s.uknots, s.umults = [-2 * hp, -hp, 0.0, hp, 2 * hp], [3, 2, 2, 2, 3]
    s.vknots = nodes
    s.vmults = [4] + [3] * (len(nodes) - 2) + [4]
    w = math.sqrt(0.5)
    # control polygon of the circle (-cos u, sin u), u = -pi .. pi
    square = [(1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0)]
    s.weights = []
    for k, (a, b) in enumerate(square):
        a, b = a * radius / rmaj, b * radius

        def pos(v, a=a, b=b):
            p, r = hx.pos(v), hx.frame(v, a, 0.0, b)
            return tuple(p[j] + r[j] for j in range(3))

        def tan(v, a=a):
            t, r = hx.tan(v), hx.frame(v, 0.0, a, 0.0)
            return tuple(t[j] + r[j] for j in range(3))

        row = _hermite(nodes, pos, tan)
        s.poles.append(row)
        s.weights.append([w if k % 2 else 1.0] * len(row))
    return s


class _Tokens(object):
    """Minimal record stand-in for resolving a nested subtype block."""

    def __init__(self, tokens):
        self.tokens = tokens


class Model(object):
    """Resolved view of one SAB stream."""

    def __init__(self, reader):
        self.header = reader.header
        self.records = reader.records
        self.subtypes = None  # [(record, open_index, close_index)]

    # -- generic helpers ---------------------------------------------------

    def rec(self, index):
        if index is None or index < 0 or index >= len(self.records):
            return None
        return self.records[index]

    @staticmethod
    def kind(rec):
        return rec.name.rsplit("-", 1)[-1]

    @staticmethod
    def ref(rec, i):
        tag, val = rec.tokens[i]
        if tag != sab.TAG_REF:
            raise Unsupported("%s: field %d is not a reference" % (rec.name, i))
        return val.index

    @staticmethod
    def flag(rec, i):
        tag = rec.tokens[i][0]
        if tag == sab.TAG_TRUE:
            return True
        if tag == sab.TAG_FALSE:
            return False
        raise Unsupported("%s: field %d is not a boolean" % (rec.name, i))

    @staticmethod
    def number(rec, i):
        tag, val = rec.tokens[i]
        if tag not in (sab.TAG_DOUBLE, sab.TAG_FLOAT, sab.TAG_LONG):
            raise Unsupported("%s: field %d is not a number" % (rec.name, i))
        return float(val)

    def chain(self, first, next_field=3):
        """Follow a 'next' linked list starting at record index first."""
        out, seen = [], set()
        idx = first
        while idx is not None and idx >= 0 and idx not in seen:
            seen.add(idx)
            r = self.rec(idx)
            if r is None:
                break
            out.append(r)
            idx = self.ref(r, next_field)
        return out

    # -- topology ----------------------------------------------------------

    def top_level(self):
        n = self.header.bodies
        if n <= 0 or n > len(self.records):
            n = len(self.records)
        return self.records[:n]

    def bodies(self):
        return [r for r in self.top_level() if self.kind(r) == "body"]

    def body_transform(self, body):
        r = self.rec(self.ref(body, 5))
        if r is None:
            return Transform()
        c = Cursor(r.tokens, 2)  # transforms lack the third common field
        rows = (c.vec(), c.vec(), c.vec())
        t = Transform(rows, c.vec(), c.f())
        try:
            c.bool()
            t.mirror = c.bool()
        except Unsupported:
            pass
        return t

    def lumps(self, body):
        return self.chain(self.ref(body, 3))

    def shells(self, lump):
        return self.chain(self.ref(lump, 4))

    def faces(self, shell):
        return self.chain(self.ref(shell, 5))

    def loops(self, face):
        return self.chain(self.ref(face, 4))

    def coedges(self, loop):
        out, seen = [], set()
        idx = self.ref(loop, 4)
        while idx >= 0 and idx not in seen:
            seen.add(idx)
            r = self.rec(idx)
            if r is None:
                break
            out.append(r)
            idx = self.ref(r, 3)
        return out

    # -- subtype table -----------------------------------------------------

    def _build_subtypes(self):
        table = []
        for rec in self.records:
            toks = rec.tokens
            stack = []
            for i, (tag, _) in enumerate(toks):
                if tag == sab.TAG_OPEN:
                    name = toks[i + 1][1] if i + 1 < len(toks) else None
                    if name == "ref":
                        stack.append(None)
                    else:
                        stack.append(len(table))
                        table.append([rec, i, len(toks)])
                elif tag == sab.TAG_CLOSE and stack:
                    k = stack.pop()
                    if k is not None:
                        table[k][2] = i
        self.subtypes = table

    def subtype_block(self, rec, open_index):
        """Resolve the '{ ... }' block starting at open_index (following refs).

        Returns (tokens, start, end, name) where start is the index just
        after the subtype name and end the index of the closing brace.
        """
        toks = rec.tokens
        for _ in range(16):
            name = toks[open_index + 1][1]
            if name != "ref":
                depth = 0
                for j in range(open_index, len(toks)):
                    tag = toks[j][0]
                    if tag == sab.TAG_OPEN:
                        depth += 1
                    elif tag == sab.TAG_CLOSE:
                        depth -= 1
                        if depth == 0:
                            return toks, open_index + 2, j, name
                return toks, open_index + 2, len(toks), name
            if self.subtypes is None:
                self._build_subtypes()
            n = toks[open_index + 2][1]
            if not isinstance(n, int) or n < 0 or n >= len(self.subtypes):
                raise Unsupported("dangling subtype reference %r" % (n,))
            rec, open_index, _ = self.subtypes[n]
            toks = rec.tokens
        raise Unsupported("subtype reference loop")

    @staticmethod
    def _find_bs(toks, start, end, surface):
        """Index of the first top-level spline of the wanted kind in a block."""
        depth = 0
        i = start
        while i < end:
            tag, val = toks[i]
            if tag == sab.TAG_OPEN:
                depth += 1
            elif tag == sab.TAG_CLOSE:
                depth -= 1
            elif depth == 0 and tag in _STR_TAGS and val in ("nubs", "nurbs") and i + 2 < end:
                if toks[i + 1][0] == sab.TAG_LONG:
                    is_surface = toks[i + 2][0] == sab.TAG_LONG
                    if is_surface == surface:
                        return i
            i += 1
        return -1

    @staticmethod
    def _close(toks, open_index):
        """Index of the brace closing the one at open_index."""
        depth = 0
        for j in range(open_index, len(toks)):
            tag = toks[j][0]
            if tag == sab.TAG_OPEN:
                depth += 1
            elif tag == sab.TAG_CLOSE:
                depth -= 1
                if depth == 0:
                    return j
        return len(toks)

    def _block_curve(self, toks, open_index):
        """B-spline of the int_cur block opening at toks[open_index]."""
        rec = _Tokens(toks)
        btoks, start, end, name = self.subtype_block(rec, open_index)
        i = self._find_bs(btoks, start, end, False)
        if i >= 0:
            return read_bs_curve(btoks, i, end)[0]
        if name == "helix_int_cur":
            return helix_curve(btoks, start, end)
        raise Unsupported("curve '%s' has no stored NURBS approximation" % name)

    def _extrusion(self, toks, start, end, reverse):
        """'cyl_spl_sur' without approximation: profile curve swept along a vector."""
        c = Cursor(toks, start, end)
        _skip_version(c)
        tag, val = c.next()
        if tag not in _STR_TAGS or val != "intcurve":
            raise Unsupported("extrusion of a '%s' profile" % (val,))
        curve_reverse = c.bool()
        if c.peek()[0] != sab.TAG_OPEN:
            raise Unsupported("malformed extrusion profile")
        bs = self._block_curve(toks, c.i)
        c.i = self._close(toks, c.i) + 1
        c.limit()
        c.limit()
        axis = c.vec()
        try:  # the part of the profile in use, in the profile's own sense
            c.vec()
            _enum(c)
            lo, hi = c.limit(), c.limit()
            if lo is not None and hi is not None:
                if curve_reverse:
                    lo, hi = -hi, -lo
                nurbs.reseam_curve(bs, lo, hi)
        except Unsupported:
            pass
        return {"bs": bs, "axis": axis, "reverse": reverse ^ curve_reverse}

    # -- geometry ----------------------------------------------------------

    def surface(self, rec):
        """Returns (kind, data) describing a surface record."""
        kind = rec.name
        c = Cursor(rec.tokens, 3)
        if kind == "plane-surface":
            return "plane", {"root": c.vec(), "normal": c.vec(), "u": c.vec()}
        if kind == "cone-surface":
            d = {"center": c.vec(), "axis": c.vec(), "major": c.vec(), "ratio": c.f()}
            c.limit()
            c.limit()
            d["sin"] = c.f()
            d["cos"] = c.f()
            return "cone", d
        if kind == "sphere-surface":
            d = {"center": c.vec(), "radius": c.f()}
            d["u"] = c.vec()
            d["pole"] = c.vec()
            return "sphere", d
        if kind == "torus-surface":
            d = {"center": c.vec(), "axis": c.vec(), "major": c.f(), "minor": c.f()}
            d["u"] = c.vec()
            return "torus", d
        if kind == "spline-surface":
            reverse = c.bool()
            if c.peek()[0] != sab.TAG_OPEN:
                raise Unsupported("spline-surface without subtype data")
            toks, start, end, name = self.subtype_block(rec, c.i)
            i = self._find_bs(toks, start, end, True)
            # the part of the surface this face uses: u lo/hi, v lo/hi
            t = Cursor(rec.tokens, self._close(rec.tokens, c.i) + 1)
            try:
                urange, vrange = (t.limit(), t.limit()), (t.limit(), t.limit())
            except Unsupported:
                urange = vrange = (None, None)
            if i >= 0:
                bs, _ = read_bs_surface(toks, i, end)
                if reverse and None not in urange:
                    urange = (-urange[1], -urange[0])
                nurbs.reseam_surface(bs, urange, vrange)
            elif name == "helix_spl_line":
                bs = helix_line_surface(toks, start, end, urange, vrange)
            elif name == "helix_spl_circ":
                bs = helix_circ_surface(toks, start, end, vrange)
            elif name == "cyl_spl_sur":
                return "extrusion", self._extrusion(toks, start, end, reverse)
            else:
                raise Unsupported("spline surface '%s' has no stored NURBS approximation" % name)
            return "spline", {"bs": bs, "reverse": reverse, "subtype": name}
        raise Unsupported("surface type '%s'" % kind)

    def curve(self, rec):
        """Returns (kind, data) describing a curve record."""
        kind = rec.name
        c = Cursor(rec.tokens, 3)
        if kind == "straight-curve":
            return "line", {"root": c.vec(), "dir": c.vec()}
        if kind == "ellipse-curve":
            return "ellipse", {"center": c.vec(), "normal": c.vec(), "major": c.vec(), "ratio": c.f()}
        if kind.endswith("intcurve-curve"):
            reverse = c.bool()
            if c.peek()[0] != sab.TAG_OPEN:
                raise Unsupported("intcurve without subtype data")
            bs = self._block_curve(rec.tokens, c.i)
            return "spline", {"bs": bs, "reverse": reverse}
        raise Unsupported("curve type '%s'" % kind)

    def point(self, vertex):
        p = self.rec(self.ref(vertex, 5))
        return Cursor(p.tokens, 3).vec()
