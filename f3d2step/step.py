"""STEP (ISO 10303-21, AP214) writer for B-rep bodies taken from a Model."""

import datetime
import math

from . import nurbs
from .model import Unsupported

CM_TO_MM = 10.0  # Fusion stores lengths in centimetres


def fmt(x):
    x = float(x)
    if x != x or x in (float("inf"), float("-inf")):
        raise Unsupported("non-finite number in geometry")
    if abs(x) < 1e-13:
        return "0."
    s = repr(x)
    if "e" in s:
        m, e = s.split("e")
        if "." not in m:
            m += "."
        return m + "E" + e
    if s.endswith(".0"):
        s = s[:-1]
    return s


def _b(v):
    return ".T." if v else ".F."


def _norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n < 1e-300:
        raise Unsupported("zero-length direction")
    return (v[0] / n, v[1] / n, v[2] / n), n


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _str(s):
    out = []
    wide = False  # inside a \X2\ ... \X0\ run
    for ch in s:
        o = ord(ch)
        if 32 <= o < 127:
            if wide:
                out.append("\\X0\\")
                wide = False
            out.append("''" if ch == "'" else "\\\\" if ch == "\\" else ch)
            continue
        if o > 0xFFFF:
            o -= 0x10000
            units = (0xD800 + (o >> 10), 0xDC00 + (o & 0x3FF))
        else:
            units = (o,)
        if not wide:
            out.append("\\X2\\")
            wide = True
        out.extend("%04X" % u for u in units)
    if wide:
        out.append("\\X0\\")
    return "'" + "".join(out) + "'"


class StepFile(object):
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.solids = []
        self.surfaces = []

    def add(self, text):
        self.lines.append(text)
        return len(self.lines)

    def ids(self, ids):
        return "(" + ",".join("#%d" % i for i in ids) + ")"

    def point(self, p):
        return self.add("CARTESIAN_POINT('',(%s,%s,%s))" % (fmt(p[0]), fmt(p[1]), fmt(p[2])))

    def direction(self, d):
        return self.add("DIRECTION('',(%s,%s,%s))" % (fmt(d[0]), fmt(d[1]), fmt(d[2])))

    def axis2(self, origin, z, x):
        z, _ = _norm(z)
        k = _dot(x, z)
        x = (x[0] - k * z[0], x[1] - k * z[1], x[2] - k * z[2])
        if _dot(x, x) < 1e-20:
            a = (1.0, 0.0, 0.0) if abs(z[0]) < 0.9 else (0.0, 1.0, 0.0)
            x = _cross(_cross(z, a), z)
        x, _ = _norm(x)
        return self.add(
            "AXIS2_PLACEMENT_3D('',#%d,#%d,#%d)" % (self.point(origin), self.direction(z), self.direction(x))
        )

    def write(self, fp):
        n = _str(self.name)
        now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        a = self.add
        ctx = a("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
        a("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2010,#%d)" % ctx)
        pctx = a("PRODUCT_CONTEXT('',#%d,'mechanical')" % ctx)
        prod = a("PRODUCT(%s,%s,'',(#%d))" % (n, n, pctx))
        a("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#%d))" % prod)
        pdf = a("PRODUCT_DEFINITION_FORMATION('','',#%d)" % prod)
        pdc = a("PRODUCT_DEFINITION_CONTEXT('part definition',#%d,'design')" % ctx)
        pd = a("PRODUCT_DEFINITION('design','',#%d,#%d)" % (pdf, pdc))
        pds = a("PRODUCT_DEFINITION_SHAPE('','',#%d)" % pd)
        lu = a("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
        au = a("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        su = a("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
        unc = a("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-05),#%d,'distance_accuracy_value','')" % lu)
        gctx = a(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#%d))"
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#%d,#%d,#%d))REPRESENTATION_CONTEXT('',''))" % (unc, lu, au, su)
        )
        origin = self.axis2((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
        reps = []
        if self.solids:
            reps.append(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('',%s,#%d)" % (self.ids([origin] + self.solids), gctx)
            )
        if self.surfaces:
            reps.append(
                "MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',%s,#%d)" % (self.ids([origin] + self.surfaces), gctx)
            )
        if len(reps) == 1:
            a("SHAPE_DEFINITION_REPRESENTATION(#%d,#%d)" % (pds, a(reps[0])))
        else:
            base = a("SHAPE_REPRESENTATION('',(#%d),#%d)" % (origin, gctx))
            a("SHAPE_DEFINITION_REPRESENTATION(#%d,#%d)" % (pds, base))
            for r in reps:
                a("SHAPE_REPRESENTATION_RELATIONSHIP('','',#%d,#%d)" % (base, a(r)))

        fp.write("ISO-10303-21;\nHEADER;\n")
        fp.write("FILE_DESCRIPTION(('converted from Autodesk Fusion .f3d by f3d2step'),'2;1');\n")
        fp.write("FILE_NAME(%s,'%s',(''),(''),'f3d2step','f3d2step','');\n" % (n, now))
        fp.write("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 3 1 1 }'));\nENDSEC;\nDATA;\n")
        for i, line in enumerate(self.lines):
            fp.write("#%d=%s;\n" % (i + 1, line))
        fp.write("ENDSEC;\nEND-ISO-10303-21;\n")


class BodyExporter(object):
    """Writes the bodies of one Model into a StepFile."""

    def __init__(self, step, model, warn):
        self.step = step
        self.m = model
        self.warn = warn
        self.xf = None
        self.scale = CM_TO_MM
        self._vertices = {}
        self._edges = {}
        self._surfaces = {}

    # -- coordinates -------------------------------------------------------

    def P(self, p):
        q = self.xf.point(p)
        return (q[0] * CM_TO_MM, q[1] * CM_TO_MM, q[2] * CM_TO_MM)

    def D(self, d):
        return self.xf.direction(d)

    def L(self, length):
        return length * self.scale

    # -- geometry ----------------------------------------------------------

    def _bs_points(self, pts):
        return [self.step.point(self.P(p)) for p in pts]

    def _bs_curve(self, bs):
        s = self.step
        pts = s.ids(self._bs_points(bs.poles))
        mults = "(" + ",".join(str(m) for m in bs.umults) + ")"
        knots = "(" + ",".join(fmt(k) for k in bs.uknots) + ")"
        if not bs.rational:
            return s.add(
                "B_SPLINE_CURVE_WITH_KNOTS('',%d,%s,.UNSPECIFIED.,.F.,.F.,%s,%s,.UNSPECIFIED.)"
                % (bs.udeg, pts, mults, knots)
            )
        w = "(" + ",".join(fmt(x) for x in bs.weights) + ")"
        return s.add(
            "(BOUNDED_CURVE()B_SPLINE_CURVE(%d,%s,.UNSPECIFIED.,.F.,.F.)"
            "B_SPLINE_CURVE_WITH_KNOTS(%s,%s,.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()"
            "RATIONAL_B_SPLINE_CURVE(%s)REPRESENTATION_ITEM(''))" % (bs.udeg, pts, mults, knots, w)
        )

    def _bs_surface(self, bs):
        s = self.step
        rows = "(" + ",".join(s.ids(self._bs_points(row)) for row in bs.poles) + ")"
        mu = "(" + ",".join(str(m) for m in bs.umults) + ")"
        mv = "(" + ",".join(str(m) for m in bs.vmults) + ")"
        ku = "(" + ",".join(fmt(k) for k in bs.uknots) + ")"
        kv = "(" + ",".join(fmt(k) for k in bs.vknots) + ")"
        if not bs.rational:
            return s.add(
                "B_SPLINE_SURFACE_WITH_KNOTS('',%d,%d,%s,.UNSPECIFIED.,.F.,.F.,.F.,%s,%s,%s,%s,.UNSPECIFIED.)"
                % (bs.udeg, bs.vdeg, rows, mu, mv, ku, kv)
            )
        w = "(" + ",".join("(" + ",".join(fmt(x) for x in row) + ")" for row in bs.weights) + ")"
        return s.add(
            "(BOUNDED_SURFACE()B_SPLINE_SURFACE(%d,%d,%s,.UNSPECIFIED.,.F.,.F.,.F.)"
            "B_SPLINE_SURFACE_WITH_KNOTS(%s,%s,%s,%s,.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()"
            "RATIONAL_B_SPLINE_SURFACE(%s)REPRESENTATION_ITEM('')SURFACE())"
            % (bs.udeg, bs.vdeg, rows, mu, mv, ku, kv, w)
        )

    def _ellipse(self, center, normal, major, ratio):
        """Circle/ellipse entity running counter-clockwise about normal."""
        s = self.step
        xdir, r = _norm(major)
        r = self.L(r)
        if abs(ratio - 1.0) < 1e-12:
            return s.add("CIRCLE('',#%d,%s)" % (s.axis2(self.P(center), self.D(normal), self.D(xdir)), fmt(r)))
        a, b = r, r * abs(ratio)
        if b > a:  # STEP expects the first semi-axis along the placement's X
            xdir = _cross(_norm(normal)[0], xdir)
            a, b = b, a
        pos = s.axis2(self.P(center), self.D(normal), self.D(xdir))
        return s.add("ELLIPSE('',#%d,%s,%s)" % (pos, fmt(a), fmt(b)))

    def curve(self, rec, lo, hi):
        """Returns (step id, reversed); lo..hi is the parameter range in use."""
        s = self.step
        kind, d = self.m.curve(rec)
        if kind == "line":
            v = s.add("VECTOR('',#%d,1.)" % s.direction(_norm(self.D(d["dir"]))[0]))
            return s.add("LINE('',#%d,#%d)" % (s.point(self.P(d["root"])), v)), False
        if kind == "ellipse":
            return self._ellipse(d["center"], d["normal"], d["major"], d["ratio"]), False
        if d["reverse"]:
            lo, hi = -hi, -lo
        nurbs.reseam_curve(d["bs"], lo, hi)
        return self._bs_curve(d["bs"]), d["reverse"]

    def _seam_dir(self, face, center, axis, default):
        """Reference direction for a surface of revolution, chosen so that its
        seam passes through a vertex of the face boundary.  Importers cope
        badly with boundary edges that cross the seam in mid-span."""
        m = self.m
        axis, _ = _norm(axis)
        for loop in m.loops(face):
            for ce in m.coedges(loop):
                e = m.rec(m.ref(ce, 6))
                if e is None or m.ref(e, 8) < 0 or m.ref(e, 3) == m.ref(e, 5):
                    continue
                v = m.rec(m.ref(e, 3))
                if v is None:
                    continue
                p = m.point(v)
                r = (p[0] - center[0], p[1] - center[1], p[2] - center[2])
                k = _dot(r, axis)
                r = (r[0] - k * axis[0], r[1] - k * axis[1], r[2] - k * axis[2])
                if _dot(r, r) > 1e-16:
                    return r
        return default

    def surface(self, rec, face):
        """Returns (step id, flipped): flipped if the STEP normal opposes the ASM one."""
        if rec.index in self._surfaces:
            return self._surfaces[rec.index]
        s = self.step
        kind, d = self.m.surface(rec)
        flip = False
        if kind == "plane":
            sid = s.add("PLANE('',#%d)" % s.axis2(self.P(d["root"]), self.D(d["normal"]), self.D(d["u"])))
        elif kind == "cone":
            xdir, r = _norm(d["major"])
            sin, cos = d["sin"], d["cos"]
            flip = cos < 0
            axis = d["axis"]
            circular = abs(d["ratio"] - 1.0) < 1e-12
            if circular:
                xdir = self._seam_dir(face, d["center"], axis, xdir)
            if abs(sin) < 1e-12:
                if circular:
                    pos = s.axis2(self.P(d["center"]), self.D(axis), self.D(xdir))
                    sid = s.add("CYLINDRICAL_SURFACE('',#%d,%s)" % (pos, fmt(self.L(r))))
                else:
                    e = self._ellipse(d["center"], axis, d["major"], d["ratio"])
                    v = s.add("VECTOR('',#%d,1.)" % s.direction(_norm(self.D(axis))[0]))
                    sid = s.add("SURFACE_OF_LINEAR_EXTRUSION('',#%d,#%d)" % (e, v))
            else:
                if not circular:
                    raise Unsupported("elliptical cone")
                if abs(cos) < 1e-12:
                    raise Unsupported("degenerate cone")
                tan = sin / cos  # radius growth per unit along the axis
                if tan < 0:
                    axis = (-axis[0], -axis[1], -axis[2])
                    tan = -tan
                pos = s.axis2(self.P(d["center"]), self.D(axis), self.D(xdir))
                sid = s.add("CONICAL_SURFACE('',#%d,%s,%s)" % (pos, fmt(self.L(r)), fmt(math.atan(tan))))
        elif kind == "sphere":
            flip = d["radius"] < 0
            xdir = self._seam_dir(face, d["center"], d["pole"], d["u"])
            pos = s.axis2(self.P(d["center"]), self.D(d["pole"]), self.D(xdir))
            sid = s.add("SPHERICAL_SURFACE('',#%d,%s)" % (pos, fmt(self.L(abs(d["radius"])))))
        elif kind == "torus":
            major, minor = d["major"], d["minor"]
            flip = minor < 0
            xdir = self._seam_dir(face, d["center"], d["axis"], d["u"])
            pos = s.axis2(self.P(d["center"]), self.D(d["axis"]), self.D(xdir))
            if major <= 0 or abs(minor) >= major:
                # Self-intersecting ("apple"/"lemon") torus.  STEP has
                # DEGENERATE_TOROIDAL_SURFACE for it, but importers handle it
                # poorly, so write the generating circle revolved about the axis.
                z, _ = _norm(d["axis"])
                k = _dot(xdir, z)
                x, _ = _norm((xdir[0] - k * z[0], xdir[1] - k * z[1], xdir[2] - k * z[2]))
                c = d["center"]
                cc = (c[0] + major * x[0], c[1] + major * x[1], c[2] + major * x[2])
                circle = s.add(
                    "CIRCLE('',#%d,%s)"
                    % (s.axis2(self.P(cc), self.D(_cross(x, z)), self.D(x)), fmt(self.L(abs(minor))))
                )
                ax = s.add("AXIS1_PLACEMENT('',#%d,#%d)" % (s.point(self.P(c)), s.direction(self.D(z))))
                sid = s.add("SURFACE_OF_REVOLUTION('',#%d,#%d)" % (circle, ax))
            else:
                sid = s.add("TOROIDAL_SURFACE('',#%d,%s,%s)" % (pos, fmt(self.L(major)), fmt(self.L(abs(minor)))))
        elif kind == "extrusion":
            v = s.add("VECTOR('',#%d,1.)" % s.direction(_norm(self.D(d["axis"]))[0]))
            sid = s.add("SURFACE_OF_LINEAR_EXTRUSION('',#%d,#%d)" % (self._bs_curve(d["bs"]), v))
            flip = d["reverse"]
        else:
            sid = self._bs_surface(d["bs"])
            flip = d["reverse"]
        self._surfaces[rec.index] = (sid, flip)
        return sid, flip

    # -- topology ----------------------------------------------------------

    def vertex(self, index):
        if index not in self._vertices:
            v = self.m.rec(index)
            if v is None:
                raise Unsupported("edge without vertex")
            p = self.step.point(self.P(self.m.point(v)))
            self._vertices[index] = self.step.add("VERTEX_POINT('',#%d)" % p)
        return self._vertices[index]

    def edge(self, e):
        """Returns the EDGE_CURVE id, or None for a degenerate (curveless) edge."""
        if e.index in self._edges:
            return self._edges[e.index]
        m = self.m
        crec = m.rec(m.ref(e, 8))
        if crec is None:
            eid = None
        else:
            ereversed = m.flag(e, 9)
            lo, hi = m.number(e, 4), m.number(e, 6)
            if ereversed:
                lo, hi = -hi, -lo
            cid, creversed = self.curve(crec, lo, hi)
            same = not (ereversed ^ creversed)
            eid = self.step.add(
                "EDGE_CURVE('',#%d,#%d,#%d,%s)"
                % (self.vertex(m.ref(e, 3)), self.vertex(m.ref(e, 5)), cid, _b(same))
            )
        self._edges[e.index] = eid
        return eid

    def face(self, f):
        m, s = self.m, self.step
        srec = m.rec(m.ref(f, 7))
        if srec is None:
            raise Unsupported("face without surface")
        bounds = []
        for loop in m.loops(f):
            oriented = []
            last_vertex = None
            for ce in m.coedges(loop):
                e = m.rec(m.ref(ce, 6))
                if e is None:
                    continue
                eid = self.edge(e)
                if eid is None:
                    last_vertex = m.ref(e, 3)
                    continue
                oriented.append(s.add("ORIENTED_EDGE('',*,*,#%d,%s)" % (eid, _b(not m.flag(ce, 7)))))
            if oriented:
                lid = s.add("EDGE_LOOP('',%s)" % s.ids(oriented))
            elif last_vertex is not None:
                lid = s.add("VERTEX_LOOP('',#%d)" % self.vertex(last_vertex))
            else:
                continue
            bounds.append(s.add("FACE_BOUND('',#%d,.T.)" % lid))
        sid, flip = self.surface(srec, f)
        same = not (m.flag(f, 8) ^ flip)
        if not bounds:
            bounds.append(s.add("FACE_BOUND('',#%d,.T.)" % self._seam_loop(srec, same)))
        return s.add("ADVANCED_FACE('',%s,#%d,%s)" % (s.ids(bounds), sid, _b(same)))

    def _seam_loop(self, srec, same):
        """Boundary for a face covering a whole sphere or torus.  ASM leaves
        such faces without loops; STEP readers want the seams as edges."""
        s = self.step
        kind, d = self.m.surface(srec)
        c = d["center"]

        def off(*terms):
            return tuple(c[j] + sum(k * v[j] for k, v in terms) for j in range(3))

        def circle(center, normal, xdir, r):
            return s.add("CIRCLE('',#%d,%s)" % (s.axis2(self.P(center), self.D(normal), self.D(xdir)), fmt(self.L(r))))

        def vertex(p):
            return s.add("VERTEX_POINT('',#%d)" % s.point(self.P(p)))

        def oriented(edge, forward):
            return s.add("ORIENTED_EDGE('',*,*,#%d,%s)" % (edge, _b(forward)))

        if kind == "sphere":
            z, _ = _norm(d["pole"])
            k = _dot(d["u"], z)
            x, _ = _norm(tuple(d["u"][j] - k * z[j] for j in range(3)))
            r = abs(d["radius"])
            south, north = vertex(off((-r, z))), vertex(off((r, z)))
            # meridian from the south pole through x to the north pole
            mer = circle(c, _cross(x, z), tuple(-v for v in z), r)
            e = s.add("EDGE_CURVE('',#%d,#%d,#%d,.T.)" % (south, north, mer))
            edges = [oriented(e, True), oriented(e, False)]
        elif kind == "torus" and d["major"] > abs(d["minor"]):
            z, _ = _norm(d["axis"])
            k = _dot(d["u"], z)
            x, _ = _norm(tuple(d["u"][j] - k * z[j] for j in range(3)))
            R, r = d["major"], abs(d["minor"])
            v = vertex(off((R + r, x)))
            big = circle(c, z, x, R + r)
            tube = circle(off((R, x)), _cross(x, z), x, r)
            ea = s.add("EDGE_CURVE('',#%d,#%d,#%d,.T.)" % (v, v, big))
            eb = s.add("EDGE_CURVE('',#%d,#%d,#%d,.T.)" % (v, v, tube))
            edges = [oriented(ea, True), oriented(eb, True), oriented(ea, False), oriented(eb, False)]
        else:
            raise Unsupported("unbounded %s face" % kind)
        if not same:
            edges.reverse()
        return s.add("EDGE_LOOP('',%s)" % s.ids(edges))

    def _is_free(self, coedge):
        """True for a coedge on an open boundary (degenerate edges do not count)."""
        m = self.m
        if m.ref(coedge, 5) >= 0:
            return False
        e = m.rec(m.ref(coedge, 6))
        return e is not None and m.ref(e, 8) >= 0

    def _shell_extent(self, shell):
        lo, hi = [1e300] * 3, [-1e300] * 3
        m = self.m
        for f in m.faces(shell):
            for loop in m.loops(f):
                for ce in m.coedges(loop):
                    e = m.rec(m.ref(ce, 6))
                    v = m.rec(m.ref(e, 3)) if e is not None else None
                    if v is None:
                        continue
                    p = m.point(v)
                    for k in range(3):
                        lo[k] = min(lo[k], p[k])
                        hi[k] = max(hi[k], p[k])
        return sum(max(0.0, hi[k] - lo[k]) for k in range(3))

    def body(self, body, name):
        """Export one body; returns (solids, sheets, faces written, faces skipped)."""
        m, s = self.m, self.step
        self.xf = m.body_transform(body)
        self.scale = CM_TO_MM * abs(self.xf.scale)
        if self.xf.mirror:
            self.warn("%s: mirrored body transform; face orientation may be wrong" % name)
        self._vertices, self._edges, self._surfaces = {}, {}, {}
        stats = [0, 0, 0, 0]
        for lump in m.lumps(body):
            closed, sheets = [], []
            for shell in m.shells(lump):
                ids = []
                is_closed = True
                for f in m.faces(shell):
                    try:
                        ids.append(self.face(f))
                        stats[2] += 1
                    except Unsupported as ex:
                        self.warn("%s: face skipped (%s)" % (name, ex))
                        stats[3] += 1
                        is_closed = False
                        continue
                    if m.flag(f, 9):  # double-sided: sheet face
                        is_closed = False
                    elif is_closed:
                        for loop in m.loops(f):
                            if any(self._is_free(ce) for ce in m.coedges(loop)):
                                is_closed = False
                                break
                if not ids:
                    continue
                if is_closed:
                    closed.append((self._shell_extent(shell), ids))
                else:
                    sheets.append(ids)
            if closed:
                # the outer shell is the largest one; any others are voids
                closed.sort(key=lambda c: -c[0])
                outer = s.add("CLOSED_SHELL('',%s)" % s.ids(closed[0][1]))
                if len(closed) == 1:
                    s.solids.append(s.add("MANIFOLD_SOLID_BREP(%s,#%d)" % (_str(name), outer)))
                else:
                    voids = []
                    for _, ids in closed[1:]:
                        cs = s.add("CLOSED_SHELL('',%s)" % s.ids(ids))
                        voids.append(s.add("ORIENTED_CLOSED_SHELL('',*,#%d,.T.)" % cs))
                    s.solids.append(s.add("BREP_WITH_VOIDS(%s,#%d,%s)" % (_str(name), outer, s.ids(voids))))
                stats[0] += 1
            for ids in sheets:
                sh = s.add("OPEN_SHELL('',%s)" % s.ids(ids))
                s.surfaces.append(s.add("SHELL_BASED_SURFACE_MODEL(%s,(#%d))" % (_str(name), sh)))
                stats[1] += 1
        return tuple(stats)
