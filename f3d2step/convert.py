"""High-level conversion: pick the B-rep bodies of a .f3d and write STEP."""

import collections
import os

from . import sab
from .container import F3dFile
from .model import Model, Unsupported
from .step import CM_TO_MM, BodyExporter, StepFile


class Result(object):
    def __init__(self):
        self.outputs = []
        self.solids = 0
        self.sheets = 0
        self.faces = 0
        self.skipped_faces = 0
        self.duplicates = 0
        self.warnings = []


class BodyInfo(object):
    """One body found in a .f3d, numbered from 1 in file order."""

    def __init__(self, number, blob, model, body):
        self.number = number
        self.blob = blob
        self.model = model
        self.body = body
        self.faces = 0
        self.kinds = collections.Counter()
        self.lo = self.hi = None  # extent of the vertices, in mm
        self.duplicate_of = None
        self.signature = None
        self.solid = True  # False for surface (sheet) bodies

    @property
    def name(self):
        return "Body%d" % self.number


def _select_blobs(names, all_blobs):
    """Designs with history keep the bodies in .smbh; the .smb next to them
    only holds sketch profiles and construction data."""
    if all_blobs:
        return names
    hist = [n for n in names if n.lower().endswith(".smbh")]
    return hist or [n for n in names if n.lower().endswith(".smb")]


def _inspect(info):
    """Fill in face statistics, extent and a geometry fingerprint."""
    model, body = info.model, info.body
    xf = model.body_transform(body)
    pts = set()
    for lump in model.lumps(body):
        for shell in model.shells(lump):
            for f in model.faces(shell):
                info.faces += 1
                srec = model.rec(model.ref(f, 7))
                kind = srec.name.split("-")[0] if srec else "none"
                info.kinds[kind] += 1
                pts.add((kind, model.flag(f, 8)))
                if model.flag(f, 9):  # double-sided face
                    info.solid = False
                for loop in model.loops(f):
                    for ce in model.coedges(loop):
                        e = model.rec(model.ref(ce, 6))
                        if model.ref(ce, 5) < 0 and e is not None and model.ref(e, 8) >= 0:
                            info.solid = False  # open boundary
                        v = model.rec(model.ref(e, 3)) if e is not None else None
                        if v is not None:
                            pts.add(tuple(round(x, 7) for x in xf.point(model.point(v))))
    coords = [p for p in pts if len(p) == 3]
    # circles often carry a single vertex; add their true extent
    for lump in model.lumps(body):
        for shell in model.shells(lump):
            for f in model.faces(shell):
                for loop in model.loops(f):
                    for ce in model.coedges(loop):
                        e = model.rec(model.ref(ce, 6))
                        c = model.rec(model.ref(e, 8)) if e is not None else None
                        if c is None or c.name != "ellipse-curve" or model.ref(e, 3) != model.ref(e, 5):
                            continue
                        _, d = model.curve(c)
                        center, n = xf.point(d["center"]), xf.direction(d["normal"])
                        r = sum(x * x for x in d["major"]) ** 0.5 * abs(xf.scale) * max(1.0, abs(d["ratio"]))
                        ext = [r * max(0.0, 1.0 - n[j] * n[j]) ** 0.5 for j in range(3)]
                        coords.append(tuple(center[j] - ext[j] for j in range(3)))
                        coords.append(tuple(center[j] + ext[j] for j in range(3)))
    if coords:
        info.lo = tuple(min(p[j] for p in coords) * CM_TO_MM for j in range(3))
        info.hi = tuple(max(p[j] for p in coords) * CM_TO_MM for j in range(3))
    info.signature = (info.faces, frozenset(pts))


def scan(path, all_blobs=False, warn=None):
    """Returns the BodyInfo list of all non-empty bodies in the file."""
    f3d = F3dFile(path)
    bodies = []
    try:
        names = _select_blobs(f3d.brep_names(), all_blobs)
        if not names:
            raise Unsupported("no B-rep data found")
        seen = {}
        for blob in names:
            model = Model(sab.parse(f3d.read(blob)))
            for body in model.bodies():
                info = BodyInfo(len(bodies) + 1, os.path.basename(blob), model, body)
                try:
                    _inspect(info)
                except Unsupported as ex:
                    if warn:
                        warn("%s: body skipped (%s)" % (info.blob, ex))
                    continue
                if info.faces == 0:
                    continue  # wire or empty body
                info.duplicate_of = seen.setdefault(info.signature, info.number)
                if info.duplicate_of == info.number:
                    info.duplicate_of = None
                bodies.append(info)
    finally:
        f3d.close()
    return bodies


def _write(path, name, infos, res):
    step = StepFile(name)
    for info in infos:
        exporter = BodyExporter(step, info.model, res.warnings.append)
        solids, sheets, faces, skipped = exporter.body(info.body, info.name)
        res.solids += solids
        res.sheets += sheets
        res.faces += faces
        res.skipped_faces += skipped
    if not (step.solids or step.surfaces):
        res.warnings.append("%s: nothing convertible, file not written" % os.path.basename(path))
        return
    with open(path, "w", encoding="ascii", newline="\n") as fp:
        step.write(fp)
    res.outputs.append(path)


def convert(src, dst, all_blobs=False, keep_duplicates=False, select=None, split=False, solids_only=False):
    """Convert src (.f3d) to dst (.step).

    select: optional set of body numbers (as shown by scan) to export.
    split:  write one file per body, named <dst stem>-BodyN.step.
    solids_only: leave out surface bodies.
    """
    res = Result()
    infos = scan(src, all_blobs, res.warnings.append)
    if select is not None:
        missing = sorted(set(select) - set(i.number for i in infos))
        if missing:
            raise Unsupported("no such body: %s" % ", ".join(str(n) for n in missing))
        infos = [i for i in infos if i.number in select]
    elif not keep_duplicates:
        res.duplicates = sum(1 for i in infos if i.duplicate_of is not None)
        infos = [i for i in infos if i.duplicate_of is None]
    if solids_only:
        infos = [i for i in infos if i.solid]
    if not infos:
        raise Unsupported("no bodies found")
    name = os.path.splitext(os.path.basename(src))[0]
    if split:
        stem, ext = os.path.splitext(dst)
        for info in infos:
            _write("%s-%s%s" % (stem, info.name, ext or ".step"), "%s %s" % (name, info.name), [info], res)
    else:
        _write(dst, name, infos, res)
    if not res.outputs:
        raise Unsupported("no convertible bodies found")
    return res
