"""High-level conversion: pick the B-rep bodies of a .f3d and write STEP."""

import collections
import os

from . import sab
from .container import F3dFile
from .model import Model, Unsupported
from .names import blob_markers, body_records
from .step import CM_TO_MM, BodyExporter, StepFile


class Result(object):
    def __init__(self):
        self.outputs = []
        self.solids = 0
        self.sheets = 0
        self.faces = 0
        self.skipped_faces = 0
        self.duplicates = 0
        self.hidden = 0
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

        self.label = None  # name given in Fusion, when it could be recovered
        self.visible = None  # True/False from the design stream, None if unknown

    @property
    def name(self):
        return "Body%d" % self.number

    @property
    def title(self):
        return self.label or self.name


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
                if srec is not None and kind in ("sphere", "torus", "cone", "plane"):
                    try:  # faces without vertices (full spheres) still need a position
                        _, d = model.surface(srec)
                        c = d.get("center") or d.get("root")
                        pts.add((kind,) + tuple(round(x, 7) for x in xf.point(c)))
                    except Unsupported:
                        pass
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


def scan(path, all_blobs=False, warn=None, default_name="Body%d"):
    """Returns the BodyInfo list of all non-empty bodies in the file.

    default_name: pattern for bodies Fusion shows with a default name; the
    stored design does not say which language the UI used."""
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
        _assign_names(f3d, bodies, default_name)
    finally:
        f3d.close()
    return bodies


def _assign_names(f3d, bodies, default_name):
    """Match design-stream body records to B-rep bodies through their tags."""
    records, markers = [], {}
    for n in f3d.names():
        parts = n.split("/")
        if parts[-1] == "BulkStream.dat" and len(parts) > 1 and "Design" in parts[-2]:
            try:
                data = f3d.read(n)
                records.extend(body_records(data))
                markers.update(blob_markers(data))
            except Exception:  # names are a nicety; never fail a conversion over them
                continue
    by_tag = {}
    for r in records:
        if r.tag is not None:
            by_tag.setdefault(r.tag, []).append(r)
    hits_of = {}
    for info in bodies:
        try:
            tags = info.model.body_tags(info.body)
        except Unsupported:
            continue
        hits = [r for t in tags for r in by_tag.get(t, [])]
        hits_of[info.number] = hits
        if not hits:
            # copied/mirrored bodies: the record carries the next tag number
            hits = [r for t in tags for r in by_tag.get(t + 1, [])]
        if hits:
            titles = set(r.title(default_name) for r in hits)
            info.label = titles.pop() if len(titles) == 1 else None
            info.visible = any(r.visible for r in hits)
    if len(bodies) == 1 and len(records) == 1 and bodies[0].label is None:
        bodies[0].label = records[0].title(default_name)  # designs without history carry no tags
        bodies[0].visible = records[0].visible
    used = set()
    for info in bodies:
        if info.label is not None:
            used.update(id(r) for r in hits_of.get(info.number, []))
    _assign_inserted(bodies, records, markers, default_name, used)
    _assign_by_order(bodies, records, default_name, used)


def _assign_by_order(bodies, records, default_name, used):
    """Tagless records in a component with tagged ones: records are stored
    in creation order, so a tagless record between records with tags a and
    b describes a body whose tag lies between a and b."""
    anchors = []  # (pos, tag) with tags increasing along the stream
    for r in records:
        if r.tag is not None and (not anchors or r.tag > anchors[-1][1]):
            anchors.append((r.pos, r.tag))
    groups = {}
    for r in records:
        if r.tag is not None or id(r) in used:
            continue
        lo = max([t for p, t in anchors if p < r.pos] or [0])
        hi = min([t for p, t in anchors if p > r.pos] or [1 << 31])
        groups.setdefault((lo, hi), []).append(r)
    pending = {}
    for info in bodies:
        if info.label is not None or info.visible is not None:
            continue
        try:
            tags = info.model.body_tags(info.body)
        except Unsupported:
            continue
        tag = min(tags)  # the body's own tag; later ones come from features
        for (lo, hi), recs in groups.items():
            if lo < tag < hi:
                pending.setdefault((lo, hi), []).append((tag, info))
    for key, cands in pending.items():
        recs = sorted(groups[key], key=lambda r: r.pos)
        if len(recs) != len(cands):
            recs = [r for r in recs if r.visible]
            if len(recs) != len(cands):
                continue
        for (_, info), r in zip(sorted(cands, key=lambda c: c[0]), recs):
            info.label = r.title(default_name)
            info.visible = r.visible
            used.add(id(r))


def _assign_inserted(bodies, records, markers, default_name, used):
    """Bodies of inserted components carry no usable tags; their records
    (tagless, visible) sit next to the mention of the component's blob."""
    by_blob = {}
    for info in bodies:
        by_blob.setdefault(info.blob, []).append(info)
    free = [r for r in records if r.tag is None and r.visible and id(r) not in used]
    for blob, infos in sorted(by_blob.items(), key=lambda kv: kv[1][0].number):
        if any(i.label is not None or i.visible is not None for i in infos):
            continue
        guid = blob[5:-5] if blob.endswith(".smbh") else blob[5:-4]
        if guid not in markers:
            continue
        at = markers[guid]
        near = sorted(free, key=lambda r: abs(r.pos - at))[: len(infos)]
        if len(near) < len(infos) or any(abs(r.pos - at) > 200000 for r in near):
            continue
        near.sort(key=lambda r: r.pos)
        for info, r in zip(infos, near):
            info.label = r.title(default_name)
            info.visible = True
            free.remove(r)
            used.add(id(r))


def _filename(text):
    """Body name made safe for use in a file name."""
    text = "".join("_" if c in '/\\:*?"<>|' or ord(c) < 32 else c for c in text)
    return text.strip(" .")


def _write(path, name, infos, res):
    step = StepFile(name)
    for info in infos:
        exporter = BodyExporter(step, info.model, res.warnings.append)
        solids, sheets, faces, skipped = exporter.body(info.body, info.title)
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


def convert(
    src,
    dst,
    all_blobs=False,
    keep_duplicates=False,
    select=None,
    split=False,
    solids_only=False,
    include_hidden=False,
    default_name="Body%d",
):
    """Convert src (.f3d) to dst (.step).

    select: optional set of body numbers (as shown by scan) to export.
    split:  write one file per body, named <dst stem>-<body name>.step.
    solids_only: leave out surface bodies.
    include_hidden: also export bodies Fusion shows as hidden.
    default_name: pattern for bodies that keep Fusion's default name.
    """
    res = Result()
    infos = scan(src, all_blobs, res.warnings.append, default_name)
    if select is not None:
        missing = sorted(set(select) - set(i.number for i in infos))
        if missing:
            raise Unsupported("no such body: %s" % ", ".join(str(n) for n in missing))
        infos = [i for i in infos if i.number in select]
    else:
        if not keep_duplicates:
            res.duplicates = sum(1 for i in infos if i.duplicate_of is not None)
            infos = [i for i in infos if i.duplicate_of is None]
        if not include_hidden:
            res.hidden = sum(1 for i in infos if i.visible is False)
            infos = [i for i in infos if i.visible is not False]
    if solids_only:
        infos = [i for i in infos if i.solid]
    if not infos:
        raise Unsupported("no bodies found")
    name = os.path.splitext(os.path.basename(src))[0]
    if split:
        stem, ext = os.path.splitext(dst)
        used = set()
        for info in infos:
            # named after the body; the number is added only to tell namesakes apart
            suffix = _filename(info.title) or info.name
            if suffix.lower() in used:
                suffix = "%s-%d" % (suffix, info.number)
            used.add(suffix.lower())
            _write("%s-%s%s" % (stem, suffix, ext or ".step"), "%s %s" % (name, info.title), [info], res)
    else:
        _write(dst, name, infos, res)
    if not res.outputs:
        raise Unsupported("no convertible bodies found")
    return res
