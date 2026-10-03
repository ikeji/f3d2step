"""Body names from the design stream.

The design stream (BulkStream.dat) is a dump of Fusion's object graph and is
not self-describing, so this does not parse it.  It only looks for the byte
pattern of body records.  Such a record holds the physical properties, the
number used for the default name ("Body<n>") and the user-given name; just
before it comes the entity tag that the B-rep body carries as an attribute.
"""

import re
import struct

# entity tag followed by an ASCII object key
_HEAD = re.compile(rb"(?s)(.{4})\x00{5}\x01\x00{5}[\x01-\x09]\x00\x00\x00\d")
# GUID string, two doubles, default-name number, colour, visibility flag
_ANCHOR = re.compile(rb"(?s)\x24\x00\x00\x00(?:[0-9a-fA-F-]\x00){36}.{16}(.{4})\x01\x00{4}.{4}\x01([\x00\x01])")
_REF7 = re.compile(rb"(?s)\x01.{4}\x00{7}")
_REF6 = re.compile(rb"(?s)\x01.{4}\x00{6}")
_FLAGS = re.compile(rb"\x01\x01\x00{9}")


def _ustr(data, pos):
    if pos + 4 > len(data):
        return None, pos
    (n,) = struct.unpack_from("<I", data, pos)
    if n > 512 or pos + 4 + 2 * n > len(data):
        return None, pos
    try:
        return data[pos + 4 : pos + 4 + 2 * n].decode("utf-16-le"), pos + 4 + 2 * n
    except UnicodeDecodeError:
        return None, pos


def _name_after(data, pos):
    """Parse the tail of a body record starting at pos; returns the stored
    name ('' if the body was never renamed) or None if it does not parse."""
    m = _FLAGS.search(data, pos, pos + 256)
    if m is None:
        return None
    pos = m.end()
    for _ in range(2):  # two GUID strings
        s, pos = _ustr(data, pos)
        if not s:
            return None
    pos += 4
    s, pos = _ustr(data, pos)  # material name
    if s is None:
        return None
    m = _REF7.match(data, pos)
    if m is None:
        return None
    pos = m.end()
    if data[pos : pos + 4] == b"\x24\x00\x00\x00":  # older files: one more GUID
        s, pos = _ustr(data, pos)
    m = _REF6.match(data, pos)
    if m is None:
        return None
    name, _ = _ustr(data, m.end())
    if name is None or any(ord(c) < 32 for c in name):
        return None
    return name


_BLOB = re.compile(rb"B\x00R\x00E\x00P\x00\.\x00((?:[0-9a-f-]\x00){36})\.\x00s\x00m\x00b\x00")


class BodyRecord(object):
    __slots__ = ("tag", "name", "index", "visible", "pos")

    def __init__(self, tag, name, index, visible, pos):
        self.tag = tag  # entity tag shared with the B-rep body, or None
        self.name = name  # '' when the body keeps its default name
        self.index = index  # number in the default name (localized "Body<n>")
        self.visible = visible
        self.pos = pos  # offset in the stream

    def title(self, default="Body%d"):
        return self.name or default % self.index


def body_records(data):
    """Returns the BodyRecords found in a design stream, in file order."""
    out = []
    last = 0
    for m in _ANCHOR.finditer(data):
        name = _name_after(data, m.end())
        if name is None:
            continue
        (index,) = struct.unpack("<I", m.group(1))
        if not name and not 0 < index < 100000:
            continue
        tag = None
        for h in _HEAD.finditer(data, last, m.start()):
            (tag,) = struct.unpack("<I", h.group(1))
        out.append(BodyRecord(tag, name, index, m.group(2) == b"\x01", m.start()))
        last = m.end()
    return out


def blob_markers(data):
    """Returns {blob GUID: offset} for the B-rep blob names mentioned in a
    design stream.  The bodies of an inserted component are described by
    tagless records near the mention of that component's blob."""
    out = {}
    for m in _BLOB.finditer(data):
        out.setdefault(m.group(1).decode("utf-16-le"), m.start())
    return out
