"""Tokenizer for ASM/ACIS binary B-rep streams (SAB) as embedded in .f3d files.

A stream is a header followed by records.  Each record is a list of tagged
tokens terminated by TAG_TERMINATOR.  Records reference each other by index.
"""

import struct

TAG_CHAR = 0x02
TAG_SHORT = 0x03
TAG_LONG = 0x04
TAG_FLOAT = 0x05
TAG_DOUBLE = 0x06
TAG_STR8 = 0x07
TAG_STR16 = 0x08
TAG_STR32 = 0x09
TAG_TRUE = 0x0A
TAG_FALSE = 0x0B
TAG_REF = 0x0C
TAG_IDENT = 0x0D
TAG_SUBIDENT = 0x0E
TAG_OPEN = 0x0F
TAG_CLOSE = 0x10
TAG_TERM = 0x11
TAG_STR32B = 0x12
TAG_POSITION = 0x13
TAG_VECTOR = 0x14
TAG_ENUM = 0x15
TAG_VECTOR2D = 0x16
TAG_INT64 = 0x17


class SabError(Exception):
    pass


class Ref(object):
    """Reference to another record by index (-1 is null)."""

    __slots__ = ("index",)

    def __init__(self, index):
        self.index = index

    def __repr__(self):
        return "$%d" % self.index


class Record(object):
    __slots__ = ("index", "name", "tokens")

    def __init__(self, index, name, tokens):
        self.index = index
        self.name = name  # e.g. "plane-surface", "body"
        self.tokens = tokens  # list of (tag, value)

    def __repr__(self):
        return "-%d %s %s" % (self.index, self.name, " ".join(_fmt(t) for t in self.tokens))


def _fmt(tok):
    tag, val = tok
    if tag == TAG_OPEN:
        return "{"
    if tag == TAG_CLOSE:
        return "}"
    if tag == TAG_TERM:
        return "#"
    if tag in (TAG_STR8, TAG_STR16, TAG_STR32, TAG_STR32B):
        return "@%d %s" % (len(val), val)
    if tag == TAG_TRUE:
        return "T"
    if tag == TAG_FALSE:
        return "F"
    if tag in (TAG_POSITION, TAG_VECTOR, TAG_VECTOR2D):
        return " ".join(repr(v) for v in val)
    if tag == TAG_ENUM:
        return "e%d" % val
    return str(val)


class Header(object):
    def __init__(self):
        self.format = ""
        self.version = 0
        self.records = 0
        self.bodies = 0
        self.flags = 0
        self.product = ""
        self.product_version = ""
        self.date = ""
        self.scale = 1.0  # millimetres per model unit
        self.resabs = 1e-6
        self.resnor = 1e-10


class Reader(object):
    def __init__(self, data):
        self.data = data
        self.pos = 0
        self.header = Header()
        self.records = []
        self.long_fmt = "<i"
        self.long_size = 4

    def _long(self):
        (v,) = struct.unpack_from(self.long_fmt, self.data, self.pos)
        self.pos += self.long_size
        return v

    def _token(self):
        d = self.data
        tag = d[self.pos]
        self.pos += 1
        p = self.pos
        if tag in (TAG_LONG, TAG_REF, TAG_ENUM):
            v = self._long()
            if tag == TAG_REF:
                v = Ref(v)
            return tag, v
        if tag == TAG_DOUBLE:
            self.pos = p + 8
            return tag, struct.unpack_from("<d", d, p)[0]
        if tag in (TAG_STR8, TAG_IDENT, TAG_SUBIDENT):
            n = d[p]
            self.pos = p + 1 + n
            return tag, d[p + 1 : p + 1 + n].decode("utf-8", "replace")
        if tag == TAG_STR16:
            (n,) = struct.unpack_from("<H", d, p)
            self.pos = p + 2 + n
            return tag, d[p + 2 : p + 2 + n].decode("utf-8", "replace")
        if tag in (TAG_STR32, TAG_STR32B):
            (n,) = struct.unpack_from("<I", d, p)
            self.pos = p + 4 + n
            return tag, d[p + 4 : p + 4 + n].decode("utf-8", "replace")
        if tag in (TAG_TRUE, TAG_FALSE, TAG_OPEN, TAG_CLOSE, TAG_TERM):
            return tag, None
        if tag in (TAG_POSITION, TAG_VECTOR):
            self.pos = p + 24
            return tag, struct.unpack_from("<3d", d, p)
        if tag == TAG_VECTOR2D:
            self.pos = p + 16
            return tag, struct.unpack_from("<2d", d, p)
        if tag == TAG_CHAR:
            self.pos = p + 1
            return tag, d[p]
        if tag == TAG_SHORT:
            self.pos = p + 2
            return tag, struct.unpack_from("<h", d, p)[0]
        if tag == TAG_FLOAT:
            self.pos = p + 4
            return tag, struct.unpack_from("<f", d, p)[0]
        if tag == TAG_INT64:
            self.pos = p + 8
            return tag, struct.unpack_from("<q", d, p)[0]
        raise SabError("unknown tag 0x%02X at offset %d" % (tag, p - 1))

    def _read_header(self):
        h = self.header
        magic = self.data[:15]
        if magic not in (b"ASM BinaryFile4", b"ASM BinaryFile8", b"ACIS BinaryFile"):
            raise SabError("not an ASM/ACIS binary stream")
        h.format = magic.decode()
        if magic.endswith(b"8"):
            self.long_fmt, self.long_size = "<q", 8
        self.pos = 15
        h.version = self._long()
        h.records = self._long()
        h.bodies = self._long()
        h.flags = self._long()
        h.product = self._token()[1]
        h.product_version = self._token()[1]
        h.date = self._token()[1]
        h.scale = abs(self._token()[1])
        h.resabs = self._token()[1]
        h.resnor = self._token()[1]

    def read(self):
        """Parse the stream; returns the list of active-model records."""
        self._read_header()
        n = len(self.data)
        index = 0
        in_history = False
        while self.pos < n:
            names = []
            tag, val = self._token()
            while tag == TAG_SUBIDENT:
                names.append(val)
                tag, val = self._token()
            if tag != TAG_IDENT:
                raise SabError("expected record name at offset %d" % self.pos)
            names.append(val)
            name = "-".join(names)
            if name.startswith("End-of-"):
                if name.endswith("-data"):
                    break
                in_history = False
                continue
            tokens = []
            while self.pos < n:
                tok = self._token()
                if tok[0] == TAG_TERM:
                    break
                tokens.append(tok)
            if name.startswith("Begin-of-"):
                in_history = True
                continue
            if in_history:
                continue
            self.records.append(Record(index, name, tokens))
            index += 1
        return self.records


def parse(data):
    r = Reader(data)
    r.read()
    return r
