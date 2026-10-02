"""Access to the members of a .f3d archive.

A .f3d file is a zip archive.  Older files use deflate; newer ones use
Zstandard (zip method 93), which the zipfile module only handles from
Python 3.14 on, so a few fallbacks are provided.
"""

import os
import shutil
import struct
import subprocess
import zipfile

ZSTD_METHOD = 93


class ContainerError(Exception):
    pass


def _zstd_decompress(data, size):
    try:
        from compression import zstd  # Python >= 3.14

        return zstd.decompress(data)
    except ImportError:
        pass
    try:
        import zstandard

        return zstandard.ZstdDecompressor().decompressobj().decompress(data)
    except ImportError:
        pass
    try:
        import pyzstd

        return pyzstd.decompress(data)
    except ImportError:
        pass
    exe = shutil.which("zstd")
    if exe:
        p = subprocess.run([exe, "-d", "-c", "-q"], input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if p.returncode == 0 or len(p.stdout) == size:
            return p.stdout
        raise ContainerError("zstd failed: %s" % p.stderr.decode("utf-8", "replace").strip())
    raise ContainerError(
        "this file is Zstandard-compressed; use Python >= 3.14, or install the "
        "'zstandard' Python package or the 'zstd' command"
    )


class F3dFile(object):
    def __init__(self, path):
        if not os.path.isfile(path):
            raise ContainerError("no such file")
        if not zipfile.is_zipfile(path):
            raise ContainerError("not a .f3d file (zip archive expected)")
        self.zip = zipfile.ZipFile(path)

    def names(self):
        return [n for n in self.zip.namelist() if not n.endswith("/")]

    def read(self, name):
        info = self.zip.getinfo(name)
        if info.compress_type != ZSTD_METHOD or ZSTD_METHOD in getattr(zipfile, "compressor_names", {}):
            try:
                return self.zip.read(name)
            except NotImplementedError:
                if info.compress_type != ZSTD_METHOD:
                    raise ContainerError("%s: unsupported zip compression method %d" % (name, info.compress_type))
        fp = self.zip.fp
        fp.seek(info.header_offset)
        hdr = fp.read(30)
        if hdr[:4] != b"PK\x03\x04":
            raise ContainerError("%s: bad local file header" % name)
        nlen, xlen = struct.unpack("<HH", hdr[26:30])
        fp.seek(info.header_offset + 30 + nlen + xlen)
        return _zstd_decompress(fp.read(info.compress_size), info.file_size)

    def brep_names(self):
        return [n for n in self.names() if "/Breps.BlobParts/" in n]

    def close(self):
        self.zip.close()
