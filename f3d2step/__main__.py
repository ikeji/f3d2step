import argparse
import os
import sys

from . import sab
from .container import ContainerError
from .convert import convert, scan
from .model import Unsupported


def _list(path, all_blobs, default_name):
    blob = None
    for info in scan(path, all_blobs, lambda w: print("warning: %s" % w, file=sys.stderr), default_name):
        if info.blob != blob:
            blob = info.blob
            h = info.model.header
            print("  %s  (%s, %s)" % (blob, h.product_version, h.date))
        detail = ", ".join("%s x%d" % kv for kv in sorted(info.kinds.items()))
        size = ""
        if info.lo is not None:
            size = "  %.2f x %.2f x %.2f mm" % tuple(info.hi[j] - info.lo[j] for j in range(3))
            size += " at (%.2f, %.2f, %.2f)" % info.lo
        dup = "  [same as Body%d]" % info.duplicate_of if info.duplicate_of else ""
        kind = "solid  " if info.solid else "surface"
        label = ' "%s"' % info.label if info.label else ""
        state = {True: "visible", False: "hidden ", None: "?      "}[info.visible]
        print("    %-8s %s %s %4d faces (%s)%s%s%s" % (info.name, kind, state, info.faces, detail, size, label, dup))


def _parse_bodies(text):
    """'1,3-5' -> {1, 3, 4, 5}"""
    out = set()
    for part in text.replace("Body", "").replace("body", "").split(","):
        part = part.strip()
        if not part:
            continue
        a, _, b = part.partition("-")
        try:
            lo, hi = int(a), int(b or a)
        except ValueError:
            raise argparse.ArgumentTypeError("bad body list: %r" % text)
        if lo < 1 or hi < lo:
            raise argparse.ArgumentTypeError("bad body range: %r" % part)
        out.update(range(lo, hi + 1))
    if not out:
        raise argparse.ArgumentTypeError("empty body list")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="f3d2step", description="Convert Autodesk Fusion .f3d files to STEP without Fusion."
    )
    ap.add_argument("input", nargs="+", help=".f3d file(s)")
    ap.add_argument("-o", "--output", help="output .step file (single input) or directory")
    ap.add_argument("-l", "--list", action="store_true", help="list the bodies found instead of converting")
    ap.add_argument(
        "-b", "--bodies", type=_parse_bodies, metavar="LIST", help="export only these bodies, e.g. 1,3-5 (see --list)"
    )
    ap.add_argument("--split", action="store_true", help="write one STEP file per body (<name>-<body name>.step)")
    ap.add_argument(
        "--default-name",
        default="Body%d",
        metavar="PATTERN",
        help='name for bodies left with Fusion\'s default name, e.g. "ボディ%%d" for a Japanese UI (default: Body%%d)',
    )
    ap.add_argument("--include-hidden", action="store_true", help="also export bodies that are hidden in Fusion")
    ap.add_argument("--solids-only", action="store_true", help="leave out surface (sheet) bodies")
    ap.add_argument("--all-blobs", action="store_true", help="also read auxiliary .smb blobs (sketch profiles etc.)")
    ap.add_argument("--keep-duplicates", action="store_true", help="do not drop bodies that are stored twice")
    ap.add_argument("-q", "--quiet", action="store_true", help="suppress warnings")
    args = ap.parse_args(argv)

    out_is_dir = args.output is not None and (os.path.isdir(args.output) or len(args.input) > 1)
    status = 0
    for src in args.input:
        try:
            if args.list:
                print("%s:" % src)
                _list(src, args.all_blobs, args.default_name)
                continue
            base = os.path.splitext(os.path.basename(src))[0] + ".step"
            if args.output is None:
                dst = os.path.join(os.path.dirname(src), base)
            elif out_is_dir:
                os.makedirs(args.output, exist_ok=True)
                dst = os.path.join(args.output, base)
            else:
                dst = args.output
            res = convert(
                src,
                dst,
                all_blobs=args.all_blobs,
                keep_duplicates=args.keep_duplicates,
                select=args.bodies,
                split=args.split,
                solids_only=args.solids_only,
                include_hidden=args.include_hidden,
                default_name=args.default_name,
            )
            if not args.quiet:
                for w in res.warnings:
                    print("warning: %s" % w, file=sys.stderr)
            target = dst if len(res.outputs) == 1 and res.outputs[0] == dst else "%d files" % len(res.outputs)
            msg = "%s -> %s: %d solid(s), %d sheet(s), %d faces" % (src, target, res.solids, res.sheets, res.faces)
            if res.skipped_faces:
                msg += ", %d faces SKIPPED" % res.skipped_faces
                status = 2
            if res.duplicates:
                msg += ", %d duplicate body(ies) dropped" % res.duplicates
            if res.hidden:
                msg += ", %d hidden body(ies) skipped" % res.hidden
            print(msg)
        except (OSError, ContainerError, Unsupported, sab.SabError) as ex:
            print("error: %s: %s" % (src, ex), file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
