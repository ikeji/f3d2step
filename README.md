# f3d2step

English | [日本語](README.ja.md)

A command-line tool that converts Autodesk Fusion (Fusion 360) `.f3d` files to STEP (AP214) without Fusion.
It needs nothing but the Python standard library (Python 3.8 or later).

A `.f3d` file is a zip archive whose geometry is stored as binary B-rep data of Autodesk ShapeManager
(ASM, an ACIS derivative). This tool reads that data directly and rewrites the topology (faces, edges,
vertices) and the surface/curve definitions as a STEP B-rep. Nothing is tessellated: planes, cylinders,
cones, spheres and tori stay analytic surfaces, and free-form surfaces come out as NURBS.

## Usage

```sh
python3 -m f3d2step part.f3d                 # writes part.step next to the input
python3 -m f3d2step part.f3d -o out.step     # choose the output file
python3 -m f3d2step *.f3d -o outdir/         # convert many files

python3 -m f3d2step -l part.f3d              # list the bodies in the file
python3 -m f3d2step part.f3d -b 2,5-7        # export only the listed body numbers
python3 -m f3d2step part.f3d --split         # one file per body: part-<body name>.step
```

After `pip install .` the tool is also available as the `f3d2step` command.

| Option | Meaning |
| --- | --- |
| `-o, --output` | Output file (single input) or directory |
| `-l, --list` | List the bodies (face count, size, position) instead of converting |
| `-b, --bodies LIST` | Body numbers to export, e.g. `1,3-5` |
| `--split` | One file per body, named `<name>-<body name>.step` (a number is added to tell namesakes apart) |
| `--include-hidden` | Also export bodies that are hidden in Fusion |
| `--default-name PATTERN` | Name for bodies that keep Fusion's default name; use `"ボディ%d"` for designs saved from a Japanese UI (default: `Body%d`) |
| `--solids-only` | Leave out surface bodies (sheets without thickness) |
| `--keep-duplicates` | Also export bodies stored twice with identical shape and position |
| `--all-blobs` | Also read the auxiliary data (`.smb`: sketch profiles and the like) |
| `-q, --quiet` | Suppress warnings |

Exit status: 0 = success, 1 = some file could not be converted, 2 = some faces could not be converted and were skipped.

Files written by recent versions of Fusion are Zstandard-compressed. Python 3.14 and later read them as is;
older Pythons need either the `zstd` command or the `zstandard` package (`pip install zstandard`).

## Limitations

**Deciding which bodies are the "final shape" works only partially.**

The geometry data in a `.f3d` also keeps bodies that were consumed or removed along the timeline, and hidden bodies.
The tool reads the visibility of body records in the design stream (`BulkStream.dat`) and matches them to the
bodies in the geometry data, but it does not decode the stream as a whole, so some bodies cannot be matched.
The `-l` listing shows the state of each body:

| State | Meaning | Default |
| --- | --- | --- |
| `visible` | Matched to a record and shown in Fusion | exported |
| `hidden` | Matched to a record and hidden | not exported (`--include-hidden` exports it) |
| `?` | No record could be matched (e.g. bodies of inserted components) | exported |

- For single-body parts, face counts and volumes match Fusion's own STEP export.
- `?` bodies may include bodies that are not visible in Fusion. Check the `-l` listing and pick bodies with `-b`.
- Component placement (occurrence transforms) and instancing are not applied. Every body is written in its
  component's local coordinates, so assemblies put together by moving or jointing components come out displaced,
  and a component placed several times appears only once. Per-occurrence visibility is not read either.

Also:

- Body names are used only when they can be recovered from the design stream (shown as `"name"` in the `-l` listing).
  Bodies that were never renamed get a default name in the language of the Fusion UI, so set `--default-name` to match.
  Bodies without a recovered name are called `Body1`, `Body2`, ... (the listing numbers). Colours and materials are not exported.
- Sketches, construction planes, mesh bodies and T-spline editing data are out of scope.
- Free-form surfaces use the NURBS approximations ASM stores. Of the surface kinds without a stored approximation,
  only helical surfaces (threads, coils) and extrusions are rebuilt analytically. Any other kind makes the tool
  skip that face with a warning (the body is then written as a surface model instead of a solid).
- Tapered helices (conical coils) and elliptical cones are not supported.
- Output units are always millimetres.

## Validation

About 90 publicly available `.f3d` files from GitHub (saved by Fusion versions from 2021 to 2026, up to more than
10,000 faces) were converted and the output was loaded into OpenCascade for a shape validity check (BRepCheck).
For parts published together with Fusion's own STEP export, face count, volume, surface area and bounding box were compared.

## How it works

| File | Role |
| --- | --- |
| `f3d2step/container.py` | Extracts the B-rep blobs from the zip (deflate / zstd) |
| `f3d2step/sab.py` | Tokenizer for the ASM binary format (SAB) |
| `f3d2step/model.py` | Interprets the records as topology and geometry; extracts NURBS approximations, builds helical shapes |
| `f3d2step/nurbs.py` | Moves the seams of periodic splines |
| `f3d2step/step.py` | STEP writer |
| `f3d2step/names.py` | Recovers body names from the design stream |
| `f3d2step/convert.py` | Body enumeration, selection and the conversion entry point |

In designs with a history, the bodies of each component live in `Breps.BlobParts/*.smbh` (`.smb` holds sketch profiles
and the like). Designs without a history, and files in the newer format, keep the bodies in `.smb`. Lengths are stored
in centimetres and converted to millimetres on output.

Tests: `python3 -m unittest discover -s tests -t .`

## Acknowledgements

The analysis of the SAB record structure in [InventorLoader](https://github.com/jmplonka/InventorLoader) (an Inventor/Fusion
importer for FreeCAD) was used as a reference while working out the format. No code was taken from it.
