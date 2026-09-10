#!/usr/bin/env python3
r"""Split inputs/A5_THEORY.md into one file per numbered section.

Split points: lines matching '^## (\d+)\. ' (heading line is kept in the file).
Text before the first heading goes to A5_sec_00_preamble.md.
Each part is a byte-exact slice of the source; a single trailing newline is
appended only if the slice did not already end with one. inputs/ is read-only.
"""
import hashlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SRC = HERE / "inputs" / "A5_THEORY.md"
OUTDIR = HERE / "inputs" / "derive"
SIZES = OUTDIR / "A5_SIZES.txt"
HEADING = re.compile(rb"(?m)^## (\d+)\. ")


def main() -> int:
    raw = SRC.read_bytes()
    matches = list(HEADING.finditer(raw))
    if not matches:
        print("FATAL: no '## N. ' headings found", file=sys.stderr)
        return 1

    bounds = []
    if matches[0].start() > 0:
        bounds.append(("A5_sec_00_preamble.md", 0, matches[0].start()))
    for i, m in enumerate(matches):
        num = int(m.group(1))
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(raw)
        bounds.append((f"A5_sec_{num:02d}.md", start, end))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, reassembled = [], bytearray()
    for name, start, end in bounds:
        chunk = raw[start:end]
        nl_added = not chunk.endswith(b"\n")
        if nl_added:
            chunk += b"\n"
        (OUTDIR / name).write_bytes(chunk)
        reassembled += raw[start:end]
        rows.append({
            "name": name, "bytes": len(chunk), "lines": len(chunk.splitlines()),
            "nl_added": nl_added, "sha256": hashlib.sha256(chunk).hexdigest(),
        })

    lines = ["# A5_THEORY.md section split - byte/line counts",
             f"# source: inputs/A5_THEORY.md bytes={len(raw)} lines={len(raw.splitlines())} "
             f"sha256={hashlib.sha256(raw).hexdigest()}",
             f"# sections written: {len(rows)}  (line counts are of the written file)",
             "# filename\tbytes\tlines\ttrailing_newline_added\tsha256"]
    for r in rows:
        lines.append(f"{r['name']}\t{r['bytes']}\t{r['lines']}\t{r['nl_added']}\t{r['sha256']}")
    lines.append(f"TOTAL\t{sum(r['bytes'] for r in rows)}\t{sum(r['lines'] for r in rows)}\t-\t-")
    SIZES.write_text("\n".join(lines) + "\n")

    ok = bytes(reassembled) == raw
    print(f"wrote {len(rows)} section files to {OUTDIR}")
    print(f"source bytes={len(raw)} lines={len(raw.splitlines())}")
    print(f"sum(part lines)={sum(r['lines'] for r in rows)}  parts_concat_equals_source={ok}")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
