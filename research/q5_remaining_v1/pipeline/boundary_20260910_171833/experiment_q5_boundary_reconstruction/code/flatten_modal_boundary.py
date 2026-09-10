#!/usr/bin/env python3
"""Flatten every leaf of inputs/modal_boundary.json to 'path = <json>' lines.

Leaves are scalars (str/int/float/bool/None). Nested dicts use '.' separators;
list elements are addressed with '[i]'. Output is sorted by path.
inputs/ is only read.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SRC = HERE / "inputs" / "modal_boundary.json"
OUT = HERE / "inputs" / "derive" / "MODAL_BOUNDARY_KEYS.txt"


def leaves(node, prefix=""):
    if isinstance(node, dict):
        if not node:
            yield prefix or "$", "{}"
            return
        for k in node:
            p = f"{prefix}.{k}" if prefix else str(k)
            yield from leaves(node[k], p)
    elif isinstance(node, list):
        if not node:
            yield f"{prefix}[]", "[]"
            return
        for i, v in enumerate(node):
            yield from leaves(v, f"{prefix}[{i}]")
    else:
        yield prefix, json.dumps(node, ensure_ascii=False)


def describe(node):
    if isinstance(node, dict):
        return f"dict(size={len(node)})"
    if isinstance(node, list):
        return f"list(size={len(node)})"
    return f"{type(node).__name__} value={json.dumps(node, ensure_ascii=False)}"


def main() -> int:
    with SRC.open("rb") as fh:
        data = json.load(fh)

    pairs = sorted(leaves(data), key=lambda kv: kv[0])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8") as fh:
        fh.write(f"# flattened leaves of inputs/modal_boundary.json  (root type: {type(data).__name__})\n")
        fh.write(f"# leaf_count: {len(pairs)}\n")
        fh.write("# top-level keys:\n")
        for k, v in data.items():
            fh.write(f"#   {k} : {describe(v)}\n")
        fh.write("# format: <leaf path> = <json value>\n")
        if isinstance(data, list):
            fh.write("# note: root is a list; elements are addressed as [0], [1], ...\n")
        for path, val in pairs:
            fh.write(f"{path} = {val}\n")

    print(f"[modal_boundary] root type: {type(data).__name__}")
    print(f"[modal_boundary] leaf count: {len(pairs)}")
    print("[modal_boundary] top-level keys:")
    for k, v in data.items():
        print(f"  {k} : {describe(v)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
