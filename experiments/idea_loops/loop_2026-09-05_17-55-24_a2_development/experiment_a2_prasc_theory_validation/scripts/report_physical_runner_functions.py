#!/usr/bin/env python3
"""Report top-level functions of a2val/physical_runner.py (AST-only, no runs).

Reads the module source directly with the stdlib ``ast`` module so no
experiment code or heavy imports are executed.  Emits:

* results/physical_runner_functions.txt  -- list + full E3-related bodies
* stdout                                  -- full report (list + bodies)

The E3 body set is the module's ``# E3 helpers on physical tangents`` section
(theorem3_physical, prior_sandwich_physical, theorem7_physical) plus
split_frequency_frames, which was explicitly requested even though the module
labels its section ``# E5 helpers on physical frequency frames``.
"""

from __future__ import annotations

import ast
import datetime as _dt
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "a2val" / "physical_runner.py"
OUT = ROOT / "results" / "physical_runner_functions.txt"

# Explicitly requested E3/related full-source functions.
E3_FULL_BODIES = (
    "theorem3_physical",
    "prior_sandwich_physical",
    "theorem7_physical",
    "split_frequency_frames",
)


def first_docstring_sentence(node: ast.FunctionDef) -> str:
    """Return the first sentence of the function docstring, or '' if absent."""
    doc = ast.get_docstring(node, clean=True)
    if not doc:
        return ""
    doc = doc.strip()
    # Stop at the first period followed by whitespace or end of the paragraph.
    m = re.search(r"\.(?=\s|$)", doc)
    if m:
        return doc[: m.end()].strip()
    # No sentence terminator: use the first paragraph.
    para = doc.split("\n\n", 1)[0].strip()
    return para


def source_signature(source: str, node: ast.FunctionDef) -> str:
    """Extract the exact ``def ... :`` header (bracket/quote aware)."""
    lines = source.splitlines()
    depth = 0
    quote: str | None = None
    out: list[str] = []
    for line_no in range(node.lineno - 1, min(node.end_lineno, len(lines))):
        line = lines[line_no]
        i = 0
        while i < len(line):
            ch = line[i]
            if quote is not None:
                out.append(ch)
                if ch == "\\" and i + 1 < len(line):
                    out.append(line[i + 1])
                    i += 2
                    continue
                if line.startswith(quote, i):
                    triple = (quote * 3) == line[i : i + 3]
                    if triple:
                        out.append(quote * 2)
                        i += 3
                        quote = None
                        continue
                    if not triple:
                        quote = None
                i += 1
                continue
            if ch in "\"'":
                if line.startswith(ch * 3, i):
                    quote = ch * 3
                    out.append(ch * 3)
                    i += 3
                else:
                    quote = ch
                    out.append(ch)
                    i += 1
                continue
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == ":" and depth == 0:
                out.append(":")
                return "".join(out).strip()
            elif ch == "#" and depth == 0:
                break
            out.append(ch)
            i += 1
        out.append("\n")
    return "".join(out).strip()


def main() -> None:
    source = MODULE.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(MODULE))
    funcs: list[tuple[ast.FunctionDef, str, str, bool]] = []
    for item in tree.body:
        if isinstance(item, ast.FunctionDef):
            sig = source_signature(source, item)
            doc = first_docstring_sentence(item)
            full = ast.get_source_segment(source, item) or ""
            funcs.append((item, sig, doc, item.name in E3_FULL_BODIES))

    lines: list[str] = []
    lines.append("physical_runner.py top-level function report")
    lines.append("file: %s" % MODULE)
    lines.append("source lines: %d" % len(source.splitlines()))
    lines.append("generated: %s" % _dt.datetime.now().astimezone().isoformat(
        timespec="seconds"
    ))
    lines.append("method: AST source introspection only; no experiment executed")
    lines.append("")
    lines.append("TOP-LEVEL FUNCTIONS (%d)" % len(funcs))
    lines.append("-" * 100)
    for i, (node, sig, doc, _) in enumerate(funcs, start=1):
        lines.append(
            "%2d. L%-4d %s" % (i, node.lineno, sig)
        )
        lines.append("     First docstring sentence: %s" % (
            doc if doc else "(no docstring)"
        ))
        lines.append("")

    lines.append("-" * 100)
    lines.append(
        "E3-RELATED FULL SOURCE BODIES "
        "(module '# E3 helpers on physical tangents' section + "
        "explicitly requested split_frequency_frames)"
    )
    lines.append("-" * 100)
    for node, sig, _doc, is_e3 in funcs:
        if not is_e3:
            continue
        body = ast.get_source_segment(source, node) or ""
        lines.append("")
        lines.append(
            "### %s  (lines %d-%d)" % (node.name, node.lineno, node.end_lineno)
        )
        lines.append("```python")
        lines.append(body)
        lines.append("```")
    lines.append("")
    lines.append(
        "NOTE: remaining top-level functions are E2 helpers "
        "(declared_E_basis, t5_adjacent_pair), module utilities, and E5 "
        "frequency-frame helpers (stack_vs_sum_physical, "
        "sequential_innovation, rank_budget_physical); only "
        "split_frequency_frames was additionally requested from that section."
    )

    report = "\n".join(lines) + "\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(report, encoding="utf-8")
    print(report, end="")
    print("[report written to %s]" % OUT)


if __name__ == "__main__":
    main()
