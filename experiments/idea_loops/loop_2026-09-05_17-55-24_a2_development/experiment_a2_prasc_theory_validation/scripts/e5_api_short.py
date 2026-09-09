#!/usr/bin/env python
"""Short API inventory for a2val/e5.py (AST-based, no code bodies).

Output (results/e5_api_short.txt):
  * one line per top-level function: line | signature | first docstring
    sentence (<= 120 chars, "(no docstring)" if absent);
  * module-level names assigned values that are not functions/classes/modules,
    with value repr truncated to 200 chars.
"""

from __future__ import annotations

import ast
import re
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
E5_PATH = ROOT / "a2val" / "e5.py"
OUT_PATH = ROOT / "results" / "e5_api_short.txt"
MAX_DOC = 120
MAX_REPR = 200


def _sig(node: ast.FunctionDef) -> str:
    prefix = "async def " if isinstance(node, ast.AsyncFunctionDef) else "def "
    params = ast.unparse(node.args)
    suffix = f" -> {ast.unparse(node.returns)}" if node.returns else ""
    return f"{prefix}{node.name}({params}){suffix}:"


def _first_sentence(doc: str | None) -> str:
    if not doc:
        return "(no docstring)"
    flat = re.sub(r"\s+", " ", doc).strip()
    m = re.search(r"[.!?](?=\s|$)", flat)
    sentence = flat[: m.end()] if m else flat
    sentence = sentence.strip()
    if len(sentence) > MAX_DOC:
        sentence = sentence[: MAX_DOC - 3].rstrip() + "..."
    return sentence


def main() -> int:
    source = E5_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(E5_PATH))

    funcs: list[tuple[int, str, str]] = []
    assigned_names: dict[str, int] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs.append(
                (
                    node.lineno,
                    _sig(node),
                    _first_sentence(ast.get_docstring(node, clean=True)),
                )
            )
            continue
        targets: list[ast.expr] = []
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target]
        for target in targets:
            for name_node in ast.walk(target):
                if isinstance(name_node, ast.Name):
                    assigned_names.setdefault(name_node.id, node.lineno)

    lines = [
        "E5 API short inventory - a2val/e5.py",
        f"Module: {E5_PATH}",
        f"Top-level defs: {len(funcs)}",
        "Fields: lineno | signature | first docstring sentence (<=120 chars)",
        "",
    ]
    for lineno, sig, doc in funcs:
        lines.append(f"{lineno:4d} | {sig} | {doc}")

    lines += [
        "",
        "Module-level non-function variables (names and repr, repr <=200 chars):",
    ]
    sys.path.insert(0, str(ROOT))
    import a2val.e5 as e5_mod  # import after sys.path setup

    found_any = False
    for name in sorted(assigned_names):
        try:
            value = getattr(e5_mod, name)
        except Exception as exc:  # pragma: no cover - defensive
            lines.append(f"  {name} = <unreadable: {exc}>")
            found_any = True
            continue
        if isinstance(
            value,
            (types.FunctionType, types.ModuleType, type, types.MethodType),
        ):
            continue
        found_any = True
        try:
            rep = repr(value)
        except Exception as exc:  # pragma: no cover - defensive
            rep = f"<repr failed: {exc}>"
        if len(rep) > MAX_REPR:
            rep = rep[: MAX_REPR - 15] + "...[truncated]"
        lines.append(f"  {name} = {rep}  (line {assigned_names[name]})")
    if not found_any:
        lines.append("  (none)")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
