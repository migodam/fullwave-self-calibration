"""Publish the prepared, scoped Gaussian portable snapshot into public_release."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reconcile-root", type=Path, required=True)
    args = parser.parse_args()

    dest_root = Path(__file__).resolve().parents[1]
    source_root = args.source.resolve()
    reconcile_root = args.reconcile_root.resolve()
    source_manifest = json.loads((source_root.parent / "MANIFEST.json").read_text())

    records = []
    for published_name, expected in sorted(source_manifest["files"].items()):
        source = source_root.parent / published_name
        if not source.is_file():
            raise SystemExit(f"missing source file: {source}")
        actual_hash = sha256(source)
        input_path = source
        origin = "Gaussian portable snapshot 2026-09-11"
        if actual_hash != expected["sha256"]:
            reconciled = reconcile_root / published_name
            if not reconciled.is_file():
                raise SystemExit(f"source hash mismatch: {published_name}")
            input_path = reconciled
            actual_hash = sha256(input_path)
            origin = "Gaussian portable snapshot 2026-09-11; runtime metadata reconciled from workspace"
        dest = dest_root / published_name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, dest)
        records.append({
            "path": published_name,
            "bytes": source.stat().st_size,
            "sha256": actual_hash,
            "origin": origin,
        })

    manifest_path = dest_root / "PUBLIC_MANIFEST.json"
    public_manifest = json.loads(manifest_path.read_text())
    public_manifest["records"] = [
        row for row in public_manifest["records"]
        if not row["path"].startswith("Gaussian/")
    ] + records
    public_manifest["gaussian_snapshot"] = {
        "source_snapshot": "Gaussian/delegated/a2_packaging/portable_20260911_234526",
        "reconciled_against": "Gaussian/ (only runtime timing metadata in 3 result files)",
        "scope": source_manifest["scope"],
        "excluded": source_manifest["excluded"],
        "files": len(records),
        "bytes": sum(row["bytes"] for row in records),
    }
    manifest_path.write_text(json.dumps(public_manifest, ensure_ascii=False, indent=2) + "\n")

    index_path = dest_root / "FILE_INDEX.md"
    index = index_path.read_text()
    marker = "\n## Gaussian\n"
    if marker in index:
        index = index.split(marker, 1)[0].rstrip() + "\n"
    index += "\n## Gaussian\n\n"
    index += "\n".join(
        f"- [{row['path']}](<{row['path']}>)" for row in records
    ) + "\n"
    index_path.write_text(index)
    print(json.dumps({
        "published_files": len(records),
        "published_bytes": sum(row["bytes"] for row in records),
        "manifest": str(manifest_path),
        "index": str(index_path),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
