"""Parent continuation: independent high-frequency convergence, not timing test."""
import sys
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT.parents[1] / "delegated" / "a3_maxwell_fft"))
from run_maxwell_fft import audit_case, SCENES, receivers

if __name__ == "__main__":
    dest = OUT / "results" / "maxwell3d_fine_forward.json"
    rows = json.loads(dest.read_text()) if dest.exists() else []
    for k, h in [(18., .015), (18., .012), (18., .01), (18., .0075), (9., .01)]:
        if any(r["k"] == k and r["spacing"] == h for r in rows):
            continue
        row, stop = audit_case(SCENES[1], k, h, receivers())
        rows.append(row)
        dest.write_text(json.dumps(rows, indent=2) + "\n")
        print(json.dumps({key: row.get(key) for key in ["k", "spacing", "status", "voxels", "relative_field_error_high", "solve_wall_seconds", "memory_estimate_mib"]}), flush=True)
        if stop:
            break
