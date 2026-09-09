"""Run NEW A4 cases in a previously nonexistent directory.

Stored frozen evidence is never overwritten. This wrapper is post-freeze
packaging, not a new scientific method or an externally preregistered test.
"""
from pathlib import Path
import argparse
import os
import json
import hashlib
import sys
from datetime import datetime, timezone
from time import perf_counter

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'src'))

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--branch-per-kind', type=int, default=1)
    p.add_argument('--calibration-scenes', type=int, default=2)
    p.add_argument('--max-cells', type=int, default=6000)
    p.add_argument('--include-dda', action='store_true')
    args = p.parse_args()
    if min(args.branch_per_kind, args.calibration_scenes, args.max_cells) < 1:
        p.error('scene counts and max-cells must be positive')
    dest = args.output.expanduser().resolve()
    if dest.exists():
        p.error('output already exists; refusing to overwrite: '+str(dest))
    if dest == ROOT or ROOT in dest.parents:
        p.error('choose a destination outside this immutable evidence directory')
    dest.mkdir(parents=True, exist_ok=False)
    protocol = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'kind': 'fresh local rerun; not original frozen evidence',
                'seed': args.seed, 'branch_per_kind': args.branch_per_kind,
                'calibration_scenes': args.calibration_scenes, 'max_cells': args.max_cells,
                'include_dda': args.include_dda,
                'source_hashes': {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                                  for f in sorted((ROOT/'src').glob('*.py'))}}
    (dest/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    import numpy as np
    from experiments import branch_experiment, calibration_experiment, spectral_experiment, dda_experiment, jsonable
    def write(name, data):
        (dest/name).write_text(json.dumps(data, default=jsonable, indent=2, allow_nan=False)+'\n')
    start = perf_counter()
    try:
        write('spectral.json', spectral_experiment(args.seed))
        br, raw = branch_experiment(args.seed+1, args.branch_per_kind, max_cells=args.max_cells)
        write('branch.json', br); np.savez_compressed(dest/'branch_data.npz', **raw)
        cal, raw = calibration_experiment(args.seed+2, args.calibration_scenes, 32)
        write('calibration.json', cal); np.savez_compressed(dest/'calibration_data.npz', **raw)
        if args.include_dda:
            write('dda.json', dda_experiment())
        write('completed.json', {'wall_seconds': perf_counter()-start, 'original_results_untouched': True})
    except Exception:
        import traceback
        (dest/'FAILED.txt').write_text(traceback.format_exc())
        raise
    print(dest)

if __name__ == '__main__':
    main()
