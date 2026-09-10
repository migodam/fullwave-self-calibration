#!/usr/bin/env bash
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
DEST=${1:-$(mktemp -d)}
if [ -e "$DEST/research/q5_window_v2" ]; then
    echo 'Destination exists; choose a new directory. Existing evidence is immutable.' >&2
    exit 1
fi
mkdir -p "$DEST/research/q5_window_v2/tests" "$DEST/research/trispace_self_calibration/a3_research"
cp "$HERE"/*.py "$DEST/research/q5_window_v2/"
cp "$HERE/tests/test_q5.py" "$DEST/research/q5_window_v2/tests/"
cp "$ROOT/research/trispace_self_calibration/a3_research/maxwell3d.py" "$DEST/research/trispace_self_calibration/a3_research/"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$DEST/research/q5_window_v2"
mkdir results
for SCRIPT in interval_certificate readout_certificate run_experiment diagnostics finite_pair finite_windows; do
    python "$SCRIPT.py" > "results/$SCRIPT.log" 2>&1
done
python -m pytest -q tests > results/tests.log 2>&1
python export_summary.py > results/export.log
cat results/tests.log
printf 'Evidence directory: %s\n' "$DEST/research/q5_window_v2/results"
