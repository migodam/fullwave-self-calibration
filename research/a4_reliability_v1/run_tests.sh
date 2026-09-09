#!/bin/sh
set -eu
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m pytest -q
