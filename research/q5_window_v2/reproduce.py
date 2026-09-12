"""Reproduce in a NEW isolated directory; never mutate stored evidence.

Usage: python research/q5_window_v2/reproduce.py --out /tmp/q5-v2-check
The destination is a small repository-shaped tree containing the unchanged
inherited DDA dependency, all new source, and regenerated evidence.
"""
import argparse, os, shutil, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def run(out):
    out=out.resolve()
    if out.exists():raise FileExistsError(f'Refusing to overwrite {out}')
    target=out/'research/q5_window_v2'
    shutil.copytree(HERE,target,ignore=shutil.ignore_patterns('results','audit','__pycache__','*.pyc'))
    (target/'results').mkdir()
    dda=ROOT/'research/trispace_self_calibration/a3_research/maxwell3d.py'
    if not dda.is_file():raise FileNotFoundError(f'Inherited DDA dependency not found: {dda}')
    dest=out/dda.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dda,dest)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    commands=[['certify_modal.py'],['certify_modal.py','--size-relative-ppm','100','--out','results/modal_size_100ppm.json'],
      ['certify_readout.py'],['certify_geometry_window.py'],['geometry_readout_budget.py'],
      ['validate_solver.py'],['recovery_v2.py'],['finite_worlds.py'],['evidence_export.py'],
      ['-m','unittest','discover','-s','tests','-v']]
    for index,args in enumerate(commands):
        command=[sys.executable,*args];print('RUN',command,flush=True)
        log=target/'results'/f'reproduction_{index:02d}.log'
        with log.open('w') as stream:
            completed=subprocess.run(command,cwd=target,env=env,stdout=stream,stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f'Command failed with status {completed.returncode}; see {log}')
    print('Evidence stored at',target/'results')
    print('This reproduces registered diagnostics, not a new final test or hardware validation.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    run(parser.parse_args().out)
