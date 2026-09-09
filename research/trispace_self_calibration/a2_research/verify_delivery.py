"""Bounded acceptance checks for delivery; never reruns frozen final studies."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RESULT=HERE/"results"

def main():
    env=os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1",OMP_NUM_THREADS="1")
    tasks=[("physical_core","research/delegated/a2_physics/test_physics.py"),
           ("solver_checks","research/delegated/a2_solver_audited/test_solver.py"),
           ("parent_solver_checks","research/delegated/a2_solver_audited/test_parent_audit.py"),
           ("passivity","research/trispace_self_calibration/a2_research/test_passivity.py"),
           ("larger_map","research/delegated/a2_highdim/test_highdim.py")]
    rep="research/delegated/a2_theory_validation/replication"
    tasks += [("theory_"+n,f"{rep}/tests/test_{n}.py") for n in ("e1","e2","e3","e5","physical_runner")]
    tasks += [("final_parent","research/trispace_self_calibration/a2_research/final_checks.py")]
    checks=[]
    for name,rel in tasks:
        t=time.perf_counter()
        p=subprocess.run([sys.executable,str(ROOT/rel)],cwd=ROOT,env=env,text=True,capture_output=True,timeout=120)
        log=RESULT/("delivery_check_"+name+".txt")
        log.write_text(p.stdout+p.stderr)
        checks.append(dict(name=name,command=[sys.executable,str(ROOT/rel)],exit_code=p.returncode,wall_seconds=time.perf_counter()-t,log=str(log)))
        print(name,p.returncode,flush=True)
    frozen_dir=ROOT/"research/delegated/a2_solver_audited"
    frozen=json.loads((frozen_dir/"frozen_parent.json").read_text())
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(frozen_dir.glob("*.py"))}
    assert hashes==frozen["code_hash"],"Frozen code drift"
    hd=[json.loads(l) for l in (ROOT/"research/delegated/a2_highdim/records.jsonl").read_text().splitlines() if l.strip()]
    test=[r for r in hd if r["phase"]=="test"]
    assert len(test)==48
    assert len({(r["seed"],r["radius_index"],r["method"]) for r in test})==48
    counts={phase:sum(r["phase"]==phase for r in hd) for phase in set(r["phase"] for r in hd)}
    manifest={}
    patterns=["Theory/Questions/A2*","research/delegated/a2_physics/*.py","research/delegated/a2_solver_audited/*.py","research/delegated/a2_solver_audited/frozen_parent.json","research/delegated/a2_solver_audited/final_*.jsonl","research/delegated/a2_highdim/*.py","research/delegated/a2_highdim/records.jsonl","research/delegated/a2_highdim/frozen.json","research/trispace_self_calibration/a2_research/*.py","research/trispace_self_calibration/a2_research/*.md","communication/A2_RESEARCH_REPORT_ZH.md"]
    for pattern in patterns:
        for p in ROOT.glob(pattern):
            if p.is_file():manifest[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    result=dict(status="pass" if all(c["exit_code"]==0 for c in checks) else "failed",checks=checks,frozen_code_unchanged=True,larger_map_counts=counts,sha256=manifest,scope="implementation and artifact integrity, not journal/scientific acceptance")
    (RESULT/"delivery_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    if result["status"]!="pass":raise SystemExit(1)
    print("PASS: bounded checks, freeze integrity, larger-map pairing",flush=True)

if __name__=="__main__": main()
