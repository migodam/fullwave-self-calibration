"""Post-run metadata repair; never represents post-hoc hashes as pre-run."""
import json,hashlib,shutil,sys,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
meas=ROOT/'runs/a2/measured/results.json';x=json.loads(meas.read_text());x['scope']='One published 2-D TM cylinder gain-profiled diagnostic. Per-frequency gain is an extra scattered-training nuisance after incident estimation; no source normalization, hardware calibration, material recovery, 3-D, or diverse-target claim.';x['executed_primary_source_hash']=None;x['postrun_metadata_repair']=True;meas.write_text(json.dumps(x,indent=2))
old=ROOT/'runs/a2/ports/ports_n24_regression.json';legacy=ROOT/'runs/a2/ports/legacy_unmanifested';legacy.mkdir(exist_ok=True)
if old.exists(): shutil.move(str(old),str(legacy/'ports_n24_regression.json'))
