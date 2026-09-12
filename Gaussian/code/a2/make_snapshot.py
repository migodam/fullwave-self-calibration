"""Post-run audit snapshot. Never backdates an execution manifest."""
from pathlib import Path
import datetime,hashlib,importlib.metadata,json,platform,sys
ROOT=Path(__file__).resolve().parents[2]
def main():
    files={}
    for folder in ['code/a2','runs/a2']:
        for p in (ROOT/folder).rglob('*'):
            if not p.is_file() or p.suffix not in {'.py','.json','.npz'}:continue
            if '__pycache__' in p.parts or p.name in {'postrun_snapshot.json','package_validation.json'}:continue
            if 'invalid_v1' in p.parts or 'portable_recheck' in p.parts:continue
            files[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'kind':'POST-RUN snapshot, not an execution-source attestation or preregistration',
         'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),
         'dependencies':{k:importlib.metadata.version(k) for k in ['numpy','scipy','matplotlib']},
         'lineage_gaps':['Measured primary executed source hash unavailable','Many earlier runs freeze runner but not all imported dependency bytes','Original Pro attachments not obtained'],
         'files_sha256':files}
    (ROOT/'runs/a2/postrun_snapshot.json').write_text(json.dumps(out,indent=2));print(json.dumps({'files':len(files),'dependencies':out['dependencies']}))
if __name__=='__main__':main()
