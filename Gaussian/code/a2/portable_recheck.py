"""Re-execute quick checks in an extracted delivery copy; retain long-run evidence."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time,zipfile
ROOT=Path(__file__).resolve().parents[2]
def main():
    archive=ROOT/'deliverables/Gaussian_A1_A1_2_A2_Pro_Package.zip'
    target=ROOT/'delegated/a2_packaging'/('portable_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    target.mkdir(parents=True,exist_ok=False)
    with zipfile.ZipFile(archive) as z:z.extractall(target)
    copied=target/'Gaussian'
    records=[]
    checks=['check_certificate.py','check_parameter_geometry.py','check_quadrature_certificate.py','check_matrix_free.py','vector_certificate_probe.py','extensions/theory_checks.py','extensions/check_manifold_native.py','extensions/check_rbf_capacity.py']
    for name in checks:
        p=copied/'code/a2'/name;t=time.perf_counter()
        run=subprocess.run([sys.executable,str(p)],cwd=target,capture_output=True,text=True,timeout=60,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
        records.append({'script':name,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'exit_code':run.returncode,'seconds':time.perf_counter()-t,'stdout':run.stdout[-2000:],'stderr':run.stderr[-2000:]})
    result={'input_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
            'extraction_directory':str(target.relative_to(ROOT)),
            'scope':'Eight small numerical checks in extracted copy; no long imaging, measured or training rerun.',
            'all_passed':all(r['exit_code']==0 for r in records),'records':records,
            'tested_python_sources':{p.relative_to(copied).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (copied/'code').rglob('*.py')}}
    out=ROOT/'runs/a2/portable_recheck';out.mkdir(parents=True,exist_ok=True)
    (out/'results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result));assert result['all_passed']
if __name__=='__main__':main()
