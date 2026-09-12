"""Validate and summarize the post-hoc grid-refinement evidence."""
from pathlib import Path
import hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3];RUN=ROOT/'runs/a2/gpu/representation_campaign/main128'
def main():
    d=json.loads((RUN/'refined_field_posthoc.json').read_text());p=RUN/'gaussian_parameters_posthoc.json';par=json.loads(p.read_text())
    assert hashlib.sha256(p.read_bytes()).hexdigest()==d['parameter_json_sha256']
    assert hashlib.sha256((ROOT/'code/a2/gpu/refined_field_check.py').read_bytes()).hexdigest()==d['source_sha256']
    assert hashlib.sha256((ROOT/'code/a2/gpu/extract_parameters.py').read_bytes()).hexdigest()==par['source_sha256']
    for name,sha in d['dependencies_sha256'].items():
        loc=ROOT/'code/a2'/name if name=='physics.py' else ROOT/'code/a2/gpu'/name
        assert hashlib.sha256(loc.read_bytes()).hexdigest()==sha
    methods=['gaussian_K16','gaussian_K64','gaussian_K144','voxel'];assert {(x['case'],x['method']) for x in d['records']}=={(i,m) for i in range(30) for m in methods};assert len(d['records'])==120
    assert d['counts']['forward_rhs']==4320 and d['counts']['adjoint_rhs']==0
    out={'scope':d['scope'],'posthoc':True,'setup_seconds':d['setup_seconds'],'total_seconds':d['total_seconds'],'groups':{},'parameter_diagnostic':{}}
    for f in range(3):
        out['groups'][str(f)]={}
        for m in methods:
            r=[x for x in d['records'] if x['case']//10==f and x['method']==m]
            for x in r:
                assert np.isfinite([x[k] for k in ['held_clean_original_N128','held_clean_N256','material_relative_N256']]).all()
                old=json.loads((RUN/f"case_{x['case']:02d}_{m}.json").read_text());assert abs(old['held_clean_field_relative']-x['held_clean_original_N128'])<1e-12
            out['groups'][str(f)][m]={k:float(np.mean([x[k] for x in r])) for k in ['held_clean_original_N128','held_clean_N256','material_relative_N256']}
    for m in methods[:-1]:
        r=[x for x in par['records'] if x['method']==m];out['parameter_diagnostic'][m]={'active_subgrid_components':sum(x['active_subgrid_sigma_count'] for x in r),'models_with_active_subgrid':sum(x['active_subgrid_sigma_count']>0 for x in r),'min_sigma_m':min(x['min_principal_sigma_m'] for x in r)}
    (RUN/'REFINEMENT_SUMMARY.json').write_text(json.dumps(out,indent=2));print(json.dumps({'all_checks_passed':True,'records':120,'seconds':d['total_seconds']}))
if __name__=='__main__':main()
