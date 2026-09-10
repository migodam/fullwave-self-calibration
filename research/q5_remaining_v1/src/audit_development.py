"""Correctly whiten diagnostics and audit immutable development outcomes."""
from recovery_cycle1 import *


def main():
    target=HERE/'results/development_audit.json'
    if target.exists():raise RuntimeError('Preserve old audit')
    data=json.loads((HERE/'results/development_cycle1.json').read_text())
    second=json.loads((HERE/'results/development_cycle2.json').read_text())
    assert data['complete'] and second['complete']
    for record in (data,second):
        for path,digest in record['source_hashes'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert len(data['rows'])==24 and len(second['rows'])==12
    struct=mx.receivers(36,.15)
    dda=mx.DipoleVIE(CENTERS,RADII,.011,K,fill_quadrature=4)
    corrected=[]
    for row in data['rows']:
        truth=row['truth']; c={'scene':row['scene'],'data_class':row['data_class'],'methods':{}}
        if row['data_class']=='same_model':ss=forward([*truth[:2],0],struct,4)
        else:
            p=dda.currents(np.array(truth[:2])+1j*LOSS)
            ss=(mx.dipole_kernel(struct,dda.points,K)@p).ravel()
        for name,m in row['methods'].items():
            t=m['theta'];stat=2*m['rss'];dof=288+(2 if name=='joint_reference' else 0)-5
            passed=bool(stat<=chi2.ppf(.99,dof))
            c['methods'][name]={'whitened_rss':stat,'heuristic_residual_pass':passed,
                'wrong_heuristic_pass':passed and not m['task_success'],
                'structural_field_error_fixed_world':float(np.linalg.norm(forward([*t[:2],0],struct)-ss)/np.linalg.norm(ss))}
        corrected.append(c)
    summary={}
    for cls in ('same_model','independent_DDA'):
        rows=[r for r in data['rows'] if r['data_class']==cls]
        for name in ('no_reference','joint_reference'):
            ms=[r['methods'][name] for r in rows]
            cs=[r['methods'][name] for r in corrected if r['data_class']==cls]
            summary[cls+'/'+name]={'n':len(ms),'success':sum(m['task_success'] for m in ms),
                'material_median_errors':np.median([m['material_errors'] for m in ms],axis=0).tolist(),
                'geometry_median_mm':float(np.median([m['geometry_error_mm'] for m in ms])),
                'wrong_heuristic_pass':sum(c['wrong_heuristic_pass'] for c in cs),
                'heuristic_pass':sum(c['heuristic_residual_pass'] for c in cs),
                'structural_median_fixed_world':float(np.median([c['structural_field_error_fixed_world'] for c in cs]))}
    for name in ('selected_EM','random_EM','fixed_EM'):
        ms=[r['methods'][name] for r in second['rows']]
        summary['cycle2/'+name]={'n':len(ms),'success':sum(m['task_success'] for m in ms),
            'material_median_errors':np.median([m['material_errors'] for m in ms],axis=0).tolist(),
            'wrong_heuristic_pass':sum(m['heuristic_residual_pass'] and not m['task_success'] for m in ms)}
    chosen=[]
    for r in second['rows']:
        chosen.append(data['rows'][r['scene']]['methods']['joint_reference'] if r['policy_action']=='reference'
                      else r['methods']['selected_EM'])
    summary['cycle2/policy']={'n':len(chosen),'success':sum(m['task_success'] for m in chosen),
        'reference_choices':sum(r['policy_action']=='reference' for r in second['rows'])}
    result={'passed':True,'summary':summary,'cycle1_corrected':corrected,
            'cycle1_seconds':data['wall_seconds'],'cycle2_seconds':second['wall_seconds'],
            'corrections':['Cycle1 residual flags withdrawn: use 2*rss for chi-square.',
                           'Cycle1 structural error superseded: fixed physical points, not pose-shifted surface.'],
            'scientific_acceptance':'unresolved','final_test_generated':False}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cycle1_corrected'},indent=2))


if __name__=='__main__':main()
