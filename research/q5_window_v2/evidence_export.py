"""Compact review tables; complete raw observations are preserved separately."""
import csv,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
RESULTS=HERE/'results'

def run():
    dest=RESULTS/'REVIEW_RESULTS.json'
    if dest.exists():raise FileExistsError(dest)
    recovery=json.loads((RESULTS/'recovery_v2/summary.json').read_text())
    modal=json.loads((RESULTS/'modal_certificate.json').read_text())
    readout=json.loads((RESULTS/'readout_certificate.json').read_text())
    geo=json.loads((RESULTS/'geometry_window_certificate.json').read_text())
    oracle={}
    scenes=[json.loads(p.read_text()) for p in sorted((RESULTS/'recovery_v2').glob('scene_*.json'))]
    for name in scenes[0]['oracle_diagnostics']:
        values=[s['oracle_diagnostics'][name] for s in scenes]
        oracle[name]={'successes':sum(v.get('material_success',False) for v in values),'attempts':len(values)}
    review={'evidence_status':'restricted_certificates_and_development_diagnostics_not_TAP_ready',
      'modal_status':modal['status'],'total_nominal_epsilon_boxes':sum(r['coverage']['count'] for r in modal['rows']),
      'modal_bounds':[{k:v['outward_decimal'] for k,v in r['bounds'].items()} for r in modal['rows']],
      'radial_geometry_nominal_kR':[.6,1.6],
      'geometry_noise_failure_upper':geo['noise_event_failure_upper']['outward_decimal'],
      'readout_material_error_bounds':[{'sphere':r['sphere'],'repetitions':r['repetitions'],
           'bound':r['quantities']['material_error_bound']['outward_decimal']} for r in readout['budget_rows']],
      'modal_path_competition_eta':readout['competing_worlds']['eta_required']['outward_decimal'],
      'independent_recovery_summary':recovery['summary'],'oracle_diagnostics':oracle,
      'raw_sha256':recovery['raw_sha256'],
      'unresolved':['Continuum Maxwell error bound for class C','Covered finite class C separation',
        'Physical probe/geometry/relative-path error budget','Complete novelty clearance','New method superiority'],
      'evidence_location':'Original complete raw data and all starts are in the delivered archive; reproduce.py regenerates the same frozen stream in a separate directory; this is not extra statistical evidence.'}
    dest.write_text(json.dumps(review,indent=2)+'\n')
    rows=[];starts=[]
    for p in sorted((RESULTS/'recovery_v2').glob('scene_*.json')):
        scene=json.loads(p.read_text())
        for group in ['methods','oracle_diagnostics']:
            for name,r in scene[group].items():
                rows.append({'scene':scene['scene'],'group':group,'method':name,
                  'epsilon1_error':r.get('material_errors',[None,None])[0],
                  'epsilon2_error':r.get('material_errors',[None,None])[1],
                  'shift_error_mm':r.get('geometry_error_mm'),'gain_error':r.get('gain_error'),
                  'success':r.get('material_success',False),'cost':r.get('cost'),
                  'evaluations':r.get('forward_evaluations'),'failure':r.get('failure','')})
                for j,s in enumerate(r['starts']):
                    starts.append({'scene':scene['scene'],'group':group,'method':name,'start':j,
                      'record_json':json.dumps(s,separators=(',',':'))})
    for name,items in [('scene_summary.csv',rows),('all_starts.csv',starts)]:
        path=RESULTS/name
        if path.exists():raise FileExistsError(path)
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(items[0]));writer.writeheader();writer.writerows(items)
    hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
      for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
      and p.name not in ['MANIFEST.sha256.json']}
    (RESULTS/'MANIFEST.sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(json.dumps(review,indent=2))
if __name__=='__main__':run()
