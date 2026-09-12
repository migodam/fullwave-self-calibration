"""Post-run fixed-gain diagnostic from the selected gain-profiled states."""
from __future__ import annotations
import hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code'))
from a2.measured import MeasuredProblem,fit
OUT=ROOT/'runs/a2/measured'
def main():
    path=OUT/'results.json';r=json.loads(path.read_text());p=MeasuredProblem();fixed={}
    for name,count,kind in [('1_gaussian',1,'gaussian'),('2_gaussian',2,'gaussian'),('4_gaussian',4,'gaussian'),('known_shape_disk_oracle',0,'disk')]:
        start=r['models'][name]['parameters'] if name in r['models'] else r[name]['parameters']
        fixed[name]=fit(p,count,start,kind,maxiter=45,profile_gain=False)
        fixed[name].pop('image',None)
    r['gain_profiled_diagnostic']=True
    r['gain_profile_interpretation']='Complex gain is profiled independently at each frequency from the scattered training data. It is an additional nuisance parameter after incident/source coefficients were estimated; it does not establish source normalization or remove material-amplitude ambiguity.'
    r['fixed_gain_one_refits_from_selected_profiled_states']=fixed
    r['paired_start_reporting_limit']='The completed run retained the ten training objectives per model and best-fit metrics only. It did not persist every paired held-out error, time, or parameter vector, so it is not a paired statistical comparison.'
    r['gaussian_boundary_proxy_against_known_shape_disk']=r.pop('gaussian_boundary_floor_against_known_shape_disk')
    r['gaussian_boundary_proxy_interpretation']='Image disagreement with a fitted known-shape softened disk proxy. It is neither a true Gaussian approximation floor nor a bound.'
    r['executed_primary_source_hash']='not captured before the long primary run; current source was subsequently amended only for post-run diagnostics and labels'
    r['postrun_fixed_gain_source_sha256']=hashlib.sha256((ROOT/'code/a2/measured.py').read_bytes()).hexdigest()
    path.write_text(json.dumps(r,indent=2));print(json.dumps({k:{'objective':v['objective'],'held':[x['held_receiver_scattered_relative'] for x in v['metrics']]} for k,v in fixed.items()},indent=2))
if __name__=='__main__':main()
