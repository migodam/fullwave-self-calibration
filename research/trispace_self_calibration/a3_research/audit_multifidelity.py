"""Recompute paired development conclusions without modifying raw results."""
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parent


def audit():
    rows=json.loads((OUT/'results/multifidelity_calibration.json').read_text())
    methods=['fine_direct','coarse_only','coarse_warm_fine','value_corrected','tangent_corrected']
    assert len(rows)==10 and len({(r['seed'],r['method']) for r in rows})==10
    comparisons=[]
    for seed in [8101,8102]:
        group={r['method']:r for r in rows if r['seed']==seed}
        assert set(group)==set(methods)
        fine,warm=group['fine_direct'],group['coarse_warm_fine']
        for method,row in group.items():
            for event in row['events']:
                if event['accepted']:
                    assert event['trial_loss']<event['anchor_loss'] and event['ratio']>.1
            objective_gap=(row['final_fine_objective']-fine['final_fine_objective'])/fine['final_fine_objective']
            if method in ['coarse_warm_fine','tangent_corrected']:
                assert abs(objective_gap)<1e-8
                assert np.linalg.norm(np.array(row['estimated'])-fine['estimated'])<1e-6
            comparisons.append(dict(seed=seed,method=method,
                relative_fine_objective_gap=objective_gap,
                time_ratio_vs_cold=row['online_seconds']/fine['online_seconds'],
                time_ratio_vs_warm=row['online_seconds']/warm['online_seconds']))
    summary=[]
    for method in methods:
        group=[r for r in rows if r['method']==method]
        summary.append(dict(method=method,mean_seconds=float(np.mean([r['online_seconds'] for r in group])),
            mean_fine_rhs=float(np.mean([r['fine_rhs'] for r in group])),
            mean_coarse_rhs=float(np.mean([r['coarse_rhs'] for r in group])),
            mean_pose_mm=1000*float(np.mean([r['pose_error_m'] for r in group])),
            mean_material_percent=100*float(np.mean([r['material_relative_error'] for r in group]))))
    out=dict(passed=True,scene_count=2,summary=summary,paired=comparisons,
        scope='two inspected DEVELOPMENT scenes, serial single runs; no confidence or novelty claim',
        decision='value/tangent correction not supported beyond simple coarse warm start; retain as negative ablation')
    (OUT/'results/multifidelity_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':audit()
