"""Invoke the installed AI Scientist development loop, not a replacement loop."""
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNTIME = Path('/Users/migodam/.codex/tools/agentic-ai-scientist')
sys.path.insert(0, str(RUNTIME))
os.environ['OPENAI_AGENTS_DISABLE_TRACING'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.chdir(ROOT)
from run_research_loop import configure_model_provider, run_development_loop

if __name__ == '__main__':
    task = ROOT / 'research/delegated/a3_algebra/TASK.md'
    idea = {
        'Name': 'a3_conditional_calibration_checks',
        'Title': 'Independent checks of fixed-chart derivatives, gain quotients and branch validation',
        'Short Hypothesis': 'A3 conditional algebra can be implemented with explicit failure controls; no performance or novelty claim follows.',
        'Related Work': 'SOM means subspace-based optimization, not self-organizing maps. Generic projection, gain closure and dual residual methods are prior art. Parent owns literature and scientific judgment.',
        'Abstract': 'Bounded independent reconstruction of A3 theorems 1–4. Read '+str(task)+'. Root is '+str(ROOT)+'. Implement the supplied task without changing the scientific topic. No ACL review, no global novelty claim, no recursive workers or pipeline. Deliver auditable code/tests/raw outputs in research/delegated/a3_algebra/.',
        'Experiments': {
            'datasets': ['deterministic complex linear systems and Gaussian validation fixtures'],
            'baselines': ['correct vs omitted residual derivative', 'rank-one vs nonseparable fields', 'covered vs missing candidate bank', 'certified vs uncontrolled resolvent'],
            'metrics': ['relative identity errors', 'empirical selection risk', 'conditional bounds', 'negative control violations'],
            'compute_estimate': 'single CPU thread, small NumPy/SciPy arrays, 10000 noise draws, existing venv',
            'steps': ['Read bounded task', 'Implement and execute checks', 'Audit controls and preserve failures', 'Summarize scope and remaining work']},
        'Risk Factors and Limitations': ['Conditional fixtures are not nonlinear calibration proof', 'No global candidate coverage', 'No SOM advantage evidence', 'Parent must independently audit']}
    loop = ROOT/'experiments/idea_loops'/time.strftime('loop_%Y-%m-%d_%H-%M-%S_a3_development')
    loop.mkdir(parents=True, exist_ok=False)
    args = SimpleNamespace(model='deepseek-v4-pro', worker='codex',
        codex_profile='deepseek-flash', max_safety_rounds=2, min_dev_rounds=1,
        final_max_turns=28, codex_timeout=2400, num_cite_rounds=1,
        writeup_retries=1, wiki_path=str(HERE/'pipeline_wiki.jsonl'))
    (HERE/'pipeline_invocation.json').write_text(json.dumps({
        'entry':'installed run_development_loop', 'loop':str(loop), 'args':vars(args),
        'task':str(task), 'initial_novelty':'parent-owned evidence QA; supplied bounded validation',
        'completion_boundary':'pipeline completion is not paper acceptance'}, indent=2)+'\n')
    (HERE/'pipeline_seed.json').write_text(json.dumps(idea, indent=2)+'\n')
    configure_model_provider(args.model)
    print('A3 installed development loop:', loop, flush=True)
    run_development_loop(idea, str(loop), args)
