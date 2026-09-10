"""Isolated actual AI-Scientist development loop, Flash-only, bounded checks."""
import os
from pathlib import Path
import sys
import json
from types import SimpleNamespace
import time
import tomllib

HERE = Path(__file__).resolve().parents[1]
RUNTIME = Path('/Users/migodam/.codex/tools/agentic-ai-scientist')
MODEL = 'deepseek-flash'
PROFILE = 'q5-deepseek41-flash'
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
os.environ['OPENAI_AGENTS_DISABLE_TRACING'] = '1'
if not os.environ.get('DEEPSEEK_API_KEY'):
    for line in Path('/Users/migodam/.codex/.env').read_text().splitlines():
        name, sep, value = line.removeprefix('export ').partition('=')
        if sep and name.strip() == 'DEEPSEEK_API_KEY':
            os.environ[name.strip()] = value.strip().strip('\"\'')
profile = tomllib.loads(Path('/Users/migodam/.codex/'+PROFILE+'.config.toml').read_text())
assert profile['model'] == MODEL
assert profile['model_providers']['deepseek']['base_url'] == 'https://api.deepseek.com/'

sys.path.insert(0, str(RUNTIME))
import httpx


async def guard(request):
    if request.method == 'POST':
        body = json.loads(request.content)
        if body.get('model') != MODEL:
            raise RuntimeError('Forbidden or omitted external model')


from agents import AsyncOpenAI, set_default_openai_client, set_default_openai_api, set_tracing_disabled
import ai_scientist.agents_common as common


def configure(model):
    if model != MODEL:
        raise RuntimeError('Only verified V4.1 Flash API alias is authorized')
    set_tracing_disabled(True)
    client = AsyncOpenAI(api_key=os.environ['DEEPSEEK_API_KEY'],
        base_url='https://api.deepseek.com', max_retries=0,
        http_client=httpx.AsyncClient(event_hooks={'request':[guard]}), timeout=180)
    set_default_openai_client(client, use_for_tracing=False)
    set_default_openai_api('responses')
    return MODEL


common.configure_model_provider = configure
import run_research_loop as loop
loop.configure_model_provider = configure


def parent_writeup(idea, report, workdir, loop_dir, args, **kwargs):
    # Deliberately stop at a documented development handoff. This is not a
    # replacement research loop, and no native paper-stage completion is claimed.
    Path(loop_dir,'PARENT_REVIEW_REQUIRED.json').write_text(json.dumps({
        'workdir':workdir,'native_development_completed':True,
        'native_writeup_executed':False,'scientific_acceptance':False}, indent=2))
    return workdir


loop.run_writeup = parent_writeup


def main():
    configure(MODEL)
    # Small real Responses call confirms endpoint/model before launching workers.
    import asyncio
    async def smoke():
        client = AsyncOpenAI(api_key=os.environ['DEEPSEEK_API_KEY'], base_url='https://api.deepseek.com',
                             max_retries=0, timeout=90)
        response = await client.responses.create(model=MODEL, input='Reply with exactly READY.', max_output_tokens=128)
        result = {'requested_model':MODEL,'returned_model':response.model,
                  'status':response.status,'ready':response.output_text.strip()=='READY',
                  'official_mapping':'https://api-docs.deepseek.com/quick_start/pricing/',
                  'usage':response.usage.model_dump() if response.usage else None}
        with (HERE/'results/provider_smoke.json').open('x') as f: json.dump(result,f,indent=2)
        if response.model != MODEL or not result['ready']:
            raise RuntimeError('Model smoke not confirmed; no fallback')
    asyncio.run(smoke())
    task = HERE/'docs/WORKER_TASK.md'
    idea = {'Name':'q5_boundary_reconstruction','Title':'Independent bounded reconstruction of supplied material/gain identities',
            'Short Hypothesis':'The supplied finite-risk constants and fixed-loss inverse are reproducible under their declared restricted assumptions.',
            'Related Work':'All general linear algebra and binary testing are prior tools. Parent owns novelty and theorem decisions.',
            'Abstract':'Implement only the bounded checks in '+str(task)+'. No topic expansion, paper, hardware, external agents, or literature calls.',
            'Experiments':{'steps':['Read bounded task','Implement independent checks','Execute tests','Report evidence limits'],
                           'metrics':['identity errors','finite-risk constants','interval coverage status'],
                           'compute_estimate':'Single CPU thread, small arrays'},
            'Risk Factors and Limitations':['Restricted known-geometry modal class, not joint two-sphere recovery','Parent review required']}
    work = HERE/'pipeline'/time.strftime('boundary_%Y%m%d_%H%M%S')
    work.mkdir(parents=True, exist_ok=False)
    args = SimpleNamespace(model=MODEL, worker='codex',codex_profile=PROFILE,
        max_safety_rounds=2,min_dev_rounds=1,final_max_turns=18,codex_timeout=900,
        wiki_path=str(work/'wiki.jsonl'),num_cite_rounds=0,writeup_retries=0)
    (work/'invocation.json').write_text(json.dumps({'args':vars(args),'entry':'installed run_development_loop',
        'runtime':str(RUNTIME),'writeup':'parent-owned; native paper stage not executed'},indent=2))
    loop.run_development_loop(idea,str(work),args)


if __name__ == '__main__': main()
