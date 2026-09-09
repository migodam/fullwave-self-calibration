"""Enter the installed AI Scientist development loop for a supplied research task.

No replacement orchestration framework. Initial novelty ideation is intentionally
not used: its first attempt confused subspace optimization with self-organizing
maps. Codex owns the separately recorded evidence audit and supplied direction.
The installed develop/evaluate/revise/writeup pipeline itself is unchanged.
"""
from pathlib import Path
from types import SimpleNamespace
import json
import os
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNTIME = Path('/Users/migodam/.codex/tools/agentic-ai-scientist')
sys.path.insert(0, str(RUNTIME))
os.environ['OPENAI_AGENTS_DISABLE_TRACING'] = '1'
os.chdir(ROOT)
from run_research_loop import configure_model_provider, run_development_loop

if __name__ == '__main__':
    idea = json.loads((HERE / 'start_idea.json').read_text())[0]
    task_path = ROOT / 'research/delegated/a2_theory_validation/TASK.md'
    idea['Abstract'] += (' Read the full bounded task at ' + str(task_path)
                        + '. All source paths are relative to ' + str(ROOT)
                        + '. SOM is subspace-based optimization, never self-organizing maps.')
    loop = ROOT / 'experiments/idea_loops' / time.strftime('loop_%Y-%m-%d_%H-%M-%S_a2_development')
    loop.mkdir(parents=True, exist_ok=False)
    args = SimpleNamespace(model='deepseek-v4-pro', worker='codex',
        codex_profile='deepseek-flash', max_safety_rounds=3, min_dev_rounds=1,
        final_max_turns=35, codex_timeout=3600, num_cite_rounds=3,
        writeup_retries=1, wiki_path=str(HERE/'pipeline_wiki.jsonl'))
    (HERE/'pipeline_invocation.json').write_text(json.dumps({
        'entry':'installed run_development_loop', 'loop':str(loop),
        'args':vars(args), 'initial_novelty':'parent-owned; automated acronym drift interrupted',
        'task':str(task_path)}, indent=2)+'\n')
    configure_model_provider(args.model)
    print('A2 bounded development loop:', loop, flush=True)
    run_development_loop(idea, str(loop), args)
