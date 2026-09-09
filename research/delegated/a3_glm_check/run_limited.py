"""One bounded supported-Codex GLM invocation; never persist a credential.

No retries, no provider switching, no reasoning-stream persistence. This is a
client launcher, not a new orchestration service. Child receives only a named
profile and an ephemeral environment token; main config remains untouched.
"""
import json
import os
import re
import signal
import subprocess
import time
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
if __name__=='__main__':
    status_path=OUT/'status.json'
    if status_path.exists():
        raise SystemExit('A previous invocation is recorded; no automatic quota-consuming retry.')
    lines=(ROOT/'User'/'API.md').read_text().splitlines()
    matches=[re.findall(r'[A-Za-z0-9_.-]{24,}',line) for line in lines if re.match(r'(?i)^\s*glm\s*[:=：]',line)]
    if len(matches)!=1 or len(matches[0])!=1:
        raise SystemExit('Credential field ambiguous; nothing sent.')
    token=matches[0][0]
    env=os.environ.copy();env['A3_GLM_API_KEY']=token
    command=['codex','exec','--ignore-user-config','--profile','a3-glm-limited',
             '--ephemeral','--skip-git-repo-check','--sandbox','read-only','-C',str(OUT),
             '--json','-o',str(OUT/'FINAL.md'),'-']
    start=time.monotonic()
    p=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                       text=True,env=env,start_new_session=True)
    state=dict(status='running',pid=p.pid,profile='a3-glm-limited',model='glm-5.3-flash',
               maximum_wall_seconds=240,automatic_retries=0,task='mechanical table verification')
    status_path.write_text(json.dumps(state,indent=2)+'\n')
    try:
        output,_=p.communicate((OUT/'TASK.md').read_text(),timeout=240)
        state.update(status='completed' if p.returncode==0 else 'failed',exit_code=p.returncode)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid,signal.SIGTERM)
        try: output,_=p.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid,signal.SIGKILL);output,_=p.communicate()
        state.update(status='wall_limit_stopped',exit_code=p.returncode)
    state['elapsed_seconds']=time.monotonic()-start
    # Keep only event type counts and final usage; discard reasoning/tool streams.
    events={};errors=[]
    for line in output.splitlines():
        try: item=json.loads(line)
        except json.JSONDecodeError: continue
        typ=item.get('type','unknown');events[typ]=events.get(typ,0)+1
        if typ=='turn.completed':state['usage']=item.get('usage')
        if typ=='error':errors.append(str(item.get('message',''))[:500].replace(token,'[redacted]'))
    state['event_counts']=events;state['errors']=errors
    status_path.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(state))
