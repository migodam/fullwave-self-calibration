"""Subscription-auth Codex stage adapter; persists only final answer, not reasoning/events."""
import argparse, datetime, hashlib, json, os, pathlib, subprocess, time
p=argparse.ArgumentParser();p.add_argument('task',type=pathlib.Path);p.add_argument('--model',choices=['gpt-5.6-sol','gpt-5.6-terra'],required=True);p.add_argument('--out',type=pathlib.Path,required=True);p.add_argument('--timeout',type=int,default=1800);p.add_argument('--schema',type=pathlib.Path);a=p.parse_args()
root=pathlib.Path(__file__).resolve().parents[1]; out=a.out.resolve();out.relative_to(root);out.mkdir(parents=True,exist_ok=True)
record={'model':a.model,'transport':'codex exec ChatGPT login','task_sha256':hashlib.sha256(a.task.read_bytes()).hexdigest(),'started':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running'}
state=out/'state.json';state.write_text(json.dumps(record,indent=2));t=time.time()
cmd=['codex','exec','--ignore-user-config','--ephemeral','--skip-git-repo-check','-C',str(root),'-m',a.model,'-c','model_reasoning_effort="high"','--sandbox','workspace-write','-o',str(out/'final.md'),'-']
if a.schema:cmd[-1:-1]=['--output-schema',str(a.schema.resolve())]
try:
    # Suppress raw event/reasoning streams and credential-bearing diagnostics.
    r=subprocess.run(cmd,input=a.task.read_text(),text=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=a.timeout)
    record.update(status='completed' if r.returncode==0 and (out/'final.md').exists() else 'failed',exit_code=r.returncode)
except subprocess.TimeoutExpired:record.update(status='timeout')
if record['status']=='completed' and a.schema:
    try:json.loads((out/'final.md').read_text());record['structured_json_valid']=True
    except (ValueError,OSError):record.update(status='invalid_output',structured_json_valid=False)
record['wall_seconds']=time.time()-t;state.write_text(json.dumps(record,indent=2));print(json.dumps(record))
