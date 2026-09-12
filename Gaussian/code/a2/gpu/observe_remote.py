"""Read-only remote execution receipt. Requires an already authenticated SSH master.
No credentials are read, accepted, or stored by this script.
"""
from pathlib import Path
import subprocess,json,datetime
ROOT=Path(__file__).resolve().parents[3]
cmd=['ssh','-S','/tmp/gaussian-xinan-ssh','jayzh@100.121.97.39','powershell -NoProfile -Command "Get-Process -Id 31320 | Select-Object Id,CPU,PeakWorkingSet64,StartTime | ConvertTo-Json"']
r=subprocess.run(cmd,capture_output=True,timeout=15)
record=dict(observed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),host='XINAN',remote_project='D:/home/research/SSH_workspace/Gaussian_A2',process_id=31320,exit_code=r.returncode)
if r.returncode==0:record['process']=json.loads(r.stdout.decode('utf-8'))
else:record['process']='Target process no longer returned by Get-Process; check completion log separately.'
p=ROOT/'runs/a2/gpu/remote_observations.json';items=json.loads(p.read_text()) if p.exists() else [];items.append(record);p.write_text(json.dumps(items,indent=2));print(json.dumps(record))
