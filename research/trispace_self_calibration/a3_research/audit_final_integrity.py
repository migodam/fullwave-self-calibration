"""Post-collection integrity audit; does not rerun or modify final decisions."""
import hashlib
import json
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
manifest=json.loads((OUT/'results/rom_final_source_manifest.json').read_text())
rows=json.loads((OUT/'results/rom_final.json').read_text())
assert len(rows)==240
for path,digest in manifest['source_hashes'].items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
assert hashlib.sha256((OUT/'ROM_FINAL_PROTOCOL.md').read_bytes()).hexdigest()==manifest['protocol_sha256']
analyzer=hashlib.sha256((OUT/'analyze_rom_final.py').read_bytes()).hexdigest()
assert analyzer=='1bb4cdf353e1a3dda7e63e3163fd62e4880cfcde06878845e044d5e0284e2cd5'
exceptions=[dict(seed=r['seed'],method=r['method'],error=r.get('error')) for r in rows if r['status']=='exception']
checks=dict(completed_rows=240,source_hashes_match=True,protocol_hash_matches=True,
    frozen_analyzer_hash_matches=True,exceptions=exceptions,
    exception_stage_ambiguity_affects_results=bool(exceptions),
    all_accepted_timestamps_before_deadline=all(e['seconds']<8 for r in rows for e in r.get('events',[])),
    scope='frozen-record/source integrity, not continuum accuracy, novelty or external review')
assert checks['all_accepted_timestamps_before_deadline']
(OUT/'results/rom_final_integrity.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
