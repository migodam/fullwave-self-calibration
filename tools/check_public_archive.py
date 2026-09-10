"""Read-only source-integrity, root-navigation and secret-pattern checks."""
import hashlib
import json
from pathlib import Path
import re
import zipfile
from urllib.parse import unquote
from build_public_archive import sensitive

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT/'PUBLIC_MANIFEST.json').read_text())
errors = []
for row in manifest['records']:
    p = ROOT/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != row['sha256']:
        errors.append('hash:'+row['path'])
for name in ('README.md','RESEARCH_NAVIGATION.md','GPT_READING_GUIDE.md','CODE_AND_REPRODUCIBILITY.md','Q5_CHATGPT_START_HERE.md'):
    for link in re.findall(r'\]\(([^)]+)\)', (ROOT/name).read_text()):
        link = unquote(link.strip('<>').split('#')[0])
        if link and '://' not in link and not (ROOT/link).exists():
            errors.append('link:'+name+':'+link)
scanned = 0
for p in ROOT.rglob('*'):
    if not p.is_file() or '.git' in p.parts or '__pycache__' in p.parts:
        continue
    scanned += 1
    if p.suffix == '.zip':
        with zipfile.ZipFile(p) as z:
            for item in z.infolist():
                if sensitive(z.read(item), Path(item.filename)):
                    errors.append('sensitive archive member:'+item.filename)
        continue
    if sensitive(p.read_bytes(), p):
        errors.append('sensitive:'+str(p.relative_to(ROOT)))
print(json.dumps(dict(passed=not errors, source_files=len(manifest['records']),
                      scanned_files=scanned, errors=errors),indent=2))
raise SystemExit(bool(errors))
