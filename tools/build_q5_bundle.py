"""Build a deterministic ChatGPT upload bundle from already screened public files."""
import hashlib
import json
from pathlib import Path
import zipfile
from build_public_archive import sensitive

root=Path(__file__).resolve().parents[1]
manifest=json.loads((root/'PUBLIC_MANIFEST.json').read_text())
paths=[]
for row in manifest['records']:
    p=row['path']
    if (p.startswith(('research/q5_remaining_v1/','research/a5_material_observability/'))
        or p.startswith('Theory/Questions/A5_')
        or p in ('Theory/Questions/Q5.md','Theory/Questions/PAPER_DRAFT_A5.md',
                 'Theory/Questions/modal_boundary.json','Theory/Questions/Q4.md',
                 'Theory/Questions/A4_2.md','communication/Q5_REMAINING_V1_RESEARCH_REPORT_ZH.md',
                 'communication/A4_COMPLETE_RESEARCH_REPORT_ZH.md',
                 'research/trispace_self_calibration/a3_research/maxwell3d.py')):
        paths.append(p)
paths+=['Q5_CHATGPT_START_HERE.md','CODE_AND_REPRODUCIBILITY.md','PUBLICATION_SCOPE.md']
records=[]
target=root/'Q5_CHATGPT_BUNDLE.zip'
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for name in sorted(set(paths)):
        data=(root/name).read_bytes()
        if sensitive(data,Path(name)):raise RuntimeError('Publication screening failed: '+name)
        info=zipfile.ZipInfo(name,(2026,9,10,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
        z.writestr(info,data)
        records.append({'path':name,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
with zipfile.ZipFile(target) as z:
    assert z.testzip() is None
    assert all(hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256'] for r in records)
report={'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'bytes':target.stat().st_size,
        'files':len(records),'records':records,'scope':'Public snapshot integrity, not scientific acceptance'}
(root/'Q5_BUNDLE_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='records'}))
