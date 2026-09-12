"""Build an explicit A2 delivery snapshot; never bundle raw third-party data."""
from pathlib import Path
import argparse, hashlib, json, os, re, zipfile
ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT/'deliverables/Gaussian_A1_A1_2_A2_Pro_Package.zip'
def selected():
    paths = set()
    for folder in ['code/a2','Theory/a2','protocols/a2','reviews/a2','runs/a2','figures/a2',
                   'delegated/a2_extensions','delegated/a2_gpu','delegated/a2_audit','delegated/a2_ports','delegated/a2_review','delegated/a2_som',
                   'Theory/a12','delegated/a12_audit','delegated/a12_literature']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and p.suffix in {'.py','.md','.json','.txt','.npz','.png','.csv'}:
                if '__pycache__' in p.parts or p.name == 'CONVERSATION.json': continue
                if 'invalid_v1' in p.parts and p.name != 'INVALID.md': continue
                if 'legacy_unmanifested' in p.parts: continue
                if 'invalid_pilot' in p.parts and p.name != 'INVALID.md': continue
                if p.name in {'package_validation.json','postrun_snapshot.json'}: continue
                paths.add(p)
    for name in ['README.md','Theory/ A1.md','Theory/A1_2.md','Theory/A2.md',
                 'code/a12_graph/core.py','code/codex_stage.py','protocols/stage_output.schema.json','data/README.md','data/a2/DATA_INDEX_ZH.md',
                 'runs/a2/gpu/main128_interrupted.log','runs/a2/gpu/main128_resume.log','runs/a2/gpu/refined_field.log',
                 'runs/a2/postrun_snapshot.json','deliverables/requirements_a2.txt','deliverables/requirements_a2_gpu.txt']:
        paths.add(ROOT/name)
    for p in (ROOT/'deliverables').glob('*.md'): paths.add(p)
    # Compact historical evidence, not a claim that all A12 training has rerun.
    for name in ['runs/a12_graph/training.json','runs/a12_graph/diagnostics.json',
                 'runs/a12_theory/results.json','reviews/FINAL_VERDICT_ZH.md']:
        paths.add(ROOT/name)
    return sorted(p for p in paths if p.is_file())
def main():
    files={}
    for p in selected():
        arc='Gaussian/'+p.relative_to(ROOT).as_posix(); data=p.read_bytes()
        preserve=p in [ROOT/'Theory/ A1.md',ROOT/'Theory/A1_2.md'] or p.parent.name=='sources'
        if p.suffix=='.md' and not preserve:
            txt=data.decode('utf-8')
            # Convert only the workspace-specific Markdown targets, not source prose/code.
            prefix=str(ROOT)+'/'
            def replace(m):
                target=m.group(1); rel=os.path.relpath(ROOT/target,p.parent)
                return '('+rel.replace(' ','%20')+')'
            txt=re.sub(r'\('+re.escape(prefix)+r'([^\n)]*)\)',replace,txt)
            data=txt.encode('utf-8')
        files[arc]=data
    files['START_HERE.md']=b'# Gaussian A1 + A1_2 + A2\n\nStart at [Chinese reading guide](Gaussian/deliverables/START_HERE.md).\n'
    manifest={'scope':'File integrity only; not scientific acceptance or preregistration.',
              'excluded':['raw Fresnel measurements','raw conversation with user uploads','virtual environments','third-party PDFs','invalid numerical v1','original Pro attachments not obtained'],
              'files':{k:{'sha256':hashlib.sha256(v).hexdigest(),'bytes':len(v)} for k,v in files.items()}}
    files['MANIFEST.json']=json.dumps(manifest,ensure_ascii=False,indent=2).encode()
    with zipfile.ZipFile(DEST,'w',zipfile.ZIP_DEFLATED) as z:
        for k,v in files.items(): z.writestr(k,v)
    sha=hashlib.sha256(DEST.read_bytes()).hexdigest()
    (ROOT/'deliverables/PACKAGE_SHA256.txt').write_text(f'{sha}  {DEST.name}\n')
    print(json.dumps({'file':str(DEST),'files':len(files),'bytes':DEST.stat().st_size,'sha256':sha}))
if __name__=='__main__': main()
