"""Check A2 archive integrity, primary-source preservation and declared links."""
from pathlib import Path
import ast,hashlib,json,re,zipfile
from urllib.parse import unquote,urlsplit
ROOT=Path(__file__).resolve().parents[2]
def validate():
    archive=ROOT/'deliverables/Gaussian_A1_A1_2_A2_Pro_Package.zip'
    errors=[]; checked=0; links=0
    with zipfile.ZipFile(archive) as z:
        names=set(z.namelist()); manifest=json.loads(z.read('MANIFEST.json'))
        if z.testzip(): errors.append('CRC failed')
        for name,meta in manifest['files'].items():
            data=z.read(name)
            if hashlib.sha256(data).hexdigest()!=meta['sha256']:errors.append('hash: '+name)
            if name.endswith('.py'):ast.parse(data,filename=name);checked+=1
        inputs=json.loads((ROOT/'runs/a2/input_manifest.json').read_text())['inputs']
        for name,meta in inputs.items():
            if name.endswith('CONVERSATION.json'):continue
            if hashlib.sha256(z.read(name)).hexdigest()!=meta['sha256']:errors.append('original altered: '+name)
        extensions=json.loads((ROOT/'runs/a2/extensions/source_manifest.json').read_text())
        for name,meta in extensions['files'].items():
            if hashlib.sha256(z.read('Gaussian/'+name)).hexdigest()!=meta['sha256']:errors.append('new source altered: '+name)
        recheck='Gaussian/runs/a2/portable_recheck/results.json'
        if recheck in names:
            checked_sources=json.loads(z.read(recheck))
            if not checked_sources['all_passed']:errors.append('portable recheck failed')
            for name,sha in checked_sources['tested_python_sources'].items():
                if hashlib.sha256(z.read('Gaussian/'+name)).hexdigest()!=sha:errors.append('changed since portable check: '+name)
        for name in names:
            if not name.endswith('.md'):continue
            if not (name.startswith('Gaussian/deliverables/') or name=='Gaussian/README.md' or name=='Gaussian/Theory/A2.md'):continue
            for target in re.findall(r'!?\[[^\]]*\]\(([^\n)]*)\)',z.read(name).decode()):
                target=unquote(target.strip('<>'))
                if urlsplit(target).scheme or target.startswith('#'):continue
                target=target.split('#')[0]
                resolved=__import__('posixpath').normpath(str(Path(name).parent/target));links+=1
                if resolved.endswith('.zip'):continue # companion delivery, not nested archive
                if resolved not in names:errors.append(f'link {name}: {target}')
    result={'zip_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'integrity_errors':errors,'python_files_parsed':checked,'local_links_checked':links,'scope':'Integrity and selected Markdown links only. Numerical rerun is separately recorded.'}
    (ROOT/'runs/a2/package_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result));assert not errors
if __name__=='__main__':validate()
