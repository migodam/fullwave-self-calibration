"""Copy a scoped research archive; never execute imported research code.

Usage: python tools/build_public_archive.py --source /path/to/Inv_SLAM [--copy]
Sources stay unchanged. Missing external resources are documented, not bundled.
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import zipfile

DEST = Path(__file__).resolve().parents[1]
ROOTS = ('Theory', 'communication', 'research', 'experiments')
SUFFIXES = {'.md', '.py', '.json', '.jsonl', '.bib', '.tex', '.png', '.svg',
            '.npz', '.npy', '.toml', '.yaml', '.yml', '.sh', '.csv', '.tsv'}
BLOCK_DIR = {'external', 'knowledge_bank', 'node_modules', '__pycache__',
             'data', 'datasets', 'downloads', 'papers', 'logs', 'sessions',
             'transcripts', 'worker_logs', 'raw', 'raw_data', 'tmp'}
BLOCK_NAMES = {'evidence.json', 'scholarqa_collect.json', 'AGENTS.md'}
SECRET = re.compile(r'(?:sk-[A-Za-z0-9_-]{20,}|gh[pousr]_[A-Za-z0-9]{20,}'
                    r'|AIza[0-9A-Za-z_-]{30,}|-----BEGIN [A-Z ]*PRIVATE KEY-----'
                    r'|Bearer\s+[A-Za-z0-9_.-]{20,})')
ASSIGNMENT = re.compile(r'''(?ix)(?:api[_-]?key|access[_-]?token|password|client[_-]?secret)
    ["']?\s*[:=]\s*["']([^"'\n]{12,})["']''')


def reason(path, size):
    parts = path.parts
    if str(path).startswith('research/q5_remaining_v1/pipeline/'):
        if not any(p in parts for p in ('code', 'inputs', 'results')):
            return 'live pipeline control state; scientific snapshots only'
    if any(p.startswith('.') or p in BLOCK_DIR or p.startswith('.venv') for p in parts):
        return 'private/cache/third-party/raw-data directory'
    name = path.name.lower()
    if name in {n.lower() for n in BLOCK_NAMES} or name.startswith(('retrieval_', 'transport_')):
        return 'raw retrieval or transport material'
    if any(x in name for x in ('credential', 'worker_log', 'transcript', 'conversation', 'reasoning')):
        return 'private or raw agent material'
    a4_numeric_log = str(path).startswith('research/a4_reliability_v1/') and path.suffix in ('.txt', '.log')
    if path.suffix not in SUFFIXES and name not in ('requirements.txt', 'checksums.txt') and not a4_numeric_log:
        return 'not an authored-document/code/result allowlisted type'
    if size > 32 * 1024 * 1024:
        return 'file exceeds 32 MiB public snapshot cap'
    if path.suffix in ('.npz', '.npy') and ('measured' in str(path).lower() or 'literature' in str(path).lower()):
        return 'measured-data redistribution not reviewed'
    return None


def sensitive(data, path):
    if path.suffix in ('.png', '.npz', '.npy'):
        return False
    text = data.decode('utf-8', errors='replace')
    if SECRET.search(text):
        return True
    for match in ASSIGNMENT.finditer(text):
        value = match.group(1)
        if not any(x in value.lower() for x in ('getenv', 'environ', 'placeholder', 'example', 'your_', 'insert_', 'redact', 'none')):
            return True
    if path.suffix == '.jsonl':
        for line in text.splitlines()[:10]:
            try:
                obj = json.loads(line)
            except ValueError:
                continue
            if isinstance(obj, dict) and (obj.get('type') in ('response_item', 'session_meta', 'turn_context')
                                         or obj.get('role') in ('assistant', 'system', 'user')):
                return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--copy', action='store_true')
    args = parser.parse_args()
    source = args.source.resolve()
    assert source != DEST and (source / 'Theory/Questions').is_dir()
    selected = []; excluded = collections.Counter(); flagged = []
    for root in ROOTS:
        for path in sorted((source / root).rglob('*')):
            if not path.is_file() or path.is_symlink():
                continue
            rel = path.relative_to(source)
            why = reason(rel, path.stat().st_size)
            if why:
                excluded[why] += 1
                continue
            data = path.read_bytes()
            if rel.suffix == '.json':
                try:
                    json.loads(data)
                except (ValueError, UnicodeError):
                    excluded['incomplete or invalid JSON snapshot'] += 1
                    continue
            if sensitive(data, rel):
                flagged.append(str(rel)); continue
            selected.append((rel, data, 'workspace'))
    archive = source / 'Theory/Questions/A4_fullwave_research_v1.zip'
    archive_hash = hashlib.sha256(archive.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive) as bundle:
        for entry in bundle.infolist():
            rel = PurePosixPath(entry.filename)
            if entry.is_dir():
                continue
            assert not rel.is_absolute() and '..' not in rel.parts
            assert str(rel).startswith('research/a4_reliability_v1/')
            why = reason(rel, entry.file_size)
            if why:
                excluded['A4: '+why] += 1; continue
            data = bundle.read(entry)
            if sensitive(data, rel):
                flagged.append(str(rel)); continue
            selected.append((Path(rel), data, 'user-supplied A4 zip'))
    report = dict(selected_files=len(selected), selected_bytes=sum(len(d) for _,d,_ in selected),
                  exclusions=dict(excluded), sensitive_file_paths=flagged)
    print(json.dumps(report, indent=2))
    if flagged:
        raise SystemExit('Review flagged files before publication; no copy performed.')
    assert len({str(p) for p,_,_ in selected}) == len(selected), 'Duplicate targets'
    if not args.copy:
        return
    records = []
    for rel, data, origin in selected:
        dest = DEST / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        # Byte-identical bulk export; source hashes are not rewritten.
        dest.write_bytes(data)
        records.append(dict(path=str(rel), bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), origin=origin))
    (DEST/'PUBLIC_MANIFEST.json').write_text(json.dumps(dict(
        date='2026-09-10', a4_archive_sha256=archive_hash, records=records,
        exclusion_counts=dict(excluded), scope='Archive integrity only; not scientific replication.'), indent=2)+'\n')
    lines = ['# Complete public file index', '',
             'Paths, source hashes and archive origin: [PUBLIC_MANIFEST.json](PUBLIC_MANIFEST.json).', '',
             'Historical claims are not automatically current conclusions. See [research navigation](RESEARCH_NAVIGATION.md).', '']
    for root in ROOTS:
        lines += ['## '+root, '']
        for rel, _, _ in selected:
            if rel.parts[0] == root:
                lines.append(f'- [{rel}](<{rel.as_posix()}>)')
        lines.append('')
    (DEST/'FILE_INDEX.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
