"""Fail-closed model discovery; never expose credentials or response bodies."""
import json
import os
from pathlib import Path
import urllib.request
import urllib.error

HERE = Path(__file__).resolve().parents[1]
ALLOWED = 'deepseek-v4.1-flash'


def main():
    token = os.environ.get('DEEPSEEK_API_KEY')
    if not token:
        p = Path('/Users/migodam/.codex/.env')
        if p.exists():
            for line in p.read_text().splitlines():
                name, sep, value = line.removeprefix('export ').partition('=')
                if sep and name.strip() == 'DEEPSEEK_API_KEY':
                    token = value.strip().strip('\"\'')
    result = {'allowed_model': ALLOWED, 'pro_forbidden': True, 'model_call_made': False}
    if not token:
        result['status'] = 'missing_credential'
    else:
        request = urllib.request.Request('https://api.deepseek.com/models',
                                         headers={'Authorization': 'Bearer '+token})
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                data = json.load(response)
            ids = [row['id'] for row in data.get('data', [])]
            result['available_model_ids'] = ids
            result['status'] = 'exact_model_found' if ALLOWED in ids else 'exact_model_not_listed'
        except urllib.error.HTTPError as error:
            result.update(status='http_error', http_status=error.code)
        except Exception as error:
            result.update(status='transport_error', error_type=type(error).__name__)
    path = HERE/'results/provider_preflight.json'
    path.parent.mkdir(exist_ok=True)
    with path.open('x') as file:
        json.dump(result, file, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
