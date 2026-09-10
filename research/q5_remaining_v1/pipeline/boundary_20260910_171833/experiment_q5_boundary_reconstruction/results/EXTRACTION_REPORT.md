# Q5 boundary-reconstruction package: section split + modal-boundary key flatten

Scope: text extraction only. No experiment/numerical code was written, no
provider calls, no installs. Everything was read from `inputs/` (read-only) and
written to `inputs/derive/`, `code/`, `logs/`, `results/`.

## Interpreter

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/trispace_self_calibration/a3_research/.venv3d/bin/python`
-> `/opt/homebrew/opt/python@3.13/bin/python3.13`, Python 3.13.13.

## Step 1 - A5_THEORY.md section split

Command: `"$PY" code/split_a5_sections.py`

- Split regex: `(?m)^## (\d+)\. ` on raw bytes; heading line retained in-file.
- Preamble (bytes before the first heading) -> `A5_sec_00_preamble.md`.
- 12 numbered sections -> `A5_sec_01.md` ... `A5_sec_12.md` (13 files total).
- Each part is a byte-exact source slice. Exactly one part, `A5_sec_12.md`, sat
  at EOF and lacked a final newline; one `\n` was appended so the file is
  self-contained. Tracked per-row as `trailing_newline_added` in `A5_SIZES.txt`.
- Source: 19015 bytes, 308 lines (`splitlines`), 307 `\n`, sha256
  `5b135f2ea486d5ffbfcd92bdc1a57c103aaecd0eeba3775da49f70052d746e49`.
- Sum of written bytes 19016 = 19015 + the single appended newline; sum of
  lines 308 = source lines 308.
- Reassembly check: concatenating the source slices (before the appended
  newline) reproduces the source byte-for-byte -> `parts_concat_equals_source=True`.
- Heading check: every `A5_sec_NN.md` starts with its own `## NN. ` line; no
  mismatches; no file is empty.

## Step 2 - modal_boundary.json leaf flatten

Command: `"$PY" code/flatten_modal_boundary.py` -> `inputs/derive/MODAL_BOUNDARY_KEYS.txt`

- Root type `dict`; 24 top-level keys.
- 866 leaves, all unique paths, sorted lexicographically by path string.
- Path convention: dict keys joined with `.`; list elements addressed as
  `[i]`; empty containers would appear as `[]`/`{}` (none occur here).
- All list elements expand down to scalars (list sizes 2 throughout).
- Leaf count 866 reconfirmed by an independent iterative-stack counter.
- Format in file: `path = <json>` with `ensure_ascii=False`, plus `#`-prefixed
  header lines carrying leaf count and top-level key summary.

## Step 3 - stdout

Command: `logs/step3_stdout.txt` was produced by cat-ing `A5_SIZES.txt` and
re-running the flatten script. Only the byte/line size table and the top-level
key summary were emitted; the 866-line leaf list was not printed to stdout.

## Integrity

`shasum -a 256` on every file under `inputs/` matches `inputs/SOURCE_HASHES.json`:

| file | sha256 | matches record |
| --- | --- | --- |
| A5_THEORY.md | 5b135f2e... | yes (source + copy) |
| modal_boundary.json | 036491f8... | yes |
| modal_boundary_pretty.json | 570b7ace... | yes (derived_artifacts) |
| WORKER_TASK.md | 8272fb3c... | yes |

mtimes of the four `inputs/` files are unchanged (16:53, 16:54, 17:18, 17:19).
Only the new directory `inputs/derive/` was added under `inputs/`.

## Uncertainties

- `modal_boundary.json` carries its own internal `source_sha256` field
  (`daad3707daf72d6a8586035b7fa5f81777049b824baa98dd5e35f5471d4ac5fc`) which
  is *not* the hash of the file itself; it refers to some upstream artifact.
  The file-level hash `036491f8...` is the one that matches `SOURCE_HASHES.json`.
- The chosen path convention (`[i]` for list indices) is an unstated detail of
  the spec; it is documented in the output header and in this report.
