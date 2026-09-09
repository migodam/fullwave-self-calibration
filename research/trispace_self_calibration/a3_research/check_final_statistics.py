"""Analytic endpoint checks; does not read any final performance outcomes."""
import json
from pathlib import Path
import math
from scipy.stats import binom
from analyze_rom_final import cp_upper,cp_lower,median_interval

out=Path(__file__).resolve().parent/'results/final_statistics_checks.json'
assert abs(cp_upper(0,40)-(1-.025**(1/40)))<1e-14
assert cp_lower(0,40)==0 and cp_upper(40,40)==1
assert abs(cp_lower(40,40)-.025**(1/40))<1e-14
assert median_interval([1.2]*40)==[1.2,1.2]
assert median_interval([])==[0.,None]
assert median_interval([1.2]*4)==[0.,None]
assert math.isclose(float(binom.cdf(0,10,.5)),2**-10,abs_tol=1e-15)
out.write_text(json.dumps(dict(passed=True,checks=['Clopper-Pearson zero/all endpoints',
    'constant and insufficient-sample median interval','exact zero-count sign tail'],
    reads_final_results=False),indent=2)+'\n')
print('Analytic statistics checks passed; no final outcomes accessed.')
