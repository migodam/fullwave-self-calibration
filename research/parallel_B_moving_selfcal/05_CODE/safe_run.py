"""Replay entrypoint with a floating-point-only bound initialization guard.

The frozen experiment is unchanged. abs(exp/log polar gain) can exceed a
bound by one rounding unit. Clamp only violations <= 1e-12; larger violations
still raise. This does not change any physical domain or initial material.
"""
import numpy as np
import experiment
_original = experiment.least_squares

def guarded_least_squares(fun, x0, *args, **kwargs):
    x = np.asarray(x0, dtype=float)
    if 'bounds' in kwargs:
        lo, hi = map(np.asarray, kwargs['bounds'])
        violation = max(float(np.max(lo-x)), float(np.max(x-hi)), 0.)
        if violation > 1e-12:
            raise ValueError(f'initial guess genuinely outside domain: {violation}')
        x = np.maximum(lo, np.minimum(hi, x))
    return _original(fun, x, *args, **kwargs)

experiment.least_squares = guarded_least_squares
if __name__ == '__main__':
    experiment.main()
