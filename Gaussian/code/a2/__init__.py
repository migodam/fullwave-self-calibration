"""Bounded Gaussian A2 physics and response-derived port infrastructure."""
from .physics import Geometry, VIE, FFTGreen
from .ports import GaussianComponent, LocalTNetwork, render, pack, unpack, parameter_scales
