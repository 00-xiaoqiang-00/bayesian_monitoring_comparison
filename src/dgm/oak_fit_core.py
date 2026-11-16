"""Deterministic OAK curve calculations retained from the original fit.

The 40-node Gauss-Legendre integration and Cox curve diagnostic are preserved
to reproduce the numerical objective used for the locked OAK DGM.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.stats import norm

EPS = 1e-12
GAUSS_NODES, GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(40)

@dataclass(frozen=True)
class MixtureParameters:
    rho: float
    pi: float
    mu_r1: float
    sigma_r1: float
    mu_r2: float
    sigma_r2: float
    mu_d: float
    sigma_d: float


def read_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",")
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(f"Expected two headerless columns in {path}")
    order = np.argsort(data[:, 0], kind="stable")
    t = data[order, 0].astype(float)
    s = np.minimum.accumulate(np.clip(data[order, 1].astype(float), 0.0, 1.0))
    return t, s


def interp_curve(curve: tuple[np.ndarray, np.ndarray], t: np.ndarray) -> np.ndarray:
    x, y = curve
    return np.interp(t, x, y, left=1.0, right=float(y[-1]))


def lognormal_survival(t: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    positive = t > 0
    out[positive] = norm.sf((np.log(t[positive]) - mu) / sigma)
    return out


def component_pfs_survival(
    t: np.ndarray,
    mu_r: float,
    sigma_r: float,
    mu_d: float,
    sigma_d: float,
    rho: float,
) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    positive = t > 0
    log_t = np.log(t[positive])
    a = (log_t - mu_r) / sigma_r
    b = (log_t - mu_d) / sigma_d
    # Integrate over u=Phi(Z_D), which is uniform.  Gauss-Legendre nodes avoid
    # both endpoints, so qnorm remains finite even for very early times.
    u0 = norm.cdf(b)
    half_width = 0.5 * (1.0 - u0)
    midpoint = 0.5 * (1.0 + u0)
    u = midpoint[:, None] + half_width[:, None] * GAUSS_NODES[None, :]
    z_d = norm.ppf(np.clip(u, EPS, 1.0 - EPS))
    conditional = norm.sf(
        (a[:, None] - rho * z_d) / math.sqrt(1.0 - rho * rho)
    )
    out[positive] = half_width * (conditional @ GAUSS_WEIGHTS)
    return np.clip(out, 0.0, 1.0)


def mixture_pfs_survival(t: np.ndarray, p: MixtureParameters) -> np.ndarray:
    s1 = component_pfs_survival(
        t, p.mu_r1, p.sigma_r1, p.mu_d, p.sigma_d, p.rho
    )
    s2 = component_pfs_survival(
        t, p.mu_r2, p.sigma_r2, p.mu_d, p.sigma_d, p.rho
    )
    return p.pi * s1 + (1.0 - p.pi) * s2


def survival_median(t: np.ndarray, s: np.ndarray) -> float:
    below = np.flatnonzero(s <= 0.5)
    if len(below) == 0:
        return math.inf
    i = int(below[0])
    if i == 0:
        return float(t[0])
    if s[i - 1] == s[i]:
        return float(t[i])
    return float(np.interp(0.5, [s[i], s[i - 1]], [t[i], t[i - 1]]))


def cox_hr_from_survival(
    times: np.ndarray, s_trt: np.ndarray, s_ctrl: np.ndarray
) -> float:
    """Population equal-allocation Cox HR from two marginal KM curves.

    The calculation uses survival drops as event masses.  It is intended as a
    deterministic curve-level validation diagnostic, not as a replacement for
    the final censored-IPD Cox analysis.
    """
    s1 = np.minimum.accumulate(np.clip(s_trt, EPS, 1.0))
    s0 = np.minimum.accumulate(np.clip(s_ctrl, EPS, 1.0))
    d1 = np.maximum(s1[:-1] - s1[1:], 0.0)
    d0 = np.maximum(s0[:-1] - s0[1:], 0.0)
    y1 = s1[:-1]
    y0 = s0[:-1]

    def score(beta: float) -> float:
        eb = math.exp(beta)
        expected_fraction = y1 * eb / np.maximum(y0 + y1 * eb, EPS)
        return float(np.sum(d1 - expected_fraction * (d0 + d1)))

    lo, hi = -10.0, 10.0
    flo, fhi = score(lo), score(hi)
    if flo * fhi > 0:
        return math.nan
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = score(mid)
        if flo * fm <= 0:
            hi = mid
            fhi = fm
        else:
            lo = mid
            flo = fm
    return math.exp(0.5 * (lo + hi))
