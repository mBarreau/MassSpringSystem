import numpy as np
from scipy.optimize import fsolve

def A0Hurwitz(model):
    m = model.m
    g = model.g
    kv = model.kv
    vs = model.vs
    muS = model.muS
    muC = model.muC

    def theta(v):
        return v * np.exp(-v**2 / vs**2) - kv * vs**2 / (2 * m * g * (muS - muC))

    try:
        vref1 = fsolve(theta, 0, xtol=1e-10)[0]
    except:
        return np.nan, np.nan

    if np.isnan(vref1):
        return np.nan, np.nan

    v0 = vs * np.sqrt(5)
    for _ in range(6):
        try:
            vref2 = fsolve(theta, v0, xtol=1e-10)[0]
        except:
            vref2 = np.nan

        if np.isnan(vref2) or vref2 <= 1.001 * vref1:
            v0 *= 1.1
        else:
            return vref1, vref2

    return vref1, vref1