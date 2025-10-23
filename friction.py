import numpy as np
def smooth_sign(x, eps=1e-3):
    return x / np.sqrt(x**2 + eps**2)

def F_nl(thetaDot, model):
    return model.m*model.g * (model.muC + (model.muS - model.muC) *
                              np.exp(-(np.abs(thetaDot)/model.vs)**2)) * smooth_sign(thetaDot)

def friction(thetaDot, model):
    return (
        model.m*model.g * (model.muC + (model.muS - model.muC) * np.exp(-(np.abs(thetaDot)/model.vs)**2)) * smooth_sign(thetaDot)
        + model.kv * thetaDot
    )