import numpy as np

def lambda_min(physics, vref, rl):
    m = physics.m
    g = physics.g
    muS = physics.muS
    muC = physics.muC
    vs = physics.vs

    def f(x):
        return 1 - np.exp((x**2 + 2*x*vref)/vs**2) + 2*x*(x+vref)/vs**2

    epsilon = 3 
    i = max(0, int(np.ceil(-np.log10(vref/10))))
    epsilon = max(epsilon, i)
    v = -vref

    while i <= epsilon:
        if f(v) * f(v + 10**(-i)) <= 0:
            if v * (v + 10**(-i)) >= 1e-5:
                i += 1
                continue
        v = v + 10**(-i)

    v_star_minus = v
    v_star_plus = v + 10**(-epsilon)

    def partial_phi(eps1):
        return -2 * m * g * (muS - muC) * (eps1 + vref) / vs**2 * np.exp(-(eps1 + vref)**2 / vs**2)

    def phi(eps1):
        return m * g * (muS - muC) * (np.exp(-((eps1 + vref)/vs)**2) - np.exp(-(vref/vs)**2))

    if max(abs(v_star_minus), abs(v_star_plus)) <= rl:
        lambda_min_val = max(-partial_phi(v_star_minus), -partial_phi(v_star_plus))
    else:
        lambda_min_val = -min(-phi(-rl)/rl, phi(rl)/rl)

    return lambda_min_val
