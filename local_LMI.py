import numpy as np
import cvxpy as cp
from lambda_min import lambda_min
from A0Hurwitz import A0Hurwitz

def local_LMI(A,B,C, model, vref):
    m, g, kv, k = model.m, model.g, model.kv, model.k
    muS, muC, vs = model.muS, model.muC, model.vs

    vref1, vref2 = A0Hurwitz(model)
    if vref1 <= vref <= vref2:
        return None  

    Pl = cp.Variable((2, 2), symmetric=True)
    tau = cp.Variable(nonneg=True)
    eta = cp.Variable(nonneg=True)

    """ A = SS.A
    B = SS.B
    C = SS.C """

    Gamma = -2*m*g*(muS - muC)*vref*np.exp(-vref**2 / vs**2)/vs**2
    A0 = A #A + B @ (Gamma * C)

    PlMin = None
    etaMin = np.inf

    rlList = np.linspace(vref, 1e-5, 5)
    for rl in rlList:
        lam = lambda_min(model, vref, rl)

        ThetaLocal = cp.bmat([
            [A0.T @ Pl + Pl @ A0 - 2*tau*Gamma*(Gamma + lam)*(C.T @ C), Pl @ B - tau*(2*Gamma + lam)*C.T],
            [(Pl @ B - tau*(2*Gamma + lam)*C.T).T, cp.reshape(-2*tau, (1,1))]
        ])

        ThetaIncl = cp.bmat([
            [Pl, C.T],
            [C, cp.reshape(rl**2, (1,1))]
        ])
        constraints = [
            Pl >> 0,
            ThetaLocal << -1e-5*np.eye(ThetaLocal.shape[0]),
            ThetaIncl >> 0,
            Pl << eta*np.eye(2)
        ]

        prob = cp.Problem(cp.Minimize(eta), constraints)
        try:
            prob.solve(solver=cp.SCS, verbose=False)
        except cp.SolverError:
            print(f"Solver failed at rl={rl:.5e}")
            continue

        if prob.status == cp.OPTIMAL:
            if eta.value < etaMin:
                PlMin = Pl.value
                etaMin = eta.value
            else:
                break

    return PlMin
