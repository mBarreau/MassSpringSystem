import numpy as np
import cvxpy as cp
from lambda_min import lambda_min
from friction import friction  
import math

def global_AS_LMI(SS, model, vref):
    globalAsympStab = 0

    FS = model.muS * model.m * model.g
    FC = model.muC * model.m * model.g

    Pg = cp.Variable((2, 2), symmetric=True)
    Pl = cp.Variable((2, 2), symmetric=True)
    tau = cp.Variable(nonneg=True)
    tau1 = cp.Variable(nonneg=True)
    tau2 = cp.Variable(nonneg=True)
    tau3 = cp.Variable(nonneg=True)
    tau4 = cp.Variable(nonneg=True)
    tau5 = cp.Variable(nonneg=True)

    A = SS.A
    B = SS.B
    C = SS.C

    Fnl = friction(vref, model) - model.kv * vref
    D = np.hstack([A, B, -Fnl * B])
    F = np.eye(2, 4)

    Gamma = -2 * model.m * model.g * (model.muS - model.muC) * vref * np.exp(-vref**2 / model.vs**2) / model.vs**2
    A0 = A + B @ (Gamma * C)

    # Theta for global LMIs
    Theta = D.T @ Pg @ F + F.T @ Pg @ D

    # Pi matrices
    pi1 = np.array([[1, 0, 0, 0]])
    pi3 = np.array([[0, 0, 1, 0]])
    pi4 = np.array([[0, 0, 0, 1]])

    Pi0 = pi4.T @ pi4 - F.T @ Pg @ F
    Pi1 = pi3.T @ pi3 - FS**2 * (pi4.T @ pi4)
    Pi2 = -pi3.T @ pi3 + FC**2 * (pi4.T @ pi4)
    Pi3 = -(pi1 + vref*pi4).T @ pi3
    Pi3 = Pi3 + Pi3.T
    Pi4 = (pi1 + vref*pi4).T @ (pi1 + vref*pi4)

    tau0List = np.linspace(0, -2*np.max(np.real(np.linalg.eigvals(A))), 30)
    rlList = np.linspace(vref, 1e-5, 30)

    for rl in rlList:
        lam = lambda_min(model, vref, rl)

        ThetaLocal = cp.bmat([
    [A0.T @ Pl + Pl @ A0 - 2*tau*Gamma*(Gamma + lam)*(C.T @ C),
     Pl @ B - tau*(2*Gamma + lam)*C.T],
    [B.T @ Pl - tau*(2*Gamma + lam)*C,
     cp.reshape(-2*tau, (1, 1))]
])


        ThetaIncl = cp.bmat([
            [Pl, C.T],
            [C, rl**2 * np.ones((1,1))]
        ])

        for tau0 in tau0List:
            ThetaBar = Theta - tau0*Pi0 - tau1*Pi1 - tau2*Pi2 - tau3*Pi3
            ThetaBar2 = Theta - tau0*Pi0 - tau5*Pi1 - tau4*Pi4

            constraints = [
                Pl >= 0,
                Pg - Pl >= 0,
                tau1 >= 0,
                tau5 >= 0,
                tau2 >= 0,
                tau3 >= 0,
                ThetaBar << -1e-5 * np.eye(ThetaBar.shape[0]),
                ThetaBar2 << -1e-5 * np.eye(ThetaBar2.shape[0]),
                tau >= 0,
                ThetaLocal <= -1e-5 * np.eye(ThetaLocal.shape[0]),
                ThetaIncl >> 0
            ]

            prob = cp.Problem(cp.Minimize(0), constraints)
            prob.solve(solver=cp.SCS, verbose=True)

            if prob.status == cp.OPTIMAL:
                globalAsympStab = 1
                return globalAsympStab

    return globalAsympStab
