from friction import friction 

import numpy as np
import cvxpy as cp
import math

def global_LMI(SS, model, vref, tol_eps=1e-6):
    FS = model.muS * model.m * model.g
    FC = model.muC * model.m * model.g

    Pg = cp.Variable((2, 2), symmetric=True)
    tau0 = None                 
    tau1 = cp.Variable(nonneg=True)
    tau2 = cp.Variable(nonneg=True)
    tau3 = cp.Variable(nonneg=True)
    tau4 = cp.Variable()        
    tau5 = cp.Variable(nonneg=True)
    eta = cp.Variable()

    A = SS.A
    B = SS.B


    Fnl = friction(vref, model) - model.kv * vref
    D = np.hstack([A, B, -Fnl * B])
    F = np.eye(2, 4)
    #F = np.hstack([np.eye(2), np.zeros((2,2))])  


    Theta = D.T @ Pg @ F + F.T @ Pg @ D
   

    pi1 = np.array([[1, 0, 0, 0]])
    pi3 = np.array([[0, 0, 1, 0]])
    pi4 = np.array([[0, 0, 0, 1]])

    Pi0 = pi4.T @ pi4 - F.T @ Pg @ F
    Pi1 = pi3.T @ pi3 - (FS**2) * (pi4.T @ pi4)
    Pi2 = (FC**2) * (pi4.T @ pi4) - pi3.T @ pi3
    tmp = (pi1 + vref * pi4).T @ pi3      
    Pi3 = -(tmp + tmp.T)
    Pi4 = (pi1 + vref * pi4).T @ (pi1 + vref * pi4)
    #Pi4 = (math.pi + vref * pi4).T @ (math.pi + vref * pi4)

    eigA = np.linalg.eigvals(A)
    tau0_max = -2 * np.max(np.real(eigA))
    tau0_list = np.linspace(0, tau0_max, 30)

    best_Pg = None
    best_eta = -1e12

    eps = max(tol_eps, 1e-10)


    for tau0 in tau0_list:
        ThetaBar = Theta - tau0 * Pi0 - tau1 * Pi1 - tau2 * Pi2 - tau3 * Pi3
        ThetaBar2 = Theta - tau0 * Pi0 - tau5 * Pi1 - tau4 * Pi4

        n = Pg.shape[0]
        constraints = [
            Pg - eta * np.eye(2) >> 0,
            tau1 >= 0,
            tau2 >= 0,
            tau3 >= 0,
            # tau4 unconstrained (do NOT add >=0)
            tau5 >= 0,
            ThetaBar << -eps * np.eye(ThetaBar.shape[0]),
            ThetaBar2 << -eps * np.eye(ThetaBar2.shape[0]),
            Pg >> 1e-8 * np.eye(2)   # ensure strictly positive definite
        ]

        prob = cp.Problem(cp.Maximize(eta), constraints)
        prob.solve(solver=cp.SCS, verbose=True)

        if prob.status in ["optimal", "optimal_inaccurate"]:
            val = eta.value if eta.value is not None else -1e12
            if val > best_eta:
                best_eta = float(val)
                best_Pg = Pg.value.copy()

    return best_Pg