import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import warnings
warnings.filterwarnings("ignore")

from utils import plot_results
from friction import friction as F, F_nl
from local_LMI import local_LMI

def equilibrium_z(v_ref):
    F_total = F(v_ref, parameters)
    return -F_total / parameters.k

def closed_loop_dynamics_linearized(SS, vref, Kp):
    A = SS.A.copy(); B = SS.B.copy(); C = SS.C.copy()
    m, g, muS, muC, vs = parameters.m, parameters.g, parameters.muS, parameters.muC, parameters.vs
    Gamma = -2.0 * m * g * (muS - muC) * (vref / (vs**2)) * np.exp(-(vref**2) / (vs**2))
    Krow = np.array([[Kp, 0.0]]) 
    Acl = A + (B * Gamma) @ C - B @ Krow     
    return Acl

def closed_loop_dynamics_nonlinear(SS, v_ref, Kp, eps):
    A = SS.A.copy(); B = SS.B.copy(); C = SS.C.copy()
    eps1, eps2 = eps
    u = -Kp * eps1
    phi = F_nl(eps1 + v_ref, parameters) - F_nl(v_ref, parameters)   
    deps = A @ eps.reshape(-1,1) + B * (phi + u)
    return deps.flatten(), u

def simulate_stage_P(Kp, v0, v_ref, z0, z_inf,next_entry):
    #apply control until state is inside next ROA
    steps = int(np.ceil(settings.tmax/settings.dt))
    t_hist, v_hist, z_hist, u_hist = [], [], [], []
    eps = np.array([v0-v_ref, z0-z_inf])

    for k in range(steps):
        t = k*settings.dt
        deps, u = closed_loop_dynamics_nonlinear(SS, v_ref, Kp, eps)
        eps = eps + deps*settings.dt
        v, z = eps[0]+v_ref, eps[1]+z_inf
        t_hist.append(t); v_hist.append(v); z_hist.append(z); u_hist.append(u)

        #if v converged, check if its inside the next ROA
        """ if next_entry is not None and max(abs(v-v_ref), abs(z-z_inf)) < settings.tol_rel:
            eps_next = np.array([v - next_entry['v_ref'], z - next_entry['z_inf']])
            if float(eps_next.T @ next_entry['P'] @ eps_next) <= 1: 
                print(f"Entered next ROA at t={t:.3f}s")
                break """
        
        if next_entry is not None:
            eps_next = np.array([v - next_entry['v_ref'], z - next_entry['z_inf']])
            if float(eps_next.T @ next_entry['P'] @ eps_next) <= 1: 
                print(f"Entered next ROA at t={t:.3f}s")
                break

    return (np.array(t_hist), np.array(v_hist), np.array(z_hist),
            np.array(u_hist), v, z)

def entries_constructor(v_start, v_target, SS, Kp):
    #add a new vref such that the last vref is contained in the region of attraction of this vref
    direction = np.sign(v_target - v_start)
    step = settings.initial_step
    v_current = v_start 
    entries = []

    #Pl must exist and Acl must be Schur
    A_cl_current = closed_loop_dynamics_linearized(SS, v_current, Kp)
    Pl_current = local_LMI(A_cl_current, SS.B, SS.C, parameters, v_current)
    if Pl_current is None or not np.all(np.real(np.linalg.eigvals(A_cl_current)) < 0): 
        return []
    entries.append({'v_ref': v_current, 'P': Pl_current, 'z_inf': equilibrium_z(v_current)})

    while direction*(v_target - v_current) > 0:
        v_candidate = np.clip(v_current + direction*step, min(v_current,v_target), max(v_current,v_target))
        z_candidate = equilibrium_z(v_candidate)

        A_cl_candidate = closed_loop_dynamics_linearized(SS, v_candidate, Kp)
        Pl_candidate = local_LMI(A_cl_candidate, SS.B, SS.C, parameters, v_candidate)
        if Pl_candidate is None or not np.all(np.real(np.linalg.eigvals(A_cl_candidate)) < 0): 
            print(Pl_candidate)
            return []
        eps = np.array([entries[-1]['v_ref'] - v_candidate, entries[-1]['z_inf'] - z_candidate])

        if float(eps.T @ Pl_candidate @ eps) <= 1:
            entries.append({'v_ref': v_candidate, 'P': Pl_candidate, 'z_inf': z_candidate})
            v_current = v_candidate
            print(f"Added stage v_ref={v_candidate:.3f}")
        else:
            step *= settings.backtrack_factor

        if np.isclose(v_current, v_target):
            break

    return entries

def controller(v_target, v0, z0, SS, Kp):
    dt = 0.001
    entries = entries_constructor(v0, v_target, SS, Kp)

    T_all, V_all, Z_all, U_all, switch_times, switch_vrefs = [], [], [], [], [], []
    v, z = v0, z0
    time_offset = 0.0

    for idx, entry in enumerate(entries):
        v_ref = entry['v_ref']
        z_inf = entry['z_inf']
        next_entry = entries[idx + 1] if idx + 1 < len(entries) else None

        t, vhist, zhist, uhist, v, z = simulate_stage_P(Kp, v, v_ref, z, z_inf, next_entry)

        T_all.extend((t + time_offset).tolist()), V_all.extend(vhist.tolist()), Z_all.extend(zhist.tolist()), U_all.extend(uhist.tolist()), switch_times.append(T_all[-1]), switch_vrefs.append(v_ref)
        time_offset = T_all[-1] + dt

    return (np.array(T_all), np.array(V_all), np.array(Z_all), np.array(U_all), switch_times, switch_vrefs, entries)


if __name__ == "__main__":
    parameters = SimpleNamespace(m=1.0, muC=0.2997, muS=0.5994, kv=1.0, vs=0.8, g=9.81, k=2.0, l0=0.0)
    
    SS = SimpleNamespace(A=np.array([[-parameters.kv/parameters.m, -parameters.k/parameters.m],
                                     [1.0, 0.0]]),
                         B=np.array([[-1.0/parameters.m],[0.0]]),
                         C=np.array([[1.0,0.0]]))
    
    settings = SimpleNamespace(dt = 0.001, initial_step=4.0, backtrack_factor=0.75, tmax=10.0, tol_rel=1e-3)
    
    v0, z0 = 12, -5
    v_target = 1.7
    Kp = 0.1
    
    result = controller(v_target, v0, z0, SS, Kp)

    T,V,Z,U,switch_times,switch_vrefs,Pl_entries = result
    print("Simulation finished, plotting results...")
    plot_results(T,V,Z,U,v_target,switch_times,switch_vrefs,Pl_entries)

    """ 
    right now it´s full state model with output feedback, we could put an observer on z
    does it need to converge fully before changing target?
    we need to show that theorem 2 (lemmas 2 and 3) apply to new closed loop A0
    proof by induction (base case and then inclusion)
    only runs for closed loop A Hurwitz but create alternatives like adapting Kp, this is lemma 2
    it is not gain scheduling, Kp is constant
    low velocities don´t work well, should work until 1.25, lower vel require lower Kp so it would make sense to do gain scheduling and lower the gain
    vale a pena discretizar?
    """