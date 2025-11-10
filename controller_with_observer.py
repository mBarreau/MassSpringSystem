import numpy as np
from types import SimpleNamespace
from scipy.signal import place_poles
import warnings
warnings.filterwarnings("ignore")

from friction import friction as F, F_nl
from local_LMI import local_LMI_cl as local_LMI
from utils import plot_results_observer, plot_results  

def equilibrium_z(v_ref):
    F_total = F(v_ref, parameters)
    return -F_total / parameters.k

def closed_loop_dynamics_linearized(vref, Kp):
    A = SS.A.copy(); B = SS.B.copy(); C = SS.C.copy()
    m, g, muS, muC, vs = parameters.m, parameters.g, parameters.muS, parameters.muC, parameters.vs
    Gamma = -2.0 * m * g * (muS - muC) * (vref / (vs**2)) * np.exp(-(vref**2) / (vs**2))
    Krow = np.array([[Kp, 0.0]]) 
    Acl = A + (B * Gamma) @ C - B @ Krow     
    return Acl

def closed_loop_dynamics_nonlinear(v_ref, Kp, eps):
    phi = F_nl(eps[0] + v_ref, parameters) - F_nl(v_ref, parameters)   
    deps = SS.A @ eps.reshape(-1,1) + SS.B * (phi - Kp * eps[0])
    return deps.flatten()

def closed_loop_dynamics_nonlinear_observer(v_ref, v, Kp, L, eps, eps_hat):
    phi = F_nl(v, parameters) - F_nl(v_ref, parameters)
    deps = SS.A @ eps.reshape(-1,1) + SS.B * (phi - Kp * eps_hat[0])
    deps_hat = SS.A @ eps_hat.reshape(-1,1) + SS.B * (phi - Kp * eps_hat[0]) + L @ SS.C @ (eps.reshape(-1,1) - eps_hat.reshape(-1,1))
    return deps.flatten(), deps_hat.flatten()

def simulate_stage_P(Kp, v0, v_ref, z0, z_inf, next_entry):
    steps = int(np.ceil(settings.tmax/settings.dt))
    t_hist, v_hist, z_hist, u_hist = [], [], [], []
    eps = np.array([v0 - v_ref, z0 - z_inf])

    for k in range(steps):
        t = k * settings.dt
        v, z = eps[0] + v_ref, eps[1] + z_inf
        deps = closed_loop_dynamics_nonlinear(v_ref, Kp, eps)
        u = -Kp * eps[0]
        eps = eps + deps * settings.dt
        
        t_hist.append(t); v_hist.append(v); z_hist.append(z); u_hist.append(u)

    return np.array(t_hist), np.array(v_hist), np.array(z_hist), np.array(u_hist), v, z

def simulate_stage_P_observer(Kp, L, v0, v_ref, z0, z_inf, next_entry):
    steps = int(np.ceil(settings.tmax / settings.dt))
    t_hist, v_hist, z_hist, zhat_hist, u_hist, utrue_hist = [], [], [], [], [], []
    eps = np.array([v0 - v_ref, z0 - z_inf])
    eps_hat = np.zeros_like(eps)

    for k in range(steps):
        t = k * settings.dt
        v = eps[0] + v_ref
        z = eps[1] + z_inf
        deps, deps_hat = closed_loop_dynamics_nonlinear_observer(v_ref, v, Kp, L, eps, eps_hat)
        utrue = -Kp * eps[0]
        eps = eps + deps * settings.dt 
        u = -Kp * eps_hat[0]
        eps_hat = eps_hat + deps_hat * settings.dt

        zhat = eps_hat[1] + z_inf
        t_hist.append(t); v_hist.append(v); z_hist.append(z); zhat_hist.append(zhat); u_hist.append(u); utrue_hist.append(utrue)

    return np.array(t_hist), np.array(v_hist), np.array(z_hist), np.array(zhat_hist), np.array(u_hist), np.array(utrue_hist), v, z

def entries_constructor(v_start, v_target, SS, Kp):
    direction = np.sign(v_target - v_start)
    step = settings.initial_step
    v_current = v_start 
    entries = []

    A_cl_current = closed_loop_dynamics_linearized(v_current, Kp)
    Pl_current = local_LMI(A_cl_current, SS.B, SS.C, parameters, v_current, Kp)
    if Pl_current is None or not np.all(np.real(np.linalg.eigvals(A_cl_current)) < 0): return []
    entries.append({'v_ref': v_current, 'P': Pl_current, 'z_inf': equilibrium_z(v_current)})

    while direction * (v_target - v_current) > 0:
        v_candidate = np.clip(v_current + direction * step, min(v_current, v_target), max(v_current, v_target))
        z_candidate = equilibrium_z(v_candidate)
        A_cl_candidate = closed_loop_dynamics_linearized(v_candidate, Kp)
        Pl_candidate = local_LMI(A_cl_candidate, SS.B, SS.C, parameters, v_candidate, Kp)
        if Pl_candidate is None or not np.all(np.real(np.linalg.eigvals(A_cl_candidate)) < 0): return []
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

def controller(v_target, v0, z0, SS, Kp, use_observer):
    entries = entries_constructor(v0, v_target, SS, Kp)
    T_all, V_all, Z_all, U_all, switch_times, switch_vrefs = [], [], [], [], [], []
    
    Zhat_all, Utrue_all = [], []
    L = None
    
    if use_observer:
        L = augmented_stability_check(SS, Kp, L, v_target)

    v, z = v0, z0
    time_offset = 0.0

    for idx, entry in enumerate(entries):
        v_ref = entry['v_ref']
        z_inf = entry['z_inf']
        next_entry = entries[idx+1] if idx+1 < len(entries) else None

        if use_observer:
            t, vhist, zhist, zhat_hist, uhist, utrue_hist, v, z = simulate_stage_P_observer(Kp, L, v, v_ref, z, z_inf, next_entry)
            Utrue_all.extend(utrue_hist.tolist()); Zhat_all.extend(zhat_hist.tolist())
        else:
            t, vhist, zhist, uhist, v, z = simulate_stage_P(Kp, v, v_ref, z, z_inf, next_entry)

        T_all.extend((t + time_offset).tolist()); V_all.extend(vhist.tolist()); Z_all.extend(zhist.tolist()); U_all.extend(uhist.tolist())
        switch_times.append(T_all[-1]); switch_vrefs.append(v_ref)
        time_offset = T_all[-1] + settings.dt

    if use_observer:
        return np.array(T_all), np.array(V_all), np.array(Z_all), np.array(Zhat_all), np.array(U_all), np.array(Utrue_all), switch_times, switch_vrefs, entries
    else:
        return np.array(T_all), np.array(V_all), np.array(Z_all), np.array(U_all), switch_times, switch_vrefs, entries

def augmented_stability_check(SS, Kp, L, v_ref):
    desired_poles = [-1+3j, -1-3j]
    A_cl = closed_loop_dynamics_linearized(v_ref, Kp)
    L = place_poles(A_cl.T, SS.C.T, desired_poles).gain_matrix.T
    A_e = A_cl - L @ SS.C
    upper = np.hstack([A_cl, -SS.B @ np.array([[Kp]]) @ SS.C])
    lower = np.hstack([np.zeros_like(A_cl), A_e])
    A_aug = np.vstack([upper, lower])
    eigvals = np.real(np.linalg.eigvals(A_aug))
    stable = np.all(eigvals < 0) 
    print(f"Augmented system Hurwitz: {stable}")   
    return L

if __name__ == "__main__":
    parameters = SimpleNamespace(m=1.0, muC=0.2997, muS=0.5994, kv=1.0, vs=0.8, g=9.81, k=2.0, l0=0.0)
    SS = SimpleNamespace(A=np.array([[-parameters.kv/parameters.m, -parameters.k/parameters.m],
                                     [1.0, 0.0]]),
                         B=np.array([[-1.0/parameters.m],[0.0]]),
                         C=np.array([[1.0,0.0]]))
    
    settings = SimpleNamespace(
        dt=0.001,
        initial_step=4.0,
        backtrack_factor=0.75,
        tmax=10.0,
        tol_rel=1e-3,
    )

    v0, z0 = 24, -5
    v_target = 1.7
    Kp = 0.1
    use_observer = False  

    print(f"Running controller (observer = {use_observer}) ...")
    
    if use_observer:
        T,V,Z,Zhat,U,Utrue,switch_times,switch_vrefs,entries = controller(v_target,v0,z0,SS,Kp,use_observer)
        plot_results_observer(T,V,Z,Zhat,U,Utrue,v_target,switch_times,switch_vrefs,entries)
    else:
        T,V,Z,U,switch_times,switch_vrefs,entries = controller(v_target,v0,z0,SS,Kp,use_observer)
        plot_results(T,V,Z,U,v_target,switch_times,switch_vrefs,entries)