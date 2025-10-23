from global_LMI import global_LMI
from global_AS_LMI import global_AS_LMI
from local_LMI import local_LMI
from utils import add_ellipse

import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from matplotlib.patches import Ellipse
import numpy.linalg as la
import warnings
warnings.filterwarnings("ignore")
from scipy.integrate import solve_ivp

def smooth_sign(x, eps=1e-3):
    return x / np.sqrt(x**2 + eps**2)

def friction(thetaDot, model):
    return (
        model.m*model.g * (model.muC + (model.muS - model.muC) * np.exp(-(np.abs(thetaDot)/model.vs)**2)) * smooth_sign(thetaDot)
        + model.kv * thetaDot
    )

def dynamics(t, x, v_ref, model):
    v, z = x
    Ff = friction(v, model)
    dv = (-Ff - model.k * z) / model.m
    dz = v - v_ref
    return [dv, dz]

def simulate(v_ref, v0, z0, model, tmax=20.0, dt=0.001):
    t_eval = np.arange(0.0, tmax + dt, dt)
    
    sol = solve_ivp(
        lambda t, x: dynamics(t, x, v_ref, model),
        (0.0, tmax),
        [v0, z0],
        t_eval=t_eval,
        method='Radau',
        rtol=1e-6,
        atol=1e-8
    )
    
    t = sol.t
    v = sol.y[0, :]
    z = sol.y[1, :]
    u = np.zeros_like(t)
    
    return t, v, z, u

def add_ellipse(ax, P, v_ref, z_inf, color, label):
    if P is None:
        return
    eigvals, eigvecs = la.eigh(P)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    vx, vy = eigvecs[:, 0]
    angle = np.degrees(np.arctan2(vy, vx))
    a = 1.0 / np.sqrt(eigvals[0])
    b = 1.0 / np.sqrt(eigvals[1])
    ell = Ellipse(
        xy=(v_ref, z_inf), width=2*a, height=2*b, angle=angle,
        fill=False, edgecolor=color, linewidth=2, label=label
    )
    ax.add_patch(ell)


def run_test(test_id, model, SS):
    tests = {
        1: {"v_ref": 0.07, "ICs": [(0.1, -2.91), (0.2, -3.05)],
            "xlim": (-0.05, 0.3), "ylim": (-3.2, -2.75),
            "title": "test_1"},
        2: {"v_ref": 1.0, "ICs": [(6.0, 0.0), (10.0, 10.0)],
            "xlim": (-15, 15), "ylim": (-15, 10),
            "title": "test_2"},
        3: {"v_ref": 1.45, "ICs": [(1.7, -2.0), (1.75, -1.7)],
            "xlim": (-0.25, 3), "ylim": (-4, -1),
            "title": "test_3"},
        4: {"v_ref": 10.0, "ICs": [(6.0, 0.0), (10.0, 10.0)],
            "xlim": (-5, 25), "ylim": (-15, 15),
            "title": "test_4"}
    }

    cfg = tests[test_id]
    v_ref = cfg["v_ref"]

    trajectories = []
    for v0, z0 in cfg["ICs"]:
        t, v, z, _ = simulate(v_ref, v0, z0, model=model, tmax=20.0)
        trajectories.append((t, v, z, f"IC: ({v0}, {z0})"))

    Pg = global_LMI(SS, model, v_ref)
    Pl = local_LMI(SS.A, SS.B, SS.C, model, v_ref)
    is_global_AS = global_AS_LMI(SS, model, v_ref)

    print("Global asymptotic stability Thm3:", is_global_AS)

    fig, ax = plt.subplots(figsize=(10,6))  
    for (_, v, z, lbl) in trajectories:
        ax.plot(v, z, label=lbl)

    F_total_vref = friction(v_ref, model)
    z_inf = -F_total_vref / model.k

    add_ellipse(ax, Pg, v_ref, z_inf, 'blue', 'Pg Thm1 (inner approx)')
    add_ellipse(ax, Pl, v_ref, z_inf, 'red', 'Pl Thm2 (outer approx)')

    ax.set_xlim(cfg["xlim"])
    ax.set_ylim(cfg["ylim"])
    ax.set_xlabel("velocity v [m/s]")
    ax.set_ylabel("z [m]")
    ax.set_title(cfg["title"] + f", v_ref={v_ref}")
    ax.grid(True)
    ax.legend()
    plt.savefig("imagens/" + cfg["title"]+ f"_v2.png", dpi=300) 
    plt.show()

if __name__ == "__main__":
    model = SimpleNamespace(
        m=1.0, muC=0.2997, muS=0.5994, kv=1.0, vs=0.8,
        g=9.81, k=2.0, l0=0.0
    )
    SS = SimpleNamespace(
        A=np.array([[-model.kv/model.m, -model.k/model.m],[1, 0]]),
        B=np.array([[-1/model.m],[0]]),
        C=np.array([[1, 0]])
    )

    run_test(2, model, SS) 