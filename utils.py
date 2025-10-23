import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def add_ellipse(ax, P, v_ref, z_inf, color='orange', linewidth=1.0, label=None):
    if P is None:
        return
    vals, vecs = np.linalg.eigh(P)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    a = 1.0 / np.sqrt(vals[0])
    b = 1.0 / np.sqrt(vals[1])
    vx, vy = vecs[:, 0]
    angle = np.degrees(np.arctan2(vy, vx))
    ell = Ellipse(xy=(v_ref, z_inf), width=2 * a, height=2 * b,
                  angle=angle, edgecolor=color, facecolor='none',
                  linewidth=linewidth, label=label)
    ax.add_patch(ell)


def plot_results(T, V, Z, U, v_target, switch_times, switch_vrefs, Pl_entries):
    tab_colors = plt.get_cmap('tab20').colors
    colors = np.array([tab_colors[i % 20] for i in range(len(Pl_entries))])
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    axs[0].plot(T, V, label="gain-scheduled v(t)")
    axs[0].axhline(v_target, linestyle='--', color='k', label='v_target')
    axs[0].scatter(switch_times, switch_vrefs, color='red', marker='o', label='stage switch')
    axs[0].set_title("Velocity vs time")
    axs[0].set_xlabel("t [s]")
    axs[0].set_ylabel("v [m/s]")
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(T, U, label="u(t)", color="tab:green")
    axs[1].set_title("Control input")
    axs[1].set_xlabel("t [s]")
    axs[1].set_ylabel("u")
    axs[1].legend()
    axs[1].grid(True)

    switch_indices = []
    for t in switch_times:
        idx_t = np.searchsorted(T, t, side='right') - 1
        switch_indices.append(max(0, idx_t))

    for i, entry in enumerate(Pl_entries):
        color = colors[i] if i < len(colors) else 'gray'
        add_ellipse(axs[2], entry['P'], entry['v_ref'], entry['z_inf'],
                    color=color, linewidth=1.2,
                    label=f"Pl v={entry['v_ref']:.3g}")
        seg_start = switch_indices[i-1] if i-1 >= 0 and i-1 < len(switch_indices) else 0
        seg_end = switch_indices[i] if i < len(switch_indices) else len(T)-1
        seg_start = max(0, seg_start)
        seg_end = min(len(V)-1, seg_end)
        if seg_end > seg_start:
            axs[2].plot(V[seg_start:seg_end+1], Z[seg_start:seg_end+1],
                        '-', color=color, linewidth=1.6)

    axs[2].set_title("Phase portrait with ROA ellipses (trajectory segments colored to match Pl)")
    axs[2].set_xlabel("v")
    axs[2].set_ylabel("z")
    axs[2].grid(True)
    axs[2].legend(loc='best', fontsize='small')
    """ axs[2].set_ylim(-5,0)
    axs[2].set_xlim(0, 5) """
    plt.show()
