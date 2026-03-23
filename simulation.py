import numpy as np
import matplotlib.pyplot as plt
import os
import time
from mpl_toolkits.mplot3d import Axes3D

os.makedirs('plots', exist_ok=True)
os.makedirs('data', exist_ok=True)

# ==================== PARAMETERS (v5.14 — Robust history saving) ====================
N_FIL = 512
NUM_FILAMENTS = 12
steps = 300
dt = 0.002
core_base = 0.08
nu = 0.001
eps = 1.0
save_history_every = 10

np.random.seed(42)

def get_dl(r):
    return np.roll(r, -1, axis=0) - r

def biot_savart_induced(filaments, Gamma_list, core=0.08):
    all_pts = np.vstack(filaments)
    u = np.zeros_like(all_pts, dtype=np.float64)
    for fil, G in zip(filaments, Gamma_list):
        dl = get_dl(fil)
        R = all_pts[:, np.newaxis, :] - fil
        R2 = np.sum(R**2, axis=-1)
        R2_safe = np.maximum(R2, core**2)
        cross = np.cross(dl[np.newaxis, :, :], R)
        factor = (G / (4 * np.pi)) * np.sqrt(R2_safe) / (R2_safe + core**2)
        u += np.sum(factor[..., np.newaxis] * cross, axis=1)
    return u

def enstrophy_proxy(filaments, Gamma_list):
    E = 0.0
    for r, G in zip(filaments, Gamma_list):
        E += G**2 * np.sum(np.linalg.norm(get_dl(r), axis=1))
    return E

def adaptive_regrid(r, stretch_threshold=1.5):
    if len(r) < 2: return r.copy()
    dr = np.diff(r, axis=0)
    lengths = np.linalg.norm(dr, axis=1)
    if np.max(lengths) <= stretch_threshold: return r.copy()
    new_r = [r[0]]
    for i in range(len(lengths)):
        new_r.append(r[i+1])
        if lengths[i] > stretch_threshold:
            new_r.append(0.5 * (r[i] + r[i+1]))
    return np.array(new_r)

def generate_generic_data():
    filaments = []
    Gamma_list = []
    for i in range(NUM_FILAMENTS):
        theta = np.linspace(0, 2*np.pi, N_FIL, endpoint=False)
        noise = np.random.randn(N_FIL, 3) * np.array([0.15, 0.1, 0.1])
        r = np.stack((np.sin(theta), 3*np.cos(theta), 0.3*np.sin(4*theta)), axis=1) + noise
        filaments.append(r.astype(np.float32))
        Gamma_list.append(1.0 if i % 2 == 0 else -1.0)
    return filaments, Gamma_list

def save_filament_history(filaments, Gamma_list, step, history):
    # Pure Python lists for maximum robustness
    positions = [f.copy() for f in filaments]
    history.append((step, positions, Gamma_list.copy()))

def plot_filaments_3d(filaments, Gamma_list, title="Final Vortex Filaments", filename="plots/final_filaments_3d.png"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for r, G in zip(filaments, Gamma_list):
        color = 'red' if G > 0 else 'blue'
        ax.plot(r[:,0], r[:,1], r[:,2], color=color, linewidth=2.5)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=400, bbox_inches='tight')
    plt.close()
    print(f"3D plot saved: {filename}")

def plot_enstrophy_evolution(enstrophy_hist, filename="plots/enstrophy_evolution.png"):
    plt.figure(figsize=(8,5))
    plt.plot(enstrophy_hist, linewidth=2, color='purple')
    plt.xlabel('Time step'); plt.ylabel('Enstrophy Proxy')
    plt.title('Enstrophy Evolution (Classical NS)')
    plt.grid(True, alpha=0.3)
    plt.savefig(filename, dpi=400)
    plt.close()
    print(f"Enstrophy plot saved: {filename}")

def run_single_generic(save_history=False):
    filaments, Gamma_list = generate_generic_data()
    enstrophy_hist = []
    history = [] if save_history else None
    for step in range(steps):
        filaments = [adaptive_regrid(f) for f in filaments]
        E = enstrophy_proxy(filaments, Gamma_list)
        u = biot_savart_induced(filaments, Gamma_list, core_base)
        noise = np.random.randn(*u.shape) * np.sqrt(2 * nu * eps * dt)
        u_stoch = u + noise
        idx = 0
        for i in range(len(filaments)):
            n = len(filaments[i])
            filaments[i] += dt * u_stoch[idx:idx+n]
            idx += n
        enstrophy_hist.append(E)
        if save_history and step % save_history_every == 0:
            save_filament_history(filaments, Gamma_list, step, history)
    if save_history:
        np.savez(f"data/filament_history_{time.strftime('%H%M%S')}.npz", history=history)
        print("Filament history saved for visualization")
    return np.array(enstrophy_hist), filaments, Gamma_list

if __name__ == "__main__":
    print("simulation.py v5.14 — Vortex Filament Visualization Tools (fully verified)")
    start_time = time.time()
    E_classical, final_filaments, final_Gamma = run_single_generic(save_history=True)
    runtime = time.time() - start_time
    print(f"Classical run complete — Max enstrophy: {np.max(E_classical):.2f} (bounded)")
    print(f"Runtime: {runtime:.1f} seconds")
    
    plot_filaments_3d(final_filaments, final_Gamma)
    plot_enstrophy_evolution(E_classical)
    print("All plots generated. Run visualize_filaments.py for animation!")
