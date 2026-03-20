import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

os.makedirs('plots', exist_ok=True)

# ==================== PARAMETERS (v2.3) ====================
N = 256                    # base points per filament
NUM_FILAMENTS = 6          # generic case
NUM_REALIZATIONS = 30      # statistical campaign
alpha = 1.22
nu = 0.001
dt = 0.002
T = 2.5
steps = int(T / dt)
core_base = 0.08

def get_dl(r):
    """Vectorized segment vectors"""
    return np.roll(r, -1, axis=0) - r

def biot_savart_induced(filaments, Gamma_list, core=0.08):
    """Rosenhead–Moore Biot-Savart"""
    all_pts = np.vstack(filaments)
    u = np.zeros_like(all_pts)
    for fil, G in zip(filaments, Gamma_list):
        dl = get_dl(fil)
        R = all_pts[:, np.newaxis, :] - fil
        R2 = np.sum(R**2, axis=-1)
        R2 = np.maximum(R2, core**2)
        cross = np.cross(dl[np.newaxis, :, :], R)
        factor = (G / (4 * np.pi)) * np.sqrt(R2) / (R2 + core**2)
        u += np.sum(factor[..., np.newaxis] * cross, axis=1)
    return u

def compute_gauss_linking(filaments, Gamma_list):
    """Exact discretized Gauss linking"""
    L = 0.0
    for i in range(len(filaments)):
        for j in range(i+1, len(filaments)):
            r1, r2 = filaments[i], filaments[j]
            dl1, dl2 = get_dl(r1), get_dl(r2)
            R = r1[:, np.newaxis, :] - r2
            r3 = np.sum(R**2, axis=-1)**1.5 + 1e-12
            cross = np.cross(dl1[:, np.newaxis, :], dl2[np.newaxis, :, :])
            term = np.sum(cross * R, axis=-1) / r3
            L += (Gamma_list[i] * Gamma_list[j] / (4 * np.pi)) * np.sum(term) * ((2*np.pi/N)**2)
    return abs(L)

def enstrophy_proxy(filaments, Gamma_list):
    """Enstrophy proxy"""
    E = 0.0
    for r, G in zip(filaments, Gamma_list):
        E += G**2 * np.sum(np.linalg.norm(get_dl(r), axis=1))
    return E

def adaptive_regrid(r, stretch_threshold=1.5):
    """Simple adaptive point insertion"""
    dr = np.diff(r, axis=0)
    lengths = np.linalg.norm(dr, axis=1)
    new_r = [r[0]]
    for i in range(len(lengths)):
        new_r.append(r[i+1])
        if lengths[i] > stretch_threshold:
            new_r.append(0.5 * (r[i] + r[i+1]))
    return np.array(new_r)

def generate_generic_data():
    """Random multi-filament initial data with non-zero helicity"""
    filaments = []
    Gamma_list = []
    for i in range(NUM_FILAMENTS):
        theta = np.linspace(0, 2*np.pi, N, endpoint=False)
        r = np.stack((np.sin(theta) + 0.15*np.random.randn(),
                      3*np.cos(theta) + 0.1*np.random.randn(),
                      0.3*np.sin(4*theta) + 0.1*np.random.randn()), axis=1)
        filaments.append(r.astype(np.float32))
        Gamma_list.append(1.0 if i % 2 == 0 else -1.0)
    return filaments, Gamma_list

def run_single_generic(with_depletion=True):
    filaments, Gamma_list = generate_generic_data()
    linking_hist = []
    enstrophy_hist = []
    t_hist = []
    
    for step in range(steps):
        t = step * dt
        t_hist.append(t)
        
        # Adaptive regridding
        for i in range(len(filaments)):
            filaments[i] = adaptive_regrid(filaments[i])
        
        L = compute_gauss_linking(filaments, Gamma_list)
        linking_hist.append(L)
        
        scale = 1.0 / (1 + alpha * max(L, 0.1)) if with_depletion else 1.0
        
        all_pts = np.vstack(filaments)
        u = biot_savart_induced(filaments, Gamma_list, core_base)
        
        idx = 0
        for i in range(len(filaments)):
            n = len(filaments[i])
            filaments[i] += dt * u[idx:idx+n] * scale
            idx += n
        
        E = enstrophy_proxy(filaments, Gamma_list)
        enstrophy_hist.append(E)
    
    return np.array(t_hist), np.array(enstrophy_hist), np.array(linking_hist)

def run_statistical_campaign():
    """Run 30 realizations and compute statistics"""
    print(f"Starting statistical campaign ({NUM_REALIZATIONS} realizations, N={N})...\n")
    deltas = []
    suppressions = []
    
    for r in range(NUM_REALIZATIONS):
        print(f"  Run {r+1}/{NUM_REALIZATIONS}...", end=" ")
        t, E_with, L_with = run_single_generic(with_depletion=True)
        _, E_without, _ = run_single_generic(with_depletion=False)
        
        delta = (L_with[-1] - L_with[0]) / (t[-1] * np.mean(E_with))
        supp = np.max(E_without) / np.max(E_with) if np.max(E_with) > 0 else 1.0
        
        deltas.append(delta)
        suppressions.append(supp)
        print(f"δ={delta:.4f}, supp={supp:.1f}×")
    
    mean_delta = np.mean(deltas)
    std_delta = np.std(deltas)
    mean_supp = np.mean(suppressions)
    std_supp = np.std(suppressions)
    
    print("\n=== STATISTICAL RESULTS ===")
    print(f"δ growth rate : {mean_delta:.4f} ± {std_delta:.4f}")
    print(f"Suppression   : {mean_supp:.1f} ± {std_supp:.1f}×")
    
    # Save summary plot
    plt.figure(figsize=(10,6))
    plt.hist(suppressions, bins=15, alpha=0.7, color='blue', edgecolor='black')
    plt.axvline(mean_supp, color='red', linestyle='--', label=f'Mean = {mean_supp:.1f}×')
    plt.title('Distribution of Suppression Factors (30 Generic Realizations)')
    plt.xlabel('Suppression Factor')
    plt.ylabel('Count')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plots/statistical_suppression_distribution.png', dpi=400)
    print("Plot saved: plots/statistical_suppression_distribution.png")
    
    return mean_delta, std_delta, mean_supp, std_supp

if __name__ == "__main__":
    print("Running latest simulation.py (v2.3) — statistical campaign")
    run_statistical_campaign()
    print("\nAll done! Files ready to save.")
