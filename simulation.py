import numpy as np
import matplotlib.pyplot as plt
import os
import time

os.makedirs('plots', exist_ok=True)

# ==================== PARAMETERS (v5.0 - High-Res + Timer + Batching) ====================
N_FIL = 512
NUM_FILAMENTS = 12
NUM_REALIZATIONS = 30
BATCH_SIZE = 10
alpha = 1.22
steps = 300
dt = 0.002
core_base = 0.08
nu = 0.001
eps = 1.0

beta_tb = 0.8
gamma_kh = 0.6
delta_sft = 0.45
epsilon_neural = 0.3
THEORETICAL_DELTA = 0.0672

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

def compute_gauss_linking(filaments, Gamma_list):
    L = 0.0
    dls = [get_dl(f) for f in filaments]
    for i in range(len(filaments)):
        for j in range(i + 1, len(filaments)):
            r1, r2 = filaments[i], filaments[j]
            dl1, dl2 = dls[i], dls[j]
            R = r1[:, np.newaxis, :] - r2
            r3 = np.sum(R**2, axis=-1)**1.5 + 1e-12
            cross = np.cross(dl1[:, np.newaxis, :], dl2[np.newaxis, :, :])
            term = np.sum(cross * R, axis=-1) / r3
            L += (Gamma_list[i] * Gamma_list[j] / (4 * np.pi)) * np.sum(term) * ((2*np.pi/N_FIL)**2)
    return abs(L)

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

def dynamic_topological_proxies(filaments, L, E):
    curv = 0.0
    for r in filaments:
        dl = get_dl(r)
        d2l = get_dl(dl)
        curv += np.sum(np.linalg.norm(np.cross(dl, d2l), axis=1))
    TB = curv / 1000.0
    Kh = L * np.log(1 + E) * 0.3
    SFT = 0.0
    for r in filaments:
        dl = get_dl(r)
        d2l = get_dl(dl)
        kappa2 = np.sum(np.linalg.norm(np.cross(dl, d2l), axis=1)**2)
        SFT += kappa2
    SFT /= 10000.0
    Q = 0.8 * np.exp(-0.01 * E) + 0.2 * np.random.randn()
    return TB, Kh, SFT, Q

def multi_topological_weight(L, TB, Kh, SFT, Q):
    return 1 / (1 + alpha * L + beta_tb*abs(TB) + gamma_kh*Kh + delta_sft*SFT + epsilon_neural*Q)

def run_single_generic(with_depletion=True):
    filaments, Gamma_list = generate_generic_data()
    linking_hist = []
    enstrophy_hist = []
    t_hist = []
    for step in range(steps):
        t = step * dt
        t_hist.append(t)
        filaments = [adaptive_regrid(f) for f in filaments]
        L = compute_gauss_linking(filaments, Gamma_list)
        E = enstrophy_proxy(filaments, Gamma_list)
        TB, Kh, SFT, Q = dynamic_topological_proxies(filaments, L, E)
        linking_hist.append(L)
        scale = multi_topological_weight(L, TB, Kh, SFT, Q) if with_depletion else 1.0
        u = biot_savart_induced(filaments, Gamma_list, core_base)
        noise = np.random.randn(*u.shape) * np.sqrt(2 * nu * eps * dt)
        u_stoch = u + noise
        idx = 0
        for i in range(len(filaments)):
            n = len(filaments[i])
            filaments[i] += dt * u_stoch[idx:idx+n] * scale
            idx += n
        enstrophy_hist.append(E)
    return np.array(t_hist), np.array(enstrophy_hist), np.array(linking_hist)

def run_statistical_campaign(start_batch=0):
    print(f"Starting high-resolution stochastic campaign v5.0 (Batches of {BATCH_SIZE})...\n")
    overall_start = time.perf_counter()
    deltas = []
    suppressions = []
    num_batches = (NUM_REALIZATIONS + BATCH_SIZE - 1) // BATCH_SIZE
    for b in range(start_batch, num_batches):
        batch_start = time.perf_counter()
        batch_realizations = min(BATCH_SIZE, NUM_REALIZATIONS - b * BATCH_SIZE)
        print(f"Batch {b+1}/{num_batches} ({batch_realizations} realizations) starting...\n")
        for r in range(batch_realizations):
            global_r = b * BATCH_SIZE + r
            print(f"  Run {global_r+1:2d}/{NUM_REALIZATIONS}...", end=" ")
            t0 = time.perf_counter()
            t, E_with, L_with = run_single_generic(with_depletion=True)
            _, E_without, _ = run_single_generic(with_depletion=False)
            delta = (L_with[-1] - L_with[0]) / (t[-1] * np.mean(E_with))
            supp = np.max(E_without) / np.max(E_with) if np.max(E_with) > 0 else 1.0
            deltas.append(delta)
            suppressions.append(supp)
            elapsed = time.perf_counter() - t0
            print(f"δ={delta:.4f}, supp={supp:.1f}×  ({elapsed:.1f}s)")
        batch_time = time.perf_counter() - batch_start
        print(f"  Batch {b+1} complete — {batch_time/60:.1f} minutes\n")
    mean_delta = np.mean(deltas)
    std_delta = np.std(deltas)
    mean_supp = np.mean(suppressions)
    std_supp = np.std(suppressions)
    total_time = time.perf_counter() - overall_start
    print("=== FINAL STATISTICAL RESULTS (v5.0 High-Res) ===")
    print(f"Observed δ growth rate : {mean_delta:.4f} ± {std_delta:.4f}")
    print(f"Suppression factor     : {mean_supp:.1f} ± {std_supp:.1f}×")
    print(f"Total runtime          : {total_time/60:.1f} minutes")
    plt.figure(figsize=(10,6))
    plt.hist(suppressions, bins=15, alpha=0.7, color='blue', edgecolor='black')
    plt.axvline(mean_supp, color='red', linestyle='--', label=f'Mean = {mean_supp:.1f}×')
    plt.title('Distribution of Suppression Factors (30 High-Res Realizations)')
    plt.xlabel('Suppression Factor')
    plt.ylabel('Count')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plots/statistical_suppression_distribution_v5.png', dpi=400)
    print("Plot saved: plots/statistical_suppression_distribution_v5.png")
    return mean_delta, std_delta, mean_supp, std_supp

if __name__ == "__main__":
    print("Running simulation.py (v5.0 — high-res stochastic with timer + batching)")
    run_statistical_campaign(start_batch=0)
    print("\nFull campaign complete — project ready for final submission!")
