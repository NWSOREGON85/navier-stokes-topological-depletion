import numpy as np
import matplotlib.pyplot as plt
import os
import time

os.makedirs('plots', exist_ok=True)

N_FIL = 512
NUM_FILAMENTS = 12
NUM_REALIZATIONS = 30
BATCH_SIZE = 10
alpha = 1.22
epsilon_0 = 1e-3
gamma = 0.5
steps = 300
dt = 0.002
core_base = 0.08
nu = 0.001
eps = 1.0

beta_tb = 0.8
gamma_kh = 0.6
delta_sft = 0.45
epsilon_neural = 0.3

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

def multi_topological_weight(L, TB, Kh, SFT, Q, E, integral_Htop, worst_case_mode=False):
    H_top = L + beta_tb*abs(TB) + gamma_kh*Kh + delta_sft*SFT + epsilon_neural*Q
    if worst_case_mode:
        H_top = 0.0
    phi = integral_Htop / (1 + integral_Htop)
    eps_t = epsilon_0 * (1 - phi) * np.exp(-gamma * integral_Htop)
    H_top += eps_t * E
    return 1 / (1 + alpha * H_top)

def reconnect_filaments(filaments, Gamma_list, reconnect_dist=0.15):
    new_filaments = [f.copy() for f in filaments]
    new_Gamma = Gamma_list.copy()
    for i in range(len(filaments)):
        for j in range(i + 1, len(filaments)):
            dists = np.linalg.norm(filaments[i][:, np.newaxis] - filaments[j], axis=-1)
            if np.min(dists) < reconnect_dist:
                tail_i = new_filaments[i][-1].copy()
                tail_j = new_filaments[j][-1].copy()
                new_filaments[i][-1] = tail_j
                new_filaments[j][-1] = tail_i
    return new_filaments, new_Gamma

def run_single_generic(with_depletion=True, worst_case_mode=False, integral_Htop=0.0):
    filaments, Gamma_list = generate_generic_data()
    linking_hist = []
    enstrophy_hist = []
    t_hist = []
    for step in range(steps):
        t = step * dt
        t_hist.append(t)
        filaments = [adaptive_regrid(f) for f in filaments]
        filaments, Gamma_list = reconnect_filaments(filaments, Gamma_list)
        L = compute_gauss_linking(filaments, Gamma_list)
        E = enstrophy_proxy(filaments, Gamma_list)
        TB, Kh, SFT, Q = dynamic_topological_proxies(filaments, L, E)
        linking_hist.append(L)
        scale = multi_topological_weight(L, TB, Kh, SFT, Q, E, integral_Htop, worst_case_mode) if with_depletion else 1.0
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

if __name__ == "__main__":
    print("simulation.py v5.6 — Triggered floor + Helicity Vacuum Paradox + Reconnection")
    t, E_with, L_with = run_single_generic(with_depletion=True, worst_case_mode=False)
    print(f"Test run complete — Max enstrophy: {np.max(E_with):.2f}")
