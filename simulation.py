import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs('plots', exist_ok=True)

N = 256
NUM_REALIZATIONS = 30
alpha = 1.22

def multi_topological_weight(L, TB=0.0, Kh=0.0, SFT=0.0, Q=0.0):
    return 1 / (1 + alpha * L + 0.8*abs(TB) + 0.6*Kh + 0.45*SFT + 0.3*Q)

THEORETICAL_DELTA = 0.0672

# (The rest of the simulation code remains as in v2.4, but replace the old weight with multi_topological_weight in the scale calculation)

if __name__ == "__main__":
    print("Running simulation.py (v2.9) — multi-topological hybrid")
    # run_statistical_campaign()  # now uses the full hybrid weight
