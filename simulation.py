import numpy as np
import matplotlib.pyplot as plt
import os

os.makedirs('plots', exist_ok=True)

N = 256
NUM_REALIZATIONS = 30
alpha = 1.22
beta = 0.8  # Thurston-Bennequin coefficient

def hybrid_weight(L, TB=0.0):
    return 1 / (1 + alpha * L + beta * abs(TB))

# (rest of the file remains as v2.4 but with hybrid_weight used in place of the old scale)
THEORETICAL_DELTA = 0.0412
