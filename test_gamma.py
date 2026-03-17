"""
Compare gamma.py (new local-window) vs gamma_original.py (full NxM)
on a synthetic profile pair. Prints max difference and pass rates.
"""
import numpy as np
import importlib.util, sys

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

base = r"c:\Users\nknutson\GitRepositories\BeamDataComparisonTools"
new_mod  = load(base + r"\gamma.py",          "gamma_new")
orig_mod = load(base + r"\gamma_original.py", "gamma_orig")

# Synthetic profiles: Gaussian with a small shift
x1 = np.arange(-15, 15, 0.25)
x2 = np.arange(-15, 15, 0.25)
d1 = np.exp(-x1**2 / (2 * 3**2))
d2 = np.exp(-(x2 - 0.1)**2 / (2 * 3**2))   # 1 mm shift

dd  = 0.005   # 0.5%
dta = 0.05    # 0.05 cm (0.5 mm)

gx_new,  gv_new  = new_mod.gamma( x1, d1, x2, d2, dd, dta, 1, 0.01, 0.01)
gx_orig, gv_orig = orig_mod.gamma(x1, d1, x2, d2, dd, dta, 1, 0.01, 0.01)

# Trim to same length for comparison
n = min(len(gv_new), len(gv_orig))
diff = np.abs(gv_new[:n] - gv_orig[:n])

print(f"Points compared : {n}")
print(f"Max |difference|: {diff.max():.2e}")
print(f"Mean difference : {diff.mean():.2e}")
print(f"Pass rate (new) : {(gv_new[:n] <= 1).mean()*100:.2f}%")
print(f"Pass rate (orig): {(gv_orig[:n] <= 1).mean()*100:.2f}%")

if diff.max() < 1e-6:
    print("\nPASS — results are numerically identical.")
else:
    print("\nWARNING — differences exceed 1e-6, investigate.")
