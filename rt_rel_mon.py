import pandas as pd
import numpy as np
import uproot
import matplotlib.pyplot as plt
from numpy.polynomial import Chebyshev

# ==========================================================
# USER SETTINGS
# ==========================================================
input_file = "hits_0.csv"
output_file = "perpendicular_hits_rt_converged_monotonic.csv"
t0 = 489.624
tmax = 729.339
r_max = 14.6  # mm
bins = 1000
deg = 20
tol = 0.00025
max_iter = 100
# ==========================================================

print(f"📂 Loading data: {input_file}")
df = pd.read_csv(input_file)

df["t_corr"] = df[" corr_time"] - t0
df = df[(df["t_corr"] > 0) & (df["t_corr"] < tmax)]

hist, edges = np.histogram(
    df["t_corr"], bins=bins, range=(0, df["t_corr"].max()), density=True
)

# Initial integrated curve
cdf = np.cumsum(hist)
cdf = cdf / np.max(cdf) * r_max

t_vals = 0.5 * (edges[:-1] + edges[1:])  # bin centers

# ----------------------------------------------------------
# Initialize r grid (same length as t_vals)
# ----------------------------------------------------------
r_grid = cdf.copy()

print("🔁 Starting Chebyshev iteration (degree=20, monotonic enforced)")

for it in range(max_iter):

    # Fit Chebyshev to current r(t)
    cheb = Chebyshev.fit(t_vals, r_grid, deg)
    coef = cheb.convert().coef

    r_new = Chebyshev(coef)(t_vals)

    # Enforce physical constraints
    r_new = np.clip(r_new, 0, r_max)
    r_new = np.maximum.accumulate(r_new)  # monotonic force
    r_new = r_new / r_new.max() * r_max

    # Interpolate back to hit radii
    df["drift_radius"] = np.interp(df["t_corr"], t_vals, r_new)

    # Check convergence on grid
    diff = np.max(np.abs(r_new - r_grid))
    print(f"Iteration {it+1}: max Δr = {diff:.5f} mm")

    r_grid = r_new.copy()

    if diff < tol:
        print("✅ Converged")
        break
else:
    print("⚠️ Did not converge in max iterations")

df.to_csv(output_file, index=False)
print(f"💾 Saved converged r-t → {output_file}")

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------
plt.figure(figsize=(8,6))

plt.plot(t_vals, cdf, "--", label="Initial integral r(t)")
plt.plot(t_vals, r_grid, linewidth=1, color='red',
         label="Final Chebyshev r(t) (monotonic)")

# Convert to numpy arrays (ensure float64 for ROOT)
t_rt = np.array(t_vals, dtype=np.float64)
r_rt = np.array(r_grid, dtype=np.float64)

with uproot.recreate("rt_relation.root") as f:
    f["rt_tree"] = {
        "time_ns": t_rt,
        "radius_mm": r_rt
    }

print("✅ Saved r–t curve as ROOT TTree → rt_relation.root")
print("   Branches: time_ns, radius_mm")

plt.xlabel("Drift Time - t₀ [ns]")
plt.ylabel("Drift Radius [mm]")
plt.title("MDT r-t Calibration — Chebyshev degree 20 (monotonic)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("rt_iter_mon.png")

