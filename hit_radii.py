import uproot
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# === Load hits CSV ===
hits_file = "hits_10.csv"
df = pd.read_csv(hits_file)

# === Load RT relation from ROOT file ===
rt_file = "rt_relation.root"
with uproot.open(rt_file) as f:
    tree = f["rt_tree"]
    time = tree["time_ns"].array(library="np")      # drift time in ns
    radius = tree["radius_mm"].array(library="np")  # drift radius in mm

# --- Interpolation function ---
# Ensure monotonic increasing
sorted_idx = np.argsort(time)
time_sorted = time[sorted_idx]
radius_sorted = radius[sorted_idx]

rt_interp = interp1d(time_sorted, radius_sorted, kind='linear',
                     bounds_error=False, fill_value=(radius_sorted[0], radius_sorted[-1]))

# --- Apply RT relation to hits ---
df["drift_radius"] = rt_interp(df[" drift_time"]-489.624)

# --- Save updated CSV ---
output_file = "hits_10_with_radius.csv"
df.to_csv(output_file, index=False)
print(f"✅ Drift radius added and saved to {output_file}")

