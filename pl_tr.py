import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================
# USER INPUT
# ============================
CSV_FILE = "tracked_out_0.csv"
TRACK_ID = 6   # <<< change ONLY this
# ============================

# ----------------------------
# Load CSV and select track
# ----------------------------
df = pd.read_csv(CSV_FILE)
df_track = df[df["track_id"] == TRACK_ID]

if df_track.empty:
    raise RuntimeError(f"No entries found for track_id = {TRACK_ID}")

# ----------------------------
# Extract hit info
# ----------------------------
x_cm = df_track["x"].values          # cm
y_cm = df_track["y"].values          # cm
r_cm = df_track["drift_radius"].values / 10.0  # mm → cm

# Line parameters: a x + b y + c = 0
a = df_track.iloc[0]["a"]
b = df_track.iloc[0]["b"]
c = df_track.iloc[0]["c"]

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(6,10))

# Drift circles (blue) and centers (black)
for x, y, r in zip(x_cm, y_cm, r_cm):
    circle = plt.Circle(
        (x, y), r,
        edgecolor="blue",
        facecolor="none",
        linewidth=2
    )
    plt.gca().add_patch(circle)
    plt.plot(x, y, "ko", markersize=4)

# ----------------------------
# Draw fitted track (robust)
# ----------------------------
if abs(b) > abs(a):
    # Solve y(x)
    x_vals = np.linspace(min(x_cm) - 20, max(x_cm) + 20, 500)
    y_vals = (-a * x_vals - c) / b
else:
    # Solve x(y)  → stable for vertical tracks
    y_vals = np.linspace(min(y_cm) - 20, max(y_cm) + 20, 500)
    x_vals = (-b * y_vals - c) / a

plt.plot(x_vals, y_vals, "r-", linewidth=3)

# ----------------------------
# Cosmetics
# ----------------------------
plt.xlabel("X [cm]")
plt.ylabel("Y [cm]")
plt.title(f"Track ID {TRACK_ID} with Drift Circles")
plt.axis("equal")
plt.grid(True)

plt.show()

