import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================
# USER INPUT
# ============================
CSV_FILE = "tracked_out_2.csv"
TRACK_ID = 0 #change to plot the track of choice
# ============================

# ----------------------------
# MDT GEOMETRY (cm → mm)
# ----------------------------
geometry_layer_cm = {
    0:(1.5,1.5), 1:(4.5,1.5), 2:(7.5,1.5), 3:(10.5,1.5),
    4:(13.5,1.5),5:(16.5,1.5),6:(19.5,1.5),7:(22.5,1.5),
    8:(3.0,4.1), 9:(6.0,4.1),10:(9.0,4.1),11:(12.0,4.1),
    12:(15.0,4.1),13:(18.0,4.1),14:(21.0,4.1),15:(24.0,4.1),
    16:(1.5,6.7),17:(4.5,6.7),18:(7.5,6.7),19:(10.5,6.7),
    20:(13.5,6.7),21:(16.5,6.7),22:(19.5,6.7),23:(22.5,6.7)
}

# Convert to mm
geometry_layer = {k:(v[0]*10.0, v[1]*10.0) for k,v in geometry_layer_cm.items()}

tdc_offsets_cm = [
    (-96,0),(-96,34.7),(-72,0),(-72,34.7),
    (-48,0),(-48,34.7),(-24,0),(-24,34.7),
    (0,0),(0,34.7),(24,0),(24,34.7),
    (48,0),(48,34.7),(72,0),(72,34.7),
    (96,0),(96,34.7)
]

tdc_offsets = [(dx*10.0, dy*10.0) for dx,dy in tdc_offsets_cm]

TUBE_RADIUS_MM = 15.0   # 1.45 cm

# ----------------------------
# Load CSV
# ----------------------------
df = pd.read_csv(CSV_FILE)
df_track = df[df["track_id"] == TRACK_ID]

if df_track.empty:
    raise RuntimeError(f"No entries found for track_id {TRACK_ID}")

# ----------------------------
# Extract hit info (mm)
# ----------------------------
x_mm = df_track["x"].values
y_mm = df_track["y"].values
r_mm = df_track["drift_radius"].values

a = df_track.iloc[0]["a"]
b = df_track.iloc[0]["b"]
c = df_track.iloc[0]["c"]

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(7,10))

# ----------------------------
# Draw MDT tubes (grey)
# ----------------------------
for dx, dy in tdc_offsets:
    for lx, ly in geometry_layer.values():
        tube = plt.Circle(
            (lx + dx, ly + dy),
            TUBE_RADIUS_MM,
            edgecolor="lightgrey",
            facecolor="none",
            linewidth=0.8,
            zorder=1
        )
        plt.gca().add_patch(tube)

# ----------------------------
# Draw drift circles + centers
# ----------------------------
for x, y, r in zip(x_mm, y_mm, r_mm):
    plt.gca().add_patch(
        plt.Circle((x, y), r,
                   edgecolor="blue",
                   facecolor="none",
                   linewidth=2,
                   zorder=3)
    )
    plt.plot(x, y, "ko", markersize=4, zorder=4)

# ----------------------------
# Fixed axes (mm)
# ----------------------------
XMIN, XMAX = -1000.0, 1250.0 # change if you need to zoom in on the track
YMIN, YMAX = 0.0, 500.0
plt.xlim(XMIN, XMAX)
plt.ylim(YMIN, YMAX)

# ----------------------------
# Draw full track: ax + by + c = 0
# ----------------------------
points = []

if abs(b) > 1e-9:
    y1 = (-a * XMIN - c) / b
    y2 = (-a * XMAX - c) / b
    if YMIN <= y1 <= YMAX: points.append((XMIN, y1))
    if YMIN <= y2 <= YMAX: points.append((XMAX, y2))

if abs(a) > 1e-9:
    x1 = (-b * YMIN - c) / a
    x2 = (-b * YMAX - c) / a
    if XMIN <= x1 <= XMAX: points.append((x1, YMIN))
    if XMIN <= x2 <= XMAX: points.append((x2, YMAX))

if len(points) >= 2:
    (x1, y1), (x2, y2) = points[:2]
    plt.plot([x1, x2], [y1, y2], "r-", linewidth=3, zorder=5)
else:
    print("WARNING: track does not intersect plotting window")

# ----------------------------
# Cosmetics
# ----------------------------
plt.xlabel("X [mm]")
plt.ylabel("Y [mm]")
plt.title(f"Track ID {TRACK_ID}")
plt.grid(True)
plt.savefig(f"{TRACK_ID}.png")


