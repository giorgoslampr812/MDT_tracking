import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Load grouped hits (long format)
df = pd.read_csv("grouped_track_hits_long0.csv")

# Load summary to get (a, b) for lines
summary = pd.read_csv("fitted_tracks_summary0.csv")

# Build lookup: track_id → (a, b)
track_params = {
    int(row.track_id): (row.a, row.b)
    for _, row in summary.iterrows()
}

# ============================================================
# Load *full geometry* from fitted_tracks_hits.csv
# so we can draw all tubes for each TDCID.
# ============================================================
geom_df = pd.read_csv("grouped_track_hits_long0.csv")
geom_df = geom_df[["TDCID", "CHNLID", "geom_x", "geom_y"]].drop_duplicates()

geometry_by_tdc = {
    int(tdc): gdf[["geom_x", "geom_y"]].to_numpy()
    for tdc, gdf in geom_df.groupby("TDCID")
}

# Unique groups
group_ids = sorted(df["group_id"].unique())
print(f"Found {len(group_ids)} groups to plot.")

# Plot settings
x_range = (-125, 125)
y_range = (0, 45)


for gid in group_ids:

    df_g = df[df["group_id"] == gid]
    if df_g.empty:
        continue

    fig, ax = plt.subplots(figsize=(10, 6))

    # ============================================================
    # Draw *all tubes* from the TDCIDs involved in this group
    # ============================================================
    tdcids = sorted(df_g["TDCID"].unique())
    for tdc in tdcids:
        if tdc in geometry_by_tdc:
            xy = geometry_by_tdc[tdc]
            ax.scatter(xy[:, 0], xy[:, 1], s=25, color="lightgray", alpha=0.5)

    # ============================================================
    # Draw the actual hit tube centers
    # ============================================================
    ax.scatter(df_g["geom_x"], df_g["geom_y"], s=60, color="gray", alpha=0.6)

    # ============================================================
    # Draw hit drift circles
    # ============================================================
    for _, h in df_g.iterrows():
        cx, cy = h["geom_x"], h["geom_y"]
        r = h["drift_radius_cm"]
        circ = Circle((cx, cy), r, fill=False,
                      edgecolor="blue", linewidth=1.6, alpha=0.9)
        ax.add_patch(circ)

    # ============================================================
    # Draw track lines (both tracks in the pair)
    # ============================================================
    track_ids = df_g["track_id"].unique()
    for tid in track_ids:
        if tid not in track_params:
            continue

        a, b = track_params[tid]
        xs = np.linspace(x_range[0], x_range[1], 400)
        ys = a * xs + b
        ax.plot(xs, ys, alpha=0.9, linewidth=2)

    # Formatting
    ax.set_title(f"Track Pair {gid} – 6 hits + full TDC geometry")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.set_xlim(*x_range)
    ax.set_ylim(*y_range)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"trackpair_{gid}.png", dpi=200)
    plt.close()
    if gid == 50:
        break
print("Done! Saved trackpair_<group_id>.png for all matched pairs.")

