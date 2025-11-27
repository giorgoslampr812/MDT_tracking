import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Try to import numba for JIT speed-up
try:
    import numba
    from numba import njit
    NUMBA_AVAILABLE = True
    print("Numba found: JIT enabled for chi2.")
except Exception:
    NUMBA_AVAILABLE = False
    print("Numba not available: running without JIT.")
# === Load hits file ===
df = pd.read_csv("hits_10_with_radius.csv")
df["drift_radius_cm"] = df["drift_radius"] / 10.0
track_hits_file = "fitted_tracks_hits10.csv"
track_sum_file = "fitted_tracks_summary10.csv"
# === MDT geometry in cm ===
geometry_layer = {
    0:(1.5,1.5),1:(4.5,1.5),2:(7.5,1.5),3:(10.5,1.5),
    4:(13.5,1.5),5:(16.5,1.5),6:(19.5,1.5),7:(22.5,1.5),
    8:(3.0,4.1),9:(6.0,4.1),10:(9.0,4.1),11:(12.0,4.1),
    12:(15.0,4.1),13:(18.0,4.1),14:(21.0,4.1),15:(24.0,4.1),
    16:(1.5,6.7),17:(4.5,6.7),18:(7.5,6.7),19:(10.5,6.7),
    20:(13.5,6.7),21:(16.5,6.7),22:(19.5,6.7),23:(22.5,6.7)
}

tdc_offsets = [
    (-96.0,0), (-96.0,34.7), (-72.0,0), (-72.0,34.7),
    (-48.0,0), (-48.0,34.7), (-24.0,0), (-24.0,34.7),
    (0.0,0.0), (0.0,34.7), (24.0,0.0), (24.0,34.7),
    (48.0,0.0), (48.0,34.7), (72.0,0.0), (72.0,34.7),
    (96.0,0.0), (96.0,34.7)
]

geometry_tdcs = {
    tdcid:{ch:(x+dx,y+dy) for ch,(x,y) in geometry_layer.items()}
    for tdcid,(dx,dy) in enumerate(tdc_offsets)
}

# === Multilayer definitions ===
top_ch_list = list(range(16, 24))
mid_ch_list = list(range(8, 16))
bottom_ch_list = list(range(0, 8))

# === Parameters ===
window_size = 2000
tolerance = 0.1         # cm tolerance for residual
angle_cut = (45, 90)    # degrees
chi2_cut = 3.0
min_hits = 3

# The three specific configurations to try per top hit
configs = [
    (-8, -16),
    (-8, -15),
    (-9, -17),
    (-9, -16)
]

# === Prepare trigger windows ===
trigger_min = int(df["triggerledge"].min())
trigger_max = int(df["triggerledge"].max())
n_windows = int(np.ceil((trigger_max - trigger_min) / window_size))
print(f"Processing {n_windows} trigger windows...")

# Output containers
track_id = 0
tracks_all = []
tracks_summary = []

# --- numba-jitted chi2 (if available) ---
if NUMBA_AVAILABLE:
    @njit
    def chi2_numba(a, b, x_arr, y_arr, r_arr, tol):
        s = 0.0
        for i in range(x_arr.shape[0]):
            val = abs(y_arr[i] - (a * x_arr[i] + b)) - r_arr[i]
            t = val / tol
            s += t * t
        return s

    # wrapper to match scipy minimize signature: takes param array
    def chi2_wrapper(params, x_arr, y_arr, r_arr, tol):
        return float(chi2_numba(params[0], params[1], x_arr, y_arr, r_arr, tol))
else:
    # pure python chi2 wrapper (fast enough for modest datasets)
    def chi2_wrapper(params, x_arr, y_arr, r_arr, tol):
        a, b = params
        val = np.abs(y_arr - (a * x_arr + b)) - r_arr
        return float(np.sum((val / tol) ** 2))

# -------------------------
# Main loop over windows
# -------------------------
for w_idx in range(n_windows):
    wmin = trigger_min + w_idx * window_size
    wmax = wmin + window_size
    df_win = df[(df["triggerledge"] >= wmin) & (df["triggerledge"] < wmax)]
    if df_win.empty:
        continue

    print(f"Window {w_idx+1}/{n_windows}: [{wmin},{wmax}) hits={len(df_win)}")

    # Group hits by TDCID for quick access
    hits_by_tdc = {tdc: grp for tdc, grp in df_win.groupby("TDCID")}

    # For each TDC
    for tdcid, geom in geometry_tdcs.items():
        if tdcid not in hits_by_tdc:
            continue
        hits_tdc = hits_by_tdc[tdcid]

        # Build quick lookup: CHNLID -> DataFrame (rows)
        hits_by_chn = {ch: grp for ch, grp in hits_tdc.groupby("CHNLID")}

        # Quick existence checks
        if not any(c in hits_by_chn for c in top_ch_list):
            continue
        if not any(c in hits_by_chn for c in mid_ch_list):
            continue
        if not any(c in hits_by_chn for c in bottom_ch_list):
            continue

        # Iterate top hits only
        top_present_chs = [ch for ch in top_ch_list if ch in hits_by_chn]
        for ch_top in top_present_chs:
            top_rows = hits_by_chn[ch_top]
            # There may be multiple top hits — iterate them
            for _, top_row in top_rows.iterrows():

                # Try the three chosen configurations
                for off_mid, off_bot in configs:
                    mid_ch_cfg = ch_top + off_mid
                    bot_ch_cfg = ch_top + off_bot

                    # Bound checks
                    if mid_ch_cfg < 0 or mid_ch_cfg > 23:
                        continue
                    if bot_ch_cfg < 0 or bot_ch_cfg > 23:
                        continue
                    if mid_ch_cfg not in hits_by_chn or bot_ch_cfg not in hits_by_chn:
                        continue

                    mid_row = hits_by_chn[mid_ch_cfg].iloc[0]
                    bot_row = hits_by_chn[bot_ch_cfg].iloc[0]

                    # Build triple
                    hits_combo = [top_row, mid_row, bot_row]

                    x_coords = np.empty(3, dtype=np.float64)
                    y_coords = np.empty(3, dtype=np.float64)
                    r_vals = np.empty(3, dtype=np.float64)
                    idxs = []

                    for i_h, h in enumerate(hits_combo):
                        ch = int(h["CHNLID"])
                        xg, yg = geom[ch]
                        x_coords[i_h] = xg
                        y_coords[i_h] = yg
                        r_vals[i_h] = float(h["drift_radius_cm"])
                        idxs.append(h.name)

                    # minimize using wrapper that calls numba-compiled inner loop (if available)
                    res = minimize(lambda p: chi2_wrapper(p, x_coords, y_coords, r_vals, tolerance),
                                   [0.0, np.mean(y_coords)], method="Nelder-Mead")
                    if not res.success:
                        continue

                    a_fit, b_fit = res.x
                    angle = np.degrees(np.arctan(a_fit))
                    chi2ndf = res.fun / (len(x_coords) - 2)

                    # Cuts
                    if abs(angle) < angle_cut[0] or abs(angle) > angle_cut[1]:
                        continue
                    if chi2ndf > chi2_cut:
                        continue

                    residuals = np.abs(y_coords - (a_fit * x_coords + b_fit)) - r_vals
                    matched_hits = np.sum(np.abs(residuals) < tolerance)
                    if matched_hits < min_hits:
                        continue

                    # Save exactly these three hits (use df_win rows)
                    for i_h, idx_h in enumerate(idxs):
                        h = df_win.loc[idx_h]
                        ch = int(h["CHNLID"])
                        xg, yg = geom[ch]
                        row = h.to_dict()
                        row.update({
                            "track_id": track_id,
                            "track_a": float(a_fit),
                            "track_b": float(b_fit),
                            "residual": float(residuals[i_h]),
                            "geom_x": float(xg),
                            "geom_y": float(yg)
                        })
                        tracks_all.append(row)

                    tracks_summary.append({
                        "track_id": track_id,
                        "TDCID": tdcid,
                        "trigger_min": wmin,
                        "trigger_max": wmax,
                        "a": float(a_fit),
                        "b": float(b_fit),
                        "angle_deg": float(angle),
                        "chi2ndf": float(chi2ndf),
                        "n_hits": int(matched_hits)
                    })

                    track_id += 1

# === Save outputs ===
df_tracks = pd.DataFrame(tracks_all)
df_summary = pd.DataFrame(tracks_summary)

# Keep only unique 3 hits per track (safety)
if not df_tracks.empty:
    df_tracks = df_tracks.drop_duplicates(subset=["track_id", "TDCID", "CHNLID"])

df_tracks.to_csv(track_hits_file, index=False)
df_summary.to_csv(track_sum_file, index=False)

print(f"\n✅ Done — found {len(df_summary)} tracks across windows")
print("Saved fitted_tracks_hits.csv and fitted_tracks_summary.csv")
"""
# === Plot geometry + drift circles + tracks ===
plt.figure(figsize=(14,8))

# draw tube centers (full geometry)
for tdcid, geom in geometry_tdcs.items():
    for ch, (x, y) in geom.items():
        plt.plot(x, y, 'ko', markersize=3, alpha=0.25)

# draw drift circles for hits used in tracks
if not df_tracks.empty:
    for _, hit in df_tracks.iterrows():
        x0 = hit["geom_x"]
        y0 = hit["geom_y"]
        r  = hit["drift_radius_cm"]
        circ = plt.Circle((x0, y0), r, fill=False, linewidth=0.9, alpha=0.6, color='C1')
        plt.gca().add_patch(circ)

# draw fitted track lines
if not df_summary.empty:
    for _, tr in df_summary.iterrows():
        a = tr["a"]; b = tr["b"]
        xs_line = np.linspace(-125, 125, 500)
        ys_line = a * xs_line + b
        plt.plot(xs_line, ys_line, linewidth=1.4, color='C0', alpha=0.9)

plt.xlim(-125, 125)
plt.ylim(0, 45)
plt.xlabel("x [cm]")
plt.ylabel("y [cm]")
plt.title(f"Tracks & drift circles (triggerledge {trigger_min}..{trigger_max})")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("tracks_plot.png", dpi=200)
plt.close()
"""
