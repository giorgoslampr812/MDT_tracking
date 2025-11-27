import pandas as pd
import numpy as np
from numba import njit

# ==========================================
# CONFIG
# ==========================================
INPUT_FILE  = "fitted_tracks_summary10.csv"
OUTPUT_FILE = "matched_track_pairs10.csv"

a_tol = 0.03
b_tol = 0.5
CHUNK_SIZE = 1600
# ==========================================


# ==============================================================
# Numba JIT function for pairing tracks
# ==============================================================

@njit
def match_even_odd(a_e, b_e, a_o, b_o, odd_available,
                   a_tol, b_tol, chunk_size):
    """
    Return an array of match indices (odd_j for each even i), -1 if no match.
    """
    Ne = len(a_e)
    No = len(a_o)
    match_idx = np.full(Ne, -1)

    for i in range(Ne):  # loop over even tracks

        best_j = -1
        best_score = 1e18  # large number

        for start in range(0, No, chunk_size):
            end = min(start + chunk_size, No)

            # Check odd availability in this chunk
            for local in range(start, end):

                if not odd_available[local]:
                    continue

                da = abs(a_e[i] - a_o[local])
                if da >= a_tol:
                    continue

                db = abs(b_e[i] - b_o[local])
                if db >= b_tol:
                    continue

                score = da + db

                if score < best_score:
                    best_score = score
                    best_j = local

        # Store match, mark odd as used
        if best_j >= 0:
            match_idx[i] = best_j
            odd_available[best_j] = False

    return match_idx


# ==============================================================
# Main Python matching logic
# ==============================================================

df = pd.read_csv(INPUT_FILE)
df = df[df["n_hits"] == 3].copy()

matched_rows = []
group_id = 0

if "trigger_min" not in df.columns:
    raise ValueError("Column 'trigger_min' is required.")

# ---------------------------------------------------------
# Group by trigger window
# ---------------------------------------------------------
for trig, group in df.groupby("n_hits"):

    ev = group[group["TDCID"] % 2 == 0].reset_index(drop=True)
    od = group[group["TDCID"] % 2 == 1].reset_index(drop=True)

    if ev.empty or od.empty:
        continue

    # Convert to NumPy for numba
    a_e = ev["a"].to_numpy()
    b_e = ev["b"].to_numpy()
    a_o = od["a"].to_numpy()
    b_o = od["b"].to_numpy()

    odd_available = np.ones(len(a_o), dtype=np.bool_)

    # ---- Numba accelerated matching ----
    match_idx = match_even_odd(
        a_e, b_e, a_o, b_o, odd_available,
        a_tol, b_tol, CHUNK_SIZE
    )

    # ---- Convert matches back to Python dict entries ----
    for i, j in enumerate(match_idx):
        if j < 0:
            continue

        matched_rows.append({
            "group_id": group_id,
            "trigger_min": trig,
            "even_track_id": int(ev.at[i, "track_id"]),
            "odd_track_id":  int(od.at[j, "track_id"]),

            "even_TDCID": int(ev.at[i, "TDCID"]),
            "odd_TDCID":  int(od.at[j, "TDCID"]),

            "a_even": float(ev.at[i, "a"]),
            "a_odd":  float(od.at[j, "a"]),
            "a_diff": abs(ev.at[i, "a"] - od.at[j, "a"]),

            "b_even": float(ev.at[i, "b"]),
            "b_odd":  float(od.at[j, "b"]),
            "b_diff": abs(ev.at[i, "b"] - od.at[j, "b"]),

            "angle_even": float(ev.at[i, "angle_deg"]),
            "angle_odd":  float(od.at[j, "angle_deg"]),

            "chi2_even": float(ev.at[i, "chi2ndf"]),
            "chi2_odd":  float(od.at[j, "chi2ndf"]),
        })

        group_id += 1


# -------------------------------------------------
# Save output
# -------------------------------------------------
df_out = pd.DataFrame(matched_rows)
df_out.to_csv(OUTPUT_FILE, index=False)

print(f"Matched {len(df_out)} odd–even track pairs (Numba accelerated).")

