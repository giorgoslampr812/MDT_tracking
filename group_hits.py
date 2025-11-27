import pandas as pd

# Load files
pairs = pd.read_csv("matched_track_pairs10.csv")
hits  = pd.read_csv("fitted_tracks_hits10.csv")

# Ensure track_id is integer
hits["track_id"] = hits["track_id"].astype(int)
pairs["even_track_id"] = pairs["even_track_id"].astype(int)
pairs["odd_track_id"] = pairs["odd_track_id"].astype(int)

grouped_rows = []
grouped_rows_long = []

for _, row in pairs.iterrows():
    gid  = row["group_id"]
    tid_e = row["even_track_id"]
    tid_o = row["odd_track_id"]

    # Select exactly the 3 hits of each track
    hits_e = hits[hits["track_id"] == tid_e]
    hits_o = hits[hits["track_id"] == tid_o]

    if len(hits_e) != 3 or len(hits_o) != 3:
        continue  # safety

    # --- LONG FORMAT (6 rows) ---
    tmp_e = hits_e.copy()
    tmp_e["group_id"] = gid
    tmp_e["track_type"] = "even"

    tmp_o = hits_o.copy()
    tmp_o["group_id"] = gid
    tmp_o["track_type"] = "odd"

    grouped_rows_long.append(tmp_e)
    grouped_rows_long.append(tmp_o)

    # --- WIDE FORMAT (single merged row with 6 hits) ---
    merged = {"group_id": gid,
              "even_track_id": tid_e,
              "odd_track_id": tid_o}



# Save outputs

df_grouped_long = pd.concat(grouped_rows_long, ignore_index=True)
df_grouped_long.to_csv("grouped_track_hits_long10.csv", index=False)

print(f"Saved:\n  grouped_track_hits.csv (wide format)\n  grouped_track_hits_long.csv (6 rows per pair)")

