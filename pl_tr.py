import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------
# Load tracked hits CSV
# ----------------------------
df = pd.read_csv("tracked_output_top_best_two_iterations.csv")

# ----------------------------
# MDT geometry in cm -> mm
# ----------------------------
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
    tid: {ch: ((x+dx)*10.0, (y+dy)*10.0) for ch,(x,y) in geometry_layer.items()}
    for tid,(dx,dy) in enumerate(tdc_offsets)
}

# ----------------------------
# Plot per track_id
# ----------------------------
for track_id in df["track_id"].unique():
    dftrack = df[df["track_id"]==track_id]

    plt.figure(figsize=(12,6))
    plt.title(f"Track ID {track_id}")
    plt.xlabel("X [mm]")
    plt.ylabel("Y [mm]")

    # Plot tubes as circles
    for idx,row in dftrack.iterrows():
        tdc = int(row["TDCID"])
        ch = int(row["CHNLID"])
        x0,y0 = geometry_tdcs[tdc][ch]
        r = row["drift_radius"]
        circle = plt.Circle((x0,y0), r, color='blue', alpha=0.3, fill=False)
        plt.gca().add_patch(circle)
        plt.plot(x0, y0, 'ko')  # tube center

        # optional: show residual as line from center to track
        # res_line_end = (x0 + r, y0)
        # plt.plot([x0,res_line_end[0]],[y0,res_line_end[1]],'r--')

    # Draw fitted line: a*x + b*y + c = 0
    a = dftrack.iloc[0]["a"]
    b = dftrack.iloc[0]["b"]
    c = dftrack.iloc[0]["c"]

    x_vals = np.linspace(dftrack["x"].min()-50, dftrack["x"].max()+50, 500)
    if abs(b) > 1e-6:
        y_vals = (-a*x_vals - c)/b
    else:
        y_vals = np.full_like(x_vals, -c/a)
    plt.plot(x_vals, y_vals, 'r-', label='Fitted line')

    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    if track_id < 30:
     plt.savefig(f"{track_id} track.png")
    plt.close()
