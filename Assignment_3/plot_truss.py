import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

print("Running the Post_Processor")

# Define paths
out_folder = "outputs/"
nodes_file = os.path.join(out_folder, "undeformed_nodes.csv")
elems_file = os.path.join(out_folder, "undeformed_elements.csv")

# Safety check
if not os.path.exists(nodes_file) or not os.path.exists(elems_file):
    print(f"Error: Could not find CSV files in {out_folder}")
    exit()

# Load Data
nodes = pd.read_csv(nodes_file, index_col="Node_ID")
elems = pd.read_csv(elems_file)

# Initialize Plot
fig, ax = plt.subplots(figsize=(14, 6))
ax.set_aspect("equal")

# Plot Elements
for _, elem in elems.iterrows():
    n_i = elem["Node_i"]
    n_j = elem["Node_j"]

    x_coords = [nodes.loc[n_i, "X"], nodes.loc[n_j, "X"]]
    y_coords = [nodes.loc[n_i, "Y"], nodes.loc[n_j, "Y"]]
    ax.plot(x_coords, y_coords, "k-", linewidth=2.5, zorder=1)

# Plot Boundary Conditions
for node_id, row in nodes.iterrows():
    if row["Fix_X"] == 1.0 and row["Fix_Y"] == 1.0:
        # Draw a right-pointing dark gray triangle
        ax.plot(
            row["X"] - 1.5,
            row["Y"],
            marker=">",
            markersize=14,
            color="dimgray",
            markeredgecolor="black",
            zorder=2,
        )

# Plot Nodes
ax.scatter(
    nodes["X"],
    nodes["Y"],
    color="white",
    edgecolor="black",
    s=80,
    linewidth=1.5,
    zorder=3,
)

# Annotate Node IDs
for node_id, row in nodes.iterrows():
    # If it's on the bottom half, put text below. Otherwise, put it above.
    if row["Y"] < 20.0:
        y_offset = -4.0
        v_align = "top"
    else:
        y_offset = 4.0
        v_align = "bottom"

    ax.text(
        row["X"],
        row["Y"] + y_offset,
        str(node_id),
        color="blue",
        fontsize=10,
        ha="center",
        va=v_align,
        fontweight="bold",
        zorder=5,
    )

# Plot Applied Loads
max_force = max(nodes["Force_X"].abs().max(), nodes["Force_Y"].abs().max())
if max_force > 0:
    for node_id, row in nodes.iterrows():
        fx, fy = row["Force_X"], row["Force_Y"]
        if fx != 0 or fy != 0:
            arrow_length = 25.0
            dx = (fx / max_force) * arrow_length
            dy = (fy / max_force) * arrow_length

            ax.annotate(
                "",
                xy=(row["X"], row["Y"]),
                xytext=(row["X"] - dx, row["Y"] - dy),
                arrowprops=dict(
                    facecolor="red",
                    edgecolor="darkred",
                    width=2.5,
                    headwidth=8,
                    headlength=10,
                ),
                zorder=4,
            )

# Formatting
ax.set_title(
    "Pre-Processing Verification: Undeformed Truss Mesh",
    fontsize=15,
    fontweight="bold",
    pad=15,
)
ax.set_xlabel("X Coordinate (mm)", fontsize=11)
ax.set_ylabel("Y Coordinate (mm)", fontsize=11)
ax.grid(True, linestyle="--", alpha=0.5)

# Buffer the axes so arrows/supports don't get clipped
plt.xlim(nodes["X"].min() - 20, nodes["X"].max() + 20)
plt.ylim(nodes["Y"].min() - 35, nodes["Y"].max() + 35)

plt.tight_layout()
plt.savefig(os.path.join(out_folder, "undeformed_truss.png"), dpi=300)
