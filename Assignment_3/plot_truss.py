import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

print("--- Truss Visualizer ---")

# Allow passing the job name via terminal, e.g., "python plot_truss.py job1"
if len(sys.argv) > 1:
    job_name = sys.argv[1]
    # Strip extension if accidentally provided
    if job_name.endswith(".inp"):
        job_name = job_name[:-4]
else:
    job_name = "job1"  # Default fallback

outfolder = "outputs"
prefix = os.path.join(outfolder, job_name + "_")
nodes_file = f"{prefix}undeformed_nodes.csv"
elems_file = f"{prefix}undeformed_elements.csv"
def_nodes_file = f"{prefix}deformed_nodes.csv"
history_file = f"{prefix}history.csv"

nodes = pd.read_csv(nodes_file, index_col="Node_ID")
elems = pd.read_csv(elems_file)

fig1, ax1 = plt.subplots(figsize=(14, 6))
ax1.set_aspect("equal")

# Plot Elements
for _, elem in elems.iterrows():
    ax1.plot(
        [nodes.loc[elem["Node_i"], "X"], nodes.loc[elem["Node_j"], "X"]],
        [nodes.loc[elem["Node_i"], "Y"], nodes.loc[elem["Node_j"], "Y"]],
        "k-",
        linewidth=2.5,
        zorder=1,
    )

# Plot Supports, Nodes, and IDs
for node_id, row in nodes.iterrows():
    if row["Fix_X"] == 1.0 and row["Fix_Y"] == 1.0:
        ax1.plot(
            row["X"] - 1.5,
            row["Y"],
            marker=">",
            markersize=14,
            color="dimgray",
            markeredgecolor="black",
            zorder=2,
        )

    y_offset = -4.0 if row["Y"] < 20.0 else 4.0
    v_align = "top" if row["Y"] < 20.0 else "bottom"
    ax1.text(
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

ax1.scatter(
    nodes["X"],
    nodes["Y"],
    color="white",
    edgecolor="black",
    s=80,
    linewidth=1.5,
    zorder=3,
)

# Plot Loads
max_force = max(nodes["Force_X"].abs().max(), nodes["Force_Y"].abs().max())
if max_force > 0:
    for node_id, row in nodes.iterrows():
        if row["Force_X"] != 0 or row["Force_Y"] != 0:
            dx = (row["Force_X"] / max_force) * 25.0
            dy = (row["Force_Y"] / max_force) * 25.0
            ax1.annotate(
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

ax1.set_title(
    "Undeformed Truss Mesh (Input Verification)", fontsize=15, fontweight="bold", pad=15
)
ax1.set_xlabel("X Coordinate (mm)")
ax1.set_ylabel("Y Coordinate (mm)")
ax1.grid(True, linestyle="--", alpha=0.5)
plt.xlim(nodes["X"].min() - 20, nodes["X"].max() + 20)
plt.ylim(nodes["Y"].min() - 35, nodes["Y"].max() + 35)

plt.tight_layout()
plt.savefig(f"{prefix}1_undeformed_truss.png", dpi=300)
print("Saved: 1_undeformed_truss.png")

# POST-PROCESSING: SOLVER RESULTS
def_nodes = pd.read_csv(def_nodes_file, index_col="Node_ID")
history = pd.read_csv(history_file)

# Plot 2: Mesh Overlay
fig2, ax2 = plt.subplots(figsize=(14, 6))
ax2.set_aspect("equal")

# Undeformed (Light Gray)
for _, elem in elems.iterrows():
    ax2.plot(
        [nodes.loc[elem["Node_i"], "X"], nodes.loc[elem["Node_j"], "X"]],
        [nodes.loc[elem["Node_i"], "Y"], nodes.loc[elem["Node_j"], "Y"]],
        color="lightgray",
        linestyle="--",
        linewidth=1.5,
        zorder=1,
    )

# Deformed (Bold Red)
for _, elem in elems.iterrows():
    ax2.plot(
        [
            def_nodes.loc[elem["Node_i"], "X_def"],
            def_nodes.loc[elem["Node_j"], "X_def"],
        ],
        [
            def_nodes.loc[elem["Node_i"], "Y_def"],
            def_nodes.loc[elem["Node_j"], "Y_def"],
        ],
        color="red",
        linestyle="-",
        linewidth=2.0,
        zorder=2,
    )

ax2.scatter(
    def_nodes["X_def"],
    def_nodes["Y_def"],
    color="white",
    edgecolor="red",
    s=50,
    zorder=3,
)
for node_id, row in nodes.iterrows():
    if row["Fix_X"] == 1.0 and row["Fix_Y"] == 1.0:
        ax2.plot(
            row["X"] - 1.5,
            row["Y"],
            marker=">",
            markersize=14,
            color="dimgray",
            markeredgecolor="black",
            zorder=4,
        )

ax2.set_title("Deformed vs. Undeformed Truss", fontsize=15, fontweight="bold", pad=15)
ax2.set_xlabel("X Coordinate (mm)")
ax2.set_ylabel("Y Coordinate (mm)")
ax2.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{prefix}2_mesh_overlay.png", dpi=300)
print("Saved: 2_mesh_overlay.png")

# Plot 3: Load-Displacement
fig3, ax3 = plt.subplots(figsize=(8, 6))


# Determine the primary applied load for the Y-axis
max_fx = history["Applied_Force_X"].abs().max()
max_fy = history["Applied_Force_Y"].abs().max()

if max_fx > max_fy:
    force = history["Applied_Force_X"].abs() / 1000.0
    load_label = "Applied Horizontal Load (kN)"
else:
    force = history["Applied_Force_Y"].abs() / 1000.0
    load_label = "Applied Vertical Load (kN)"

# Extract absolute displacements
disp_x = history["Tip_Disp_X"].abs()
disp_y = history["Tip_Disp_Y"].abs()

ax3.plot(
    disp_x,
    force,
    color="blue",
    marker="o",
    linestyle="-",
    linewidth=2,
    markersize=4,
    label="Horizontal Disp. (X)",
)


ax3.set_title(
    "Load vs. Tip Horizontal Displacement", fontsize=14, fontweight="bold", pad=15
)
ax3.set_xlabel("Tip Displacement Magnitude (mm)")
ax3.set_ylabel(load_label)
ax3.grid(True, linestyle="--", alpha=0.5)
# ax3.legend(loc="best", fontsize=11)
plt.tight_layout()
plt.savefig(f"{prefix}3_load_horizontal_displacement.png", dpi=300)
print("Saved: 3_load_horizontal_displacement.png")

fig4, ax4 = plt.subplots(figsize=(8, 6))

ax4.plot(
    disp_y,
    force,
    color="red",
    marker="s",
    linestyle="--",
    linewidth=2,
    markersize=4,
    label="Vertical Disp. (Y)",
)

ax4.set_title(
    "Load vs. Tip Vertical Displacement", fontsize=14, fontweight="bold", pad=15
)
ax4.set_xlabel("Tip Displacement Magnitude (mm)")
ax4.set_ylabel(load_label)
ax4.grid(True, linestyle="--", alpha=0.5)
# ax4.legend(loc="best", fontsize=11)

plt.tight_layout()
plt.savefig(f"{prefix}4_load_vertical_displacement.png", dpi=300)
print("Saved: 4_load_vertical_displacement.png")

fig5, ax5 = plt.subplots(figsize=(8, 6))

ax5.plot(
    np.sqrt(disp_x**2 + disp_y**2),
    force,
    color="red",
    marker="s",
    linestyle="--",
    linewidth=2,
    markersize=4,
    label="Total Disp.",
)

ax5.set_title("Load vs. Tip Total Displacement", fontsize=14, fontweight="bold", pad=15)
ax5.set_xlabel("Tip Displacement Magnitude (mm)")
ax5.set_ylabel(load_label)
ax5.grid(True, linestyle="--", alpha=0.5)
# ax5.legend(loc="best", fontsize=11)

plt.tight_layout()
plt.savefig(f"{prefix}5_load_total_displacement.png", dpi=300)
print("Saved: 5_load_total_displacement.png")

print("All plotting complete!")
