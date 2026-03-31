import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

print("--- Truss Visualizer ---")

out_folder = "outputs/"
nodes_file = os.path.join(out_folder, "undeformed_nodes.csv")
elems_file = os.path.join(out_folder, "undeformed_elements.csv")

# PRE-PROCESSING: UNDEFORMED MESH (ALWAYS RUNS)
if not os.path.exists(nodes_file) or not os.path.exists(elems_file):
    print(f"FATAL ERROR: Could not find base mesh files in {out_folder}")
    exit()

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
plt.savefig(os.path.join(out_folder, "1_undeformed_truss.png"), dpi=300)
print("Saved: 1_undeformed_truss.png")

# POST-PROCESSING: SOLVER RESULTS
def_nodes_file = os.path.join(out_folder, "deformed_nodes.csv")
history_file = os.path.join(out_folder, "history.csv")

if not os.path.exists(def_nodes_file) or not os.path.exists(history_file):
    print(
        "Notice: Solver results not found. Stopping after pre-processing plot. (Did the solver crash?)"
    )
    exit()

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
plt.savefig(os.path.join(out_folder, "2_mesh_overlay.png"), dpi=300)
print("Saved: 2_mesh_overlay.png")

# Plot 3: Load-Displacement
fig3, ax3 = plt.subplots(figsize=(8, 6))

# Determine the primary direction of the applied load
max_fx = history["Applied_Force_X"].abs().max()
max_fy = history["Applied_Force_Y"].abs().max()

if max_fx > max_fy:
    # Load is primarily horizontal
    disp = history["Tip_Disp_X"].abs()
    force = history["Applied_Force_X"].abs() / 1000.0
    x_label = "Tip Horizontal Displacement (mm)"
else:
    # Load is primarily vertical
    disp = history["Tip_Disp_Y"].abs()
    force = history["Applied_Force_Y"].abs() / 1000.0
    x_label = "Tip Vertical Displacement (mm)"

ax3.plot(
    disp, force, color="blue", marker="o", linestyle="-", linewidth=2, markersize=4
)

ax3.set_title("Load vs. Reaction Displacement", fontsize=14, fontweight="bold", pad=15)
ax3.set_xlabel(x_label)
ax3.set_ylabel("Applied Load (kN)")
ax3.grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig(os.path.join(out_folder, "3_load_displacement.png"), dpi=300)
print("Saved: 3_load_displacement.png")

print("All plotting complete!")
