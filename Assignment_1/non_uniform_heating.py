import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os

# Parameters
l = 1.0  # Length of the plate
w = 1.0  # Width of the plate
nx_elems = 4  # Number of elements in x direction
ny_elems = 4  # Number of elements in y direction
nx_nodes = nx_elems + 1
ny_nodes = ny_elems + 1

# Generate non-uniform node coordinates
x_coords = np.linspace(0, l, nx_elems + 1)
y_coords = np.linspace(0, w, ny_elems + 1)
X, Y = np.meshgrid(x_coords, y_coords)
X_nodes = np.column_stack((X.flatten(), Y.flatten()))

# Generate arrays of element indices in the x and y directions
i_idx, j_idx = np.meshgrid(np.arange(nx_elems), np.arange(ny_elems))

# Connectivity matrix
n1 = j_idx * (nx_nodes) + i_idx
n2 = n1 + 1
n3 = n2 + nx_nodes
n4 = n1 + nx_nodes

conn = np.column_stack((n1.ravel(), n2.ravel(), n3.ravel(), n4.ravel()))


# locate the current directory of this script
current_dir = (
    os.path.dirname(os.path.abspath(__file__))
    if "__file__" in locals()
    else os.getcwd()
)
file_path = os.path.join(current_dir, "corrected_deformed_coords.csv")

print(f"Looking for deformed coordinates at: {file_path}")

try:
    # Lowercase x for current configuration
    x_nodes = np.loadtxt(file_path, delimiter=",")
except FileNotFoundError:
    print("\nERROR: 'corrected_nodes.csv' not found!")
    print(
        "Please ensure your corrected coordinate file is saved in the directory shown above."
    )
    exit()

if len(x_nodes) != len(X_nodes):
    print(
        "\nERROR: Node count mismatch! Undeformed has {} nodes, but Deformed CSV has {}.".format(
            len(X_nodes), len(x_nodes)
        )
    )
    exit()


fig, ax = plt.subplots(figsize=(8, 5))


def compute_kinematics(X_elem, Y_elem, x_elem, y_elem, xi=0.0, eta=0.0):
    xi_n = np.array([-1, 1, 1, -1])
    eta_n = np.array([-1, -1, 1, 1])

    dN_dxi = 0.25 * xi_n * (1 + eta_n * eta)
    dN_deta = 0.25 * eta_n * (1 + xi_n * xi)

    J = np.array(
        [
            [np.dot(dN_dxi, X_elem), np.dot(dN_dxi, Y_elem)],
            [np.dot(dN_deta, X_elem), np.dot(dN_deta, Y_elem)],
        ]
    )

    J_inv = np.linalg.inv(J)

    dN_dnat = np.vstack((dN_dxi, dN_deta))
    dN_dphys = np.dot(J_inv, dN_dnat)

    dN_dX = dN_dphys[0, :]
    dN_dY = dN_dphys[1, :]

    F = np.array(
        [
            [np.dot(dN_dX, x_elem), np.dot(dN_dY, x_elem)],
            [np.dot(dN_dX, y_elem), np.dot(dN_dY, y_elem)],
        ]
    )

    E = 0.5 * (np.dot(F.T, F) - np.eye(2))

    return F, E


F_tensors = []
E_tensors = []

for i, element in enumerate(conn):
    # --- 1. THE MATH (Strictly 4 nodes) ---
    X_elem_math = X_nodes[element, 0]
    Y_elem_math = X_nodes[element, 1]

    x_elem_math = x_nodes[element, 0]
    y_elem_math = x_nodes[element, 1]

    # Compute kinematics using the 4 nodes
    F, E = compute_kinematics(X_elem_math, Y_elem_math, x_elem_math, y_elem_math)
    F_tensors.append(F)
    E_tensors.append(E)

    # --- 2. THE PLOTTING (5 nodes to close the visual loop) ---
    elem_idx = np.append(element, element[0])

    X_elem_plot = X_nodes[elem_idx, 0]
    Y_elem_plot = X_nodes[elem_idx, 1]

    x_elem_def_plot = x_nodes[elem_idx, 0]
    y_elem_def_plot = x_nodes[elem_idx, 1]

    # Add legends only on the first loop iteration
    if i == 0:
        ax.plot(X_elem_plot, Y_elem_plot, "b-", linewidth=1.5, label="undeformed")
        ax.plot(
            x_elem_def_plot, y_elem_def_plot, "r--", linewidth=1.5, label="deformed"
        )
    else:
        ax.plot(X_elem_plot, Y_elem_plot, "b-", linewidth=1.5)
        ax.plot(x_elem_def_plot, y_elem_def_plot, "r--", linewidth=1.5)
ax.plot(X_nodes[:, 0], X_nodes[:, 1], "bo", markerfacecolor="none", markersize=6)

ax.plot(x_nodes[:, 0], x_nodes[:, 1], "ro", markerfacecolor="none", markersize=6)

# Formatting to match Figure 1 of the assignment
ax.set_xlabel("$X$ (m)", fontsize=14)
ax.set_ylabel("$Y$ (m)", fontsize=14)
ax.set_xlim(0, 1.05)
ax.set_ylim(-0.25, 1.25)
ax.set_xticks(np.arange(0, 1.2, 0.2))
ax.set_yticks(np.arange(-0.2, 1.4, 0.2))
ax.legend(loc="upper left", framealpha=1, edgecolor="black")

plt.tight_layout()
plt.show()

# Extract strain components from the tensors
E11 = [E[0, 0] for E in E_tensors]
E22 = [E[1, 1] for E in E_tensors]
E12 = [E[0, 1] for E in E_tensors]

# Extract element vertices (using undeformed reference configuration)
verts = [X_nodes[element] for element in conn]

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
strains = [E11, E22, E12]
titles = ["$E_{11}$", "$E_{22}$", "$E_{12}$"]

for ax, strain_data, title in zip(axes, strains, titles):
    poly = PolyCollection(verts, cmap="jet", edgecolor="black", linewidth=0.5)
    poly.set_array(np.array(strain_data))

    ax.add_collection(poly)
    fig.colorbar(poly, ax=ax, fraction=0.046, pad=0.04)

    ax.set_title(title, fontsize=14)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("$X$ (m)", fontsize=12)
    ax.set_ylabel("$Y$ (m)", fontsize=12)

plt.tight_layout()
plt.show()
