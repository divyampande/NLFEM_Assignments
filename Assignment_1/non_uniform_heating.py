import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os
from scipy.interpolate import griddata

# Choose "element" for cell-centered constant values
# Choose "nodal" for smoothed nodal averaged values
# Choose "gauss_surface" for plots suggested in the assignment
COMPUTE_MODE = "gauss_surface"

# Setup the geometry and mesh
l = 1.0
w = 1.0
nx_elems = 4
ny_elems = 4
nx_nodes = nx_elems + 1
ny_nodes = ny_elems + 1

# Generate undeformed node coordinates
x_coords = np.linspace(0, l, nx_nodes)
y_coords = np.linspace(0, w, ny_nodes)
X, Y = np.meshgrid(x_coords, y_coords)
X_nodes = np.column_stack((X.flatten(), Y.flatten()))

# Connectivity matrix for 4-node quadrilaterals
i_idx, j_idx = np.meshgrid(np.arange(nx_elems), np.arange(ny_elems))
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


import numpy as np


def compute_kinematics_vectorized(X_elems, x_elems, xi=0.0, eta=0.0):
    xi_n = np.array([-1, 1, 1, -1])
    eta_n = np.array([-1, -1, 1, 1])

    dN_dxi = 0.25 * xi_n * (1 + eta_n * eta)
    dN_deta = 0.25 * eta_n * (1 + xi_n * xi)

    # Shape: (2, 4)
    dN_dnat = np.vstack((dN_dxi, dN_deta))

    # Jacobian. Shape (N, 2, 2)
    J = np.einsum("ij, Njk -> Nik", dN_dnat, X_elems)
    J_inv = np.linalg.inv(J)

    # Physical derivatives. Shape: (N, 2, 4)
    dN_dphys = np.einsum("Nij, jk -> Nik", J_inv, dN_dnat)

    # Deformation Gradient. Shape: (N, 2, 2)
    F = np.einsum("Nik, Nkj -> Nij", dN_dphys, x_elems)

    # Transposes for batched arrays
    F_T = np.transpose(F, axes=(0, 2, 1))
    I = np.eye(2)

    # Strain Tensors for ALL elements
    E = 0.5 * (np.matmul(F_T, F) - I)

    F_inv = np.linalg.inv(F)
    F_inv_T = np.transpose(F_inv, axes=(0, 2, 1))
    e = 0.5 * (I - np.matmul(F_inv_T, F_inv))

    eps = 0.5 * (F + F_T) - I

    return F, E, e, eps


# Computation
if COMPUTE_MODE == "element":
    # Shape: (N_elems, 4, 2)
    X_elems = X_nodes[conn]
    x_elems = x_nodes[conn]

    # Evaluate at element center for the entire mesh
    F_results, E_results, e_results, eps_results = compute_kinematics_vectorized(
        X_elems, x_elems, xi=0.0, eta=0.0
    )

elif COMPUTE_MODE == "nodal":
    num_nodes = len(X_nodes)

    F_sum = np.zeros((num_nodes, 2, 2))
    E_sum = np.zeros((num_nodes, 2, 2))
    e_sum = np.zeros((num_nodes, 2, 2))
    eps_sum = np.zeros((num_nodes, 2, 2))
    node_counts = np.zeros(num_nodes)

    node_nat_coords = [(-1, -1), (1, -1), (1, 1), (-1, 1)]

    X_elems = X_nodes[conn]
    x_elems = x_nodes[conn]

    for j, (xi_val, eta_val) in enumerate(node_nat_coords):
        F, E, e, eps = compute_kinematics_vectorized(
            X_elems, x_elems, xi=xi_val, eta=eta_val
        )

        # Get the global node IDs for this specific corner (j) across all elements
        global_nodes_j = conn[:, j]

        # Assembly (Scatter and Add)
        np.add.at(F_sum, global_nodes_j, F)
        np.add.at(E_sum, global_nodes_j, E)
        np.add.at(e_sum, global_nodes_j, e)
        np.add.at(eps_sum, global_nodes_j, eps)
        np.add.at(node_counts, global_nodes_j, 1)

    # Average the tensors
    counts_reshaped = node_counts[:, None, None]
    F_results = F_sum / counts_reshaped
    E_results = E_sum / counts_reshaped
    e_results = e_sum / counts_reshaped
    eps_results = eps_sum / counts_reshaped


elif COMPUTE_MODE == "gauss_surface":
    gauss_val = 1.0 / np.sqrt(3.0)
    gauss_points = [
        (-gauss_val, -gauss_val),
        (gauss_val, -gauss_val),
        (gauss_val, gauss_val),
        (-gauss_val, gauss_val),
    ]

    num_elems = len(conn)
    num_gp = len(gauss_points)

    gp_coords = np.zeros((num_elems, num_gp, 2))
    F_gp_array = np.zeros((num_elems, num_gp, 2, 2))
    E_gp_array = np.zeros((num_elems, num_gp, 2, 2))
    e_gp_array = np.zeros((num_elems, num_gp, 2, 2))
    eps_gp_array = np.zeros((num_elems, num_gp, 2, 2))

    X_elems = X_nodes[conn]
    x_elems = x_nodes[conn]

    # Loop over the 4 Gauss points
    for gp_idx, (xi_val, eta_val) in enumerate(gauss_points):
        F, E, e, eps = compute_kinematics_vectorized(
            X_elems, x_elems, xi=xi_val, eta=eta_val
        )

        N1 = 0.25 * (1 - xi_val) * (1 - eta_val)
        N2 = 0.25 * (1 + xi_val) * (1 - eta_val)
        N3 = 0.25 * (1 + xi_val) * (1 + eta_val)
        N4 = 0.25 * (1 - xi_val) * (1 + eta_val)
        N = np.array([N1, N2, N3, N4])

        # Vectorized calculation of physical GP coordinates for all elements
        gp_coords[:, gp_idx, :] = np.einsum("j, Njk -> Nk", N, x_elems)

        # Store results in the pre-allocated block
        F_gp_array[:, gp_idx, :, :] = F
        E_gp_array[:, gp_idx, :, :] = E
        e_gp_array[:, gp_idx, :, :] = e
        eps_gp_array[:, gp_idx, :, :] = eps

    # Flatten the arrays so each Gauss point acts as a separate entry
    gp_coords = gp_coords.reshape(-1, 2)
    F_results = F_gp_array.reshape(-1, 2, 2)
    E_results = E_gp_array.reshape(-1, 2, 2)
    e_results = e_gp_array.reshape(-1, 2, 2)
    eps_results = eps_gp_array.reshape(-1, 2, 2)


# Create a results subdirectory
results_dir = os.path.join(current_dir, "results")
os.makedirs(results_dir, exist_ok=True)

# Save Tensors as .npy files inside the results folder
np.save(
    os.path.join(results_dir, f"Deformation_Gradient_F_{COMPUTE_MODE}.npy"), F_results
)
np.save(
    os.path.join(results_dir, f"Green_Lagrange_Strain_E_{COMPUTE_MODE}.npy"), E_results
)
np.save(
    os.path.join(results_dir, f"Eulerian_Strain_e_{COMPUTE_MODE}.npy"),
    e_results,
)
np.save(
    os.path.join(results_dir, f"Engineering_Strain_eps_{COMPUTE_MODE}.npy"), eps_results
)
print(f"Successfully saved tensor data to {results_dir}")


# Plotting
fig, ax = plt.subplots(figsize=(8, 5))
for i in range(ny_nodes):
    start = i * nx_nodes
    end = start + nx_nodes

    # Add legends only on the first row
    if i == 0:
        ax.plot(
            X_nodes[start:end, 0],
            X_nodes[start:end, 1],
            "b-",
            linewidth=1.5,
            label="undeformed",
        )
        ax.plot(
            x_nodes[start:end, 0],
            x_nodes[start:end, 1],
            "r--",
            linewidth=1.5,
            label="deformed",
        )
    else:
        ax.plot(X_nodes[start:end, 0], X_nodes[start:end, 1], "b-", linewidth=1.5)
        ax.plot(x_nodes[start:end, 0], x_nodes[start:end, 1], "r--", linewidth=1.5)

# Plot vertical grid lines (columns)
for j in range(nx_nodes):
    ax.plot(X_nodes[j::nx_nodes, 0], X_nodes[j::nx_nodes, 1], "b-", linewidth=1.5)
    ax.plot(x_nodes[j::nx_nodes, 0], x_nodes[j::nx_nodes, 1], "r--", linewidth=1.5)

ax.plot(X_nodes[:, 0], X_nodes[:, 1], "bo", markerfacecolor="none", markersize=6)
ax.plot(x_nodes[:, 0], x_nodes[:, 1], "ro", markerfacecolor="none", markersize=6)

ax.set_xlabel("X (m)", fontsize=14)
ax.set_ylabel("Y (m)", fontsize=14)
ax.set_xlim(0, 1.05)
ax.set_ylim(-0.25, 1.25)
ax.set_xticks(np.arange(0, 1.2, 0.2))
ax.set_yticks(np.arange(-0.2, 1.4, 0.2))
ax.legend(loc="upper left", framealpha=1, edgecolor="black")

plt.tight_layout()
deform_plot_path = os.path.join(results_dir, "Deformation_Plot.png")
plt.savefig(deform_plot_path, dpi=300)
plt.show()


# Plotting Strain Components
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
titles = ["$E_{11}$", "$E_{22}$", "$E_{12}$"]
verts = [X_nodes[element] for element in conn]

if COMPUTE_MODE == "element":
    # Extract components
    E11 = E_results[:, 0, 0]
    E22 = E_results[:, 1, 1]
    E12 = E_results[:, 0, 1]
    strains = [E11, E22, E12]
    verts = [X_nodes[element] for element in conn]

    for ax, strain_data, title in zip(axes, strains, titles):
        poly = PolyCollection(verts, cmap="jet", edgecolor="black", linewidth=0.5)
        poly.set_array(np.array(strain_data))
        ax.add_collection(poly)
        fig.colorbar(poly, ax=ax, fraction=0.046, pad=0.04)

elif COMPUTE_MODE == "nodal":
    # Extract components and reshape to grid for contour plotting
    E11 = E_results[:, 0, 0].reshape((ny_nodes, nx_nodes))
    E22 = E_results[:, 1, 1].reshape((ny_nodes, nx_nodes))
    E12 = E_results[:, 0, 1].reshape((ny_nodes, nx_nodes))
    strains = [E11, E22, E12]

    for ax, strain_data, title in zip(axes, strains, titles):
        contour = ax.contourf(X, Y, strain_data, levels=50, cmap="jet")
        fig.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)

elif COMPUTE_MODE == "gauss_surface":
    E11 = E_results[:, 0, 0]
    E22 = E_results[:, 1, 1]
    E12 = E_results[:, 0, 1]
    strains = [E11, E22, E12]

    # Create a dense uniform meshgrid over the bounds of the deformed geometry
    grid_x, grid_y = np.mgrid[
        min(x_nodes[:, 0]) : max(x_nodes[:, 0]) : 200j,
        min(x_nodes[:, 1]) : max(x_nodes[:, 1]) : 200j,
    ]

    for ax, strain_data, title in zip(axes, strains, titles):
        # Interpolate scattered Gauss point strains onto the dense grid
        grid_strain = griddata(gp_coords, strain_data, (grid_x, grid_y), method="cubic")

        # Plot using contourf
        contour = ax.contourf(grid_x, grid_y, grid_strain, levels=50, cmap="jet")
        fig.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)

# Format the strain plots
for ax, title in zip(axes, titles):
    ax.set_title(f"Green Strain: {title} ({COMPUTE_MODE})", fontsize=14)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("X (m)", fontsize=12)
    ax.set_ylabel("Y (m)", fontsize=12)

plt.tight_layout()
strain_plot_path = os.path.join(results_dir, f"Strain_Plot_{COMPUTE_MODE}.png")
plt.savefig(strain_plot_path, dpi=300)
plt.show()
