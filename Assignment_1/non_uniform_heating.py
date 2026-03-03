import numpy as np
import matplotlib.pyplot as plt
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

# Loop through the elements to draw the outlines
for i, element in enumerate(conn):
    elem_idx = np.append(element, element[0])

    X_elem = X_nodes[elem_idx, 0]
    Y_elem = X_nodes[elem_idx, 1]

    x_elem_def = x_nodes[elem_idx, 0]
    y_elem_def = x_nodes[elem_idx, 1]

    # Add legends only on the first loop iteration
    if i == 0:
        ax.plot(X_elem, Y_elem, "b-", linewidth=1.5, label="undeformed")
        ax.plot(x_elem_def, y_elem_def, "r--", linewidth=1.5, label="deformed")
    else:
        ax.plot(X_elem, Y_elem, "b-", linewidth=1.5)
        ax.plot(x_elem_def, y_elem_def, "r--", linewidth=1.5)

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
