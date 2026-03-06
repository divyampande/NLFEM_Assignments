import numpy as np
import matplotlib.pyplot as plt
import os


def correct_and_plot_nodes(
    nx_elem=4,
    ny_elem=4,
    Lx=1.0,
    Ly=1.0,
):

    current_dir = (
        os.path.dirname(os.path.abspath(__file__))
        if "__file__" in globals()
        else os.getcwd()
    )

    csv_filename = os.path.join(current_dir, "Deformed_Coords.csv")
    output_filename = os.path.join(current_dir, "corrected_deformed_coords.csv")

    nx_nodes = nx_elem + 1
    ny_nodes = ny_elem + 1
    total_nodes = nx_nodes * ny_nodes
    # Load manual data
    try:
        original = np.loadtxt(csv_filename, delimiter=",")
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return

    if len(original) != total_nodes:
        print(
            f"Error: Expected {total_nodes} nodes for a {nx_elem}x{ny_elem} mesh, but found {len(original)}."
        )
        return

    corrected = original.copy()
    center_y = Ly / 2.0

    # Enforce Centerline strictly at Y = 0.5
    if ny_elem % 2 == 0:
        mid_row = ny_nodes // 2
        for col in range(nx_nodes):
            center_node_idx = mid_row * nx_nodes + col
            corrected[center_node_idx, 1] = center_y

    # Enforce Horizontal Symmetry (about Y = 0.5)
    for row in range(ny_nodes // 2):
        top_row = ny_nodes - 1 - row

        for col in range(nx_nodes):
            idx_bottom = row * nx_nodes + col
            idx_top = top_row * nx_nodes + col

            # Average the X coordinates
            avg_x = (original[idx_bottom, 0] + original[idx_top, 0]) / 2.0
            corrected[idx_bottom, 0] = avg_x
            corrected[idx_top, 0] = avg_x

            # Average the absolute distance from the physical center_y
            dist_bottom = center_y - original[idx_bottom, 1]
            dist_top = original[idx_top, 1] - center_y
            avg_dist = (dist_bottom + dist_top) / 2.0

            # Apply new symmetric Y coordinates
            corrected[idx_bottom, 1] = center_y - avg_dist
            corrected[idx_top, 1] = center_y + avg_dist

    # Enforce Boundary Conditions (X = 0 and X = 1)
    for row in range(ny_nodes):
        left_edge_idx = row * nx_nodes + 0
        right_edge_idx = row * nx_nodes + nx_elem

        corrected[left_edge_idx, 0] = 0.0
        corrected[right_edge_idx, 0] = Lx

    # Save the corrected data
    np.savetxt(output_filename, corrected, delimiter=",", fmt="%.6f")
    print(f"Corrected coordinates saved to {output_filename}")

    # Plotting Original vs. Corrected
    plt.figure(figsize=(10, 8))

    # Plot original nodes
    plt.scatter(
        original[:, 0],
        original[:, 1],
        color="blue",
        alpha=0.3,
        s=100,
        label="Original (Manual)",
    )

    # Plot corrected nodes
    plt.scatter(
        corrected[:, 0],
        corrected[:, 1],
        color="red",
        s=50,
        marker="x",
        label="Corrected (Symmetric/Constrained)",
    )

    # Draw dynamic grid lines
    for i in range(ny_nodes):
        # Rows
        start = i * nx_nodes
        end = start + nx_nodes
        plt.plot(corrected[start:end, 0], corrected[start:end, 1], "r--", alpha=0.6)

    for j in range(nx_nodes):
        # Columns
        plt.plot(corrected[j::nx_nodes, 0], corrected[j::nx_nodes, 1], "r--", alpha=0.6)

    plt.title(f"Deformed Shape Correction ({nx_elem}x{ny_elem} Mesh)")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.grid(True, linestyle=":", alpha=0.7)
    plt.legend()
    plt.axis("equal")
    plt.show()


correct_and_plot_nodes()
