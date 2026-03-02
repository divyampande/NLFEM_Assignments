import numpy as np
import matplotlib.pyplot as plt


def correct_and_plot_nodes(
    csv_filename, output_filename="NLFEM_Assignments\Assignment_1\corrected_nodes.csv"
):
    # 1. Load original manual data
    try:
        original = np.loadtxt(csv_filename, delimiter=",")
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return

    if len(original) != 25:
        print("Error: Expected exactly 25 nodes in the CSV.")
        return

    corrected = original.copy()

    # 2. Enforce Centerline strictly at Y = 0.5
    # The middle row consists of indices 10, 11, 12, 13, 14
    center_idx = [10, 11, 12, 13, 14]
    corrected[center_idx, 1] = 0.5

    # 3. Enforce Horizontal Symmetry (about Y = 0.5)
    # Pair Row 0 (bottom) with Row 4 (top)
    # Pair Row 1 (mid-bottom) with Row 3 (mid-top)
    row_pairs = [(0, 4), (1, 3)]

    for bottom_row, top_row in row_pairs:
        for col in range(5):
            idx_bottom = bottom_row * 5 + col
            idx_top = top_row * 5 + col

            # Average the X coordinates
            avg_x = (original[idx_bottom, 0] + original[idx_top, 0]) / 2.0
            corrected[idx_bottom, 0] = avg_x
            corrected[idx_top, 0] = avg_x

            # Average the absolute distance from Y=0.5
            dist_bottom = 0.5 - original[idx_bottom, 1]
            dist_top = original[idx_top, 1] - 0.5
            avg_dist = (dist_bottom + dist_top) / 2.0

            # Apply new symmetric Y coordinates
            corrected[idx_bottom, 1] = 0.5 - avg_dist
            corrected[idx_top, 1] = 0.5 + avg_dist

    # 4. Enforce Boundary Conditions (X = 0 and X = 1)
    # Done last to ensure symmetry averaging didn't pull the boundaries off the wall
    left_edge_idx = [0, 5, 10, 15, 20]
    right_edge_idx = [4, 9, 14, 19, 24]

    corrected[left_edge_idx, 0] = 0.0
    corrected[right_edge_idx, 0] = 1.0

    # 5. Save the corrected data
    np.savetxt(output_filename, corrected, delimiter=",", fmt="%.6f")
    print(f"Corrected coordinates saved to {output_filename}")

    # 6. Plotting Original vs. Corrected
    plt.figure(figsize=(10, 8))

    # Plot original nodes (Faint blue)
    plt.scatter(
        original[:, 0],
        original[:, 1],
        color="blue",
        alpha=0.3,
        s=100,
        label="Original (Manual)",
    )

    # Plot corrected nodes (Solid red)
    plt.scatter(
        corrected[:, 0],
        corrected[:, 1],
        color="red",
        s=50,
        marker="x",
        label="Corrected (Symmetric/Constrained)",
    )

    # Draw grid lines for the corrected mesh
    for i in range(5):
        # Rows
        start = i * 5
        end = start + 5
        plt.plot(corrected[start:end, 0], corrected[start:end, 1], "r--", alpha=0.6)
        # Columns
        plt.plot(corrected[i::5, 0], corrected[i::5, 1], "r--", alpha=0.6)

    plt.title("Deformed Shape: Manual Extraction vs. Corrected Data")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.grid(True, linestyle=":", alpha=0.7)
    plt.legend()
    plt.axis("equal")
    plt.show()


# Run the function (replace with your actual WebPlotDigitizer filename)
correct_and_plot_nodes("NLFEM_Assignments\Assignment_1\Deformed_Coords.csv")
