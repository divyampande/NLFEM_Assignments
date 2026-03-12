import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from numba import njit
import time


# Numba compiles this function to C-speed machine code just-in-time
@njit
def assemble_system(u, x, E, m, A):
    n_node = len(u)
    n_elem = n_node - 1

    # Pre-allocate arrays for sparse matrix creation (4 entries per 1D element)
    I = np.zeros(4 * n_elem, dtype=np.int32)
    J = np.zeros(4 * n_elem, dtype=np.int32)
    V = np.zeros(4 * n_elem, dtype=np.float64)
    f_int = np.zeros(n_node, dtype=np.float64)

    # 2-point Gauss Quadrature
    gauss_pts = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])
    gauss_wts = np.array([1.0, 1.0])

    idx = 0
    for el in range(n_elem):
        n1, n2 = el, el + 1
        le = x[n2] - x[n1]
        detJ = le / 2.0
        invJ = 1.0 / detJ

        u1, u2 = u[n1], u[n2]

        ke00, ke01, ke10, ke11 = 0.0, 0.0, 0.0, 0.0
        f1, f2 = 0.0, 0.0

        for gp in range(2):
            gw = gauss_wts[gp]
            B1 = -0.5 * invJ
            B2 = 0.5 * invJ

            strain = B1 * u1 + B2 * u2
            stress = E * np.arctan(m * strain)
            Et = (E * m) / (1.0 + (m * strain) ** 2)

            # Element stiffness and internal force integration
            factor = Et * A * detJ * gw
            ke00 += B1 * B1 * factor
            ke01 += B1 * B2 * factor
            ke10 += B2 * B1 * factor
            ke11 += B2 * B2 * factor

            f_factor = stress * A * detJ * gw
            f1 += B1 * f_factor
            f2 += B2 * f_factor

        # Store data for sparse COO format
        I[idx] = n1
        J[idx] = n1
        V[idx] = ke00
        idx += 1
        I[idx] = n1
        J[idx] = n2
        V[idx] = ke01
        idx += 1
        I[idx] = n2
        J[idx] = n1
        V[idx] = ke10
        idx += 1
        I[idx] = n2
        J[idx] = n2
        V[idx] = ke11
        idx += 1

        f_int[n1] += f1
        f_int[n2] += f2

    return I, J, V, f_int


def solve_nonlinear_bar_fast(L, A, E, m, F_total, n_elem, n_steps):
    n_node = n_elem + 1
    x = np.linspace(0, L, n_node)
    force_steps = np.linspace(0, F_total, n_steps + 1)

    u = np.zeros(n_node)
    u_history, f_history = [0.0], [0.0]
    strain_history, stress_history = [0.0], [0.0]

    for step in range(1, n_steps + 1):
        F_ext = np.zeros(n_node)
        F_ext[-1] = force_steps[step]

        # Newton-Raphson Loop
        for _ in range(50):
            I, J, V, f_int = assemble_system(u, x, E, m, A)
            R = F_ext - f_int

            # Build Sparse Matrix (COO to CSR format for fast solving)
            K_global = coo_matrix((V, (I, J)), shape=(n_node, n_node)).tocsr()

            # Slice matrix to apply BCs (Node 0 is fixed, solve only for free DOFs)
            K_free = K_global[1:, 1:]
            R_free = R[1:]

            # Fast sparse linear solve
            du_free = spsolve(K_free, R_free)
            u[1:] += du_free

            if np.linalg.norm(du_free) < 1e-6:
                break

        u_history.append(u[-1])
        f_history.append(force_steps[step])

        # Post-process stress/strain in the last element [cite: 27]
        le_last = x[-1] - x[-2]
        strain_last = (-u[-2] + u[-1]) / le_last
        stress_last = E * np.arctan(m * strain_last)

        strain_history.append(strain_last)
        stress_history.append(stress_last)

    return u_history, f_history, strain_history, stress_history


if __name__ == "__main__":
    # If you run this, you'll see a slight delay on the first execution
    # as Numba compiles the C code, and then it will execute instantly.
    L, A, E, m, F_total = 1.0, 1e-4, 100e6, 40.0, 10000.0
    n_elem, n_steps = 10000, 50

    print("Starting benchmark...")

    # 1. Start the stopwatch
    start_time = time.perf_counter()
    u_hist, f_hist, eps_hist, sig_hist = solve_nonlinear_bar_fast(
        L, A, E, m, F_total, n_elem, n_steps
    )
    # 3. Stop the stopwatch
    end_time = time.perf_counter()

    # 4. Calculate the difference
    execution_time = end_time - start_time

    print(f"Solver execution time: {execution_time:.6f} seconds")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(eps_hist, sig_hist, "b-o", linewidth=2)
    ax1.set_xlabel("Strain ($\\epsilon$)")
    ax1.set_ylabel("Stress ($\\sigma$) [Pa]")
    ax1.set_title("Stress-Strain in Last Element")
    ax1.grid(True)

    ax2.plot(u_hist, f_hist, "r-s", linewidth=2)
    ax2.set_xlabel("Displacement at $x=L$ [m]")
    ax2.set_ylabel("Applied Force [N]")
    ax2.set_title("Force vs Reaction Displacement")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()
