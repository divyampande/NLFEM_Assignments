import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def assemble_system_vectorized(u, L, A, E, m, n_elem):
    n_node = n_elem + 1
    le = L / n_elem
    detJ = le / 2.0
    invJ = 1.0 / detJ

    # Vectorized element node indices
    n1 = np.arange(n_elem)
    n2 = n1 + 1

    # Vectorized element displacements
    u1 = u[n1]
    u2 = u[n2]

    # B-matrix components (constant for C0 linear 1D elements)
    B1 = -0.5 * invJ
    B2 = 0.5 * invJ

    # Strain is calculated for all elements simultaneously
    strain = B1 * u1 + B2 * u2

    # 2-point Gauss Quadrature
    # (Since B is constant in a linear element, the loop just aggregates the weights)
    gauss_wts = np.array([1.0, 1.0])

    ke_00, ke_01, ke_10, ke_11 = (
        np.zeros(n_elem),
        np.zeros(n_elem),
        np.zeros(n_elem),
        np.zeros(n_elem),
    )
    fe_0, fe_1 = np.zeros(n_elem), np.zeros(n_elem)

    for gw in gauss_wts:
        stress = E * np.arctan(m * strain)
        Et = (E * m) / (1.0 + (m * strain) ** 2)

        # Vectorized stiffness integration
        factor = Et * A * detJ * gw
        ke_00 += B1 * B1 * factor
        ke_01 += B1 * B2 * factor
        ke_10 += B2 * B1 * factor
        ke_11 += B2 * B2 * factor

        # Vectorized internal force integration
        f_factor = stress * A * detJ * gw
        fe_0 += B1 * f_factor
        fe_1 += B2 * f_factor

    # Flatten indices and values for ultra-fast sparse COO assembly
    I = np.column_stack([n1, n1, n2, n2]).ravel()
    J = np.column_stack([n1, n2, n1, n2]).ravel()
    V = np.column_stack([ke_00, ke_01, ke_10, ke_11]).ravel()

    K_global = coo_matrix((V, (I, J)), shape=(n_node, n_node)).tocsr()

    # Fast assembly of global internal force vector using unbuffered addition
    f_int = np.zeros(n_node)
    np.add.at(f_int, n1, fe_0)
    np.add.at(f_int, n2, fe_1)

    return K_global, f_int


def solve_nonlinear_bar_vectorized(L, A, E, m, F_total, n_elem, n_steps):
    n_node = n_elem + 1
    force_steps = np.linspace(0, F_total, n_steps + 1)

    u = np.zeros(n_node)
    u_history, f_history = [0.0], [0.0]
    strain_history, stress_history = [0.0], [0.0]

    for step in range(1, n_steps + 1):
        F_ext = np.zeros(n_node)
        F_ext[-1] = force_steps[step]

        for _ in range(50):
            K_global, f_int = assemble_system_vectorized(u, L, A, E, m, n_elem)
            R = F_ext - f_int

            # Apply boundary condition (fixed at x=0) by slicing
            K_free = K_global[1:, 1:]
            R_free = R[1:]

            du_free = spsolve(K_free, R_free)
            u[1:] += du_free

            if np.linalg.norm(du_free) < 1e-6:
                break

        u_history.append(u[-1])
        f_history.append(force_steps[step])

        # Post-process stress/strain in the last element to satisfy assignment requirements
        le_last = L / n_elem
        strain_last = (-u[-2] + u[-1]) / le_last
        stress_last = E * np.arctan(m * strain_last)

        strain_history.append(strain_last)
        stress_history.append(stress_last)

    return u_history, f_history, strain_history, stress_history


if __name__ == "__main__":
    L, A, E, m, F_total = 1.0, 1e-4, 100e6, 40.0, 10000.0
    n_elem, n_steps = 1000, 20

    u_hist, f_hist, eps_hist, sig_hist = solve_nonlinear_bar_vectorized(
        L, A, E, m, F_total, n_elem, n_steps
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Raw strings (r"") prevent Python escape sequence warnings
    ax1.plot(eps_hist, sig_hist, "b-o", linewidth=2)
    ax1.set_xlabel(r"Strain ($\epsilon$)")
    ax1.set_ylabel(r"Stress ($\sigma$) [Pa]")
    ax1.set_title("Stress-Strain in Last Element")
    ax1.grid(True)

    ax2.plot(u_hist, f_hist, "r-s", linewidth=2)
    ax2.set_xlabel(r"Displacement at $x=L$ [m]")
    ax2.set_ylabel("Applied Force [N]")
    ax2.set_title("Force vs Reaction Displacement")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()
