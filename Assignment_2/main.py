import numpy as np
import matplotlib.pyplot as plt
import time


def solve_nonlinear_bar(L, A, E, m, F_total, n_elem, n_steps):
    n_node = n_elem + 1
    x = np.linspace(0, L, n_node)

    gauss_pts = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3)])
    gauss_wts = np.array([1.0, 1.0])

    force_steps = np.linspace(0, F_total, n_steps + 1)

    u = np.zeros(n_node)
    u_history, f_history = [0.0], [0.0]
    strain_history, stress_history = [0.0], [0.0]

    for step in range(1, n_steps + 1):
        F_ext = force_steps[step]

        for _ in range(50):
            K = np.zeros((n_node, n_node))
            R = np.zeros(n_node)
            R[-1] = F_ext

            for el in range(n_elem):
                n1, n2 = el, el + 1
                le = x[n2] - x[n1]
                detJ = le / 2.0
                invJ = 1.0 / detJ

                ue = u[[n1, n2]]
                ke = np.zeros((2, 2))
                f_int = np.zeros(2)

                for gp, gw in zip(gauss_pts, gauss_wts):
                    dN = np.array([-0.5, 0.5])
                    B = dN * invJ

                    strain = np.dot(B, ue)
                    stress = E * np.arctan(m * strain)
                    Et = (E * m) / (1.0 + (m * strain) ** 2)

                    ke += np.outer(B, B) * Et * A * detJ * gw
                    f_int += B * stress * A * detJ * gw

                K[np.ix_([n1, n2], [n1, n2])] += ke
                R[[n1, n2]] -= f_int

            K[0, :] = 0
            K[:, 0] = 0
            K[0, 0] = 1.0
            R[0] = 0.0

            du = np.linalg.solve(K, R)
            u += du

            if np.linalg.norm(du) < 1e-6:
                break

        u_history.append(u[-1])
        f_history.append(F_ext)

        le_last = x[-1] - x[-2]
        B_last = np.array([-0.5, 0.5]) / (le_last / 2.0)
        strain_last = np.dot(B_last, u[[-2, -1]])
        stress_last = E * np.arctan(m * strain_last)

        strain_history.append(strain_last)
        stress_history.append(stress_last)

    return u_history, f_history, strain_history, stress_history


if __name__ == "__main__":
    L = 1.0
    A = 1e-4
    E = 100e6
    m = 40.0
    F_total = 10000.0

    n_elem = 10000
    n_steps = 20
    print("Starting benchmark...")

    # 1. Start the stopwatch
    start_time = time.perf_counter()
    u_hist, f_hist, eps_hist, sig_hist = solve_nonlinear_bar(
        L, A, E, m, F_total, n_elem, n_steps
    )
    # 3. Stop the stopwatch
    end_time = time.perf_counter()

    # 4. Calculate the difference
    execution_time = end_time - start_time

    print(f"Solver execution time: {execution_time:.6f} seconds")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(eps_hist, sig_hist, "b-o", linewidth=2)
    ax1.set_xlabel("Strain ($\epsilon$)")
    ax1.set_ylabel("Stress ($\sigma$) [Pa]")
    ax1.set_title("Stress-Strain in Last Element")
    ax1.grid(True)

    ax2.plot(u_hist, f_hist, "r-s", linewidth=2)
    ax2.set_xlabel("Displacement at $x=L$ [m]")
    ax2.set_ylabel("Applied Force [N]")
    ax2.set_title("Force vs Reaction Displacement")
    ax2.grid(True)

    plt.tight_layout()
    plt.show()
