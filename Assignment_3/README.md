# Assignment 3: 2D Co-rotational Truss FEA Solver

This directory contains a custom-built, high-performance Finite Element Analysis (FEA) engine written in **Fortran 2008**. It is designed to solve 2D truss structures undergoing large displacements (geometric nonlinearity) using a **Co-rotational kinematic formulation** and a fully implicit **Newton-Raphson** iterative solver.

## Features

* **Geometric Nonlinearity**: Captures large rotations and buckling behaviors using updated Lagrangian coordinates.
* **Dynamic Memory**: Fully dynamic allocation allows for arbitrary mesh sizes and load cases without recompiling the source code.
* **Dense Linear Algebra**: Utilizes LAPACK's `dgesv` for robust, high-speed matrix inversion.
* **Automated Post-Processing**: Python script automatically generates publication-ready plots of the mesh (undeformed vs. deformed) and the nonlinear load-displacement structural response.

---

## Prerequisites

To compile and run the solver, you need a Fortran compiler and standard linear algebra libraries.

* **Linux (Ubuntu/Debian)**: `sudo apt-get install gfortran libblas-dev liblapack-dev`
* **macOS**: `brew install gcc lapack`
* **Python Packages**: `pip install pandas matplotlib numpy`

---

## Input File Format (.inp)

The solver is driven by plain-text input files. The engine reads parameters in a strict top-to-bottom order. 

**Example `job1.inp`:**

```plaintext
6               ! n_bays (Number of truss bays)
240.0           ! L_total (Total length in mm)
40.0            ! H_total (Total height in mm)
65.0            ! A_cross (Cross-sectional area in mm^2)
200000.0        ! E_mod (Young's Modulus in MPa)
90              ! n_steps (Number of Newton-Raphson load increments)
2               ! num_supports (Number of fixed boundary nodes)
0.0   0.0   1.0  1.0   ! Support 1: [X_coord, Y_coord, Fix_X, Fix_Y] (1.0 = Fixed)
0.0   40.0  1.0  1.0   ! Support 2: [X_coord, Y_coord, Fix_X, Fix_Y]
1               ! num_loads (Number of applied point loads)
240.0 40.0  0.0 -90000.0 ! Load 1: [X_coord, Y_coord, Force_X, Force_Y] (Forces in N)
```

---

## Build and Execution Instructions

### 1. Compile the Solver
Navigate to the `Assignment_3` directory and compile the Fortran modules into a single executable, linking the LAPACK and BLAS libraries:

```bash
gfortran -O3 -Wall kinds_mod.f90 material_mod.f90 solver_mod.f90 main.f90 -o fea_solver -llapack -lblas
```

### 2. Create Output Directory and Run the Engine
Before running the solver, you must create an `outputs/` directory to store the generated data. Then, execute the solver and pass your input file as an argument:

```bash
mkdir outputs
./fea_solver job1.inp
```

### 3. Generate Visualizations
Run the Python post-processor, passing the exact name of your job. The script will dynamically read the CSV outputs and generate the plots.

```bash
python plot_truss.py job1
```

---

## Output Artifacts

For a given job named `job1`, the pipeline generates the following files in the `outputs/` directory:

* **CSV Data**:
    * `job1_undeformed_nodes.csv` & `job1_undeformed_elements.csv`: Structural topology verification.
    * `job1_deformed_nodes.csv`: Final global coordinates of the converged structure.
    * `job1_history.csv`: Step-by-step logging of tip displacements and applied loads.
* **Plots**:
    * `1_undeformed_truss.png`: Verification plot of the input deck, loads, and BCs.
    * `2_mesh_overlay.png`: Visual overlay of the original (gray) and deformed (red) configurations.
    * `3_load_displacement.png`: The nonlinear structural response curves.
