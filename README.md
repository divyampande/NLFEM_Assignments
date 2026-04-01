# Nonlinear Finite Element Methods (NLFEM) Coursework

This repository contains the programming assignments, computational solvers, and project reports for the NLFEM coursework at **IIT Madras**.

**Author:** Divyam Pandey

---

## Repository Structure

* **Assignment_1/**: Finite Deformation Kinematics of a Non-Uniformly Heated Composite Structure. Includes a custom Python solver for computing the Deformation Gradient and Green-Lagrange strain tensors from extracted graphical data.
* **Assignment_2/**: Modifications and extensions to an existing MATLAB codebase for nonlinear finite element analysis. The core scripts are located within the `matlab_source/` sub-directory.
* **Assignment_3/**: 2D Co-rotational Truss Finite Element Solver. A custom-built Fortran engine that handles geometric nonlinearity using a Newton-Raphson iterative scheme. Includes a Python-based post-processing pipeline for automated visualization of load-displacement curves and deformed configurations.
* *(Future assignments and projects will be added here)*

---

## Dependencies

Because this repository covers various computational approaches, it relies on a mixed technology stack:

### Python (Post-processing & Kinematics)
* `numpy`
* `pandas`
* `matplotlib`

### Fortran (High-Performance Solvers)
* GNU Fortran Compiler (`gfortran`)
* `LAPACK` & `BLAS` libraries for dense linear algebra

### MATLAB
* Standard MATLAB environment (for Assignment 2)

---

## License

This project is for educational and academic purposes.
