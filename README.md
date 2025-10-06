# Data-Driven Model Order Reduction Techniques for Dynamic Systems

**Author**: Charukhesh B.R

## 📖 Project Overview

This repository provides a MATLAB implementation and comparative analysis of several state-of-the-art **Model Order Reduction** techniques. The project explores the application of these methods to a variety of linear and nonlinear dynamic systems, from simple LTI controllers to complex fluid dynamics and structural vibration models.

The primary goal is to demonstrate the process, advantages, and limitations of each technique, covering:
- **Balanced Truncation** for linear time-invariant (LTI) systems
- **Dynamic Mode Decomposition (DMD)** for data-driven modeling of complex systems
- **Total Least Squares DMD (TLS-DMD)** for robust system identification from noisy data
- **Proper Orthogonal Decomposition (POD) with Galerkin Projection** for reducing the complexity of nonlinear partial differential equations (PDEs)

## Techniques (Algorithms) used in this project

The project provides a hands-on implementation of the following foundational MOR methods:

1.  **Balanced Truncation:** Reduces LTI systems by finding a state-space realization where the controllability ($W_c$) and observability ($W_o$) Gramians are equal and diagonal. States corresponding to small Hankel Singular Values ($\sigma_i$) are truncated, preserving stability and providing an *a priori* error bound

2.  **Dynamic Mode Decomposition (DMD):** A data-driven method that finds a best-fit linear operator to approximate the evolution of a system from snapshot data. It is highly effective for systems with dominant oscillatory or periodic behavior, such as fluid flows. The core assumption is $x_{k+1} \approx Ax_k$

3.  **Total Least Squares DMD (TLS-DMD):** A robust variant of DMD designed to handle measurement noise. It assumes noise is present in all measurements, leading to a more accurate and unbiased estimation of the system's dynamics compared to standard DMD

4.  **POD-Galerkin ROM for Nonlinear PDEs:** An intrusive method where the solution to a PDE is approximated in a low-dimensional basis derived from Proper Orthogonal Decomposition ($u(t) \approx \Phi\hat{u}(t)$). The governing equations are then projected onto this basis (Galerkin projection) to create a small, computationally efficient system of ODEs

## 📁 Repository Structure

```
ModelReductionProject/
│
├── main.m
├── plotCylinder.m
├── mat files/
│ ├── mass_spring_damper_matrices.mat
│ └── cylinder_flow_snapshots.mat
│
├── report.pdf
└── README.md
```

- **main.m**: The main interactive MATLAB script that runs all analyses
- **plotCylinder.m**: A helper function for visualizing the cylinder flow fields in the DMD analysis (Section 6)
- **mat files/**: A directory containing the datasets required for specific analyses
- **report.pdf**: A detailed report covering the theoretical background, implementation details, and results for each technique

## How to Run and Use This Project

1.  **Clone the Repository**
    ```bash
    git clone https://github.com/Charukhesh/DataDrivenMORs.git
    cd DataDrivenMORs
    ```

2.  **Open and Run the Main File in MATLAB**
    - Open `main.m` in MATLAB
    - Run the script. You will be prompted in the command window to select which analysis you want to perform by entering a section number

### Guide to Project Sections (by Question Number)

The `main.m` script is organized into sections corresponding to the sections in the project report.

---
#### **Sections 1-4: Balanced Truncation for an LTI System**
-   **Description**: A four-part deep dive into balanced truncation.
-   **Input**:
    - `Section No.`: **1** - Computes and compares Gramians using `lyap` and `gram`
    - `Section No.`: **2** - Analyzes Gramian eigenvalues and the effect of state transformations
    - `Section No.`: **3** - Compares a manual implementation of balanced realization against MATLAB's `balreal` command
    - `Section No.`: **4** - Performs balanced truncation, compares step responses, and validates the theoretical error bounds

---
#### **Section 5: Balanced Truncation for a Structural Vibration Model**
-   **Description**: Applies balanced truncation to a 12th-order mass-spring-damper model. This case is special because the system is **not controllable**, demonstrating the need to first find a minimal realization (`minreal`) before balancing.
-   **Input**: `Section No.`: **5**

---
#### **Section 6: Dynamic Mode Decomposition (DMD) for Cylinder Flow**
-   **Description**: Uses DMD to create a low-rank linear model of a complex fluid flow (Karman vortex street). Determines the model rank based on energy criteria and visualizes the high-fidelity reconstruction.
-   **Input**: `Section No.`: **6**

---
#### **Section 7: Total Least Squares DMD (TLS-DMD)**
-   **Description**: A Monte Carlo simulation comparing the robustness of Standard DMD vs. TLS-DMD for a system with measurement noise. The results clearly show that TLS-DMD is an unbiased estimator, while Standard DMD is biased.
-   **Input**: `Section No.`: **7**

---
#### **Section 8: POD-Galerkin ROM for Burgers' Equation**
-   **Description**: Implements an intrusive POD-based ROM for the nonlinear 1D viscous Burgers' equation. It demonstrates the full workflow: FOM simulation, POD basis extraction, Galerkin projection to build the ROM, and error analysis.
-   **Input**: `Section No.`: **8**

## Key Analyses

1.  **Balanced Truncation Requires Minimal Realization**: The structural vibration model (Section 5) highlights that balanced truncation is only applicable to controllable and observable systems. The analysis shows that using `minreal` to first remove the 10 uncontrollable states results in an exact 2nd-order model, demonstrating a perfect application of the theory.

2.  **DMD for Complex Systems**: DMD successfully extracts the dominant dynamics of the cylinder flow with just 9 modes (out of a possible 150), capturing over 99.9% of the system's energy. This demonstrates its power in creating simple, linear models for high-dimensional, nonlinear phenomena.

3.  **TLS-DMD is Superior for Noisy Data**: The Monte Carlo simulation (Section 7) provides definitive proof that TLS-DMD is a more robust and accurate method for system identification from noisy data. It acts as an **unbiased estimator** with lower variance, whereas Standard DMD is systematically biased by noise.

4.  **POD-Galerkin for Nonlinear PDEs**: The analysis of Burgers' equation (Section 8) shows that its smooth dynamics are inherently low-dimensional. An optimal ROM size of **r=3 or 4** (out of 101 states) is quantitatively justified through an "accuracy vs. efficiency" analysis, achieving a massive computational speedup with negligible loss in accuracy.