# Constrained minmax density transportation for parabolic PDEs: A direct optimal control perspective
**Authors: Vaibhav Upadhyay, Siddhartha Ganguly, Kenji Kashima, Debasish Chatterjee**

This repository contains a Julia implementation for solving a **robust optimal control problem** for a **1D parabolic partial differential equation (PDE)**. The setup models a boundary-controlled diffusion system with bounded uncertainties, optimized over a finite time horizon using a bilevel optimization structure.

---

## Repository Structure

This repository contains three different problem setups, each organized into a separate folder:

- `diff_dist_control/`  
  Solves a **robust optimal control** problem for the **heat equation**, transferring a **non-zero initial distribution** to a **prescribed non-zero terminal distribution**.

- `dist_null_control/`  
  Solves a **robust null control** problem for the **heat equation**, transferring a **non-zero initial distribution** to the **zero (null) distribution**.

- `reac_advec_diff_control/`  
  Solves a **robust optimal control** problem for the **reaction-advection-diffusion equation**, transferring a **non-zero initial density** to the **null distribution** under more general PDE dynamics.

Each folder contains a `main.jl` along with its supporting files (system matrices, cost, minimization/maximization logic, and visualization scripts).


---

## Problem Structure

The problem is formulated as a **min-max optimal control** problem:
- **Minimization (inner loop):** Optimal control minimizing a cost functional (using Ipopt via JuMP).
- **Maximization (outer loop):** Worst-case disturbance selected via global optimization (Differential Evolution using `BlackBoxOptim.jl`).
- **Uncertainty Set Constraint:** Projected onto a constrained set using a custom projection routine.

---

## File Overview

| File                        | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| `main.jl`                  | Entry point that runs the full optimization. Performs bilevel optimization. |
| `system.jl`                | Defines the **state** and **control matrices** from PDE discretization.     |
| `objective_function.jl`    | Contains the cost functional evaluated during inner minimization.           |
| `constraints.jl`           | Defines constraints for the inner control optimization.                     |
| `internal_minimization.jl`| Solves the **inner minimization** using Ipopt for a fixed uncertainty.       |
| `uncert_projection.jl`     | Projects uncertainty samples onto a constrained set.                        |
| `cost_matrix_generator.jl` | Generates the **Q matrix** for continuous-time control cost integration.     |
| `visual_control.jl`        | Plots the boundary control trajectory over time.                            |
| `visual_state.jl`          | Plots the spatiotemporal state density evolution.                           |

---

## Tunable Parameters

You can modify the following parameters in 'main.jl' to experiment with the problem setup:
| Parameter                                            | Description                                                                   |
| ---------------------------------------------------- | ----------------------------------------------------------------------------- |
| `N`                                                  | Number of spatial grid points for PDE discretization.                         |
| `terminal_time`                                      | Time horizon for control; determines the time interval $[0, T]$.              |
| `initial_modes`                                      | Number of spatial modes in the initial state; adjusts smoothness/oscillation. |
| `initial_state`, `terminal_state`                    | Custom initial and target state profiles for the PDE.                         |
| `lambda`                                             | Weight on the control effort in the cost functional.                          |
| `gamma`                                              | Weight on the disturbance impact in the cost functional.                      |
| `control_bounds`                                     | Lower and upper bounds on the boundary control input.                         |
| `num_params`                                         | Number of basis functions to parameterize the disturbance.                    |
| `disturbance_upper_bound`, `disturbance_lower_bound` | Range of allowed disturbance magnitude.                                       |
| `MaxSteps` (in `bboptimize`)                         | Number of iterations for the differential evolution maximization.             |
| `time_period`                                        | Period for the disturbance parameterization; affects frequency resolution.    |

---

### Install Dependencies

Ensure the following Julia packages are installed:

```julia
using Pkg
Pkg.add(["JuMP", "Ipopt", "CairoMakie", "GLMakie", "LinearAlgebra", "StaticArrays", "Dates", "Serialization", "Optim", "Random", "BlackBoxOptim", "Colors", "QuadGK"])


