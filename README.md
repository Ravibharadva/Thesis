# Physics-Informed Neural Networks for Inverse Navier–Stokes Problems

**Spatial Wake Analysis and High-Reynolds Number Extension**

M.Sc. thesis — Mathematical Modeling, Simulation, and Optimization, Universität Koblenz
Author: **Nandkishor Virani**
Supervisors: Prof. Dr. Martin Siebenborn, Dr. Christian Kahle

---

## Overview

This repository contains the code for my master's thesis. The thesis applies Physics-Informed Neural Networks (PINNs) to the inverse problem of identifying the unknown coefficients of the incompressible Navier–Stokes equations from flow data behind a circular cylinder. It builds on the cylinder-wake example of Raissi, Perdikaris and Karniadakis (2019) and extends it in two directions:

1. **Spatial wake analysis** — instead of training on the whole wake, the network is trained on smaller spatial windows (mid, far and off-centre regions) to study how the location of the measurement data affects parameter recovery.
2. **Higher Reynolds number** — the reference simulation is extended from Re = 100 to Re = 500, and the same inverse problem is solved on the new data.

This README covers the PINN part. The other two parts of the workflow have their own READMEs:

- [`nektar++/README.md`](nektar++/README.md) — mesh generation and the Nektar++ simulation
- [`extract/README.md`](extract/README.md) — turning the simulation output into `.mat` training files

## Workflow

```
Gmsh (.geo → mesh) → Nektar++ (simulation) → FieldConvert (.vtu) → extract.py (.mat) → PINN training
       └──────────── nektar++/ ────────────┘                   └ extract/ ┘    └ this README ┘
```

## Problem formulation

The 2D incompressible Navier–Stokes equations are written with two unknown parameters:

$$
u_t + \lambda_1 (u u_x + v u_y) = -p_x + \lambda_2 (u_{xx} + u_{yy})
$$

$$
v_t + \lambda_1 (u v_x + v v_y) = -p_y + \lambda_2 (v_{xx} + v_{yy})
$$

- $u, v$ — velocity components in $x$ and $y$
- $p$ — pressure (not given as training data; it is recovered by the network up to a constant)
- $\lambda_1$ — convection coefficient (true value 1.0)
- $\lambda_2$ — kinematic viscosity, equal to $1/\mathrm{Re}$ in non-dimensional form (0.01 at Re = 100, 0.002 at Re = 500)

The network outputs a stream function $\psi$ and the pressure $p$. The velocities are taken as $u = \psi_y$ and $v = -\psi_x$, so the continuity equation $u_x + v_y = 0$ is satisfied automatically. The loss combines the mismatch to the measured velocities with the residuals of the two momentum equations. $\lambda_1$ and $\lambda_2$ are trainable variables and are learned together with the network weights, first with Adam and then with L-BFGS-B.

## Experiments included in the thesis

| Experiment | Description |
|---|---|
| TF1 → TF2 port | Raissi's original TensorFlow 1 code rewritten in TensorFlow 2 and checked against the published Re = 100 results |
| Wake regions | Training on mid, far and off-centre windows of the wake instead of the full domain |
| Re = 500 extension | New Nektar++ simulation at Re = 500 and inverse identification of $\lambda_1$, $\lambda_2$ |
| Noisy data | Training with 1 % noise added to the velocity data (Adam only) |

The parameters are always initialised without hints and without anchor terms in the loss, so that the recovered values come only from the data and the physics.

**Main results:** from the mid and far wake windows, $\lambda_1$ was recovered with an error below 0.1 % and $\lambda_2$ with an error of about 2–4 %.

## Repository structure

```
.
├── nektar++/        # Gmsh .geo file, mesh and Nektar++ session file (see nektar++/README.md)
├── extract/         # extract.py: Nektar++ output → .mat training files (see extract/README.md)
├── Original/        # TF2 port of Raissi et al.'s cylinder-wake PINN (Re = 100, full domain)
├── example_mid/     # Training on the mid-wake region
├── example_far/     # Training on the far-wake region
├── adjust_new/      # Adjusted / off-centre window experiments
├── higher_Re_ex/    # Re = 500 extension
└── README.md
```

## Training data

The training scripts read `.mat` files produced by `extract/extract.py`:

| Key | Shape | Meaning |
|---|---|---|
| `x_star`, `y_star` | N × 1 | spatial coordinates (in cylinder diameters, cylinder at the origin) |
| `t_star` | T × 1 | non-dimensional time |
| `U_star`, `V_star` | N × T | velocity components |
| `P_star` | N × T | pressure (used only for evaluation, not for training) |

The original Raissi dataset (`cylinder_nektar_wake.mat`, used in `Original/`) uses the combined keys `X_star` (N × 2) and `U_star` (N × 2 × T) instead.

## Requirements

- Python 3.x
- TensorFlow 2.x
- NumPy, SciPy, Matplotlib

```bash
pip install tensorflow numpy scipy matplotlib
```

## Usage

1. Put the `.mat` file in the folder of the experiment.
2. Check the file name and the true value of $\lambda_2$ in the training script (`0.01` for Re = 100, `0.002` for Re = 500).
3. Run the script:

```bash
cd example_mid
python <training_script>.py
```

Each script trains the PINN (Adam followed by L-BFGS-B), prints the identified $\lambda_1$ and $\lambda_2$ with their relative errors, and saves plots of exact versus predicted fields.

## Reference

M. Raissi, P. Perdikaris, G. E. Karniadakis. *Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations.* Journal of Computational Physics 378 (2019), 686–707. https://doi.org/10.1016/j.jcp.2018.10.045

## Acknowledgements

I thank Prof. Dr. Martin Siebenborn and Dr. Christian Kahle for supervising this thesis. The reference simulations were produced with the open-source spectral/hp element framework Nektar++.

## Contact

Nandkishor Virani — GitHub: [@Ravibharadva](https://github.com/Ravibharadva)
