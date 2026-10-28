# extract — Nektar++ output to PINN training data

`extract.py` turns the Nektar++ simulation output (VTU files) into the `.mat` files used for PINN training. It is the link between the CFD solver and the neural network: the training scripts never read Nektar++ files directly, only the `.mat` files produced here.

It is part of the thesis workflow:

```
nektar++/ (simulation + VTU files) → extract.py → cylinder_nektar_mid.mat / cylinder_nektar_far.mat → PINN training
```

## What the script does

1. **Reads the checkpoints.** All files matching `cyl_chk_*.vtu` in the same folder as the script are sorted by their number. A block of `CKPT_COUNT` snapshots is taken, starting at `CKPT_START` (to skip the start-up transient) with step `CKPT_STRIDE`. For every snapshot the fields `u`, `v` and `p` are read at all mesh nodes.
2. **Finds the free-stream velocity.** If `U_INF = None`, the script measures $U_\infty$ as the time-averaged velocity on the leftmost 2 % of the domain (the inlet strip).
3. **Non-dimensionalises the data** with one consistent convective scaling:

   $$x^* = \frac{x - x_\text{ref}}{D}, \quad y^* = \frac{y - y_\text{ref}}{D}, \quad u^* = \frac{u}{U_\infty}, \quad v^* = \frac{v}{U_\infty}, \quad p^* = \frac{p}{\rho U_\infty^2}, \quad t^* = \frac{t\,U_\infty}{D}$$

   Here $(x_\text{ref}, y_\text{ref})$ is the cylinder centre, so the cylinder sits at the origin of the new coordinates. The time scaling $t^*$ is important: if it is skipped, the time derivatives in the PDE residual have the wrong size and $\lambda_1$, $\lambda_2$ are identified wrongly.
4. **Cuts out the wake windows.** The windows are given in cylinder diameters from the cylinder centre:

   | Window | $x^*$ range | $y^*$ range |
   |---|---|---|
   | near | 1 – 8 | −2.5 – −0.5 |
   | mid | 8 – 15 | −2 – 2 |
   | far | 15 – 25 | −2 – 2 |

   By default (`TARGET_N = None`) all mesh nodes inside a window are kept. If `TARGET_N` is set to a number, that many nodes are picked at random (seed `RNG_SEED`).
5. **Saves** the `mid` and `far` windows as `cylinder_nektar_mid.mat` and `cylinder_nektar_far.mat` in the same folder.
6. **Checks the result** before training:
   - **[A] Mean velocity** — the mean of $u^*$ must be above 0.3. A value near zero means the mean flow was subtracted by mistake.
   - **[B] Convection speed** — the speed at which vortices travel downstream is estimated by cross-correlating $v^*$ at pairs of points on the centreline. It must be between 0.3 and 3 times the mean $u^*$; otherwise the time and space scales do not match.
   - **[D] Time resolution** — the shedding period is found from the spectrum of $v^*$, and there must be at least 8 snapshots per period.
   - **[C] Reynolds number** — prints the Reynolds number implied by $U_\infty$, $D$ and $\nu$, and the value of $\lambda_2 = \nu / (U_\infty D)$ that the PINN should recover.

   If all checks pass, the script prints `READY TO TRAIN`.

## Output format

| Key | Shape | Meaning |
|---|---|---|
| `x_star`, `y_star` | N × 1 | spatial coordinates |
| `t_star` | T × 1 | non-dimensional time |
| `U_star`, `V_star` | N × T | velocity components |
| `P_star` | N × T | pressure |

N is the number of mesh nodes in the window, T the number of snapshots.

## How to use it

1. **Create the VTU files** as described in [`../nektar++/README.md`](../nektar++/README.md). They must be named `cyl_chk_<number>.vtu`.
2. **Copy `extract.py` into the folder with the VTU files.** The script reads its input from, and writes its output to, its own folder.
3. **Set the configuration block** (section 1 of the script) for the Reynolds number you are extracting (see the table below).
4. **Run it:**

   ```bash
   python extract.py
   ```

5. **Read the output.** Check the line `using ... checkpoints` and make sure every check says `PASS` and the last line says `READY TO TRAIN`.
6. **Use the file for training.** Copy `cylinder_nektar_mid.mat` or `cylinder_nektar_far.mat` into the experiment folder.

## Settings for Re = 100 and Re = 500

The values in the configuration block must match the Nektar++ session file of the simulation they belong to. The time values are in the `<PARAMETERS>` section of the session XML, and the viscosity is the `Kinvis` parameter.

| Setting | Re = 100 | Re = 500 | Notes |
|---|---|---|---|
| `RE` | `100.0` | `500.0` | Only used for the printed target in check [C] |
| `NU` | `0.01` | `0.002` | Must equal `Kinvis` in the session file. This sets the expected $\lambda_2$ |
| `TIMESTEP` | `0.005` | from session file | `TimeStep` in `<PARAMETERS>` |
| `IO_CHECKSTEPS` | `20` | from session file | `IO_CheckSteps` in `<PARAMETERS>` |
| `X_REF`, `Y_REF` | `0.15`, `0.16` | cylinder centre of the Re = 500 mesh | Change only if the Re = 500 mesh is different |
| `CKPT_START` | `300` | after the transient | The wake must be fully periodic from this snapshot on |
| `CKPT_COUNT` | `500` | ≤ number of files − `CKPT_START` | See the note below |
| `D`, `RHO`, `U_INF` | `1.0`, `1.0`, `None` | same | Unchanged |
| `WINDOWS` | as above | same | Given in diameters, so they do not depend on Re |

Points that need care for Re = 500:

- **Number of checkpoints.** The script does not stop if `CKPT_START + CKPT_COUNT × CKPT_STRIDE` is larger than the number of VTU files; it simply uses fewer snapshots. For example, with 501 files and `CKPT_START = 300`, only 201 snapshots are used. Lower `CKPT_START` or `CKPT_COUNT` so they fit.
- **Training script.** The true value of $\lambda_2$ in the training script must also be changed from `0.01` to `0.002`, otherwise the reported error of $\lambda_2$ is wrong.
- **Separate folders.** The output names are always `cylinder_nektar_mid.mat` and `cylinder_nektar_far.mat`, so keep the Re = 100 and Re = 500 runs in separate folders (or rename the files after extraction) so that one does not overwrite the other.

Check [C] confirms the setting: for Re = 500 it should print an implied Re close to 500 and an expected $\lambda_2$ of `0.0020`.

## Requirements

```bash
pip install numpy scipy meshio
```
