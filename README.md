# Meta-LED Optimization

**Designing a nanostructured LED surface that emits light in a chosen direction and polarization — not in all directions at once.**

![Python](https://img.shields.io/badge/python-3.12-blue.svg)
![RCWA](https://img.shields.io/badge/Maxwell%20solver-S4%20%2F%20Lumerical-informational)
![Optimization](https://img.shields.io/badge/optimization-Ax%20%2F%20BoTorch-orange)

---

## Overview

A standard LED radiates spontaneous emission into (almost) every angle, which wastes light
that never reaches a detector, fiber, or viewer at the angle you actually care about. This
project designs a **meta-LED**: a GaN-based LED whose top surface is etched with a
subwavelength pattern of parallel **nanoribbons**, each with a few rectangular **notches** cut
into its sides, so that emission is redirected into a target angle (e.g. straight up, or
steered off-normal) and/or a target linear polarization.

The ribbon/notch geometry — positions, widths, and depths, all in the tens-to-hundreds-of-nm
range — is treated as a continuous design space and searched over automatically. Every
candidate design is fabrication-aware: geometric constraints enforce realistic minimum
trench/mesa widths (~50-100 nm) so the optimizer can't propose something impossible to etch.

## How it works

**1. Directivity via reciprocity.** RCWA solvers like S4 can't simulate a spontaneously
emitting dipole directly, so the figure of merit is computed via the reciprocity theorem
instead: a plane wave is launched into the structure from a given in-plane wavevector
`(kx, ky)` and polarization, and the resulting field intensity at the quantum-well location
is proportional to how much power a dipole there would radiate into that outgoing direction.
Sweeping `(kx, ky)` over the collection numerical aperture and normalizing the intensity at
the target direction against the average over the full aperture gives the **directivity, D**.

**2. A dual-polarization figure of merit.** Directivity is computed separately for s- and
p-polarization (`Ds`, `Dp`) and combined as:

```
D = Ds + Dp - |Ds - Dp|
```

which rewards designs that are directional in *both* polarizations simultaneously, rather than
rewarding a design that's extremely directional in one polarization while staying diffuse in
the other.

**3. Bayesian optimization over the geometry.** [Ax](https://ax.dev/) drives the search with a
two-phase strategy: an initial Sobol (quasi-random, space-filling) sweep to explore the design
space broadly, followed by BoTorch (Gaussian-process-based Bayesian optimization) trials that
exploit the accumulated data to propose increasingly promising geometries. Each candidate is
built and solved in [S4](https://github.com/victorliu/S4) (a Fourier Modal Method / RCWA Maxwell
solver), and its directivity is fed back to Ax as the trial outcome to maximize.

Each evaluation runs in its own subprocess purely to force memory cleanup between S4 calls
(S4's Python bindings otherwise leak memory over long optimization runs) — it's a workaround,
not a parallelization strategy.

## Repository structure

| Path | Contents |
| --- | --- |
| `S4/` | The active RCWA implementation: structure builder (`makeS4structure.py`), reciprocity/figure-of-merit code (`directivity.py`), the Ax optimization driver (`dOpt.py`), design-variable and constraint definitions (`vars_and_constraints.py`), and supporting plotting/analysis scripts. |
| `LumRCWA/` | A parallel implementation of the same optimization loop using Lumerical's commercial RCWA/FDTD solver instead of S4, as an alternate physics backend. |
| `results/` | Dated output folders from optimization runs (trial history, convergence plots, emission-pattern heatmaps) — including validation against an actual fabricated device. Local run output, not tracked in git. |
| `with-ribbon-placement-error/` | Earlier/broader project history: first working reciprocity pipeline, fabrication-tolerance sensitivity studies, and off-normal beam-steering exploration. Also local run output, not tracked in git. |

## Setup

Install the pure-Python dependencies:

```bash
pip install -r requirements.txt
```

**S4 must be installed separately** — it's a compiled Maxwell solver (the Stanford Stratified
Structure Solver), not a normal PyPI package. Build it (or a maintained fork) from source and
make sure its Python bindings are importable as `import S4` before running anything here.

## Usage

There is no command-line entry point — `S4/dOpt.py` is written as an interactive script (in the
Spyder/Jupyter "run cell by cell" style), not a script you invoke with `python dOpt.py`. The two
functions you actually want to call, from an interactive session with `S4/` on your path:

```python
import dOpt

# Launch a fresh Bayesian-optimization run over the geometry defined in
# vars_and_constraints.py:
data, client = dOpt.D_opt(dim, var_2d, params_2d, constraints,
                          (training_count, learning_count))
```

```python
import directivity, makeS4structure

# Visualize the emission pattern of one specific design:
directivity.plot_at_QWz_2d(
    makeS4structure.struct_2d(params_2d), params_2d, QW_depth, correct_apod=True
)
```

`plotting_optimization_results.py` and `plot_total_emission_across_spectral_FWHM.py` are
similarly interactive: they expect a `data.npy` trial history (or an already-computed
intensity array) to already be loaded in the workspace.
