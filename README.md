# PendulumClocks

Simulations of the quantum optomechanical pendulum clock.

This repository accompanies our paper (arXiv): https://arxiv.org/abs/2506.10666  
It is intended for peer physicists who may want to (i) reproduce the results, and/or (ii) adapt the codes to investigate new directions.

---

## Table of Contents
- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Quick start (typical workflow)](#quick-start-typical-workflow)
- [Reproducing specific figures / analyses](#reproducing-specific-figures--analyses)
- [Editing parameters](#editing-parameters)
- [Outputs and data locations](#outputs-and-data-locations)
- [Citation](#citation)
- [License](#license)

---

## Overview

The simulations are organized into folders by physical scenario (single-atom / multi-atom, zero temperature / finite temperature, etc.).  
Each folder (and in some cases subfolders) is designed to be runnable largely independently.

**Main interface:** MATLAB `run_*.m` scripts.  
**Engine:** Julia (`.jl`) is used as the main engine, but users typically do *not* run Julia directly unless they want to modify engine-level details.

---

## Repository structure

The main folders are:

- `Single_atom_zero_T/`  
  Single-atom simulations at zero temperature. Reproduces Figs. 5–7 of the paper.
  - `Single_atom_zero_T/Phase_space_limit_cycle/`  
    Focuses on Fig. 4 (phase-space limit cycle) and related analysis.

- `Single_atom_Thermodynamics/`  
  Single-atom simulations at finite temperatures. Reproduces Fig. 8, the stared data in Fig. 12, and heat currents in the appendix.
  
- `Two_atoms_zero_T/` 2-atoms simulations at zero temperature. Reproduces Fig. 10 of the paper.
- `Two_atoms_Thermodynamics/` 2-atoms simulations at finite temperature. Reproduces a stared data in Fig.12 

- `Multiple_atoms/`  
  Multi-atom scenario; produces Figs. 11–12 of the paper.


Within each scenario folder you will typically find:
- `run_factorisation_and_resample.m` (main “driver” script for that folder)
- `Factorisation.m` (core model setup/propagation parameters)
- `Detector_Filter_saturation.m` (tick detection / detector current processing)
- analysis helpers such as `Allan*.m`, `size_resample.m`, and figure/plot scripts.

---

## Requirements

You need both installed on the same machine:

- **MATLAB** (main interface and analysis)
- **Julia** (engine; typically called indirectly)

Notes:
- On a first run you may encounter missing-path errors if MATLAB cannot “see” helper functions. Usually this is resolved by running from within the relevant folder, or adding required folders to the MATLAB path.

---

## Quick start (typical workflow)

### 1) Choose a target folder
Pick the folder corresponding to the figure/physics scenario you want to reproduce (e.g. `Single_atom_zero_T/`, `Multiple_atoms/`, etc.), and set MATLAB’s **Current Folder** to that directory.

### 2) Generate trajectories and saved data (required per folder)
Open and run [the first section of]:
- `run_factorisation_and_resample.m` (in that folder)

This script is typically organized in sections:
- **First section:** Monte-Carlo simulation of clock trajectories and saving results (usually `.mat`) into the *same directory*.  
  This step can be expensive; for the paper’s runs it can take from ~day to ~week depending on parameters (number of trajectories/iterations, number of atoms, total trajectory time, etc.).
- **Subsequent sections:** Load the previously generated data and perform analysis/plotting to produce figures.

Practical workflow:
- Run the first (data-generation) section once.
- Later, you can re-open MATLAB and continue from the analysis sections without re-generating data, as long as the saved `.mat` files are still present.

---

## Reproducing specific figures / analyses

After you have generated the trajectory data by running the **first section** of `run_factorisation_and_resample.m` in the relevant folder:

- **Autocorrelation of ticks (Fig. 14):**  
  Run `run_autocorrelations.m` in `Single_atom_zero_T/Phase_space_limit_cycle/`.

- **Tick-statistics histograms vs temperature (Fig. 8):**  
  Run `run_histogram_overlay.m` inside `Single_atom_Thermodynamics/`.

- **Tick-statistics histograms and detector current vs atom number M (Fig. 11 middle/bottom):**  
  Run `run_histogram_overlay.m` inside `Multiple_atoms/`.

(Other figures are often produced by later sections of `run_factorisation_and_resample.m`, but not all.)

---

## Editing parameters

### `run_factorisation_and_resample.m`
This is the main place to start editing parameters.

Typical editable parameters include bath temperatures (hot/cold), and “sweep” controls such as:
- `imax`: how many trajectories (iterations) to simulate at fixed parameters
- `iTmax`: number of different temperatures to sweep (folder-dependent)
- `iMmax`: number of atom counts to sweep (folder-dependent)

This script typically runs:
- **unconditional evolution** (no unraveling), then
- **conditional (unravelled) evolution** which is stochastic.

Common switch:
- `ur = 0`: no unraveling (unconditional)
- `ur = 1`: unraveling (conditional, stochastic)

### `Factorisation.m`
Usually called by other scripts (not run directly).  
Contains important system parameters. Changing defaults may remove the limit cycle regime.

Important parameter:
- `tmax`: set via an `if` clause depending on `ur`.
  - For good statistics/convergence with `ur == 1`, `tmax` ideally corresponds to **thousands of mechanical periods**.
  - For faster tests, much smaller values (e.g. ~100 periods) can already show meaningful behavior.

### `Detector_Filter_saturation.m`
Called by other scripts (not run directly).  
You can modify detector/tick extraction settings, e.g.:
- tick threshold current `I*` via `threshholddown`
- smoothing/control parameter `n` (integer > 3): larger means slower but smoother detector-current processing; results should not materially change beyond providing a smoother current trace.

### `Allan.m`, `Allan_from_scratch.m`, `size_resample.m`
Helpers called by other scripts.  
In typical usage you do not need to modify them; there are no key user-facing parameters.

---

## Outputs and data locations

- Outputs are typically `.mat` files.
- They are usually saved in the **current folder** (the folder you are running in), i.e. alongside the code for that scenario.
- Analysis sections/scripts will load these `.mat` files and generate plots.

Tip: keep separate folders for separate parameter sweeps/runs (or rename/move `.mat` outputs) to avoid mixing datasets.

---

## Citation

If you use this repository, please cite:

- https://arxiv.org/abs/2506.10666

(You may paste the official BibTeX entry from the arXiv page into this README if desired.)

---

## License

MIT License (see `LICENSE`).
