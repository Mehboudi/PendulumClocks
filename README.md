# PendulumClocks

Simulations of the Quantum Optomechanical Pendulum Clock.

## Table of Contents
- [Introduction](#introduction)
- [Structure](#structure)
- [Setup](#setup)
- [Running the Simulations](#running-the-simulations)
- [Editing Parameters](#editing-parameters)

## Introduction
This repository contains MATLAB simulations for the Quantum Optomechanical Pendulum Clock. The simulations are organized into different subfolders, each dealing with specific aspects of the project.

## Structure

- `Single_atom_zero_T`: Contains the code for simulation of a single atom case. It reproduces figs. 5-8 of the paper.
- `Single_atom_zero_T/Phase_space_limit_cycle`:  Contains the code for simulation of a single atom case. It focuses on fig. 4 of the paper.
- - `Single_atom_Thermodynamics`: Contains the code for simulation of a single atom case at finite temperatures. It reproduces figs. 9 and the heat currents in the appendix.
- `Multiple_atoms`: Focuses on the multi-atom scenario, produces Figs 9-11 of the paper.

## Setup
It is assumed that you have some basic knowledge about MATLAB. In particular, on your first try, you may get errors such as file x was not found, unknown variable, etc, which you should think about how to resolve. But in principle everything should work fine.

1. Clone the repository:
    ```sh
    git clone https://github.com/Mehboudi/PendulumClocks.git
    ```
2. Navigate to the project directory:
    ```sh
    cd PendulumClocks
    ```
3. Ensure you have MATLAB installed on your system.

## Running the Simulations
You always have to simulate some data first (should be done for each folder/subfolder separately. You cannot count on the data of one folder for another folder; they focus on different aspect of the project). This is done by opening the code called run_factorisation_and_resample.m in the respective folder. You can either run the full code in one go or run it section by section. I explain the section by section style: start by running the first section of the code, which (monte-carlo) simulates the trajectories of the clock and saves them in the same directory. This can take a while depending on the parameters, namely the number of iterations (trajectories), number of atoms, total trajectory time etc; that's why the data is saved automatically so that if you close matlab, you won't lose the effort. For the graphs in the paper, this is of the order of a day-week. But for tests you can run it for shorter times (see below).
Then you can move to the next section of the code and run the section, and move forward section by section. These sections, depending on the folder, will load the data that you generated before and analyze it. They will produce the results as graphs. 
If you want to run the code again tomorrow, you don't have to run the first section again, as your data is already saved. Just move on to the analysis sections. 

Most data analysis and figures will be produced within the subsequent sections of the code run_factorisation_and_resample.m. However, not all.

Assuming that you have already run the first section of run_factorisation_and_resample.m in the folder that you are interested in:
If you want to get the autocorrelation of the ticks (Fig. 7), run run_autocorrelations.m in the phase_space_limit_cycle folder.
If you want the tick statistics histograms for various values of temperature (fig. 9), run run_histogram_overlay.m inside the single_atom_Thermodynamics.
If you want the ick statistics histograms and detector current for various values of M (middle and bottom panels of Fig. 10), run run_histogram_overlay.m inside the Multiple_atoms.

## Editing Parameters

### run_factorisation_and_resample.m

In this code you can edit some of the system parameters, namely temperatures of the cold/hot baths. But the main editable oparameters are imax iMmax and iTmax (depending on the folder) which decide how many times you simulate a trajectory with the same parameters (imax), or for how many different temperatures (iTmax) or for how many atoms (from 1 to iMmax). 

This code, first simulates the dynamics without unraveling (i.e., unconditional evolution). Then it simulates the unravelled dynamics (conditional) which is stochastic. The parameter ur=0 means no unravelling, ur=1 means unravelling. 

### factorisation.m

This code is always called from other codes, not directly.
You can edit main system parameters here. These are subtle, changing the default values may lead to a situation with no limit cycle. 

Most important parameter is tmax. It's value is determined in an if clause, either depending on ur=0 or ur=1. For having a good statistics and convergence of the results tmax for ur==1 is better to be thousands of periods of the mechanical oscillator. However, even with 100 periods you can see some meaningful results. 

### Detector_Filter_saturation

This code is also called from other codes, you don't have to open it. However, you can in principle modify some of the parameters, including the thrshhold for deciding when a tick happens (I^*); that would be threshholddown.
The parameter n should be an integer above 3. The bigger, the slower the code. But the result doesn't change. It only makes sense to increase it if you want to have a smoother code for the detector current.

### Allan and Allan_from_scratch and size_resample

These code is also called from other codes.
You don't need to touch these, there are really no parameters to modify


