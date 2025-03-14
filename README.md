# PendulumClocks

Simulations of the Quantum Optomechanical Pendulum Clock.

## Table of Contents
- [Introduction](#introduction)
- [Structure](#structure)
- [Setup](#setup)
- [Running the Simulations](#running-the-simulations)
- [Editing Parameters](#editing-parameters)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This repository contains MATLAB simulations for the Quantum Optomechanical Pendulum Clock. The simulations are organized into different subfolders, each dealing with specific aspects of the project.

## Structure
Provide a brief overview of each subfolder here. For example:

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
### Single_atom_zero_T/phase_space_limit_cycle

1. Navigate to the `src` directory:
    ```sh
    cd src
    ```
2. Run the main simulation script:
    ```matlab
    main_simulation.m
    ```

## Editing Parameters
Guide users on where and how to edit parameters. For example:

1. Open the `config.m` file in the `src` directory.
2. Edit the parameters as needed. For example:
    ```matlab
    parameter1 = 10;
    parameter2 = 5;
    ```

## Contributing
Explain how others can contribute to the project. For example:

1. Fork the repository.
2. Create a new branch for your changes.
3. Submit a pull request.

## License
State the license under which the project is distributed.






# Pendulom_clocks
This code accompanies our draft entitled "". 
It characterises the performance of a pendulum clock. Namely, it simulates its evolution (i) conditional and (ii) unconditional, under the factorisation assumption. In the conditional case, this is a monte-carlo simulations. 
The codes can illustrate the phase-space evolution of the mechanical oscillator, the populations of the atome, or the number operator of the cavity. All in both cases (i) and (ii). For (i) we also can depict the tick times, and its statistics (histogram). We will calculate the accuracy, and the Allan Variance. We then use a filter that represents the death time of the detectors. Using this filter, our statistics can improve.

# How to run
You can run the code from run_factorisation_and_resample.m in this code, imax decides how many trajectories you simulate.
The code Factorisation.mlx contains the core of the simulation. Set the parameter values, including the couplings etc, and the time-length of the evolution here. The parameter dt, cannot be too big, the simulations may break down.
The data that comes out of Factorisation is too big. size_resample shrinks the size. Good for saving and data analysis. The code also saves the relevant information for each trajectory (after resizing).

To get the data analysis (clock performance) one can run the remaining parts of run_factorisation_and_resample.m namely, using a detector, or not using a detector. If you run the full code from the command line, these will be included automatically.

# Conditional or unconditional? 
In Factorisation.mlx set ur=1 if you want to unravel (conditional evolution). set ur=0 if not.

# If you want to see the dynamics
The code run_factorisation_and_resample.m only depicts the clock performance. You cannot see the dynamics of different observables. For that, you can manually plot them (also, some codes are commented within the run_factorisation_and_resample.m code commented as load and analyse, check there too).
