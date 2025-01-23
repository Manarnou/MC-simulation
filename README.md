# MC-simulation code

This project is an ongoing development of an MC simulation recreating the Molecular Dynamics (MD) simulation presented in the paper "Correlations in the Motion of Atoms in Liquid Argon" by A. Rahman(1964). 
Main code: New-simulation-MC-relaxation.f90, it is an update to the old-simulation.f90 

## System and Simulation Parameters

- Number of molecules (N): 864 ( initially N=100 is used to for the current testing process)
- Box size: Cubic box with periodic boundary conditions
- Density: 1.374 g/cm³
- Temperature: 94.4 K
- Potential: Lennard-Jones potential for interatomic interactions

## Current Status

The code currently outputs only the total energy of each new configuration. Future work will include additional analyses, such as pair-correlation functions, to replicate the results from Rahman's paper.
