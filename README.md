# simulation_sodium_nmr

## Matlab code for simulation of sodium NMR - v02.01. 

Simulation of sodium NMR sequences using density operators, superoperators, evolution under Liouvillians, and irreducible tensors decomposition.

This simulation works for multiple tissues at a time, with option of Gaussian distribution of chemical shifts and residual quadrupolar interactions, and with option of chemical exchange (between 2 tissues only). Any RF pulse sequence can be designed.

See inside the code for details and comments.

## Run the simulation
Type 'run_simulation_sodium_nmr_0201' in Matlab to run simulation.

- Change simulation parameters in 'prepare_parameters_0201'. 
- Change relaxation parameters in 'prepare_tissue_relaxation_times_0201'.
- Change RF pulse sequence parameters in 'prepare_sequence_parameters_0201'.

To run and publish an html report, type in Matlab: publish('run_simulation_sodium_nmr_0201');

## References:
- Madelin G et al. Prog NMR Spectr 79, 14-47, 2014 (sodium NMR and MRI)
- Lee JS et al, J Chem Phys 131, 174501, 2009 (density operators and superoperators for simulation of spin 3/2 dynamics)
- Gilles A et al. Sci Rep 7, 17435, 2017 (example of use of the simulation code)
