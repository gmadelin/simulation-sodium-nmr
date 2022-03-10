=================================================
	TITLE
=================================================
Matlab code for simulation of sodium NMR - v02.01


=================================================
	PURPOSE
=================================================
Simulation of sodium NMR sequences using density operators, superoperators, evolution under Liouvillians, and irreducible tensors decomposition


=================================================
	COMMENTS
=================================================
This simulation works for multiple tissues at a time, with option of Gaussian distribution of chemical shifts and residual quadrupolar interactions, and with option of chemical exchange (between 2 tissues only)
See inside the code for details and comments


=================================================
	MATLAB SCRIPT TO RUN
=================================================
run_simulation_sodium_nmr_0201


=================================================
	PARAMETERS                                  
=================================================
Change simulation parameters in 'prepare_parameters_0201'
Change relaxation parameters in 'prepare_tissue_relaxation_times_0201'
Change RF pulse sequence parameters in 'prepare_sequence_parameters_0201' 


=================================================
	REFERENCES
=================================================
Madelin G et al. Prog NMR Spectr 79, 14-47, 2014 (sodium NMR and MRI)
Lee JS et al, J Chem Phys 131, 174501, 2009 (density operators and superoperators for simulation of spin 3/2 dynamics)
Gilles A et al. Sci Rep 7, 17435, 2017 (example of use of the simulation code)

=================================================
	NOTE
=================================================
To run and publish an html report, type in Matlab: publish('run_simulation_sodium_nmr_0201');
If you want to run the simulation with white gaussian noise (wgn) added to the signal acquisition (sim.param.acq_add_gaussian_noise = 1 in 'prepare_parameters_0201'): to use 'wgn', the following product must be licensed, installed, and enabled: Communications System Toolbox


=================================================
	DATE-AUTHOR
=================================================
2018/05 - guillaume.madelin@nyumc.org, in collaboration with alexej.jerschow@nyu.edu and jaeseung.lee@nyumc.org


=================================================
	MATLAB VERSION
=================================================
Matlab 2016a and later (code was not tested on previous versions)


=================================================
	FUNDING
=================================================
This code was developed with the support of NIH grants R01NS097494, R21CA213169 and P41EB017183