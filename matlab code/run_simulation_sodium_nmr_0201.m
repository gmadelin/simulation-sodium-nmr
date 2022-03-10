%% simulation of sodium nmr - v02.01

%% description
% * script       : run_simulation_sodium_nmr_0201
% * purpose      : simulation of sodium nmr sequences using density superoperator, evolution under liouvillians, and irreducible tensors decomposition
% * comment      : works for n tissues at a time with gaussian distribution of chemical shift and residual quadrupolar interaction, and chemical exchange between 2 compartments
% * parameters   : change simulation parameters in 'prepare_parameters_0201', change relaxation parameters in 'prepare_tissue_relaxation_times_0201', change rf pulse sequence parameters in 'prepare_sequence_parameters_0201' 
% * reference    : lee js et al, j chem phys 131, 174501, 2009; madelin g et al. prog nmr spectr 79, 14-47, 2014
% * note         : to run and publish an html report, type: publish('run_simulation_sodium_nmr_0201');
% * date-author  : 2018/05 - guillaume.madelin@nyumc.org 

%% clean up matlab workspace and figures 

clear; close all; clc;
disp(' '); disp(' simulation sodium nmr - v02.01');  

%% 1 - prepare simulation 
% * prepare tissues to simulate and parameters of the simulation 
% * prepare tissue relaxation times 
% * prepare redfield relaxation matrix
% * prepare frequency offset (w)
% * prepare residual quadrupolar interaction (wq)
% * prepare superoperators for sodium spin 3/2
% * prepare rf pulse sequence parameters (number and types of pulses, delays, phases, ...)
% * prepare rf pulse sequence
% * functions used: 'prepare_[name]_0201'

disp(' prepare simulation...'); tic;  

% 1.1 - prepare parameters for the simulation                                   => you can change simulation parameters in this file (tissues, choice of sequence, dwell time, offset frequency, etc.)                                     
[sim] = prepare_parameters_0201();                                               

% 1.2 - prepare simulation                                                          
[sim] = prepare_simulation_0201(sim);
      
% 1.3 - prepare tissue relaxation times + optimal tau (for mqf)                 => you can change relaxation times in this file                   
[sim] = prepare_tissue_relaxation_times_0201(sim);

% 1.4 - prepare redfield relaxation superoperator                                
[sim] = prepare_redfield_relaxation_0201(sim);

% 1.5 - prepare frequency offset (chemical shift) -> variable w                    
[sim] = prepare_w_0201(sim); 

% 1.6 - prepare residual quadrupolar interaction (rqi) -> variable wq                 
[sim] = prepare_wq_0201(sim); 
    
% 1.7 - prepare superoperators Ix, Iy, Iz (for spin 3/2) and initial vector    
[sim] = prepare_superoperator_0201(sim);
        
% 1.8 - prepare rf pulse sequence parameters (pulses and delays)                => you can change specific sequence parameters in this file                     
[sim] = prepare_sequence_parameters_0201(sim);

% 1.9 - prepare rf pulse sequence                                                  
[sim] = prepare_sequence_0201(sim);
    
%% 2 - simulation of superoperator evolution 
% * evolution of superoperators during the rf pulse sequence
% * decomposition in spherical tensors
% * functions used: 'simulate_[name]_0201' and 'calculate_[name]_0201' (within 'simulate_[name]_0201')

disp( ' run simulation...');  tic;

% 2.1 - superoperator evolution                                                     
[sim] = simulate_evolution_superoperator_0201(sim);

% 2.2 - spherical tensor decomposition                                              
[sim] = simulate_spherical_tensor_decomposition_0201(sim);

disp(['   timing [s] = ' num2str(toc)]);

%% 3 - data acquisition and spectra 
% * acquire nmr signal and apply fft (as in a spectrometer)
% * functions used: 'acquire_[name]_0201'

disp( ' acquire signal...');  

% 3.1 - acquire signal and apply fft (spectrum)                             
[sim] = acquire_signal_spectrum_0201(sim); 
    
%% 4 - figures
% * plot w and wq distributions for each tissue
% * plot evolution spherical tensors for each tissue
% * plot evolution magnetization for each tissue
% * plot acquired signal and spectrum
% * plot rf pulse sequence
% * functions used: 'plot_[name]_0201'

disp( ' plot figures...');  

% 4.1 - plot w and wq distributions for each tissue (individual figure for each tissue)
[sim] = plot_distributions_w_wq_0201(sim);

% 4.2 - plot evolution spherical tensors for each tissue (individual figure for each tissue)
[sim] = plot_evolution_spherical_tensor_0201(sim);  

% 4.3 - plot evolution magnetization for all tissues (all tissues in same figure)
[sim] = plot_evolution_magnetization_0201(sim);    

% 4.4 - plot acquired signal (T1-1) and spectrum (fft)
[sim] = plot_acquired_signal_spectrum_0201(sim);

% 4.5 - plot rf pulse sequence  
[sim] = plot_rf_pulse_sequence_0201(sim);    

%% 5 - display simulation information
% * display information on the simulation parameters and results
% * functions used: 'display_[name]_0201'

disp( ' display simulation info...');  

% 5.1 - display info  
[sim] = display_simulation_info_0201(sim);

% 5.2 - the end
disp(' voila!'); disp(' ');

%%
