% =================================================================================================
% function    : [sim] = prepare_parameters_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare parameters for simulation of sodium nmr
% input       : -
% output      : sim (struct) 
% comment     : -   
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_parameters_0201()

    % ---- prepare structure sim
    sim = struct;
       
    % ---- get date and time of simulation  
    sim.param.date_time = prepare_get_date_time_0201();                                     % get date and time of simulation
    
    % ---- tissues and chemical exchange rates
    sim.param.tissues                           = {'csf','ic'};                             % [{'tissue name 1','tissue name 2', ...}] -> choose tissues to simulate (will be numbered with variable 'tis' in the code): 'csf', 'parenchyma', 'ec' (extracellular), 'ic' (intracellular), 'gel2%' ,'other', 'liquid', 'cartilage', 'plasma', ...
    sim.param.n_tissues                         = length(sim.param.tissues);                % [] -> number of tissues to simulate
    sim.param.fraction                          = ones(1,length(sim.param.tissues));        % [] -> fractions of each tissue in the sample = amplitude of initial magnetization of each tissue, default = 1 for all
    sim.param.chemical_exchange                 = 0;                                        % [0,1]=[no,yes] -> apply chemical exchange (works only between 2 tissues)
    sim.param.k_ex                              = [100 100];                                % [hz] -> chemical exchange rates [k12 k21], works only for 2 tissues
    
    % ---- simulation parameters
    sim.param.dt_us                             = 100;                                      % [us] -> dwell time for the simulation 
    sim.param.acq_param_choice                  = 'dt_duration';                            % ['dt_duration','bw_pixel'] -> select how you want to choose parameters to acquire signal: 'dwell time + adc duration' or 'bandwidth of acquisition per pixel + number of pixels (data points)'
    sim.param.acq_duration_us                   = 80000;                                    % [us] -> if acq_param_choice = dt_duration: data acquisition duration (similar to adc duration) 
    sim.param.acq_dt_factor                     = 2;                                        % [] -> multiplication factor for dwell time of data acquisition (dt_acquisition = factor*dt_simulation) 
    sim.param.acq_bw                            = 10000;                                    % [hz/pixel] -> if acq_param_choice = bw_pixel: bandwidth of data acquisition during adc 
    sim.param.acq_n_pixels                      = 512;                                      % [] -> if acq_param_choice = bw_pixel: number of pixels (data points) for data acquisition during adc 
    sim.param.acq_delay_us                      = 200;                                      % [us] -> data acquisition delay after rf pulse (equivalent to echo time TE)  
    sim.param.acq_add_gaussian_noise            = 1;                                        % [0,1]=[no,yes] -> add white gaussian noise
    sim.param.acq_noise_power_dbw               = -80;                                      % [dbw] -> select noise power when gaussian noise is added: -23dB <-> mean(noise)=0.05
    sim.param.acq_after_each_pulse              = 0;                                        % [0,1]=[no,yes] -> 1: acquire data after each pulse (for multipulse sequences such as 'bssfp' or 'np'); 0: acquire data only after the last pulse

    % ---- choose sequence -> change specific sequence parameters in 'prepare_sequence_parameters_0201' 
    sim.param.seq_choice                        = '1p';                                     % ['name'] -> choice of sequences: 1p, 2p, 3p, 4p, np, ir, dir, dqf, tqf, bssfp, 12p_rand
    sim.param.pulse_choice                      = 1;                                        % [] -> see sequences in 'prepare_sequence_parameters_0201' to see options
        
    % ---- parameters for mqf
    sim.param.n_rf_pulses_mqf                   = 3;                                        % [3,4] -> number of pulses for dqf and tqf sequences (if n=4, an inversion pulse is inserted between the first 2 pulses)
    sim.param.dqf_ma                            = 0;                                        % [0,1]=[no,yes] -> apply dqf with magic angle (dqf-ma)
    
    % ---- frequency offset for each tissue -> w (note about notation: f is in Hz, w=2*pi*f will be in rad/s)
    for tis=1:sim.param.n_tissues
        % -- choose frequency offset parameters for each tissue -> w                            
        sim.tissue(tis).w.index                 =     0*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [0,1] -> 0 = single w (per tissue), 1 = gaussian distribution (per tissue) [you can add more tissues: + value*(tis==4), etc.]
        % -- one single w                                                                        
        sim.tissue(tis).w.f.f0                  =     0*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [hz] -> single offset frequency f0 in hz per tissue [you can add more tissues: + value*(tis==4), etc.]
        % -- gaussian distribution of w
        sim.tissue(tis).w.f.min                 =     0*(tis==1) +  -200*(tis==2) +     0*(tis==3);     % [hz] -> min f                                    
        sim.tissue(tis).w.f.df                  =    25*(tis==1) +    50*(tis==2) +     0*(tis==3);     % [hz] -> step df
        sim.tissue(tis).w.f.max                 =   300*(tis==1) +   200*(tis==2) +     0*(tis==3);     % [hz] -> max f                                    
        sim.tissue(tis).w.f.mean                =   200*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [hz] -> mean of distribution of f                                   
        sim.tissue(tis).w.f.std                 =   250*(tis==1) +   250*(tis==2) +     0*(tis==3);     % [hz] -> std of distribution of f                                    
    end
    
    % ---- residual quadrupolar interaction (rqi) for each tissue -> wq (note about notation: fq is in hz, wq=2*pi*fq will be in rad/s)
    for tis=1:sim.param.n_tissues
        % -- choose rqi parameters -> wq                                            
        sim.tissue(tis).wq.index                =     0*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [0,1] -> 0 = single wq (per tissue), 1 = gaussian distribution (per tissue) [you can add more tissues: + value*(tis==4), etc.]
        % -- one single wq                                                          
        sim.tissue(tis).wq.fq.fq0               =     0*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [hz] -> single rqi frequency (wq) per tissue [you can add more tissues: + value*(tis==4), etc.]                                    
        % -- gaussian distribution of wq
        sim.tissue(tis).wq.fq.min               =     0*(tis==1) +     0*(tis==2) +     0*(tis==3);     % [hz] -> min fq                                    
        sim.tissue(tis).wq.fq.dfq               =   100*(tis==1) +   100*(tis==2) +     0*(tis==3);     % [hz] -> step dfq                                    
        sim.tissue(tis).wq.fq.max               =   500*(tis==1) +   300*(tis==2) +     0*(tis==3);     % [hz] -> max fq                                   
        sim.tissue(tis).wq.fq.mean              =   300*(tis==1) +   100*(tis==2) +     0*(tis==3);     % [hz] -> mean of distribution of fq                                   
        sim.tissue(tis).wq.fq.std               =   250*(tis==1) +   250*(tis==2) +     0*(tis==3);     % [hz] -> std of distribution of fq
    end
    
    % ---- display info
    sim.display.info.list_tissue                = 0;                                % [0,1]=[no,yes] -> in prepare_simulation_0201
    sim.display.info.relaxation_times           = 0;                                % [0,1]=[no,yes] -> in prepare_tissue_relaxation_times_0201
    sim.display.info.spectral_densities         = 0;                                % [0,1]=[no,yes] -> in prepare_redfield_relaxation_0201
    sim.display.info.frequency_offset           = 0;                                % [0,1]=[no,yes] -> in prepare_w_0201
    sim.display.info.rqi                        = 0;                                % [0,1]=[no,yes] -> in prepare_wq_0201
    sim.display.info.sequence                   = 0;                                % [0,1]=[no,yes] -> in prepare_sequence_parameters_0201
    sim.display.info.acquisition                = 0;                                % [0,1]=[no,yes] -> in plot_acquired_signal_spectrum_0201
    sim.display.info.simulation                 = 1;                                % [0,1]=[no,yes] -> in plot_acquired_signal_spectrum_0201
    
    % ---- display figures
    sim.display.figure.sequence_check           = 0;                                % [0,1]=[no,yes] -> in prepare_sequence_0201
    sim.display.figure.distribution_w_wq        = 1;                                % [0,1]=[no,yes] -> in plot_distributions_w_wq_0201
    sim.display.figure.evolution_sph_tensor     = 1;                                % [0,1]=[no,yes] -> in plot_distributions_w_wq_0201
    sim.display.figure.evolution_magnetization  = 1;                                % [0,1]=[no,yes] -> in plot_evolution_magnetization_0201
    sim.display.figure.signal_spectra_all_acq   = 1;                                % [0,1]=[no,yes] -> in plot_acquired_signal_spectrum_0201
    sim.display.figure.noise_spectra_all_acq    = 0;                                % [0,1]=[no,yes] -> in plot_acquired_signal_spectrum_0201
    sim.display.figure.signal_sequence          = 1;
    sim.display.figure.sequence                 = 1;                                % [0,1]=[no,yes] -> in plot_rf_pulse_sequence_0201
    sim.display.figure.tissues                  = 'all';                            % ['all', or integer] -> select tissues to plot, in plot_evolution_spherical_tensor_0201, plot_evolution_magnetization_0201, and plot_acquired_signal_spectrum_0201
  
    % ---- set figure position and size
    scrsz = get(0,'screensize');
    set(0,'defaultfigureposition',[0.15*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(3) 0.7*scrsz(4)]); 
    set(0,'defaultfigurecolor','w');  
    warning('off','images:initsize:adjustingmag');
            
end
% =================================================================================================
