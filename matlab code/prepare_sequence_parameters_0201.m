% =================================================================================================
% function    : prepare_sequence_parameters_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare sequence pulse and delay parameters for simulation 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : pulse choice = rect, sinc, gauss, hsecn, wurst, fermi, arb  
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_sequence_parameters_0201(sim)
    
    % ---- rf pulse sequence parameters
    for tis=1:sim.param.n_tissues
        
        % ---- simulation dwell time
        sim.seq.param.dt = sim.param.dt_us * 1e-6;                       % [us] -> [s] 

        % ---- choose sequence -> change specific sequence paramaters in 'prepare_sequence_parameters_0201.m' 
        sim.seq.param.seq_choice      = sim.param.seq_choice;            % choice of sequences: 1p, 2p, 3p, 4p, np, ir, dir, dqf, tqf, bssfp, 12p_rand
        sim.seq.param.pulse_choice    = sim.param.pulse_choice;          % [] -> see sequences below to see options

        % ---- parameters for mqf
        sim.seq.param.n_rf_pulses_mqf = sim.param.n_rf_pulses_mqf;       % [3,4] -> number of pulses for dqf and tqf sequences (if n=4, an inversion pulse is inserted between the first 2 pulses)
        sim.seq.param.dqf_ma          = sim.param.dqf_ma;                % [0,1]=[no,yes] -> apply dqf with magic angle (dqf-ma)
        
        % ---- sequence parameters
        switch sim.seq.param.seq_choice    

            % ---- 1-pulse (1p) -> [ alpha-delay ]
            case '1p'

                switch sim.seq.param.pulse_choice 
                    case 1
                        sim.seq.param.rf.type      = {'rect'};           % []
                        sim.seq.param.rf.alpha     = {90*(pi/180)};      % [rad]
                        sim.seq.param.rf.amplitude = {0};                % [hz]
                        sim.seq.param.rf.duration  = {1000};             % [us]
                        sim.seq.param.rf.bandwidth = {0};                % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {100000};           % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 2
                        sim.seq.param.rf.type      = {'hsecn'};          % []
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {1000};             % [hz]
                        sim.seq.param.rf.duration  = {10000};            % [us]
                        sim.seq.param.rf.bandwidth = {2000};             % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {5.3};              % []
                        sim.seq.param.rf.n1        = {2};                % [], exponent n for hsec^n pulse design
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 3
                        sim.seq.param.rf.type      = {'wurst'};          % []
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {1000};             % [hz]
                        sim.seq.param.rf.duration  = {10000};            % [us]
                        sim.seq.param.rf.bandwidth = {2000};             % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {20};               % [], factor n in wurst pulse design
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 4
                        sim.seq.param.rf.type      = {'gauss'};          % [], note: duration*bandwidth ~ 3.36/sigma_factor  
                        sim.seq.param.rf.alpha     = {90 *(pi/180)};     % [rad]
                        sim.seq.param.rf.duration  = {5000};             % [us]
                        sim.seq.param.rf.sigma_factor = {1};             % [], scaling factor for sigma (bandwidth) of the gaussian pulse
                        sim.seq.param.rf.sigma     = {cell2mat(sim.seq.param.rf.sigma_factor)*0.2690};   % []                        
                        sim.seq.param.rf.bandwidth = {3.36./(cell2mat(sim.seq.param.rf.sigma_factor)*cell2mat(sim.seq.param.rf.duration).*1e-6)};
                        sim.seq.param.rf.amplitude = {cell2mat(sim.seq.param.rf.alpha) ./( cell2mat(sim.seq.param.rf.sigma).*cell2mat(sim.seq.param.rf.duration).*1e-6.*pi.*sqrt(2*pi) )};
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 5
                        sim.seq.param.rf.type      = {'sinc'};           % [] 
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {410};              % [hz]
                        sim.seq.param.rf.duration  = {5000};             % [us]
                        sim.seq.param.rf.bandwidth = {0};                % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0.46};             % []
                        sim.seq.param.rf.n1        = {4};                % [], number of zero crossing on the left side of the pulse
                        sim.seq.param.rf.n2        = {4};                % [], number of zero crossing on the right side of the pulse
                        sim.seq.param.rf.phase     = {0};                % [rad] 
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 6
                        sim.seq.param.rf.type      = {'arb'};            % [] 
                        sim.seq.param.rf.duration  = {2000};             % [us]
                        sim.seq.param.rf.shape.n_points = round(sim.seq.param.rf.duration{1}./(sim.seq.param.dt*1e6));          % []
                        sim.seq.param.rf.shape.t   = (0:1:sim.seq.param.rf.shape.n_points-1)*sim.seq.param.dt*1e6;              % [us]
                        sim.seq.param.rf.shape.tau = 2*sim.seq.param.rf.shape.t/sim.seq.param.rf.duration{1}-1;                 % []
                        sim.seq.param.rf.shape.amp = exp(-(sim.seq.param.rf.shape.tau.^2)./0.15);                               % [hz]; ex: shape = gaussian
                        sim.seq.param.rf.shape.phi = pi*exp(-(sim.seq.param.rf.shape.tau.^2)./0.15);                            % [rad]; ex: shape = gaussian
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {557*sim.seq.param.rf.shape.amp};                                          % [hz]
                        sim.seq.param.rf.bandwidth = {0};                % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {sim.seq.param.rf.shape.phi};                                              % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 7
                        sim.seq.param.rf.type      = {'bir4'};           % []
                        sim.seq.param.rf.alpha     = {pi/2};             % [rad]
                        sim.seq.param.rf.amplitude = {550};              % [hz]
                        sim.seq.param.rf.duration  = {6000};             % [us], multiple of 4*dt
                        sim.seq.param.rf.bandwidth = {2000};             % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {10};               % [], default = 10
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {1.4711};           % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    case 8
                        sim.seq.param.rf.type      = {'ahp_sincos'};     % []
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {1000};             % [hz]
                        sim.seq.param.rf.duration  = {5000};             % [us]
                        sim.seq.param.rf.bandwidth = {2000};             % [hz]
                        sim.seq.param.rf.sigma     = {0};                % []
                        sim.seq.param.rf.t0        = {0};                % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]                   
                    case 9
                        sim.seq.param.rf.type      = {'fermi'};          % []  
                        sim.seq.param.rf.alpha     = {0};                % [rad]
                        sim.seq.param.rf.amplitude = {250};              % [hz]
                        sim.seq.param.rf.duration  = {2000};             % [us]
                        sim.seq.param.rf.bandwidth = {0};                % [hz]
                        sim.seq.param.rf.sigma     = {0};                % [], default = 0.1345
                        sim.seq.param.rf.t0        = {0.5};              % []
                        sim.seq.param.rf.beta      = {0};                % []
                        sim.seq.param.rf.n1        = {0};                % []
                        sim.seq.param.rf.n2        = {0};                % []
                        sim.seq.param.rf.phase     = {0};                % [rad]
                        sim.seq.param.rf.freq      = {0};                % [hz]
                        sim.seq.param.rf.delay     = {20000};            % [us]
                        sim.seq.param.adc.phase    = {0};                % [rad]
                    otherwise
                        error('choose pulses');
                end

            % ---- 2-pulses (2p) -> [ alpha01-delay-alpha02-delay ]
            case '2p'

                switch sim.seq.param.pulse_choice 
                    case 1
                        sim.seq.param.rf.type      = {    'rect'      'rect' };               % []
                        sim.seq.param.rf.alpha     = {     pi/2         pi/2 };               % [rad]
                        sim.seq.param.rf.amplitude = {        0            0 };               % [hz]
                        sim.seq.param.rf.duration  = {      500          500 };               % [us]
                        sim.seq.param.rf.bandwidth = {        0            0 };               % [hz]
                        sim.seq.param.rf.sigma     = {        0            0 };               % []
                        sim.seq.param.rf.t0        = {        0            0 };               % []
                        sim.seq.param.rf.beta      = {        0            0 };               % []
                        sim.seq.param.rf.n1        = {        0            0 };               % []
                        sim.seq.param.rf.n2        = {        0            0 };               % []
                        sim.seq.param.rf.phase     = {        0            0 ; ...
                                                           pi/2         pi/2 };               % [rad]
                        sim.seq.param.rf.freq      = {        0            0 };               % [hz]
                        sim.seq.param.rf.delay     = {    10000        30000 } ;              % [us]
                        sim.seq.param.adc.phase    = {0; 0};                                  % [rad]
                    case 2                    
                        sim.seq.param.rf.type      = {  'hsecn'       'rect' };               % []
                        sim.seq.param.rf.alpha     = {        0         pi/2 };               % [rad]
                        sim.seq.param.rf.amplitude = {      500            0 };               % [hz]
                        sim.seq.param.rf.duration  = {     5000         1000 };               % [us]
                        sim.seq.param.rf.bandwidth = {     2000            0 };               % [hz]
                        sim.seq.param.rf.sigma     = {        0            0 };               % []
                        sim.seq.param.rf.t0        = {        0            0 };               % []
                        sim.seq.param.rf.beta      = {      5.3            0 };               % []
                        sim.seq.param.rf.n1        = {        1            0 };               % []
                        sim.seq.param.rf.n2        = {        0            0 };               % []
                        sim.seq.param.rf.phase     = {        0         pi/2 ; ...
                                                           pi/4            0 ; ...
                                                         3*pi/4       3*pi/4 ; ...                  
                                                             pi            0 };               % [rad]
                        sim.seq.param.rf.freq      = {        0            0 };               % [hz]
                        sim.seq.param.rf.delay     = {    10000        20000 };               % [us]
                        sim.seq.param.adc.phase    = {0; 0; 0; 0};                            % [rad]            
                    otherwise
                        error('choose pulses');
                end

            % ---- 3-pulses (3p) -> [ alpha01-delay-alpha02-delay-alpha03-delay ]
            case '3p'

                switch sim.seq.param.pulse_choice 
                    case 1
                        sim.seq.param.rf.type      = {   'rect'       'rect'       'rect' };     % []
                        sim.seq.param.rf.alpha     = {     pi/2         pi/3         pi/4 };     % [rad]
                        sim.seq.param.rf.amplitude = {        0            0            0 };     % [hz]
                        sim.seq.param.rf.duration  = {     1000         1000         1000 };     % [us]
                        sim.seq.param.rf.bandwidth = {        0            0            0 };     % [hz]
                        sim.seq.param.rf.sigma     = {        0            0            0 };     % []
                        sim.seq.param.rf.t0        = {        0            0            0 };     % []
                        sim.seq.param.rf.beta      = {        0            0            0 };     % []
                        sim.seq.param.rf.n1        = {        0            0            0 };     % []
                        sim.seq.param.rf.n2        = {        0            0            0 };     % []
                        sim.seq.param.rf.phase     = {        0            0            0 };     % [rad]
                        sim.seq.param.rf.freq      = {        0            0            0 };     % [hz]
                        sim.seq.param.rf.delay     = {     10000        7000        40000 };     % [us]
                        sim.seq.param.adc.phase    = {0};                                        % [rad]                    
                    otherwise
                        error('choose pulses');
                end

            % ---- 4-pulses (4p) -> [ alpha01-delay-alpha02-delay-alpha03-delay-alpha04-delay ]
            case '4p'

                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.rf.type      = {   'rect'       'rect'       'rect'       'rect' };    % []
                        sim.seq.param.rf.alpha     = {     pi/2       1*pi/3       1*pi/4         pi/2 };    % [rad]
                        sim.seq.param.rf.amplitude = {        0            0            0            0 };    % [hz]
                        sim.seq.param.rf.duration  = {     1000         1000         1000         1000 };    % [us]
                        sim.seq.param.rf.bandwidth = {        0            0            0            0 };    % [hz]
                        sim.seq.param.rf.sigma     = {        0            0            0            0 };    % []
                        sim.seq.param.rf.t0        = {        0            0            0            0 };    % []
                        sim.seq.param.rf.beta      = {        0            0            0            0 };    % []
                        sim.seq.param.rf.n1        = {        0            0            0            0 };    % []
                        sim.seq.param.rf.n2        = {        0            0            0            0 };    % []
                        sim.seq.param.rf.phase     = {       pi         pi/2       2*pi/3            0 };    % [rad]
                        sim.seq.param.rf.freq      = {        0            0            0            0 };    % [hz]
                        sim.seq.param.rf.delay     = {     5000        10000         7500        20000 };    % [us]
                        sim.seq.param.adc.phase    = {0};                                                    % [rad]
                    otherwise
                        error('choose pulses');
                end
                
            % ---- 12-pulses random (12p_rand) -> [ alpha01-delay-alpha02-delay-alpha03-delay-alpha04-delay-alpha05-delay-alpha06-delay-alpha07-delay-alpha08-delay-alpha09-delay-alpha10-delay-alpha11-delay-alpha12-delay ]n 
            case '12p_rand'

                switch sim.seq.param.pulse_choice  
                    
                    case 1
                        
                        rf_n_pulses = 12;                               % []       
                        rand_1 = rand(1,rf_n_pulses);                   % []
                        rf_alpha_rad = rand_1*pi;                       % [rad]  
                        
                        rand_2 = rand(1,rf_n_pulses);                   % []
                        rf_phase_rad = rand_2*pi;                       % [rad]
                        
                        rf_duration_us = 500.*ones(1,rf_n_pulses);      % [us]
                        rf_delay_us    = 5000.*ones(1,rf_n_pulses);     % [us]
                        
                        sim.seq.param.n_repetition = 3;
                        sim.seq.param.rf.type      = repmat({            'rect'             'rect'             'rect'             'rect'             'rect'             'rect'             'rect'             'rect'             'rect'              'rect'              'rect'              'rect' }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.alpha     = repmat({   rf_alpha_rad(1)    rf_alpha_rad(2)    rf_alpha_rad(3)    rf_alpha_rad(4)    rf_alpha_rad(5)    rf_alpha_rad(6)    rf_alpha_rad(7)    rf_alpha_rad(8)    rf_alpha_rad(9)    rf_alpha_rad(10)    rf_alpha_rad(11)    rf_alpha_rad(12) }, 1, sim.seq.param.n_repetition);    % [rad]
                        sim.seq.param.rf.amplitude = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % [hz]
                        sim.seq.param.rf.duration  = repmat({ rf_duration_us(1)  rf_duration_us(2)  rf_duration_us(3)  rf_duration_us(4)  rf_duration_us(5)  rf_duration_us(6)  rf_duration_us(7)  rf_duration_us(8)  rf_duration_us(9)  rf_duration_us(10)  rf_duration_us(11)  rf_duration_us(12) }, 1, sim.seq.param.n_repetition);    % [us]
                        sim.seq.param.rf.bandwidth = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % [hz]
                        sim.seq.param.rf.sigma     = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.t0        = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.beta      = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.n1        = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.n2        = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % []
                        sim.seq.param.rf.phase     = repmat({   rf_phase_rad(1)     rf_phase_rad(2)   rf_phase_rad(3)    rf_phase_rad(4)    rf_phase_rad(5)    rf_phase_rad(6)    rf_phase_rad(7)    rf_phase_rad(8)    rf_phase_rad(9)    rf_phase_rad(10)    rf_phase_rad(11)    rf_phase_rad(12) }, 1, sim.seq.param.n_repetition);    % [rad]
                        sim.seq.param.rf.freq      = repmat({                 0                  0                  0                  0                  0                  0                  0                  0                  0                   0                   0                   0 }, 1, sim.seq.param.n_repetition);    % [hz]
                        sim.seq.param.rf.delay     = repmat({    rf_delay_us(1)     rf_delay_us(2)     rf_delay_us(3)     rf_delay_us(4)     rf_delay_us(5)     rf_delay_us(6)     rf_delay_us(7)     rf_delay_us(8)     rf_delay_us(9)     rf_delay_us(10)     rf_delay_us(11)     rf_delay_us(12) }, 1, sim.seq.param.n_repetition);    % [us]
                        sim.seq.param.adc.phase    = {0};  % [rad]
                       
                    otherwise
                        error('choose pulses');
                end
                
            % ---- n-pulses (np) -> [ alpha-tr ]n 
            case 'np'

                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.n_repetition = 3;                                                                         % []
                        sim.seq.param.rf.alpha0    = 60*(pi/180);                                                               % [rad]
                        sim.seq.param.rf.tr0       = 30000;                                                                     % [us]
                        sim.seq.param.rf.type      = repmat(                     {'rect'}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.alpha     = repmat(    {sim.seq.param.rf.alpha0}, 1, sim.seq.param.n_repetition);      % [rad]
                        sim.seq.param.rf.amplitude = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.duration  = repmat(                        {800}, 1, sim.seq.param.n_repetition);      % [us]
                        sim.seq.param.rf.bandwidth = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.sigma     = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.t0        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.beta      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n1        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n2        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.phase     = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [rad]
                        sim.seq.param.rf.freq      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.tr0       = sim.seq.param.rf.tr0 - sim.seq.param.rf.duration{1};                       % [us]
                        sim.seq.param.rf.delay     = repmat(       {sim.seq.param.rf.tr0}, 1, sim.seq.param.n_repetition);      % [us]
                        sim.seq.param.adc.phase    = {0};                                                                       % [rad]
                    case 2
                        sim.seq.param.n_repetition = 10;                                                                        % []
                        sim.seq.param.rf.tr0       = 10000;                                                                     % [us]
                        sim.seq.param.rf.type      = repmat(                    {'gauss'}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.alpha     = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [rad]
                        sim.seq.param.rf.amplitude = repmat(                        {250}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.duration  = repmat(                       {5000}, 1, sim.seq.param.n_repetition);      % [us]
                        sim.seq.param.rf.bandwidth = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.sigma     = repmat(                   {2*0.1345}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.t0        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.beta      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n1        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n2        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.phase     = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [rad]
                        sim.seq.param.rf.freq      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.delay     = repmat(       {sim.seq.param.rf.tr0}, 1, sim.seq.param.n_repetition);      % [us]
                        sim.seq.param.adc.phase    = {0};                                                                       % [rad]
                    case 3
                        sim.seq.param.n_repetition = 5;                                                                         % []
                        sim.seq.param.rf.type      = repmat(                         {'arb'}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.duration  = repmat(                          {2000}, 1, sim.seq.param.n_repetition);   % [us]
                        sim.seq.param.rf.delay     = repmat(                         {10000}, 1, sim.seq.param.n_repetition);   % [us]
                        sim.seq.param.rf.shape.n_points = round(sim.seq.param.rf.duration{1}./(sim.seq.param.dt*1e6));          % []
                        sim.seq.param.rf.shape.t   = (0:1:sim.seq.param.rf.shape.n_points-1)*sim.seq.param.dt*1e6;              % [us]
                        sim.seq.param.rf.shape.tau = 2*sim.seq.param.rf.shape.t/sim.seq.param.rf.duration{1}-1;                 % []
                        sim.seq.param.rf.shape.amp = exp(-(sim.seq.param.rf.shape.tau.^2)./0.15);                               % [hz], gaussian
                        sim.seq.param.rf.shape.phi = pi*exp(-(sim.seq.param.rf.shape.tau.^2)./0.7);                             % [rad]
                        sim.seq.param.rf.alpha     = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % [rad]
                        sim.seq.param.rf.amplitude = repmat({500*sim.seq.param.rf.shape.amp}, 1, sim.seq.param.n_repetition);   % [hz]
                        sim.seq.param.rf.bandwidth = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % [hz]
                        sim.seq.param.rf.sigma     = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.t0        = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.beta      = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.n1        = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.n2        = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % []
                        sim.seq.param.rf.phase     = repmat(    {sim.seq.param.rf.shape.phi}, 1, sim.seq.param.n_repetition);   % [rad] 
                        sim.seq.param.rf.freq      = repmat(                             {0}, 1, sim.seq.param.n_repetition);   % [hz]
                        sim.seq.param.adc.phase    = {0};                                                                       % [rad]
                    case 4
                        sim.seq.param.n_repetition = 10;                                                                        % []
                        if (tis==1)
                            sim.seq.param.rf.rand1 = 1*rand(1,sim.seq.param.n_repetition)-0;                                    % rand number in range [0 1]
                            sim.seq.param.rf.rand2 = 1*rand(1,sim.seq.param.n_repetition)-0;                                    % rand number in range [0 1] 
                            sim.seq.param.rf.rand3 = 1*rand(1,sim.seq.param.n_repetition)+0.5;                                  % rand number in range [0.25 1.25] 
                        else
                            sim.seq.param.rf.rand1 = p(1).seq.param.rf.rand1;
                            sim.seq.param.rf.rand2 = p(1).seq.param.rf.rand2;
                            sim.seq.param.rf.rand3 = p(1).seq.param.rf.rand3;
                        end
                        sim.seq.param.rf.alpha0    =  90*(pi/180);                                                              % [rad]
                        sim.seq.param.rf.phase0    =   0*(pi/180);                                                              % [rad]
                        sim.seq.param.rf.tr0       = 20000;                                                                     % [us]
                        sim.seq.param.rf.type      = repmat(                     {'rect'}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.alpha     = num2cell(sim.seq.param.rf.alpha0 .* sim.seq.param.rf.rand1);               % [rad]
                        sim.seq.param.rf.amplitude = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.duration  = repmat(                        {500}, 1, sim.seq.param.n_repetition);      % [us]
                        sim.seq.param.rf.bandwidth = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.sigma     = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.t0        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.beta      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n1        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.n2        = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % []
                        sim.seq.param.rf.phase     = num2cell(sim.seq.param.rf.phase0 .* sim.seq.param.rf.rand2);               % [rad]
                        sim.seq.param.rf.freq      = repmat(                          {0}, 1, sim.seq.param.n_repetition);      % [hz]
                        sim.seq.param.rf.tr0       = sim.seq.param.rf.tr0 - sim.seq.param.rf.duration{1};                       % [us]
                        sim.seq.param.rf.delay     = num2cell(sim.seq.param.rf.tr0 .* sim.seq.param.rf.rand3);                  % [us]
                        sim.seq.param.adc.phase    = {0};                                                                       % [rad]
                    otherwise
                        error('choose pulses');
                end
                
            % ---- inversion recovery (ir) -> [ pi-ti-alpha-acq ]n 
            case 'ir'
                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.n_repetition = 2;                                                                                                 % []
                        sim.seq.param.rf.ti        = 10000;                                                                                             % [us]
                        sim.seq.param.rf.tr0       = 50000;                                                                                             % [us]
                        sim.seq.param.rf.type      = repmat({                 'rect'                    'rect' }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.alpha     = repmat({                     pi                      pi/2 }, 1, sim.seq.param.n_repetition);       % [rad]
                        sim.seq.param.rf.amplitude = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.duration  = repmat({                    500                       500 }, 1, sim.seq.param.n_repetition);       % [us]
                        sim.seq.param.rf.bandwidth = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.sigma     = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.t0        = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.beta      = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.n1        = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.n2        = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.phase     = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [rad]
                        sim.seq.param.rf.freq      = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.tr0       = sim.seq.param.rf.tr0 - sim.seq.param.rf.ti - sim.seq.param.rf.duration{1} - sim.seq.param.rf.duration{2};      % [us]
                        sim.seq.param.rf.delay     = repmat({       sim.seq.param.rf.ti   sim.seq.param.rf.tr0 }, 1, sim.seq.param.n_repetition);       % [us]
                        sim.seq.param.adc.phase    = {0};                                                                                               % [rad]
                    case 2
                        sim.seq.param.n_repetition = 3;                                                                                                 % []
                        sim.seq.param.rf.ti        = 19000;                                                                                             % [us]
                        sim.seq.param.rf.tr0       = 140000;                                                                                            % [us]
                        sim.seq.param.rf.type      = repmat({                'wurst'                    'rect' }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.alpha     = repmat({                      0               (90/180)*pi }, 1, sim.seq.param.n_repetition);       % [rad]
                        sim.seq.param.rf.amplitude = repmat({                    240                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.duration  = repmat({                  10000                       600 }, 1, sim.seq.param.n_repetition);       % [us]
                        sim.seq.param.rf.bandwidth = repmat({                   2000                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.sigma     = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.t0        = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.beta      = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.n1        = repmat({                     20                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.n2        = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % []
                        sim.seq.param.rf.phase     = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [rad]
                        sim.seq.param.rf.freq      = repmat({                      0                         0 }, 1, sim.seq.param.n_repetition);       % [hz]
                        sim.seq.param.rf.tr0       = sim.seq.param.rf.tr0 - sim.seq.param.rf.ti - sim.seq.param.rf.duration{1} - sim.seq.param.rf.duration{2};      % [us]
                        sim.seq.param.rf.delay     = repmat({ sim.seq.param.rf.ti   sim.seq.param.rf.tr0 }, 1, sim.seq.param.n_repetition);             % [us]
                        sim.seq.param.adc.phase    = {0};                                                                                               % [rad]
                    otherwise
                        error('choose pulses');
                end

            % ---- double inversion recovery (dir) -> [ pi-ti1-pi-ti2-alpha-acq ]n 
            case 'dir'
                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.n_repetition = 3;                                                                                                                         % []
                        sim.seq.param.rf.ti1       = 20000;                                                                                                                     % [us]
                        sim.seq.param.rf.ti2       = 10000;                                                                                                                     % [us]
                        sim.seq.param.rf.tr0       = 40000;                                                                                                                     % [us]
                        sim.seq.param.rf.type      = repmat({                   'rect'                   'rect'                    'rect'}, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.alpha     = repmat({                       pi                       pi              (80/180)*pi }, 1, sim.seq.param.n_repetition);     % [rad]
                        sim.seq.param.rf.amplitude = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % [hz]
                        sim.seq.param.rf.duration  = repmat({                     1000                     1000                     1000 }, 1, sim.seq.param.n_repetition);     % [us]
                        sim.seq.param.rf.bandwidth = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % [hz]
                        sim.seq.param.rf.sigma     = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.t0        = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.beta      = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.n1        = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.n2        = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % []
                        sim.seq.param.rf.phase     = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % [rad]
                        sim.seq.param.rf.freq      = repmat({                        0                        0                        0 }, 1, sim.seq.param.n_repetition);     % [hz]
                        sim.seq.param.rf.delay     = repmat({     sim.seq.param.rf.ti1      sim.seq.param.rf.ti2    sim.seq.param.rf.tr0 }, 1, sim.seq.param.n_repetition);     % [us]
                        sim.seq.param.adc.phase    = {0};                                                                                                                       % [rad]
                    otherwise
                        error('choose pulses');
                end

            % ---- double quantum filter (dqf) w/ or w/out magic angle -> 3 pulses: [ pi/2-tau-alpha-delta-alpha-acq ] or 4 pulses: [pi/2-tau/2-pi-tau/2-alpha-delta-alpha-acq ] with alpha = pi/2 or magic angle
            case 'dqf'

                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.ma = acos(1/sqrt(3));                                                     % [rad], magic angle = 54.7347 deg = 0.9553 rad
                        sim.seq.param.rf.type = repmat({'rect'},1,sim.seq.param.n_rf_pulses_mqf);           
                        switch sim.seq.param.n_rf_pulses_mqf                
                            case 3
                                if (sim.seq.param.dqf_ma),  sim.seq.param.rf.alpha = {   pi/2    sim.seq.param.ma    sim.seq.param.ma };            % [rad]
                                else,                       sim.seq.param.rf.alpha = {   pi/2                pi/2                pi/2 };  end       % [rad]
                                sim.seq.param.rf.amplitude = {        0            0            0 };            % [hz]
                                sim.seq.param.rf.duration  = {      500          500          500 };            % [us]
                                sim.seq.param.rf.bandwidth = {        0            0            0 };            % [hz]
                                sim.seq.param.rf.sigma     = {        0            0            0 };            % []
                                sim.seq.param.rf.t0        = {        0            0            0 };            % []
                                sim.seq.param.rf.beta      = {        0            0            0 };            % []
                                sim.seq.param.rf.n1        = {        0            0            0 };            % []
                                sim.seq.param.rf.n2        = {        0            0            0 };            % []
                                sim.seq.param.rf.phase     = {        0            0            0 ; ...
                                                                 1*pi/2       1*pi/2            0 ; ...
                                                                 2*pi/2       2*pi/2            0 ; ...
                                                                 3*pi/2       3*pi/2            0 };            % [rad]
                                sim.seq.param.rf.freq      = {        0            0            0 };            % [hz]
                                sim.seq.param.rf.delay     = {     4500          200        20000 };            % [us]                    
                                sim.seq.param.adc.phase    = {0; pi; 0; pi};                                    % [rad]               
                            case 4 
                                if (sim.seq.param.dqf_ma),  sim.seq.param.rf.alpha = {   pi/2     pi    sim.seq.param.ma    sim.seq.param.ma };             % [rad]
                                else,                       sim.seq.param.rf.alpha = {   pi/2     pi                pi/2                pi/2 };  end        % [rad]
                                sim.seq.param.rf.amplitude = {        0              0                0             0 };            % [hz]
                                sim.seq.param.rf.duration  = {      500            500              500           500 };            % [us]
                                sim.seq.param.rf.bandwidth = {        0              0                0             0 };            % [hz]
                                sim.seq.param.rf.sigma     = {        0              0                0             0 };            % []
                                sim.seq.param.rf.t0        = {        0              0                0             0 };            % []
                                sim.seq.param.rf.beta      = {        0              0                0             0 };            % []
                                sim.seq.param.rf.n1        = {        0              0                0             0 };            % []
                                sim.seq.param.rf.n2        = {        0              0                0             0 };            % []
                                sim.seq.param.rf.phase     = {        0              0                0             0 ; ...
                                                                 1*pi/2         1*pi/2           1*pi/2             0 ; ...
                                                                 2*pi/2         2*pi/2           2*pi/2             0 ; ...
                                                                 3*pi/2         3*pi/2           3*pi/2             0 };            % [rad]
                                sim.seq.param.rf.freq      = {        0              0                0             0 };            % [hz]
                                sim.seq.param.rf.delay     = {     4000           4000              200         20000 };            % [us]
                                sim.seq.param.adc.phase    = {0; pi; 0; pi};                                                        % [rad]               
                            otherwise
                                error('n_pulses = 3 or 4');
                        end
                        if (sim.seq.param.dqf_ma),  sim.seq.param.seq_choice = [sim.seq.param.seq_choice 'ma - ' num2str(sim.seq.param.n_rf_pulses_mqf) ' pulses'];
                        else,                       sim.seq.param.seq_choice = [sim.seq.param.seq_choice   ' - ' num2str(sim.seq.param.n_rf_pulses_mqf) ' pulses'];  end
                    otherwise
                        error('choose pulses');
                end
                
            % ---- triple quantum filter (tqf) -> 3 pulses: [ pi/2 - tau - pi/2 - delta - pi/2 - acq ] or 4 pulses: [ pi/2 - tau/2 - pi - tau/2 - pi/2 - delta - pi/2 - acq ]
            case 'tqf'

                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.rf.type = repmat({'rect'},1,sim.seq.param.n_rf_pulses_mqf);
                        switch sim.seq.param.n_rf_pulses_mqf                
                            case 3
                                sim.seq.param.rf.alpha     = {     pi/2           pi/2         pi/2 };           % [rad]
                                sim.seq.param.rf.amplitude = {        0              0            0 };           % [hz]
                                sim.seq.param.rf.duration  = {      500            500          500 };           % [us]
                                sim.seq.param.rf.bandwidth = {        0              0            0 };           % [hz]
                                sim.seq.param.rf.sigma     = {        0              0            0 };           % []
                                sim.seq.param.rf.t0        = {        0              0            0 };           % []
                                sim.seq.param.rf.beta      = {        0              0            0 };           % []
                                sim.seq.param.rf.n1        = {        0              0            0 };           % []
                                sim.seq.param.rf.n2        = {        0              0            0 };           % []
                                sim.seq.param.rf.phase     = {   1*pi/6    1*pi/6+pi/2            0 ; ...
                                                                 3*pi/6    3*pi/6+pi/2            0 ; ...
                                                                 5*pi/6    5*pi/6+pi/2            0 ; ...
                                                                 7*pi/6    7*pi/6+pi/2            0 ; ...
                                                                 9*pi/6    9*pi/6+pi/2            0 ; ...
                                                                11*pi/6   11*pi/6+pi/2            0 };           % [rad]                   
                                sim.seq.param.rf.freq      = {        0              0            0 };           % [hz]
                                sim.seq.param.rf.delay     = {     4500            200       120000 };           % [us]                    
                                sim.seq.param.adc.phase    = {0; pi; 0; pi; 0; pi};                              % [rad]               
                            case 4 
                                sim.seq.param.rf.alpha     = {     pi/2             pi             pi/2          pi/2 };         % [rad]
                                sim.seq.param.rf.amplitude = {        0              0                0             0 };         % [hz]
                                sim.seq.param.rf.duration  = {      500            500              500           500 };         % [us]
                                sim.seq.param.rf.bandwidth = {        0              0                0             0 };         % [hz]
                                sim.seq.param.rf.sigma     = {        0              0                0             0 };         % []
                                sim.seq.param.rf.t0        = {        0              0                0             0 };         % []
                                sim.seq.param.rf.beta      = {        0              0                0             0 };         % []
                                sim.seq.param.rf.n1        = {        0              0                0             0 };         % []
                                sim.seq.param.rf.n2        = {        0              0                0             0 };         % []
                                sim.seq.param.rf.phase     = {   1*pi/6         1*pi/6      1*pi/6+pi/2             0 ; ...
                                                                 3*pi/6         3*pi/6      3*pi/6+pi/2             0 ; ...
                                                                 5*pi/6         5*pi/6      5*pi/6+pi/2             0 ; ...
                                                                 7*pi/6         7*pi/6      7*pi/6+pi/2             0 ; ...
                                                                 9*pi/6         9*pi/6      9*pi/6+pi/2             0 ; ...
                                                                11*pi/6        11*pi/6     11*pi/6+pi/2             0 };         % [rad]    
                                sim.seq.param.rf.freq      = {        0              0                0             0 };         % [hz]
                                sim.seq.param.rf.delay     = {     2000           2000              200        100000 };         % [us]
                                sim.seq.param.adc.phase    = {0; pi; 0; pi; 0; pi};                                              % [rad]              
                            otherwise
                                error('n_pulses = 3 or 4');
                        end            
                        sim.seq.param.seq_choice = [sim.seq.param.seq_choice ' - ' num2str(sim.seq.param.n_rf_pulses_mqf) ' pulses']; 
                    otherwise
                        error('choose pulses');
                end
                
            % ---- balanced steady-state free precession (bssfp) -> [ alpha/2-tr/2-{(-alpha)-tr-(+alpha)-tr}n ]  
            case 'bssfp'   

                switch sim.seq.param.pulse_choice                
                    case 1
                        sim.seq.param.n_repetition = 6;                                                                                                                                         % []
                        sim.seq.param.rf.alpha0    = pi/3;                                                                                                                                      % [rad]
                        sim.seq.param.rf.tr0       = 10000;                                                                                                                                     % [us]
                        sim.seq.param.rf.type      = [                       {'rect'}  repmat({                      'rect'                      'rect' },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.alpha     = [    {sim.seq.param.rf.alpha0/2}  repmat({    -sim.seq.param.rf.alpha0     sim.seq.param.rf.alpha0 },1,sim.seq.param.n_repetition) ];      % [rad]
                        sim.seq.param.rf.amplitude = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % [hz]
                        sim.seq.param.rf.duration  = [                          {500}  repmat({                         500                         500 },1,sim.seq.param.n_repetition) ];      % [us]
                        sim.seq.param.rf.bandwidth = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % [hz]
                        sim.seq.param.rf.sigma     = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.t0        = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.beta      = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.n1        = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.n2        = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % []
                        sim.seq.param.rf.phase     = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % [rad]                                   
                        sim.seq.param.rf.freq      = [                            {0}  repmat({                           0                           0 },1,sim.seq.param.n_repetition) ];      % [hz]
                        sim.seq.param.rf.delay     = [       {sim.seq.param.rf.tr0/2}  repmat({        sim.seq.param.rf.tr0        sim.seq.param.rf.tr0 },1,sim.seq.param.n_repetition) ];      % [us]
                        sim.seq.param.adc.phase    = {0};                                                                                                                                       % [rad]
                    otherwise
                        error('choose pulses');
                end                                

        end
    end
    
    % ---- display info
    if (sim.display.info.sequence)
        disp(' sequence parameters');
        for tis=1:sim.param.n_tissues
            disp(['   tissue ' num2str(tis)]);
            disp(['     sequence = ' sim.seq.param.seq_choice]);
            disp(['     pulse choice option = ' num2str(sim.seq.param.pulse_choice)]);
        end
    end

end
% =================================================================================================
