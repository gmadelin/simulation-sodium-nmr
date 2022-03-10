% =================================================================================================
% function    : acquire_signal_spectrum_0201
% -------------------------------------------------------------------------------------------------
% purpose     : acquire signal and apply fft (spectrum) 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = acquire_signal_spectrum_0201(sim)
    
    % ---- data acquisition parameters
    for tis=1:sim.param.n_tissues
        sim.tissue(tis).acq.acq_after_each_pulse = sim.param.acq_after_each_pulse;      % [0,1]
        sim.tissue(tis).acq.add_gaussian_noise   = sim.param.acq_add_gaussian_noise;    % [0,1]
        sim.tissue(tis).acq.noise_power          = sim.param.acq_noise_power_dbw;       % [dbw]
        sim.tissue(tis).acq.acq_param_choice     = sim.param.acq_param_choice;          % ['bw_pixel','dt_duration']
    end
    
    % ---- white gaussian noise (same for all tissues, amplitude calculated on signal from first tissue)
    if (sim.param.acq_add_gaussian_noise)
        sim.tissue(1).acq.signal_t11 = sqrt(2)*sim.tissue(1).sphtens.tensor(6,:);
        white_gaussian_noise = wgn(1,length(sim.tissue(1).acq.signal_t11),sim.tissue(1).acq.noise_power,'dbw','complex');
    end
    
    % ---- data acquisition and fft
    for tis=1:sim.param.n_tissues

        % ---- data acquisition parameters: select "dt + acq duration" or "bandwidth/pixel + pixel number" for data acquisition (here: 1D pixel = data point)
        switch sim.param.acq_param_choice
            case 'dt_duration'
                for i=1:sim.seq.n_pulse
                    sim.tissue(tis).acq.dt{i}                  = sim.param.acq_dt_factor * sim.seq.dt;                                      % [s], acquisition dwell time = multiple of simulation dt -> j
                    sim.tissue(tis).acq.duration{i}            = sim.param.acq_duration_us * 1e-6;                                          % [s], acquisition duration -> k
                    sim.tissue(tis).acq.delay_before_acq{i}    = sim.param.acq_delay_us * 1e-6;                                             % [s], delay before acquisition -> m
                    sim.tissue(tis).acq.bw{i}                  = 1./sim.tissue(tis).acq.dt{i};                                              % [hz/pixel], acquisition bandwidth 
                    sim.tissue(tis).acq.n_pixels{i}            = round(sim.tissue(tis).acq.duration{i}./sim.tissue(tis).acq.dt{i});         % [], number of points (or pixels) of acquisitions
                    sim.tissue(tis).acq.nn_dt{i}               = length(sim.tissue(tis).acq.dt{i});                                         % []
                    sim.tissue(tis).acq.nn_duration{i}         = length(sim.tissue(tis).acq.duration{i});                                   % [] 
                    sim.tissue(tis).acq.nn_delay_before_acq{i} = length(sim.tissue(tis).acq.delay_before_acq{i});                           % []                    
                end        
            case 'bw_pixel'                                                                                                                 % note: works ok when all variables = single value (not vector)     
                for i=1:sim.seq.n_pulse
                    sim.tissue(tis).acq.bw{i}                  = sim.param.acq_bw;                                                          % [hz/pixel], acquisition bandwidth 
                    sim.tissue(tis).acq.n_pixels{i}            = sim.param.acq_n_pixels;                                                    % [], number of points (or pixels) of acquisitions
                    sim.tissue(tis).acq.delay_before_acq{i}    = sim.param.acq_delay_us * 1e-6;                                             % [s], delay before acquisition -> m
                    sim.tissue(tis).acq.dt{i}                  = roundn(1./sim.tissue(tis).acq.bw{i},-6);                                   % [s], acquisition dwell time = multiple of simulation dt -> j
                    sim.tissue(tis).acq.duration{i}            = roundn(sim.tissue(tis).acq.n_pixels{i}.*sim.tissue(tis).acq.dt{i},-6);     % [s], acquisition duration -> k
                    sim.tissue(tis).acq.nn_dt{i}               = length(sim.tissue(tis).acq.dt{i});                                         % []
                    sim.tissue(tis).acq.nn_duration{i}         = length(sim.tissue(tis).acq.duration{i});                                   % [] 
                    sim.tissue(tis).acq.nn_delay_before_acq{i} = length(sim.tissue(tis).acq.delay_before_acq{i});                           % [] 
                end
        end
        
        % ---- check acquisition parameters
        for i=1:sim.seq.n_pulse
            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}
                if ( sim.tissue(tis).acq.delay_before_acq{i}(m) < sim.seq.dt ),  sim.tissue(tis).acq.delay_before_acq{i}(m) = 0;  end
            end
            for j=1:sim.tissue(tis).acq.nn_dt{i} 
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    if ( sim.tissue(tis).acq.duration{i}(k) < sim.tissue(tis).acq.dt{i}(j) ),  sim.tissue(tis).acq.duration{i}(k) = sim.tissue(tis).acq.dt{i}(j);  end
                end
            end
            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}
                if ( sim.tissue(tis).acq.delay_before_acq{i}(m) >= sim.seq.delay{i}.duration ),  sim.tissue(tis).acq.delay_before_acq{i}(m) = 0;  end
            end
            for k=1:sim.tissue(tis).acq.nn_duration{i}
                for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}
                    if ( (sim.tissue(tis).acq.duration{i}(k)+sim.tissue(tis).acq.delay_before_acq{i}(m)) > sim.seq.delay{i}.duration )
                        sim.tissue(tis).acq.duration{i}(k) = sim.seq.delay{i}.duration - sim.tissue(tis).acq.delay_before_acq{i}(m);
                    end
                end
            end
        end

        % ----re-actualise number of pixels
        for i=1:sim.seq.n_pulse
            sim.tissue(tis).acq.n_pixels{i} = round(sim.tissue(tis).acq.duration{i}./sim.tissue(tis).acq.dt{i});                % [], number of points (or pixels) of acquisitions
        end        

        % ---- acquisition delay and duration
        for i=1:sim.seq.n_pulse
            sim.tissue(tis).acq.n_last_delay(i) = sim.seq.delay{end}.n_points;                                                  % []
            for j=1:sim.tissue(tis).acq.nn_dt{i}
                sim.tissue(tis).acq.n_dt_step(i,j) = round(sim.tissue(tis).acq.dt{i}(j)./sim.seq.dt);                           % [], acquisition dwell time step for reading signal
            end
            for k=1:sim.tissue(tis).acq.nn_duration{i}
                sim.tissue(tis).acq.n_duration(i,k) = round(sim.tissue(tis).acq.duration{i}(k)./sim.seq.dt);                    % []
            end
            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}
                sim.tissue(tis).acq.n_delay_before_acq(i,m) = round(sim.tissue(tis).acq.delay_before_acq{i}(m)./sim.seq.dt);    % []
            end
        end

        % ---- acquisition of signal after each pulse
        
        % ---- signal = sqrt(2)*T1-1
        sim.tissue(tis).acq.signal_t11 = sqrt(2)*sim.tissue(tis).sphtens.tensor(6,:); 
        
        % ---- signal + noise
        if (sim.tissue(tis).acq.add_gaussian_noise)
            % -- add white_gaussian_noise  
            sim.tissue(tis).acq.noise0  = white_gaussian_noise;  % same noise for all tissues
            sim.tissue(tis).acq.signal0 = sim.tissue(tis).acq.signal_t11 + sim.tissue(tis).acq.noise0;
        else
            % -- no noise
            sim.tissue(tis).acq.noise0  = 0.*sim.tissue(tis).acq.signal_t11;
            sim.tissue(tis).acq.signal0 = 1.*sim.tissue(tis).acq.signal_t11;
        end
        
        % ---- 1st pulse
        sim.tissue(tis).acq.n_delay_no_acq(1) = 1 + sim.seq.pulse{1}.n_points;    % number of delay points before first acquisition 
        for j=1:sim.tissue(tis).acq.nn_dt{1}
            for k=1:sim.tissue(tis).acq.nn_duration{1}
                for m=1:sim.tissue(tis).acq.nn_delay_before_acq{1}  
                    sim.tissue(tis).acq.signal1{1,j,k,m} = sim.tissue(tis).acq.signal0( ( sim.tissue(tis).acq.n_delay_no_acq(1) + sim.tissue(tis).acq.n_delay_before_acq(1,m) + 1 ) : sim.tissue(tis).acq.n_dt_step(1,j) : ...
                                                           ( sim.tissue(tis).acq.n_delay_no_acq(1) + sim.tissue(tis).acq.n_delay_before_acq(1,m) + sim.tissue(tis).acq.n_duration(1,k) ) );
                    sim.tissue(tis).acq.noise1{1,j,k,m}  = sim.tissue(tis).acq.noise0(  ( sim.tissue(tis).acq.n_delay_no_acq(1) + sim.tissue(tis).acq.n_delay_before_acq(1,m) + 1 ) : sim.tissue(tis).acq.n_dt_step(1,j) : ...
                                                           ( sim.tissue(tis).acq.n_delay_no_acq(1) + sim.tissue(tis).acq.n_delay_before_acq(1,m) + sim.tissue(tis).acq.n_duration(1,k) ) );
                    sim.tissue(tis).acq.n_acq{1,j,k,m}   = length(sim.tissue(tis).acq.signal1{1,j,k,m});
                end
            end
        end
        
        % ---- next pulses
        if (sim.seq.n_pulse>1)
            for i=2:sim.seq.n_pulse
                sim.tissue(tis).acq.n_delay_no_acq(i) = sim.tissue(tis).acq.n_delay_no_acq(i-1) + sim.seq.delay{i-1}.n_points + sim.seq.pulse{i}.n_points;
                for j=1:sim.tissue(tis).acq.nn_dt{i}
                    for k=1:sim.tissue(tis).acq.nn_duration{i}
                        for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}                 
                            sim.tissue(tis).acq.signal1{i,j,k,m} = sim.tissue(tis).acq.signal0( ( sim.tissue(tis).acq.n_delay_no_acq(i) + sim.tissue(tis).acq.n_delay_before_acq(i,m) + 1 ) : sim.tissue(tis).acq.n_dt_step(i,j) : ...
                                                                   ( sim.tissue(tis).acq.n_delay_no_acq(i) + sim.tissue(tis).acq.n_delay_before_acq(i,m) + sim.tissue(tis).acq.n_duration(i,k) ) );  
                            sim.tissue(tis).acq.noise1{i,j,k,m}  = sim.tissue(tis).acq.noise0(  ( sim.tissue(tis).acq.n_delay_no_acq(i) + sim.tissue(tis).acq.n_delay_before_acq(i,m) + 1 ) : sim.tissue(tis).acq.n_dt_step(i,j) : ...
                                                                   ( sim.tissue(tis).acq.n_delay_no_acq(i) + sim.tissue(tis).acq.n_delay_before_acq(i,m) + sim.tissue(tis).acq.n_duration(i,k) ) );  
                            sim.tissue(tis).acq.n_acq{i,j,k,m}   = length(sim.tissue(tis).acq.signal1{i,j,k,m});
                        end
                    end
                end
            end
        end
                    
        % ---- choose signals to process
        if sim.tissue(tis).acq.acq_after_each_pulse             
            % -- data acquired after each pulse in the sequence
            sim.tissue(tis).acq.signal = sim.tissue(tis).acq.signal1; 
            sim.tissue(tis).acq.noise  = sim.tissue(tis).acq.noise1; 
        else
            % -- acquire data only after the last pulse of the sequence
            for j=1:sim.tissue(tis).acq.nn_dt{end}
                for k=1:sim.tissue(tis).acq.nn_duration{end}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{end}  
                        sim.tissue(tis).acq.signal{1,j,k,m} = sim.tissue(tis).acq.signal1{end,j,k,m};
                        sim.tissue(tis).acq.noise{1,j,k,m}  = sim.tissue(tis).acq.noise1{end,j,k,m};
                    end
                end
            end
        end
        
        % ---- number of spectra to calculate and plot
        if (sim.tissue(tis).acq.acq_after_each_pulse)
            sim.tissue(tis).acq.n_spectra = sim.seq.n_pulse;
        else
            sim.tissue(tis).acq.n_spectra = 1;
        end
        
        % ---- get max of each acquired signal
        for i=1:sim.tissue(tis).acq.n_spectra
            for j=1:sim.tissue(tis).acq.nn_dt{i}
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i} 
                        sim.tissue(tis).acq.signal_abs_max{i,j,k,m} = max(abs(sim.tissue(tis).acq.signal{i,j,k,m}));                        
                        if (1), ind_abs_max = find( abs(sim.tissue(tis).acq.signal{i,j,k,m}) == sim.tissue(tis).acq.signal_abs_max{i,j,k,m} );
                        else    ind_abs_max = ind_abs_max(1);   end   %#ok    % test: get 1st point of acquired signal as max signal      
                        sim.tissue(tis).acq.signal_max{i,j,k,m} = sim.tissue(tis).acq.signal{i,j,k,m}(ind_abs_max);  %#ok
                    end
                end
            end
        end        
        
        % ---- fft of signal = spectrum 
        for i=1:sim.tissue(tis).acq.n_spectra
            for j=1:sim.tissue(tis).acq.nn_dt{i}
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}  
                        sim.tissue(tis).acq.fs(i,j,k,m)         = 1/sim.tissue(tis).acq.dt{i}(j);                                                               % sampling frequency: fs = 1/dw
                        sim.tissue(tis).acq.n_fft(i,j,k,m)      = length(sim.tissue(tis).acq.signal{i,j,k,m});                                                  % number of points for fft 
                        sim.tissue(tis).acq.freq{i,j,k,m}       = -(sim.tissue(tis).acq.fs(i,j,k,m)-1)/2:sim.tissue(tis).acq.fs(i,j,k,m)/sim.tissue(tis).acq.n_fft(i,j,k,m):sim.tissue(tis).acq.fs(i,j,k,m)/2;   % frequency range of the spectrum : -(fs-1)/2 to fs/2
                        sim.tissue(tis).acq.signal_fft{i,j,k,m} = fftshift(fft(sim.tissue(tis).acq.signal{i,j,k,m})./sim.tissue(tis).acq.n_fft(i,j,k,m));       % fft
                        sim.tissue(tis).acq.noise_fft{i,j,k,m}  = fftshift(fft(sim.tissue(tis).acq.noise{i,j,k,m})./sim.tissue(tis).acq.n_fft(i,j,k,m));        % fft
                    end
                end
            end
        end
        
        % ---- snr
        for i=1:sim.tissue(tis).acq.n_spectra
            for j=1:sim.tissue(tis).acq.nn_dt{i}
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}  
                        sim.tissue(tis).acq.snr_signal{i,j,k,m}     = max( abs(sim.tissue(tis).acq.signal{i,j,k,m}) ) ./ std( abs(sim.tissue(tis).acq.noise{i,j,k,m}) );
                        sim.tissue(tis).acq.snr_signal_fft{i,j,k,m} = max( abs(sim.tissue(tis).acq.signal_fft{i,j,k,m}) ) ./ std( abs(sim.tissue(tis).acq.noise_fft{i,j,k,m}) );
                    end
                end
            end
        end
        
        % ---- timing
        for i=1:sim.tissue(tis).acq.n_spectra
            for j=1:sim.tissue(tis).acq.nn_dt{i}
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}  
                        sim.tissue(tis).acq.n_points(i,j,k,m) = length(sim.tissue(tis).acq.signal{i,j,k,m});                                                % []
                        sim.tissue(tis).acq.t{i,j,k,m}        = (0:1:sim.tissue(tis).acq.n_points(i,j,k,m)-1).*sim.tissue(tis).acq.dt{i}(j) .* 1e6;         % [us]
                    end
                end
            end
        end

        % ---- acquisition timing (for display in figures) [use sim.seq.dt from simulation, not sim.tissue(tis).acq.dt from acquisition]
        for k=1:sim.tissue(tis).acq.nn_duration{1}
            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{1}
                sim.tissue(tis).acq.acq{k,m} = zeros(1,sim.seq.n_t);
            end
        end
        if (sim.tissue(tis).acq.acq_after_each_pulse)
            for i=1:sim.tissue(tis).acq.n_spectra
                for k=1:sim.tissue(tis).acq.nn_duration{i}
                    for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}  
                        sim.tissue(tis).acq.acq{k,m}( sim.tissue(tis).acq.n_delay_no_acq(i)+sim.tissue(tis).acq.n_delay_before_acq(i,m)+1 : sim.tissue(tis).acq.n_delay_no_acq(i)+sim.tissue(tis).acq.n_delay_before_acq(i,m)+sim.tissue(tis).acq.n_duration(i,k) ) = 1;
                    end
                end
            end
        else
            for k=1:sim.tissue(tis).acq.nn_duration{1}
                for m=1:sim.tissue(tis).acq.nn_delay_before_acq{1}  
                    sim.tissue(tis).acq.acq{k,m}( sim.tissue(tis).acq.n_delay_no_acq(end)+sim.tissue(tis).acq.n_delay_before_acq(end,m)+1 : sim.tissue(tis).acq.n_delay_no_acq(end)+sim.tissue(tis).acq.n_delay_before_acq(end,m)+sim.tissue(tis).acq.n_duration(end,k) ) = 1;
                end
            end
        end
        
    end

end
% =================================================================================================
