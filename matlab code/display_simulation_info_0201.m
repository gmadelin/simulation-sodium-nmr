% =================================================================================================
% function    : display_simulation_info_0201
% -------------------------------------------------------------------------------------------------
% purpose     : display info from simulation
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -   
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = display_simulation_info_0201(sim)

    % ---- display info
    if (sim.display.info.simulation)

        disp( ' simulation info'); 
        
        if (sim.param.chemical_exchange), disp(['   number of tissues = ' num2str(sim.param.n_tissues-1)]);
        else disp(['   number of tissues = ' num2str(sim.param.n_tissues)]);  end
        
        disp(['   chemical exchange = ' num2str(sim.param.chemical_exchange)]);
        
        for tis=1:sim.param.n_tissues
            
            disp(['   tissue ' num2str(tis) ' = ' sim.tissue(tis).tissue.name]);
            disp(['     fraction  = ' num2str(sim.tissue(tis).tissue.fraction)]);
            
            if (sim.param.chemical_exchange),  disp(['     k_ex [hz] = ' num2str(sim.param.k_ex)]);  end
            
            if ~( sim.param.chemical_exchange && tis==3)
                
                disp( '     relaxation times [ms]');
                disp(['       t1short = ' num2str(1000*sim.tissue(tis).relax.t1short) ', t1long = ' num2str(1000*sim.tissue(tis).relax.t1long)]);
                disp(['       t2short = ' num2str(1000*sim.tissue(tis).relax.t2short) ', t2long = ' num2str(1000*sim.tissue(tis).relax.t2long)]);
                disp( '     spectral densities [hz]');
                disp(['       j0 = ' num2str(sim.tissue(tis).relax.j(1)) ', j1 = ' num2str(sim.tissue(tis).relax.j(2)) ', j2 = ' num2str(sim.tissue(tis).relax.j(3))]);
                disp( '     frequency offset (w) parameters');
                disp(['       gaussian distribution = ' num2str(sim.tissue(tis).w.index)]);
                
                if (sim.tissue(tis).w.index),  disp(['       distribution min/df/max/mean/std [hz] = ' num2str(sim.tissue(tis).w.f.min) '/' num2str(sim.tissue(tis).w.f.df) '/' ...
                                          num2str(sim.tissue(tis).w.f.max) '/' num2str(sim.tissue(tis).w.f.mean) '/' num2str(sim.tissue(tis).w.f.std)]);
                else disp(['       frequency offset [hz] = ' num2str(sim.tissue(tis).w.f.f0)]); 
                end
                
                disp( '     residual quadrupolar interaction (wq) parameters');
                disp(['       gaussian distribution = ' num2str(sim.tissue(tis).wq.index)]);
                
                if (sim.tissue(tis).wq.index),  disp(['       distribution min/df/max/mean/std [hz] = ' num2str(sim.tissue(tis).wq.fq.min) '/' num2str(sim.tissue(tis).wq.fq.dfq) '/' ...
                                         num2str(sim.tissue(tis).wq.fq.max) '/' num2str(sim.tissue(tis).wq.fq.mean) '/' num2str(sim.tissue(tis).wq.fq.std)]);
                else disp(['       rqi frequency [hz]    = ' num2str(sim.tissue(tis).wq.fq.fq0)]);
                end
                
            end    
         
        end
        
        disp( '   rf pulse sequence parameters');        
        disp(['     simulation dwell time [us] = ' num2str(sim.seq.dt*1e6)]);
        disp(['     sequence                   = ' sim.seq.param.seq_choice]);
        disp(['     pulse choice option        = ' num2str(sim.seq.param.pulse_choice)]);    
        
        disp( '   data acquisition parameters');
        disp(['     acq after each pulse       = ' num2str(sim.param.acq_after_each_pulse)]);
        disp(['     add gaussian noise         = ' num2str(sim.param.acq_add_gaussian_noise)]);
       
        if (sim.param.acq_add_gaussian_noise),  disp(['     noise power          [dbw] = ' num2str(sim.tissue(tis).acq.noise_power)]); end
        
        disp(['     choice of acq parameters   = ' sim.tissue(tis).acq.acq_param_choice]);
        
        for tis=1:sim.param.n_tissues
            
            disp(['     tissue ' num2str(tis) ' = ' sim.tissue(tis).tissue.name]);
            
            if (sim.param.acq_after_each_pulse)
                for i=1:sim.tissue(tis).acq.n_spectra
                    disp(['       spectrum # ' num2str(i)]);
                    disp(['         acq dwell time       [us] = ' num2str(sim.tissue(tis).acq.dt{i}.*1e6)]);
                    disp(['         acq bandwidth  [hz/pixel] = ' num2str(sim.tissue(tis).acq.bw{i})]);
                    disp(['         acq duration         [us] = ' num2str(sim.tissue(tis).acq.duration{i}.*1e6)]);      
                    disp(['         acq pixel number          = ' num2str(sim.tissue(tis).acq.n_pixels{i})]);
                    disp(['         delay before acq     [us] = ' num2str((sim.tissue(tis).acq.delay_before_acq{i}).*1e6)]);
                    disp( '         snr info: i=n_spectra, j=n_dt, k=n_duration, m=n_delay_before_acq');
                    for j=1:sim.tissue(tis).acq.nn_dt{i}
                        for k=1:sim.tissue(tis).acq.nn_duration{i}
                            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}           
                                disp(['         snr abs(signal) / abs(fft) [ijkm: ' num2str(i) ' ' num2str(j) ' ' num2str(k) ' ' num2str(m) '] = ' ...
                                       num2str(roundn((sim.tissue(tis).acq.snr_signal{i,j,k,m}),-1)) '  /  ' num2str(roundn((sim.tissue(tis).acq.snr_signal_fft{i,j,k,m}),-1))]);
                            end
                        end
                    end
                end
            else                
                for i=1:sim.tissue(tis).acq.n_spectra   % n_spectra = 1 (after last pulse)
                    disp(['       spectrum # ' num2str(i)]);
                    disp(['         acq dwell time       [us] = ' num2str(sim.tissue(tis).acq.dt{end}.*1e6)]);
                    disp(['         acq bandwidth  [hz/pixel] = ' num2str(sim.tissue(tis).acq.bw{end})]);
                    disp(['         acq duration         [us] = ' num2str(sim.tissue(tis).acq.duration{end}.*1e6)]);      
                    disp(['         acq pixel number          = ' num2str(sim.tissue(tis).acq.n_pixels{end})]);
                    disp(['         delay before acq     [us] = ' num2str((sim.tissue(tis).acq.delay_before_acq{end}).*1e6)]);
                    disp( '         snr info: i=n_spectra, j=n_dt, k=n_duration, m=n_delay_before_acq');
                    for j=1:sim.tissue(tis).acq.nn_dt{end}
                        for k=1:sim.tissue(tis).acq.nn_duration{end}
                            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{end}           
                                disp(['         snr abs(signal) / abs(fft) [ijkm: ' num2str(i) ' ' num2str(j) ' ' num2str(k) ' ' num2str(m) '] = ' ...
                                       num2str(roundn((sim.tissue(tis).acq.snr_signal{i,j,k,m}),-1)) '  /  ' num2str(roundn((sim.tissue(tis).acq.snr_signal_fft{i,j,k,m}),-1))]);
                            end
                        end
                    end
                end
            end
            
        end
         
    end

end
% =================================================================================================
