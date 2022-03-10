% =================================================================================================
% function    : plot_acquired_signal_spectrum_0201
% -------------------------------------------------------------------------------------------------
% purpose     : plot acquired signal and fft spectrum 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -   
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = plot_acquired_signal_spectrum_0201(sim)

    % ---- select tissues to plot
    switch sim.display.figure.tissues
        case 'all',  tis1 = 1;                           tisn = sim.param.n_tissues;
        otherwise,   tis1 = sim.display.figure.tissues;  tisn = tis1;  
    end
    
    % ---- plot selected figures
    for tis=tis1:tisn

        % ---- figure 1: plot signal and fft spectra of all acquisitions with different dwell times, delays before acq, and acq durations (n figures)
        if (sim.display.figure.signal_spectra_all_acq)
            
            for i=1:sim.tissue(tis).acq.n_spectra
                for j=1:sim.tissue(tis).acq.nn_dt{i}
                    for k=1:sim.tissue(tis).acq.nn_duration{i}
                        for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i} 
                            
                            if (sim.param.acq_after_each_pulse)
                                name_fig1 = ['signal and spectrum - ' sim.tissue(tis).tissue.name ' - acq ' num2str(i) ' - acq dt [us] = ' num2str(sim.tissue(tis).acq.dt{i}(j)*1e6) ...
                                                ' - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{i}(k)*1e6) ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{i}(m)*1e6)]; 
                            else
                                name_fig1 = ['signal and spectrum - ' sim.tissue(tis).tissue.name ' - acq ' num2str(i) ' - acq dt [us] = ' num2str(sim.tissue(tis).acq.dt{end}(j)*1e6) ...
                                                ' - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{end}(k)*1e6) ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{end}(m)*1e6)]; 
                            end
                                        
                            if (sim.param.chemical_exchange)
                                name_fig1 = [name_fig1 ' - k_{ex} [hz] = ' num2str(sim.param.k_ex(1)) ' / ' num2str(sim.param.k_ex(2))];  %#ok
                            end
                            
                            fig1 = figure('name',name_fig1);
                                    subplot(2,2,1)
                                        hold on;
                                            plot(sim.tissue(tis).acq.t{i,j,k,m},real(sim.tissue(tis).acq.signal{i,j,k,m}),'-b','linewidth',1);
                                            plot(sim.tissue(tis).acq.t{i,j,k,m},imag(sim.tissue(tis).acq.signal{i,j,k,m}),'-r','linewidth',1);
                                            plot([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                        hold off; box on;
                                        if (sim.tissue(tis).acq.t{i,j,k,m}(end)>sim.tissue(tis).acq.t{i,j,k,m}(1))
                                            xlim([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)]); 
                                        end
                                        xlabel('time [\mus]'); ylabel('signal (real+imag)');
                                        title([sim.tissue(tis).tissue.name ' - signal T_{1-1} - real and imag'],'fontweight','bold');
                                    subplot(2,2,2)
                                        hold on;
                                            plot(sim.tissue(tis).acq.freq{i,j,k,m},real(sim.tissue(tis).acq.signal_fft{i,j,k,m}),'-b','linewidth',1);
                                            plot(sim.tissue(tis).acq.freq{i,j,k,m},imag(sim.tissue(tis).acq.signal_fft{i,j,k,m}),'-r','linewidth',1);
                                            plot([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                        hold off; box on;
                                        if sim.tissue(tis).acq.freq{i,j,k,m}(end)>sim.tissue(tis).acq.freq{i,j,k,m}(1)
                                            xlim([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)]);
                                        end
                                        xlabel('frequency [hz]'); ylabel('signal fft (real+imag)');
                                        title([sim.tissue(tis).tissue.name ' - fft of signal T_{1-1} - real and imag'],'fontweight','bold');
                                    subplot(2,2,3)
                                        hold on;
                                            plot(sim.tissue(tis).acq.t{i,j,k,m},abs(sim.tissue(tis).acq.signal{i,j,k,m}),'-k','linewidth',2);
                                            plot([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                        hold off; box on;
                                        if sim.tissue(tis).acq.t{i,j,k,m}(end)>sim.tissue(tis).acq.t{i,j,k,m}(1)
                                            xlim([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)]); 
                                        end
                                        xlabel('time [\mus]');  ylabel('signal (abs)');
                                        title([sim.tissue(tis).tissue.name ' - signal T_{1-1} - abs'],'fontweight','bold');
                                    subplot(2,2,4)
                                        hold on;
                                            plot(sim.tissue(tis).acq.freq{i,j,k,m},abs(sim.tissue(tis).acq.signal_fft{i,j,k,m}),'-k','linewidth',2);
                                            plot([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                        hold off; box on;
                                        if sim.tissue(tis).acq.freq{i,j,k,m}(end)>sim.tissue(tis).acq.freq{i,j,k,m}(1)
                                            xlim([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)]); 
                                        end
                                        xlabel('frequency [hz]'); ylabel('signal fft (abs)');
                                        title([sim.tissue(tis).tissue.name ' - fft of signal T_{1-1} - abs'],'fontweight','bold');
                                    annotation(fig1,'textbox',[0.06 0.95 0.90 0.05],'string',name_fig1,'horizontalalignment','center', ...
                                                    'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % figure title

                        end
                    end
                end
            end
            
        end
        
        % ---- figure 2: plot noise and fft spectra of all acquisitions of all acquisitions with different dwell times, delays before acq, and acq durations (n figures)
        if (sim.display.figure.noise_spectra_all_acq && sim.tissue(tis).acq.add_gaussian_noise)
            
            for i=1:sim.tissue(tis).acq.n_spectra
                for j=1:sim.tissue(tis).acq.nn_dt{i}
                    for k=1:sim.tissue(tis).acq.nn_duration{i}
                        for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i} 
                                                        
                            if (sim.param.acq_after_each_pulse)
                                name_fig2 = ['noise and spectrum - ' sim.tissue(tis).tissue.name ' - acq ' num2str(i) ' - acq dt [us] = ' num2str(sim.tissue(tis).acq.dt{i}(j)*1e6) ...
                                                ' - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{i}(k)*1e6) ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{i}(m)*1e6)]; 
                            else
                                name_fig2 = ['noise and spectrum - ' sim.tissue(tis).tissue.name ' - acq ' num2str(i) ' - acq dt [us] = ' num2str(sim.tissue(tis).acq.dt{end}(j)*1e6) ...
                                                ' - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{end}(k)*1e6) ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{end}(m)*1e6)]; 
                            end
                                                                    
                            fig2 = figure('name',name_fig2);
                                subplot(2,2,1)
                                    hold on;
                                        plot(sim.tissue(tis).acq.t{i,j,k,m},real(sim.tissue(tis).acq.noise{i,j,k,m}),'-b','linewidth',1);
                                        plot(sim.tissue(tis).acq.t{i,j,k,m},imag(sim.tissue(tis).acq.noise{i,j,k,m}),'-r','linewidth',1);
                                        plot([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                    hold off; box on;
                                    if sim.tissue(tis).acq.t{i,j,k,m}(end)>sim.tissue(tis).acq.t{i,j,k,m}(1)
                                        xlim([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)]); 
                                    end
                                    xlabel('time [\mus]'); ylabel('noise (real+imag)');
                                    title([sim.tissue(tis).tissue.name ' - noise - real and imag'],'fontweight','bold');
                                subplot(2,2,2)
                                    hold on;
                                        plot(sim.tissue(tis).acq.freq{i,j,k,m},real(sim.tissue(tis).acq.noise_fft{i,j,k,m}),'-b','linewidth',1);
                                        plot(sim.tissue(tis).acq.freq{i,j,k,m},imag(sim.tissue(tis).acq.noise_fft{i,j,k,m}),'-r','linewidth',1);
                                        plot([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                    hold off; box on;
                                    if sim.tissue(tis).acq.freq{i,j,k,m}(end)>sim.tissue(tis).acq.freq{i,j,k,m}(1)
                                        xlim([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)]);
                                    end
                                    xlabel('frequency [hz]'); ylabel('noise fft (real+imag)');
                                    title([sim.tissue(tis).tissue.name ' - fft of noise - real and imag'],'fontweight','bold');
                                subplot(2,2,3)
                                    hold on;
                                        plot(sim.tissue(tis).acq.t{i,j,k,m},abs(sim.tissue(tis).acq.noise{i,j,k,m}),'-k','linewidth',2);
                                        plot([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                    hold off; box on;
                                    if sim.tissue(tis).acq.t{i,j,k,m}(end)>sim.tissue(tis).acq.t{i,j,k,m}(1)
                                        xlim([sim.tissue(tis).acq.t{i,j,k,m}(1) sim.tissue(tis).acq.t{i,j,k,m}(end)]); 
                                    end
                                    xlabel('time [\mus]'); ylabel('noise (abs)');
                                    title([sim.tissue(tis).tissue.name ' - noise - abs'],'fontweight','bold');
                                subplot(2,2,4)
                                    hold on;
                                        plot(sim.tissue(tis).acq.freq{i,j,k,m},abs(sim.tissue(tis).acq.noise_fft{i,j,k,m}),'-k','linewidth',2);
                                        plot([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)],[0 0],':k','linewidth',1);
                                    hold off; box on;
                                    if sim.tissue(tis).acq.freq{i,j,k,m}(end)>sim.tissue(tis).acq.freq{i,j,k,m}(1)
                                        xlim([sim.tissue(tis).acq.freq{i,j,k,m}(1) sim.tissue(tis).acq.freq{i,j,k,m}(end)]); 
                                    end
                                    xlabel('frequency [hz]'); ylabel('noise fft (abs)');
                                    title([sim.tissue(tis).tissue.name ' - fft of noise - abs'],'fontweight','bold');                              
                                annotation(fig2,'textbox',[0.06 0.95 0.90 0.05],'string',name_fig2,'horizontalalignment','center', ...
                                                    'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % figure title
                                       
                        end
                    end
                end
            end 
            
        end        
                
        % ---- figure 3: plot signal acquired after each pulse (1 figure) [with acq duration and delay before acq #1]
        if (sim.display.figure.signal_sequence)
            
            color_acq = [0.75 0.75 0.75];      % color for acquisition time
            tt = sim.seq.t * 1e3;              % [us]
            
            for k=1:sim.tissue(tis).acq.nn_duration{1}
                for m=1:sim.tissue(tis).acq.nn_delay_before_acq{1} 
                    
                    name_fig3 = ['signal acquired during sequence - ' sim.tissue(tis).tissue.name ...
                                  ' - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{end}(k)*1e6) ...
                                  ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{end}(m)*1e6) ];
                                      
                    fig3 = figure('name',name_fig3);
                        subplot(3,1,1)
                            hold on;
                                plot(tt,real(sim.tissue(tis).acq.signal0).*sim.tissue(tis).acq.acq{k,m},'b','linewidth',1);
                                plot(tt,imag(sim.tissue(tis).acq.signal0).*sim.tissue(tis).acq.acq{k,m},'r','linewidth',1);
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);                
                                plot(tt, abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':k');  % pulses
                                plot(tt,-abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':k'); % pulses
                                plot(tt,sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':','color',color_acq);               % acq
                                plot(tt,-sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':','color',color_acq);              % acq
                            hold off; box on;
                            xlim([tt(1) tt(end)]); xlabel('time [ms]'); ylabel('signal (real+imag)');
                            title([sim.tissue(tis).tissue.name ' - acq signal T_{1-1} - real and imag'],'fontweight','bold');
                        subplot(3,1,2)
                            hold on;
                                plot(tt,abs(sim.tissue(tis).acq.signal0).*sim.tissue(tis).acq.acq{k,m},'k','linewidth',2);
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);                
                                plot(tt, abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':k');  % pulses
                                plot(tt,-abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':k'); % pulses
                                plot(tt,sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':','color',color_acq);               % acq
                                plot(tt,-sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.tissue(tis).acq.signal0)).*1.00,':','color',color_acq);              % acq
                            hold off; box on;
                            xlim([tt(1) tt(end)]); xlabel('time [ms]'); ylabel('signal (abs)');
                            title([sim.tissue(tis).tissue.name ' - acq signal T_{1-1} - abs'],'fontweight','bold');
                        subplot(3,1,3)
                            hold on;
                                plot(tt,mod(angle(sim.tissue(tis).acq.signal0),2*pi).*sim.tissue(tis).acq.acq{k,m},'--k','linewidth',2);
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);                
                                plot(tt, abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(mod(angle(sim.tissue(tis).acq.signal0),2*pi))).*1.00,':k');  % pulses
                                plot(tt,-abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(mod(angle(sim.tissue(tis).acq.signal0),2*pi))).*1.00,':k'); % pulses
                                plot(tt,sim.tissue(tis).acq.acq{k,m}.*max(abs(mod(angle(sim.tissue(tis).acq.signal0),2*pi))).*1.00,':','color',color_acq);               % acq
                                plot(tt,-sim.tissue(tis).acq.acq{k,m}.*max(abs(mod(angle(sim.tissue(tis).acq.signal0),2*pi))).*1.00,':','color',color_acq);              % acq
                            hold off; box on;
                            xlim([tt(1) tt(end)]); xlabel('time [ms]'); ylabel('signal (phase)');
                            title([sim.tissue(tis).tissue.name ' - acq signal T_{1-1} - phase'],'fontweight','bold');
                        annotation(fig3,'textbox',[0.06 0.95 0.90 0.05],'string',name_fig3,'horizontalalignment','center', ...
                                                    'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % figure title
 
                end
            end
            
        end
        
        % ---- display info
        if (sim.display.info.acquisition)
            
            if (sim.param.acq_after_each_pulse)
                disp(' acquisition info');
                disp(['   tissue ' num2str(tis) ' - ' sim.tissue(tis).tissue.name ' - acquisition parameters']);
                for i=1:sim.tissue(tis).acq.n_spectra
                    disp(['     spectrum # ' num2str(i)]);
                    disp(['       acq dwell time       [us] = ' num2str(sim.tissue(tis).acq.dt{i}.*1e6)]);
                    disp(['       acq bandwidth  [Hz/pixel] = ' num2str(sim.tissue(tis).acq.bw{i})]);
                    disp(['       acq duration         [us] = ' num2str(sim.tissue(tis).acq.duration{i}.*1e6)]);
                    disp(['       acq pixel number          = ' num2str(sim.tissue(tis).acq.n_pixels{i})]);
                    disp(['       delay before acq     [us] = ' num2str((sim.tissue(tis).acq.delay_before_acq{i}).*1e6)]);
                    disp( '       snr info: i=n_spectra, j=n_dt, k=n_duration, m=n_delay_before_acq');
                    for j=1:sim.tissue(tis).acq.nn_dt{i}
                        for k=1:sim.tissue(tis).acq.nn_duration{i}
                            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}           
                                disp(['       snr abs(signal) / abs(fft) [ijkm: ' num2str(i) ' ' num2str(j) ' ' num2str(k) ' ' num2str(m) '] = ' ...
                                     num2str(roundn((sim.tissue(tis).acq.snr_signal{i,j,k,m}),-1)) '  /  ' num2str(roundn((sim.tissue(tis).acq.snr_signal_fft{i,j,k,m}),-1))]);
                                disp(['       max abs(signal)            = ' num2str(sim.tissue(tis).acq.signal_abs_max{i,j,k,m})]);                       
                            end
                        end
                    end
                end
            else  
                disp(' acquisition info');
                disp(['   tissue ' num2str(tis) ' - ' sim.tissue(tis).tissue.name ' - acquisition parameters']);
                for i=1:sim.tissue(tis).acq.n_spectra   % n_spectra = 1 (after last pulse)
                    disp(['     spectrum # ' num2str(i)]);
                    disp(['       acq dwell time       [us] = ' num2str(sim.tissue(tis).acq.dt{end}.*1e6)]);
                    disp(['       acq bandwidth  [Hz/pixel] = ' num2str(sim.tissue(tis).acq.bw{end})]);
                    disp(['       acq duration         [us] = ' num2str(sim.tissue(tis).acq.duration{end}.*1e6)]);
                    disp(['       acq pixel number          = ' num2str(sim.tissue(tis).acq.n_pixels{end})]);
                    disp(['       delay before acq     [us] = ' num2str((sim.tissue(tis).acq.delay_before_acq{end}).*1e6)]);
                    disp( '       snr info: i=n_spectra, j=n_dt, k=n_duration, m=n_delay_before_acq');
                    for j=1:sim.tissue(tis).acq.nn_dt{i}
                        for k=1:sim.tissue(tis).acq.nn_duration{i}
                            for m=1:sim.tissue(tis).acq.nn_delay_before_acq{i}           
                                disp(['       snr abs(signal) / abs(fft) [ijkm: ' num2str(i) ' ' num2str(j) ' ' num2str(k) ' ' num2str(m) '] = ' ...
                                     num2str(roundn((sim.tissue(tis).acq.snr_signal{i,j,k,m}),-1)) '  /  ' num2str(roundn((sim.tissue(tis).acq.snr_signal_fft{i,j,k,m}),-1))]);
                                disp(['       max abs(signal)            = ' num2str(sim.tissue(tis).acq.signal_abs_max{i,j,k,m})]);                       
                            end
                        end
                    end
                end
            end
            
        end
        
    end

end
% =================================================================================================

