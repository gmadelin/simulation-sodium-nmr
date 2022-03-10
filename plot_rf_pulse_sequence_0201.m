% =================================================================================================
% function    : plot_rf_pulse_sequence_0201
% -------------------------------------------------------------------------------------------------
% purpose     : plot the rf pulse sequence + acquisition window 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -   
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = plot_rf_pulse_sequence_0201(sim)
      
    % ---- plot only 1 figure (same sequence for all tissues) 
    if (sim.display.figure.sequence)
        
        for tis=1:1    
       
            color_acq   = [0.75 0.75 0.75];          % color for acquisition time
            linecolors  = 1-hot(sim.seq.n_phase);    % choice: hsv, jet, prism, gray, hot, cool, bone, copper, flag, winter, pink, spring, summer, autumn, lines
            tt          = sim.seq.t*1e3;             % [ms]
            
            for k=1:sim.tissue(tis).acq.nn_duration{1}
                for m=1:sim.tissue(tis).acq.nn_delay_before_acq{1}   
                    
                    name_fig1 = ['rf pulse sequence - acq duration [us] = ' num2str(sim.tissue(tis).acq.duration{end}(k)*1e6) ' - delay before acq [us] = ' num2str(sim.tissue(tis).acq.delay_before_acq{end}(m)*1e6)];
                               
                    fig1 = figure('name',name_fig1); 
                        subplot(2,2,1);
                            hold on;
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                                plot(tt,sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.seq.w1(1,:))./(2*pi)).*1.00,':','color',color_acq);     % acq
                                plot(tt,-sim.tissue(tis).acq.acq{k,m}.*max(abs(sim.seq.w1(1,:))./(2*pi)).*1.00,':','color',color_acq);    % acq
                                for j=1:sim.seq.n_phase
                                    plot(tt,abs(sim.seq.w1(j,:))./(2*pi),'-','linewidth',1,'color',linecolors(j,:));
                                end
                            hold off; box on;
                            xlabel('time [ms]'); ylabel('amplitude [hz]');
                            title('rf pulses - amplitude','fontweight','bold'); 
                            ylim([0 1.05*max(abs(sim.seq.w1(1,:))./(2*pi))]);
                        subplot(2,2,3);
                            hold on;
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                                for j=1:sim.seq.n_phase
                                    plot(tt,mod(sim.seq.phi(j,:),2*pi),'-','linewidth',1,'color',linecolors(j,:)); 
                                end
                            hold off; box on;
                            xlabel('time [ms]'); ylabel('phase [rad]');
                            title('rf pulses - phase','fontweight','bold');
                        subplot(2,2,2);
                            hold on;
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                                for j=1:sim.seq.n_phase
                                    plot(tt,real(sim.seq.w1(j,:).*exp(1i*sim.seq.phi(j,:)))./(2*pi),'-','linewidth',1,'color',linecolors(j,:));           
                                end
                            hold off; box on;
                            xlabel('time [ms]'); ylabel('real [hz]');
                            title('rf pulses - real','fontweight','bold'); 
                            ylim(1.05*[-max(abs(sim.seq.w1(1,:))) max(abs(sim.seq.w1(1,:)))]./(2*pi));
                        subplot(2,2,4);
                            hold on;
                                plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                                for j=1:sim.seq.n_phase
                                    plot(tt,imag(sim.seq.w1(j,:).*exp(1i*sim.seq.phi(j,:)))./(2*pi),'-','linewidth',1,'color',linecolors(j,:));
                                end
                            hold off; box on;
                            xlabel('time [ms]'); ylabel('imag [hz]');
                            title('rf pulses - imag','fontweight','bold');
                            ylim(1.05*[-max(abs(sim.seq.w1(1,:))) max(abs(sim.seq.w1(1,:)))]./(2*pi));
                        annotation(fig1,'textbox',[0.43 0.95 0.16 0.043],'string',name_fig1,'horizontalalignment','center',...
                                        'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % figure title

                end
            end            

        end
        
    end

end
% =================================================================================================
