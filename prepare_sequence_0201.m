% =================================================================================================
% function    : prepare_sequence_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare rf pulse sequence: concatenate pulses and delays 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -  
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_sequence_0201(sim)

    % ---- rf pulse sequence
    for tis=1:sim.param.n_tissues
        
        % ---- initialize structures                                                   
        for i=1:length(sim.seq.param.rf.type)
            for j=1:size(sim.seq.param.rf.phase,1)
                sim.seq.pulse{i,j} = struct;  
                sim.seq.delay{i,j} = struct;
            end
        end   
        
        % ---- prepare parameters of pulses + delays from sequence parameters           
        for i=1:length(sim.seq.param.rf.type)
            for j=1:size(sim.seq.param.rf.phase,1)
                
                % -- initialize
                sim.seq.pulse{i,j} = struct;
                sim.seq.delay{i,j} = struct;
                
                % -- pulses
                sim.seq.pulse{i,j}.dt        = sim.seq.param.dt;                        % [s] -> pulse dwell time
                sim.seq.pulse{i,j}.type      = char(sim.seq.param.rf.type{i});          % [char] -> pulse type
                sim.seq.pulse{i,j}.alpha     = sim.seq.param.rf.alpha{i};               % [rad] -> flip angle, for rect pulse only
                sim.seq.pulse{i,j}.amplitude = sim.seq.param.rf.amplitude{i};           % [hz] -> pulse power, for adiabatic, gaussian, fermi, sinc pulses
                sim.seq.pulse{i,j}.duration  = sim.seq.param.rf.duration{i}.*1e-6;      % [us]->[s] -> pulse duration
                sim.seq.pulse{i,j}.bandwidth = sim.seq.param.rf.bandwidth{i};           % [hz] -> for adiabatic pulses
                sim.seq.pulse{i,j}.sigma     = sim.seq.param.rf.sigma{i};               % [hz] -> for gaussian pulses
                sim.seq.pulse{i,j}.t0        = sim.seq.param.rf.t0{i};                  % [] -> parameter for sinc, fermi pulses
                sim.seq.pulse{i,j}.beta      = sim.seq.param.rf.beta{i};                % [] -> parameter for sinc, hsecn pulses
                sim.seq.pulse{i,j}.n1        = sim.seq.param.rf.n1{i};                  % [] -> parameter for hsecn, wurst, sinc pulses
                sim.seq.pulse{i,j}.n2        = sim.seq.param.rf.n2{i};                  % [] -> parameter for hsecn, wurst, sinc pulses
                sim.seq.pulse{i,j}.phase     = sim.seq.param.rf.phase{j,i};             % [rad] -> pulse phase
                sim.seq.pulse{i,j}.freq      = sim.seq.param.rf.freq{i};                % [hz] -> pulse off-resonance frequency
                
                % -- delays after pulses
                sim.seq.delay{i,j}.dt        = sim.seq.param.dt;                        % [s]
                sim.seq.delay{i,j}.duration  = sim.seq.param.rf.delay{i}.*1e-6;         % [s]
                
            end
        end
        
        % ---- prepare pulses + delays
        for i=1:length(sim.seq.param.rf.type)
            for j=1:size(sim.seq.param.rf.phase,1)
                sim.seq.pulse{i,j} = prepare_rf_pulse_0201(sim.seq.pulse{i,j},0);       % run function to prepare rf pulse
                sim.seq.delay{i,j} = prepare_delay_0201(sim.seq.delay{i,j});            % run function to prepare delay 
            end
        end

        % ---- initialize variables
        sim.seq.n_pulse = size(sim.seq.pulse,1);                                        % number of pulses
        sim.seq.n_phase = size(sim.seq.pulse,2);                                        % number of phase cycles

        % ---- timing
        sim.seq.dt = sim.seq.param.dt;
        sim.seq.t = 0;                                                                  % first point = 0 (no rf)
        for i=1:sim.seq.n_pulse
            sim.seq.t = horzcat( sim.seq.t, sim.seq.t(end) + sim.seq.param.dt + sim.seq.pulse{i,1}.t );
            sim.seq.t = horzcat( sim.seq.t, sim.seq.t(end) + sim.seq.param.dt + sim.seq.delay{i,1}.t );
        end
        sim.seq.n_t = length(sim.seq.t);
       

        % ---- sequence = pulses (w1, phi, wrf) + delays (first point = 0: no rf)
        for j=1:sim.seq.n_phase
            sim.seq.w10{j}  = 0;                                                        % [rad/s] -> rf amplitude, 
            sim.seq.wrf0{j} = 0;                                                        % [rad/s] -> rf frequency sweep 
            sim.seq.phi0{j} = 0;                                                        % [rad] -> rf phase 
        end
        for j=1:sim.seq.n_phase
            for i=1:sim.seq.n_pulse    
                % -- w1 + delay
                sim.seq.w10{j}  = horzcat( sim.seq.w10{j},  sim.seq.pulse{i,j}.w1 );    
                sim.seq.w10{j}  = horzcat( sim.seq.w10{j},  sim.seq.delay{i,j}.delay );
                % -- phi + delay
                sim.seq.phi0{j} = horzcat( sim.seq.phi0{j}, sim.seq.pulse{i,j}.phi );    
                sim.seq.phi0{j} = horzcat( sim.seq.phi0{j}, sim.seq.delay{i,j}.delay );
                % -- wrf + delay
                sim.seq.wrf0{j} = horzcat( sim.seq.wrf0{j}, sim.seq.pulse{i,j}.wrf );    
                sim.seq.wrf0{j} = horzcat( sim.seq.wrf0{j}, sim.seq.delay{i,j}.delay );
            end
        end

        % ---- transfer data from cells to matrices
        sim.seq.w1  = [];
        sim.seq.phi = [];
        sim.seq.wrf = [];
        for j=1:sim.seq.n_phase
            sim.seq.w1  = vertcat(sim.seq.w1,  sim.seq.w10{j});
            sim.seq.phi = vertcat(sim.seq.phi, sim.seq.phi0{j});
            sim.seq.wrf = vertcat(sim.seq.wrf, sim.seq.wrf0{j});
        end
        
    end
    
    % ---- figure
    if (sim.display.figure.sequence_check) 
        linecolors = 1-hot(sim.seq.n_phase);            % choice: hsv, jet, prism, gray, hot, cool, bone, copper, flag, winter, pink, spring, summer, autumn, lines 
        fig1 = figure('name','rf pulse sequence');
            subplot(3,1,1);        
                hold on;
                    for j=1:sim.seq.n_phase,  plot(sim.seq.t,sim.seq.w1(j,:)./(2*pi),'color',linecolors(j,:));     end
                hold off; box on; 
                title('rf amplitude','fontweight','bold'); xlabel('time [s]'); ylabel('amplitude [hz]');
            subplot(3,1,2)
                hold on;
                    for j=1:sim.seq.n_phase,  plot(sim.seq.t,mod(sim.seq.phi(j,:),2*pi),'color',linecolors(j,:));  end
                hold off; box on;
                title('rf phase','fontweight','bold'); xlabel('time [s]'); ylabel('phase [rad]');
            subplot(3,1,3)
                hold on;
                    for j=1:sim.seq.n_phase,  plot(sim.seq.t,sim.seq.wrf(j,:)./(2*pi),'color',linecolors(j,:));    end
                hold off; box on;
                title('rf frequency','fontweight','bold'); xlabel('time [s]'); ylabel('frequency [hz]');
            annotation(fig1,'textbox',[0.44 0.95 0.15 0.045],'string',{'rf pulse sequence'},'horizontalalignment','center',...
                            'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % figure title

    end

end
% =================================================================================================
