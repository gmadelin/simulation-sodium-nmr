% =================================================================================================
% function    : prepare_wq_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare residual quadrupolar interaction (rqi) distribution -> wq
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_wq_0201(sim)

    % ---- residual quadrupolar interaction (rqi)
    for tis=1:sim.param.n_tissues
        
        % ---- initialize
        if (sim.tissue(tis).wq.index)
            sim.tissue(tis).wq.fq.range = sim.tissue(tis).wq.fq.min:sim.tissue(tis).wq.fq.dfq:sim.tissue(tis).wq.fq.max;        % [Hz]
        else
            sim.tissue(tis).wq.fq.range = sim.tissue(tis).wq.fq.fq0;                                                            % [Hz]
            sim.tissue(tis).wq.fq.mean  = sim.tissue(tis).wq.fq.range;                                                          % [Hz]
            sim.tissue(tis).wq.fq.std   = 1;                                                                                    % [Hz]
        end
        sim.tissue(tis).wq.wq = 2*pi*sim.tissue(tis).wq.fq.range;                                                               % [rad/s]
        sim.tissue(tis).wq.n  = length(sim.tissue(tis).wq.wq);                                                                  % []

        % ---- gaussian distribution 
        sim.tissue(tis).wq.mean     = 2*pi*sim.tissue(tis).wq.fq.mean;                                                          % [rad/s]
        sim.tissue(tis).wq.std      = 2*pi*sim.tissue(tis).wq.fq.std;                                                           % [rad/s]
        sim.tissue(tis).wq.gaussian = (1/(sim.tissue(tis).wq.std.*sqrt(2*pi))).*exp(-0.5*((sim.tissue(tis).wq.wq-sim.tissue(tis).wq.mean)/sim.tissue(tis).wq.std).^2);      % []
        sim.tissue(tis).wq.ones     = ones(1,sim.tissue(tis).wq.n);                                                             % []
        sim.tissue(tis).wq.dist     = 100.* sim.tissue(tis).wq.gaussian ./ sum(sim.tissue(tis).wq.gaussian);                    % [%]   
    end
        
    % ---- display info
    if (sim.display.info.rqi)
        disp(' residual quadrupolar interaction (wq) parameters');
        for tis=1:sim.param.n_tissues
            % ---- display info
            disp(['   tissue ' num2str(tis)]);
            disp(['     gaussian distribution = ' num2str(sim.tissue(tis).wq.index)]);
            if sim.tissue(tis).wq.index,  disp(['     distribution min/df/max/mean/std [hz] = ' num2str(sim.tissue(tis).wq.fq.min) '/' num2str(sim.tissue(tis).wq.fq.dfq) '/' num2str(sim.tissue(tis).wq.fq.max) '/' num2str(sim.tissue(tis).wq.fq.mean) '/' num2str(sim.tissue(tis).wq.fq.std)]);
            else                          disp(['     rqi frequency [hz] = ' num2str(sim.tissue(tis).wq.fq.fq0)]);
            end
        end     
    end

end
% =================================================================================================
