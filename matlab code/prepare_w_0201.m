% =================================================================================================
% function    : prepare_w_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare frequency offset distribution -> w
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_w_0201(sim)

    % ---- frequency offset (w)
    for tis=1:sim.param.n_tissues

        % ---- initialize
        if (sim.tissue(tis).w.index)
            sim.tissue(tis).w.f.range = sim.tissue(tis).w.f.min:sim.tissue(tis).w.f.df:sim.tissue(tis).w.f.max;                 % [Hz]
        else
            sim.tissue(tis).w.f.range = sim.tissue(tis).w.f.f0;                                                                 % [Hz]
            sim.tissue(tis).w.f.mean  = sim.tissue(tis).w.f.range;                                                              % [Hz]
            sim.tissue(tis).w.f.std   = 1;                                                                                      % [Hz]
        end
        sim.tissue(tis).w.w = 2*pi*sim.tissue(tis).w.f.range;                                                                   % [rad/s]
        sim.tissue(tis).w.n = length(sim.tissue(tis).w.w);                                                                      % []

        % ---- gaussian distribution 
        sim.tissue(tis).w.mean     = 2*pi*sim.tissue(tis).w.f.mean;                                                             % [rad/s]
        sim.tissue(tis).w.std      = 2*pi*sim.tissue(tis).w.f.std;                                                              % [rad/s]
        sim.tissue(tis).w.gaussian = (1/(sim.tissue(tis).w.std.*sqrt(2*pi))).*exp(-0.5*((sim.tissue(tis).w.w-sim.tissue(tis).w.mean)/sim.tissue(tis).w.std).^2);    % []
        sim.tissue(tis).w.ones     = ones(1,sim.tissue(tis).w.n);                                                               % []
        sim.tissue(tis).w.dist     = 100.* sim.tissue(tis).w.gaussian ./ sum(sim.tissue(tis).w.gaussian);                       % [%]  
        
    end
     
    % ---- display info
    if (sim.display.info.frequency_offset)
        disp(' frequency offset (w) parameters');
        for tis=1:sim.param.n_tissues
            disp(['   tissue ' num2str(tis)]);
            disp(['     gaussian distribution = ' num2str(sim.tissue(tis).w.index)]);
            if sim.tissue(tis).w.index,  disp(['     distribution min/df/max/mean/std [hz] = ' num2str(sim.tissue(tis).w.f.min) '/' num2str(sim.tissue(tis).w.f.df) '/' num2str(sim.tissue(tis).w.f.max) '/' num2str(sim.tissue(tis).w.f.mean) '/' num2str(sim.tissue(tis).w.f.std)]);
            else,                        disp(['     frequency offset [hz] = ' num2str(sim.tissue(tis).w.f.f0)]);
            end    
        end
    end
  
end
% =================================================================================================
