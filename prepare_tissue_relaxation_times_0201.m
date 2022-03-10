% =================================================================================================
% function    : prepare_tissue_relaxation_times_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare relaxation times of different tissues
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_tissue_relaxation_times_0201(sim)

    % ---- prepare relaxation times
    for tis=1:sim.param.n_tissues
        
        switch sim.tissue(tis).tissue.name
            
            case 'csf'                                              % cerebrospinal fluid
                sim.tissue(tis).relax.t1short   = 64 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 64 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 56 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 56 * 1e-3;            % [s]
                
            case 'parenchyma'                                       % brain parenchyma (average)
                sim.tissue(tis).relax.t1short   = 35 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 35 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  5 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 25 * 1e-3;            % [s]
                
            case 'ec'                                               % extracellular
                sim.tissue(tis).relax.t1short   = 46 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 46 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  3 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 30 * 1e-3;            % [s]
                
            case 'ic'                                               % intracellular
                sim.tissue(tis).relax.t1short   = 24 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 24 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  2 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 14 * 1e-3;            % [s]
                
            case 'plasma'                                           % blood plasma
                sim.tissue(tis).relax.t1short   = 45 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 45 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 15 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 25 * 1e-3;            % [s]
                
            case 'gel0%'                                            % 140 mM + water
                sim.tissue(tis).relax.t1short   = 52 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 52 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 52 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 52 * 1e-3;            % [s]

            case 'gel1%'                                            % 1% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 47 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 47 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 22 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 40 * 1e-3;            % [s]

            case 'gel2%'                                            % 2% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 43 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 43 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 11 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 37 * 1e-3;            % [s]
            
            case 'gel3%'                                            % 3% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 36 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 36 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  8 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 28 * 1e-3;            % [s]
            
            case 'gel4%'                                            % 4% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 35 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 35 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  7 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 27* 1e-3;             % [s]
            
            case 'gel5%'                                            % 5% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 33 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 33 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  6 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 26 * 1e-3;            % [s]
           
            case 'gel6%'                                            % 6% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 31 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 31 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  5 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 23 * 1e-3;            % [s]

            case 'gel7%'                                            % 7% agar + 140 mM NaCl
                sim.tissue(tis).relax.t1short   = 28 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 28 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  4 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 27 * 1e-3;            % [s]
                                
            case 'cartilage'                                        % cartilage
                sim.tissue(tis).relax.t1short   = 20 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 20 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  1 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 13 * 1e-3;            % [s]

            case 'liquid'                                           % liquid 140 mM NaCl (ideal/theoretical values with t1=t2)
                sim.tissue(tis).relax.t1short   = 60 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 60 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 60 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 60 * 1e-3;            % [s]
                
            case 'synovial'                                         % synovial fluid in knee
                sim.tissue(tis).relax.t1short   = 45 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 45 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   = 40 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 40 * 1e-3;            % [s]
                     
            case 'pf1'                                              % pf1 bacteriophage
                sim.tissue(tis).relax.t1short   = 12 * 1e-3;            % [s]
                sim.tissue(tis).relax.t1long    = 20 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2short   =  1 * 1e-3;            % [s]
                sim.tissue(tis).relax.t2long    = 10 * 1e-3;            % [s]

            case 'other'                                            % other tissue...
                sim.tissue(tis).relax.t1short   = 10000 * 1e-3;         % [s]
                sim.tissue(tis).relax.t1long    = 10000 * 1e-3;         % [s]
                sim.tissue(tis).relax.t2short   = 10000 * 1e-3;         % [s]
                sim.tissue(tis).relax.t2long    = 10000 * 1e-3;         % [s]
                
            otherwise                                               % used for virtual tissue 3 in case of chemical exchange between tissues 1 and 2
                sim.tissue(tis).relax.t1short   = NaN;                  % [s]
                sim.tissue(tis).relax.t1long    = NaN;                  % [s]
                sim.tissue(tis).relax.t2short   = NaN;                  % [s]
                sim.tissue(tis).relax.t2long    = NaN;                  % [s]           
                
        end
                
        % ---- optimal tau for mqf sequences
        sim.tissue(tis).relax.tau_opt = (sim.tissue(tis).relax.t2short.*sim.tissue(tis).relax.t2long./(sim.tissue(tis).relax.t2long-sim.tissue(tis).relax.t2short)).*log(sim.tissue(tis).relax.t2long./sim.tissue(tis).relax.t2short);             
    
        % ---- sodium gamma
        sim.tissue(tis).gamma = 7.0808493 * 1e7;                        % [rad/s/T] -> sodium gyromagnetic ratio
        
    end
    
    % ---- display info
    if (sim.display.info.relaxation_times)
        disp(' relaxation times [ms]');
        for tis=1:sim.param.n_tissues
            disp(['   tissue ' num2str(tis)]);
            disp(['     t1short = ' num2str(1000*sim.tissue(tis).relax.t1short) ', t1long = ' num2str(1000*sim.tissue(tis).relax.t1long)]);
            disp(['     t2short = ' num2str(1000*sim.tissue(tis).relax.t2short) ', t2long = ' num2str(1000*sim.tissue(tis).relax.t2long)]);
        end
    end
    
end
% =================================================================================================
