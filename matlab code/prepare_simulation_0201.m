% =================================================================================================
% function    : prepare_simulation_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare simulation parameters: tissues and data structure
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_simulation_0201(sim)

    % ---- chemical exchange only when 2 tissues
    if (sim.param.n_tissues~=2), sim.param.chemical_exchange = 0;  end

    % ---- set w and wq for tissue 3 = 0 when chemical exchange between tissues 1 and 2 
    if (sim.param.chemical_exchange)       
        % -- add 1 tissue = combined tissue 1 + 2
        sim.param.n_tissues = 3;
        % -- w
        sim.tissue(3).w.index    = 0;                                                
        sim.tissue(3).w.f.f0     = 0;
        sim.tissue(3).w.f.min    = 0; 
        sim.tissue(3).w.f.df     = 0;
        sim.tissue(3).w.f.max    = 0; 
        sim.tissue(3).w.f.mean   = 0;
        sim.tissue(3).w.f.std    = 0;
        % -- wq
        sim.tissue(3).wq.index   = 0;                                                
        sim.tissue(3).wq.fq.fq0  = 0;
        sim.tissue(3).wq.fq.min  = 0; 
        sim.tissue(3).wq.fq.dfq  = 0;
        sim.tissue(3).wq.fq.max  = 0; 
        sim.tissue(3).wq.fq.mean = 0;
        sim.tissue(3).wq.fq.std  = 0;        
    end    
        
    % ---- simulation parameters for each tissue
    if (sim.param.chemical_exchange)
        for tis=1:2  
            sim.tissue(tis).tissue.name       = char(sim.param.tissues(tis));
            sim.tissue(tis).tissue.fraction   = sim.param.fraction(tis);
            sim.tissue(tis).chemical_exchange = sim.param.chemical_exchange;   
            sim.tissue(tis).k_ex              = sim.param.k_ex;        
        end
        % -- initialize virtual new combined tissue ( = tissue 1 + tissue 2) in case of chemical exchange between two tissues
        sim.tissue(3).tissue.name             = [sim.tissue(1).tissue.name '+' sim.tissue(2).tissue.name];
        sim.tissue(3).tissue.fraction         = sim.tissue(1).tissue.fraction + sim.tissue(2).tissue.fraction; 
        sim.tissue(3).chemical_exchange       = sim.param.chemical_exchange;       
        sim.tissue(3).k_ex                    = sim.param.k_ex;  
    else 
        for tis=1:sim.param.n_tissues
            sim.tissue(tis).tissue.name       = char(sim.param.tissues(tis));
            sim.tissue(tis).tissue.fraction   = sim.param.fraction(tis);
            sim.tissue(tis).chemical_exchange = sim.param.chemical_exchange;   
            sim.tissue(tis).k_ex              = sim.param.k_ex;        
        end
    end    
            
    % ---- display info
    if (sim.display.info.list_tissue)
        disp(' list of tissues');
        for tis=1:sim.param.n_tissues
            disp(['   tissue ' num2str(tis) ' = ' sim.tissue(tis).tissue.name]);   
        end
    end
    
end
% =================================================================================================
