% =================================================================================================
% function    : prepare_redfield_relaxation_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare redfield relaxation superoperator for spin 3/2
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : Lee JS, Regatte RR, Jerschow A, J Chem Phys 131, 174501, 2009 
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - jaeseung.lee@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_redfield_relaxation_0201(sim)

    % ---- redfield relaxation superoperator
    for tis=1:sim.param.n_tissues

        % ---- calculate spectral densities j0, j1 and j2 from relaxation times 
        sim.tissue(tis).relax.relax_rates = 1./[sim.tissue(tis).relax.t1short  sim.tissue(tis).relax.t1long  sim.tissue(tis).relax.t2short  sim.tissue(tis).relax.t2long]';     % [r1short r1long r2short r2long]
        sim.tissue(tis).relax.m = [0 6 0; 0 0 6; 3 3 0; 0 3 3];                                     % [r1short r1long r2short r2long]' = m * [j0 j1 j2]' --> [j0 j1 j2]' = m \ [r1short r1long r2short r2long]' 
        sim.tissue(tis).relax.j = sim.tissue(tis).relax.m \ sim.tissue(tis).relax.relax_rates;      % notation: j0 = sim.tissue(tis).relax.j(1), j1 = sim.tissue(tis).relax.j(2), j2 = sim.tissue(tis).relax.j(3)         

        % ---- variables
        a = 3*sim.tissue(tis).relax.j(1);    % a = 3*J0
        b = 3*sim.tissue(tis).relax.j(2);    % b = 3*J1
        c = 3*sim.tissue(tis).relax.j(3);    % c = 3*J2  
        d = a + b + c;
        e = b + c;

        % ---- relaxation superoperator
        sim.tissue(tis).relax.R0 = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
        sim.tissue(tis).relax.R1 = [e -b -c 0; -b e 0 -c; -c 0 e -b; 0 -c -b e];
        sim.tissue(tis).relax.R2 = [d 0 -c; 0 e 0; -c 0 d];
        sim.tissue(tis).relax.R3 = sim.tissue(tis).relax.R2;
        sim.tissue(tis).relax.R4 = [d b; b d];
        sim.tissue(tis).relax.R5 = sim.tissue(tis).relax.R4;
        sim.tissue(tis).relax.R6 = [e 0; 0 e];
        sim.tissue(tis).relax.R  = blkdiag(sim.tissue(tis).relax.R0,sim.tissue(tis).relax.R1,sim.tissue(tis).relax.R2,sim.tissue(tis).relax.R3,sim.tissue(tis).relax.R4,sim.tissue(tis).relax.R5,sim.tissue(tis).relax.R6);
        sim.tissue(tis).relax.R(5,1) = -e;   sim.tissue(tis).relax.R(5,2) =  b;   sim.tissue(tis).relax.R(5,3) =  c;
        sim.tissue(tis).relax.R(6,1) =  b;   sim.tissue(tis).relax.R(6,2) = -e;   sim.tissue(tis).relax.R(6,4) =  c;
        sim.tissue(tis).relax.R(7,1) =  c;   sim.tissue(tis).relax.R(7,3) = -e;   sim.tissue(tis).relax.R(7,4) =  b;
        sim.tissue(tis).relax.R(8,2) =  c;   sim.tissue(tis).relax.R(8,3) =  b;   sim.tissue(tis).relax.R(8,4) = -e;
        
    end
    
    % ---- display info
    if (sim.display.info.spectral_densities)
        disp(' spectral densities [hz]');
        for tis=1:sim.param.n_tissues
            disp(['   tissue ' num2str(tis)]);
            disp(['     j0 = ' num2str(sim.tissue(tis).relax.j(1)) ', j1 = ' num2str(sim.tissue(tis).relax.j(2)) ', j2 = ' num2str(sim.tissue(tis).relax.j(3))]);
        end
    end

end
% =================================================================================================
