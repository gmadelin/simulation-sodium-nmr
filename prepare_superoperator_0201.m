% =================================================================================================
% function    : prepare_superoperator_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare superoperators Ix, Iy, Iz for spin 3/2, and initial vector 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : - 
% reference   : Lee JS, Regatte RR, Jerschow A, J Chem Phys 131, 174501, 2009 
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - jaeseung.lee@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = prepare_superoperator_0201(sim)

    % ---- prepare superoperators for simulation
    for tis=1:sim.param.n_tissues

        % ---- sodium spin
        sim.tissue(tis).superop.spin = 3/2;

        % ---- constants
        n_supop = 20;
        s = sqrt(3)/2;

        % ---- Ix superoperator
        sim.tissue(tis).superop.Ix = zeros(n_supop,n_supop);
        sim.tissue(tis).superop.Ix( 5, 9) = -s;  sim.tissue(tis).superop.Ix( 5,12) =  s;  sim.tissue(tis).superop.Ix( 6, 9) =  s;  sim.tissue(tis).superop.Ix( 6,10) = -1;  sim.tissue(tis).superop.Ix( 6,12) = -s;  sim.tissue(tis).superop.Ix( 6,13) =  1;
        sim.tissue(tis).superop.Ix( 7,10) =  1;  sim.tissue(tis).superop.Ix( 7,11) = -s;  sim.tissue(tis).superop.Ix( 7,13) = -1;  sim.tissue(tis).superop.Ix( 7,14) =  s;
        sim.tissue(tis).superop.Ix( 8,11) =  s;  sim.tissue(tis).superop.Ix( 8,14) = -s;  sim.tissue(tis).superop.Ix( 9,15) = -1;  sim.tissue(tis).superop.Ix(10,15) =  s;  sim.tissue(tis).superop.Ix(10,16) = -s;  
        sim.tissue(tis).superop.Ix(11,16) =  1;  sim.tissue(tis).superop.Ix(12,17) =  1;  sim.tissue(tis).superop.Ix(13,17) = -s;  sim.tissue(tis).superop.Ix(13,18) =  s;  
        sim.tissue(tis).superop.Ix(14,18) = -1;  sim.tissue(tis).superop.Ix(15,19) = -s;  sim.tissue(tis).superop.Ix(16,19) =  s;  sim.tissue(tis).superop.Ix(17,20) =  s;  sim.tissue(tis).superop.Ix(18,20) = -s;  
        sim.tissue(tis).superop.Ix = sim.tissue(tis).superop.Ix + sim.tissue(tis).superop.Ix';

        % ---- Iy superoperator
        sim.tissue(tis).superop.Iy = zeros(n_supop,n_supop); 
        sim.tissue(tis).superop.Iy( 5, 9) = -1i*s;  sim.tissue(tis).superop.Iy( 5,12) =  -1i*s;   sim.tissue(tis).superop.Iy( 6, 9) =   1i*s;  sim.tissue(tis).superop.Iy( 6,10) =    -1i;  sim.tissue(tis).superop.Iy( 6,12) =   1i*s;  sim.tissue(tis).superop.Iy( 6,13) = -1i;
        sim.tissue(tis).superop.Iy( 7,10) =    1i;  sim.tissue(tis).superop.Iy( 7,11) =  -1i*s;   sim.tissue(tis).superop.Iy( 7,13) =     1i;  sim.tissue(tis).superop.Iy( 7,14) =  -1i*s;
        sim.tissue(tis).superop.Iy( 8,11) =  1i*s;  sim.tissue(tis).superop.Iy( 8,14) =   1i*s;   sim.tissue(tis).superop.Iy( 9,15) =    -1i;  sim.tissue(tis).superop.Iy(10,15) =   1i*s;  sim.tissue(tis).superop.Iy(10,16) =  -1i*s;  
        sim.tissue(tis).superop.Iy(11,16) =    1i;  sim.tissue(tis).superop.Iy(12,17) =    -1i;   sim.tissue(tis).superop.Iy(13,17) =   1i*s;  sim.tissue(tis).superop.Iy(13,18) =  -1i*s;  
        sim.tissue(tis).superop.Iy(14,18) =    1i;  sim.tissue(tis).superop.Iy(15,19) =  -1i*s;   sim.tissue(tis).superop.Iy(16,19) =   1i*s;  sim.tissue(tis).superop.Iy(17,20) =  -1i*s;  sim.tissue(tis).superop.Iy(18,20) =   1i*s;  
        sim.tissue(tis).superop.Iy = sim.tissue(tis).superop.Iy + sim.tissue(tis).superop.Iy';

        % ---- Iz superoperator
        sim.tissue(tis).superop.Iz = diag([ 0  0  0  0   0  0  0  0   1  1  1   -1 -1 -1   2  2   -2 -2   3   -3]);

        % ---- Ix, Iy, Iz vectors (for info)
        sim.tissue(tis).superop.Ix_v = [  0 0 0 0   0 0 0 0       sqrt(3)/2   1     sqrt(3)/2      sqrt(3)/2  1    sqrt(3)/2    0 0   0 0   0   0];
        sim.tissue(tis).superop.Iy_v = [  0 0 0 0   0 0 0 0   -1i*sqrt(3)/2 -1i -1i*sqrt(3)/2   1i*sqrt(3)/2 1i 1i*sqrt(3)/2    0 0   0 0   0   0];
        sim.tissue(tis).superop.Iz_v = [  0 0 0 0   3/2 1/2 -1/2 -3/2   0 0 0   0 0 0   0 0   0 0   0   0];

        % ---- initial state = Iz
        sim.tissue(tis).superop.rho0 = [3/2 1/2 -1/2 -3/2   3/2 1/2 -1/2 -3/2   0 0 0   0 0 0   0 0   0 0   0   0]'; 

    end
        
end
% =================================================================================================
