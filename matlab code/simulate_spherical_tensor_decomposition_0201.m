% =================================================================================================
% function    : simulate_spherical_tensor_decomposition_0201
% -------------------------------------------------------------------------------------------------
% purpose     : spherical tensor decomposition of superoperator density rho
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - jaeseung.lee@nyu.edu, alexej.jerschow@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = simulate_spherical_tensor_decomposition_0201(sim)

    for tis=1:sim.param.n_tissues

        % ---- tensor names
        sim.tissue(tis).sphtens.name = {'T00' 'T10' 'T20' 'T30' 'T1+1' 'T1-1' 'T2+1' 'T2+2' 'T3+1' 'T3+2' 'T3+3' 'T2-1' 'T2-2' 'T3-1' 'T3-2' 'T3-3'};

        % ---- density superoperator rho to decompose
        sim.tissue(tis).rho.rho_check = sim.tissue(tis).rho.rho_avg;    % rho averaged over w, wq and phase cycling

        % ---- spherical tensor decomposition
        for i=1:sim.seq.n_t        
            
            % -- temp matrix 4x4
            rho_temp = zeros(4);
            rho_temp(1,1) = sim.tissue(tis).rho.rho_check( 5,i);  rho_temp(2,2) = sim.tissue(tis).rho.rho_check( 6,i);   
            rho_temp(3,3) = sim.tissue(tis).rho.rho_check( 7,i);  rho_temp(4,4) = sim.tissue(tis).rho.rho_check( 8,i);
            rho_temp(1,2) = sim.tissue(tis).rho.rho_check( 9,i);  rho_temp(2,3) = sim.tissue(tis).rho.rho_check(10,i);   rho_temp(3,4) = sim.tissue(tis).rho.rho_check(11,i); 
            rho_temp(2,1) = sim.tissue(tis).rho.rho_check(12,i);  rho_temp(3,2) = sim.tissue(tis).rho.rho_check(13,i);   rho_temp(4,3) = sim.tissue(tis).rho.rho_check(14,i); 
            rho_temp(1,3) = sim.tissue(tis).rho.rho_check(15,i);  rho_temp(2,4) = sim.tissue(tis).rho.rho_check(16,i); 
            rho_temp(3,1) = sim.tissue(tis).rho.rho_check(17,i);  rho_temp(4,2) = sim.tissue(tis).rho.rho_check(18,i);
            rho_temp(1,4) = sim.tissue(tis).rho.rho_check(19,i);  rho_temp(4,1) = sim.tissue(tis).rho.rho_check(20,i);
            
            % -- calculate spherical tensor coefficients for each evolution point
            sim.tissue(tis).sphtens.tensor( 1,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 0,  0, sim.tissue(tis).superop.spin);    % T00         
            sim.tissue(tis).sphtens.tensor( 2,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 1,  0, sim.tissue(tis).superop.spin);    % T10
            sim.tissue(tis).sphtens.tensor( 3,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 2,  0, sim.tissue(tis).superop.spin);    % T20
            sim.tissue(tis).sphtens.tensor( 4,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3,  0, sim.tissue(tis).superop.spin);    % T30
            sim.tissue(tis).sphtens.tensor( 5,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 1, +1, sim.tissue(tis).superop.spin);    % T1+1
            sim.tissue(tis).sphtens.tensor( 6,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 1, -1, sim.tissue(tis).superop.spin);    % T1-1    
            sim.tissue(tis).sphtens.tensor( 7,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 2, +1, sim.tissue(tis).superop.spin);    % T2+1 
            sim.tissue(tis).sphtens.tensor( 8,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 2, +2, sim.tissue(tis).superop.spin);    % T2+2
            sim.tissue(tis).sphtens.tensor( 9,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, +1, sim.tissue(tis).superop.spin);    % T3+1
            sim.tissue(tis).sphtens.tensor(10,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, +2, sim.tissue(tis).superop.spin);    % T3+2
            sim.tissue(tis).sphtens.tensor(11,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, +3, sim.tissue(tis).superop.spin);    % T3+3 
            sim.tissue(tis).sphtens.tensor(12,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 2, -1, sim.tissue(tis).superop.spin);    % T2-1
            sim.tissue(tis).sphtens.tensor(13,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 2, -2, sim.tissue(tis).superop.spin);    % T2-2
            sim.tissue(tis).sphtens.tensor(14,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, -1, sim.tissue(tis).superop.spin);    % T3-1
            sim.tissue(tis).sphtens.tensor(15,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, -2, sim.tissue(tis).superop.spin);    % T3-2
            sim.tissue(tis).sphtens.tensor(16,i) = calculate_spherical_tensor_coefficient_0201( rho_temp, 3, -3, sim.tissue(tis).superop.spin);    % T3-3
            
        end

        % ---- put all tensor components < 1e-10 to zero
        sim.tissue(tis).sphtens.tensor(abs(sim.tissue(tis).sphtens.tensor)<1e-10) = 0;

    end

end
% =================================================================================================
