% =================================================================================================
% function    : calculate_spherical_tensor_coefficient_0201
% -------------------------------------------------------------------------------------------------
% purpose     : calculate spherical tensor coefficient from rho
% input       : rho (vector), l, m, s (scalars), norm (string)
% output      : sph_tens_coeff (array (2s+1)*(2s+1)) 
% comment     : -  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - alexej.jerschow@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sph_tens_coeff] = calculate_spherical_tensor_coefficient_0201(rho,l,m,s,norm)

    % ---- input variables
    if (nargin<5),  norm = 'h';  end

    % ---- calculate spherical tensors coefficients
    sph_tens = calculate_spherical_tensor_0201(l,m,s,norm);
    sph_tens_coeff = trace(sph_tens'*rho)/trace(sph_tens'*sph_tens);

end
% =================================================================================================
