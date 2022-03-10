% =================================================================================================
% function    : calculate_spherical_tensor_0201
% -------------------------------------------------------------------------------------------------
% purpose     : calculate spherical tensor/matrix Tlm for spin s
% input       : l, m, s (scalars), norm (string)
% output      : sph_tens (array (2s+1)*(2s+1)) 
% comment     : -  
% reference   : Bowden GJ, Hutchison WD, J Magn Reson 67, 403-414, 1986
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - alexej.jerschow@nyu.edu, jaeseung.lee@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sph_tens] = calculate_spherical_tensor_0201(l,m,s,norm)

    % ---- input variables
    if (nargin<4),  norm = 'h';  end

    % ---- calculate sph_tens
    sph_tens = zeros(2*s+1);    % 4x4 matrix for spin s=3/2
    for m1=-s:s
        for m2=-s:s
            sph_tens(-m2+s+1,-m1+s+1) = calculate_clebsch_gordan_coefficient_0201(s,m1,l,m,s,m2);     % from AJ + JSL           
            % sph_tens(-m2+s+1,-m1+s+1) = calculate_clebsch_gordan_coefficient_0201(s,m1,s,m2,l,m);   % calculate_clebsch_gordan_coefficient_0201(j1,m1,j,m,j2,m2)?            
        end
    end

    % ---- normalization of sph_tens
    switch norm
        case 'h'    % after Bowden and Hutchinson
            if(2*s-l>=0),  sph_tens = sph_tens*factorial(l)*sqrt(factorial(2*s+l+1)/(2^l*factorial(2*l)*factorial(2*s-l)*(2*s+1)));  end
        case 'n'
            sph_tens = sph_tens*sqrt((2*l+1)/(2*s+1));
        otherwise
            return
    end

end
% =================================================================================================
