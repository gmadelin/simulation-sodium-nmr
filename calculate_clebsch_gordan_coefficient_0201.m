% =================================================================================================
% function    : calculate_clebsch_gordan_coefficient_0201
% -------------------------------------------------------------------------------------------------
% purpose     : calculate clebsch-gordan coefficients for spin tensor
% input       : l, m, spin, norm (scalars)
% output      : cg_coeff
% comment     : -  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - alexej.jerschow@nyu.edu, jaeseung.lee@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [cg_coeff] = calculate_clebsch_gordan_coefficient_0201(j1,m1,j2,m2,j,m)

    % ---- cg_coeff=0 in some conditions
    if ( (m1+m2~=m) || (j1+j2<j) || (j2-j1>j) || (j1-j2>j) )
        cg_coeff = 0;
        return
    end

    % ---- delta
    delta = sqrt( factorial(j1+j2-j)*factorial(j1-j2+j)*factorial(-j1+j2+j) / factorial(j1+j2+j+1) );
    delta = delta * sqrt( (2*j+1)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j+m)*factorial(j-m) );

    % ---- cg_coeff0
    cg_coeff0 = 0;
    for n=max([(-m1+j2-j),(m2+j1-j),0]):min([(j1-m1),(j2+m2),(j1+j2-j)])
        if ( (j1-m1-n>=0) && (j-j2+m1+n>=0) && (j2+m2-n>=0) && (j-j1-m2+n>=0) && (j1+j2-j-n>=0) )
            cg_coeff0 = cg_coeff0 + (-1)^n / ( factorial(j1-m1-n)*factorial(j-j2+m1+n)*factorial(j2+m2-n)*factorial(j-j1-m2+n)*factorial(n)*factorial(j1+j2-j-n) );
        end
    end

    % ---- final cg_coeff
    cg_coeff = cg_coeff0*delta;

end
% =================================================================================================
