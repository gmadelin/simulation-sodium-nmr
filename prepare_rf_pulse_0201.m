% =================================================================================================
% function    : prepare_rf_pulse_0201
% -------------------------------------------------------------------------------------------------
% purpose     : prepare rf pulses 
% input       : p (struct), plot_fig ([0,1])
% output      : p (struct) 
% comment     : -  
% reference   : - 
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [p] = prepare_rf_pulse_0201(p,plot_fig)

    % ---- input variables
    if (nargin<2),  plot_fig = 0;  end

    % ---- timing
    p.n_points  = round(p.duration/p.dt);     % []
    p.t         = (0:1:p.n_points-1)*p.dt;    % [s]
    p.tau       = 2*p.t/p.duration-1;         % [-1 1]
    p.dtau      = 2*p.dt/p.duration;          % []

    % ---- rf pulse design
    switch p.type

        case 'rect'
            % -- pulse parameters
            p.w1max     = p.alpha./p.duration;                                  % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s]
            % -- pulse design
            p.fn        = ones(1,p.n_points);                                   % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points);                              % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]
            
        case 'sinc' 
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s]
            p.t0        = 1./(2+p.n1+p.n2);                                     % []
            p.nmax      = max(p.n1,p.n2);                                       % []
            p.tau1      = linspace(-p.n1.*p.t0,p.n2.*p.t0,p.n_points);          % []
            % -- pulse design
            p.fn        = sinc(p.tau1/p.t0);                                    % []
            p.window    = (1-p.beta)+p.beta.*cos(pi*p.tau1/(p.nmax.*p.t0));     % []
            p.w1        = p.w1max.*p.window.*p.fn;                              % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points);                              % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]

        case 'hsecn'
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s]
            p.a         = 2*pi*p.bandwidth/2;                                   % [rad/s]
            % -- pulse design
            p.fn        = sech(p.beta.*p.tau.^p.n1);                            % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points) + 2*p.a.*( ((cumtrapz(p.fn.^2)).*p.dtau) ./ ((sum(p.fn.^2)).*p.dtau) - 1/2 );     % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]         

        case 'ahp_sincos'
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s]
            p.a         = 2*pi*p.bandwidth/2;                                   % [rad/s]
            % -- pulse design
            p.fn        = sin((pi/2).*((p.tau+1)/2));                           % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points) + p.w1max.*cos((pi/2).*((p.tau+1)/2));     % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]         
            
        case 'bir4'         % design from staewen rs et al. inv radiol 25, 559-567, 1990
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s]
            p.a         = 2*pi*p.bandwidth/2;                                   % [rad/s]
            % -- timing
            p.tau1      = linspace(0,1,length(p.t)/4);                          % []
            p.tau2      = linspace(1,2,length(p.t)/4);                          % []
            p.tau3      = linspace(2,3,length(p.t)/4);                          % []
            p.tau4      = linspace(3,4,length(p.t)/4);                          % []
            % -- pulse design - w1 
            p.fn1       = tanh(p.beta.*(1-p.tau1));                             % []
            p.fn2       = tanh(p.beta.*(p.tau2-1));                             % []
            p.fn3       = tanh(p.beta.*(3-p.tau3));                             % []
            p.fn4       = tanh(p.beta.*(p.tau4-3));                             % []
            p.fn        = [p.fn1 p.fn2 p.fn3 p.fn4];                            % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            % -- pulse design - wrf
            p.wrf1      = p.w.*ones(1,length(p.fn1)) + p.a.*(tan(p.n2*(p.tau1-0)))./tan(p.n2);  % [rad/s]
            p.wrf2      = p.w.*ones(1,length(p.fn2)) + p.a.*(tan(p.n2*(p.tau2-2)))./tan(p.n2);  % [rad/s]
            p.wrf3      = p.w.*ones(1,length(p.fn3)) + p.a.*(tan(p.n2*(p.tau3-2)))./tan(p.n2);  % [rad/s]
            p.wrf4      = p.w.*ones(1,length(p.fn4)) + p.a.*(tan(p.n2*(p.tau4-4)))./tan(p.n2);  % [rad/s]
            p.wrf       = [p.wrf1 p.wrf2 p.wrf3 p.wrf4];
            % -- pulse design - phi
            p.phi1      = cumtrapz(p.wrf1).*p.dt + p.phase.*ones(1,length(p.fn1));                                  % [rad]         
            p.phi2      = cumtrapz(p.wrf2).*p.dt + p.phase.*ones(1,length(p.fn2)) + p.phi1(end) + p.alpha/2 + pi;   % [rad]         
            p.phi3      = cumtrapz(p.wrf3).*p.dt + p.phase.*ones(1,length(p.fn3)) + p.phi2(end);                    % [rad]         
            p.phi4      = cumtrapz(p.wrf4).*p.dt + p.phase.*ones(1,length(p.fn4)) + p.phi3(end) - p.alpha/2 - pi;   % [rad]         
            p.phi       = [p.phi1 p.phi2 p.phi3 p.phi4];
            
        case 'wurst'
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s], freq offset
            p.a         = 2*pi*p.bandwidth/1;                                   % [rad/s]
            p.k         = 2*pi*p.bandwidth/p.duration;                          % [rad/s/s]
            % -- pulse design
            p.fn        = 1-abs(sin((pi/2)*p.tau)).^p.n1;                       % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points) + p.a.*p.tau;                 % [rad/s] 
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]

        case 'gauss'
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s], freq offset
            % -- pulse design
            p.fn        = exp((-p.tau.^2)./(2.*p.sigma.^2));                    % [], sigma = 0.1345 or 0.2690?, as duration = 7.434*sigma (with tau=[-1 1])
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]     
            p.wrf       = p.w.*ones(1,p.n_points);                              % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]

        case 'fermi'  
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s], freq offset
            p.a         = (2-p.t0)./13.81;                                      % []
            % -- pulse design
            p.fn        = 1./(1+exp((abs(p.tau)-p.t0)./p.a));                   % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]     
            p.wrf       = p.w.*ones(1,p.n_points);                              % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]

        case 'arb'
            % -- pulse parameters
            p.w1max     = 2*pi*p.amplitude;                                     % [rad/s]
            p.w         = 2*pi*p.freq;                                          % [rad/s], freq offset
            % -- pulse design
            p.fn        = ones(1,p.n_points);                                   % []
            p.w1        = p.w1max.*p.fn;                                        % [rad/s]
            p.wrf       = p.w.*ones(1,p.n_points);                              % [rad/s]
            p.phi       = cumtrapz(p.wrf).*p.dt + p.phase.*ones(1,p.n_points);  % [rad]
       
        otherwise
            error(' choose a rf pulse type');
    end  

    % ---- figure
    if (plot_fig)  
        
        figure;
            subplot(3,1,1)
                hold on;
                    plot(p.t,p.w1./(2*pi),'Linewidth',1);
                    plot([p.t(1) p.t(end)],[0 0],':k');
                hold off; box on;
                title('amplitude','FontWeight','bold'); xlabel('time [s]'); ylabel('amplitude [Hz]');
                if strcmp(p.type,'sinc'),  ylim(1.1*[min(p.w1/(2*pi)) max(p.w1/(2*pi))]);
                else                       ylim(1.1*[-10 max(p.w1/(2*pi))]);  end
            subplot(3,1,2)
                hold on;
                    plot(p.t,mod(p.phi,2*pi),'Linewidth',1);
                hold off; box on;
                title('phase','FontWeight','bold'); xlabel('time [s]'); ylabel('phase [rad]');
            subplot(3,1,3)
                hold on;
                    plot(p.t,p.wrf/(2*pi),'Linewidth',1); 
                    plot([p.t(1) p.t(end)],[p.freq p.freq],':k');
                hold off; box on;
                title('frequency sweep','FontWeight','bold'); xlabel('time [s]'); ylabel('frequency [Hz]');
                ylim(1.1*[-max(p.wrf)-0.1 max(p.wrf)+0.1]./(2*pi) + p.freq); 
                
    end

end
% =================================================================================================
