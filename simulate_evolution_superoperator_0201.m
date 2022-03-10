% =================================================================================================
% function    : simulate_evolution_superoperator_0201
% -------------------------------------------------------------------------------------------------
% purpose     : sodium spin density superoperator rho evolution under rf, relaxation, chemical shift and residual quadrupolar interaction 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : -  
% reference   : Lee JS, Regatte RR, Jerschow A, J Chem Phys 131, 174501, 2009 
% -------------------------------------------------------------------------------------------------
% date-author : 2012/08 - jaeseung.lee@nyu.edu
%               2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = simulate_evolution_superoperator_0201(sim)
    
    % ---- with chemical exchange (between 2 compartments only, with w and wq of same size = 1 only)
    if (sim.param.chemical_exchange)
   
        % ---- initialize density operator rho of the 2 tissues in chemical exchange
        for tis=1:2            
            sim.tissue(tis).rho.rho1(:,1) = sim.tissue(tis).tissue.fraction.*sim.tissue(tis).superop.rho0;   % start with Iz = T10 (1st point at t=0 with rf pulse = 0)            
        end
        
        % ---- superoperator evolution                 
        for k=1:sim.seq.n_phase                 % evolution over phase cycling            
            for i=1:sim.tissue(1).w.n           % evolution over rqi distribution wq                
                for j=1:sim.tissue(1).wq.n      % evolution over frequency offset distribution w
                    
                    % -- rqi (wq)
                    sim.tissue(1).wq.hq = (1/2)*sim.tissue(1).wq.wq(j).*diag([0 0 0 0 0 0 0 0 2 0 -2 -2 0 2 2 -2 -2 2 0 0]);  
                    sim.tissue(2).wq.hq = (1/2)*sim.tissue(2).wq.wq(j).*diag([0 0 0 0 0 0 0 0 2 0 -2 -2 0 2 2 -2 -2 2 0 0]);  
                    
                    % -- evolution step by step over rf pulse sequence 
                    for n=2:sim.seq.n_t 
                        
                        % -- rf pulse -> w1
                        sim.tissue(1).rho.t_w1 = expm(-1i*sim.seq.phi(k,n)*sim.tissue(1).superop.Iz) ...
                                                 * expm(-1i*sim.seq.w1(1,n)*sim.seq.dt*sim.tissue(1).superop.Ix) ...
                                                 * expm(1i*sim.seq.phi(k,n)*sim.tissue(1).superop.Iz);    
                        sim.tissue(2).rho.t_w1 = expm(-1i*sim.seq.phi(k,n)*sim.tissue(2).superop.Iz) ...
                                                 * expm(-1i*sim.seq.w1(1,n)*sim.seq.dt*sim.tissue(2).superop.Ix) ...
                                                 * expm(1i*sim.seq.phi(k,n)*sim.tissue(2).superop.Iz);    
                        
                        % -- chemical shift -> w
                        sim.tissue(1).rho.t_w = expm(-1i*sim.tissue(1).w.w(i)*sim.seq.dt*sim.tissue(1).superop.Iz);      
                        sim.tissue(2).rho.t_w = expm(-1i*sim.tissue(2).w.w(i)*sim.seq.dt*sim.tissue(2).superop.Iz); 
                        
                        % -- residual quadrupolar interaction -> wq
                        sim.tissue(1).rho.t_wq = expm(-1i*sim.seq.dt*sim.tissue(1).wq.hq);                    
                        sim.tissue(2).rho.t_wq = expm(-1i*sim.seq.dt*sim.tissue(2).wq.hq); 
                        
                        % -- redfield relaxation -> r
                        sim.tissue(1).rho.t_relax = expm(-sim.seq.dt*sim.tissue(1).relax.R);       
                        sim.tissue(2).rho.t_relax = expm(-sim.seq.dt*sim.tissue(2).relax.R);
                        
                        % -- rho evolution
                        sim.tissue(1).rho.rho1(:,n) = sim.tissue(1).rho.t_relax * sim.tissue(1).rho.t_w * sim.tissue(1).rho.t_wq ...
                                                      * sim.tissue(1).rho.t_w1 * sim.tissue(1).rho.rho1(:,n-1) * (1-sim.tissue(1).k_ex(1).*sim.seq.dt) ...
                                                      + sim.tissue(2).rho.t_relax * sim.tissue(2).rho.t_w * sim.tissue(2).rho.t_wq ...
                                                      * sim.tissue(2).rho.t_w1 * sim.tissue(2).rho.rho1(:,n-1) * ( +sim.tissue(2).k_ex(2).*sim.seq.dt);                    
                        sim.tissue(2).rho.rho1(:,n) = sim.tissue(1).rho.t_relax * sim.tissue(1).rho.t_w * sim.tissue(1).rho.t_wq ...
                                                      * sim.tissue(1).rho.t_w1 * sim.tissue(1).rho.rho1(:,n-1) * ( +sim.tissue(1).k_ex(2).*sim.seq.dt) ...
                                                      + sim.tissue(2).rho.t_relax * sim.tissue(2).rho.t_w * sim.tissue(2).rho.t_wq ...
                                                      * sim.tissue(2).rho.t_w1 * sim.tissue(2).rho.rho1(:,n-1) * (1-sim.tissue(2).k_ex(1).*sim.seq.dt);                    
                    
                    end
                    
                    % ---- final rho (taking into account amplitude of w and wq distributions (in %), and adc phase)
                    sim.tissue(1).rho.rho(:,:,k,i,j) = sim.tissue(1).rho.rho1(:,:).*sim.tissue(1).w.dist(i).*sim.tissue(1).wq.dist(j).*exp(1i.*sim.seq.param.adc.phase{k});
                    sim.tissue(2).rho.rho(:,:,k,i,j) = sim.tissue(2).rho.rho1(:,:).*sim.tissue(2).w.dist(i).*sim.tissue(2).wq.dist(j).*exp(1i.*sim.seq.param.adc.phase{k});
                    sim.tissue(3).rho.rho(:,:,k,i,j) = sim.tissue(1).rho.rho(:,:,k,i,j) + sim.tissue(2).rho.rho(:,:,k,i,j);     % add 2 superoperators from 2 tissues
                    
                end
            end  
        end         
        
        % ---- average rho over phase cycles + w and wq distributions
        for tis=1:3
            sim.tissue(tis).rho.rho_total1 = squeeze(sum(sim.tissue(tis).rho.rho,5));           % sum over phase cycles (k)
            sim.tissue(tis).rho.rho_total2 = squeeze(sum(sim.tissue(tis).rho.rho_total1,4));    % sum over w (i)
            sim.tissue(tis).rho.rho_total3 = squeeze(sum(sim.tissue(tis).rho.rho_total2,3));    % sum over wq (j)
            sim.tissue(tis).rho.rho_avg = sim.tissue(tis).rho.rho_total3./(sum(sim.tissue(1).w.dist).*sum(sim.tissue(1).wq.dist).*sim.seq.n_phase);  % average over phase cycles + w + wq
        end   
                    
    % ---- without chemical exchange
    else           
        
        for tis=1:sim.param.n_tissues

            % ---- initialize density operator rho for each tissue
            sim.tissue(tis).rho.rho1(:,1) = sim.tissue(tis).tissue.fraction.*sim.tissue(tis).superop.rho0;   % start with Iz = T10 (1st point at t=0 with rf pulse = 0)
                        
            % ---- superoperator evolution
            for k=1:sim.seq.n_phase                 % evolution over phase cycling
                for i=1:sim.tissue(tis).w.n         % evolution over rqi distribution wq
                    for j=1:sim.tissue(tis).wq.n    % evolution over frequency offset distribution w
                        
                        % -- rqi (wq)
                        sim.tissue(tis).wq.hq = (1/2)*sim.tissue(tis).wq.wq(j).*diag([0 0 0 0 0 0 0 0 2 0 -2 -2 0 2 2 -2 -2 2 0 0]);  
                        
                        % -- evolution step by step over rf pulse sequence 
                        for n=2:sim.seq.n_t
                            sim.tissue(tis).rho.t_w1      = expm(-1i*sim.seq.phi(k,n)*sim.tissue(tis).superop.Iz) ...
                                                            * expm(-1i*sim.seq.w1(1,n)*sim.seq.dt*sim.tissue(tis).superop.Ix) ...
                                                            * expm(1i*sim.seq.phi(k,n)*sim.tissue(tis).superop.Iz);                     % rf pulse -> w1
                            sim.tissue(tis).rho.t_w       = expm(-1i*sim.tissue(tis).w.w(i)*sim.seq.dt*sim.tissue(tis).superop.Iz);     % chemical shift -> w
                            sim.tissue(tis).rho.t_wq      = expm(-1i*sim.seq.dt*sim.tissue(tis).wq.hq);                                 % residual quadrupolar interaction -> wq
                            sim.tissue(tis).rho.t_relax   = expm(-sim.seq.dt*sim.tissue(tis).relax.R);                                  % redfield relaxation -> r
                            sim.tissue(tis).rho.rho1(:,n) = sim.tissue(tis).rho.t_relax * sim.tissue(tis).rho.t_w * sim.tissue(tis).rho.t_wq ...
                                                            * sim.tissue(tis).rho.t_w1 * sim.tissue(tis).rho.rho1(:,n-1);               % rho evolution                             
                        end
                        
                        % ---- final rho (taking into account amplitude of w and wq distributions (in %), and adc phase)
                        sim.tissue(tis).rho.rho(:,:,k,i,j) = sim.tissue(tis).rho.rho1(:,:).*sim.tissue(tis).w.dist(i).*sim.tissue(tis).wq.dist(j).*exp(1i.*sim.seq.param.adc.phase{k});
                        
                    end
                end  
            end

            % ---- average rho over phase cycles + w and wq distributions
            sim.tissue(tis).rho.rho_total1 = squeeze(sum(sim.tissue(tis).rho.rho,5));           % sum over phase cycles (k)
            sim.tissue(tis).rho.rho_total2 = squeeze(sum(sim.tissue(tis).rho.rho_total1,4));    % sum over w (i)
            sim.tissue(tis).rho.rho_total3 = squeeze(sum(sim.tissue(tis).rho.rho_total2,3));    % sum over wq (j)
            sim.tissue(tis).rho.rho_avg    = sim.tissue(tis).rho.rho_total3./(sum(sim.tissue(tis).w.dist).*sum(sim.tissue(tis).wq.dist).*sim.seq.n_phase);  % average over phase cycles + w + wq

        end

    end
    
end
% =================================================================================================
