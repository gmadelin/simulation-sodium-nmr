% =================================================================================================
% function    : plot_evolution_spherical_tensor_0201
% -------------------------------------------------------------------------------------------------
% purpose     : plot evolution of spherical tensors for each tissue during the sequence 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : individual figure for each tissue  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = plot_evolution_spherical_tensor_0201(sim)
    
    % ---- plot evolution spherical tensor
    if (sim.display.figure.evolution_sph_tensor)
        
        % ---- select tissues to plot
        switch sim.display.figure.tissues
            case 'all',  tis1 = 1;                           tisn = sim.param.n_tissues;
            otherwise,   tis1 = sim.display.figure.tissues;  tisn = tis1;  
        end
    
        % ---- plot figure for each selected tissue
        for tis=tis1:tisn        
            tt = sim.seq.t*1e3;          % [ms]
            fig1 = figure('name',['evolution spherical tensor - ' sim.tissue(tis).tissue.name]); 
                for i=1:15
                    subplot(3,5,i);
                    hold on;
                        plot(tt,real(sim.tissue(tis).sphtens.tensor(i+1,:)),'b','linewidth',1);
                        plot(tt,imag(sim.tissue(tis).sphtens.tensor(i+1,:)),'r','linewidth',1);
                        plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                        plot(tt,abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).sphtens.tensor(i+1,:))).*1.00,':k');     % pulses
                        plot(tt,-abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis).sphtens.tensor(i+1,:))).*1.00,':k');    % pulses
                    hold off; box on;
                    title(sim.tissue(tis).sphtens.name(i+1),'FontWeight','bold');
                    xlim([tt(1) tt(end)]); xlabel('time [ms]');
                    annotation(fig1,'textbox',[0.40 0.94 0.25 0.05],'string',{['evolution spherical tensor - ' sim.tissue(tis).tissue.name]},'horizontalalignment','center',...
                                    'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);  % overall figure title
                end 
            
        end
        
    end
       
end
% =================================================================================================
