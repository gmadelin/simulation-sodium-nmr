% =================================================================================================
% function    : plot_evolution_magnetization_0201
% -------------------------------------------------------------------------------------------------
% purpose     : plot evolution of magnetization of all the tissues during the sequence 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : all tissues in same figure   
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = plot_evolution_magnetization_0201(sim)

    % ---- plot figure evolution of magnetization 
    if (sim.display.figure.evolution_magnetization) 
        
        % ---- select tissues to plot
        switch sim.display.figure.tissues
            case 'all',  tis1 = 1;                           tisn = sim.param.n_tissues;
            otherwise,   tis1 = sim.display.figure.tissues;  tisn = tis1;  
        end

        % ---- calculate magnetization from spherical tensors
        for tis=tis1:tisn
            sim.tissue(tis).mag(1,:)  = sim.tissue(tis).sphtens.tensor(2,:);                                                        % Iz
            sim.tissue(tis).mag(2,:)  = (1i/sqrt(2))*(sim.tissue(tis).sphtens.tensor(5,:)+sim.tissue(tis).sphtens.tensor(6,:));     % Iy
            sim.tissue(tis).mag(3,:)  = (1/sqrt(2)) *(sim.tissue(tis).sphtens.tensor(5,:)-sim.tissue(tis).sphtens.tensor(6,:));     % Ix
        end
                 
        % ---- tissue names for display in legend
        switch sim.display.figure.tissues
            case 'all'
                tissue_name = cell(1,sim.param.n_tissues);  
                for tis=tis1:tisn,  tissue_name{tis} = [' ' sim.tissue(tis).tissue.name];  end
                tis_w1 = 1;
            otherwise   
                tissue_name = cell(1,1);                    
                for tis=tis1:tisn,  tissue_name{1}   = [' ' sim.tissue(tis).tissue.name];  end
                tis_w1 = tis1;
        end
        
        % ---- prepare colors for magnetization plots
        color_mag = [0 0 1 ; 1 0 0 ; 0 1 0 ; 1 0 1 ; 0.5 0 0 ; 0 0 0.5 ; 0.5 0.5 0.5 ; 0.5 0 0.5];  % [blue, red, green, magenta, maroon, navy, gray, purple]

        % ---- figure
        tt = sim.seq.t*1e3;                % [ms] 
        name_mag = {'M_z','M_y','M_x'};    % magnetization names
        fig1 = figure('name','evolution of magnetization'); 
            for i=1:3
                subplot(1,3,i);
                hold on
                    for tis = tis1:tisn,  plot(tt,real(sim.tissue(tis).mag(i,:)),'linewidth',2,'color',color_mag(tis,:));  end          % magnetizations
                    plot([tt(1) tt(end)],[0 0],':k','linewidth',1);
                    plot(tt,abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis_w1).mag(i,:))).*1.00,':k');         % pulses
                    plot(tt,-abs(sim.seq.w1(1,:))./max(abs(sim.seq.w1(1,:))).*max(abs(sim.tissue(tis_w1).mag(i,:))).*1.00,':k');        % pulses
                hold off; box on;
                xlim([tt(1) tt(end)]); xlabel('time [ms]'); ylabel('magnetization'); 
                title(name_mag(i),'fontweight','bold');
            end
            legend(tissue_name');        
            annotation(fig1,'textbox',[0.365 0.94 0.30 0.05],'string',{'evolution of magnetization'},'horizontalalignment','center',...
                            'fontweight','bold','fontsize',12,'fitboxtotext','on','linestyle','none','color',[0.45 0.45 0.45]);        % overall figure title

    end
        
end
% =================================================================================================
