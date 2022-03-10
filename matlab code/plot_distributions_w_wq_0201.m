% =================================================================================================
% function    : plot_distributions_w_wq_0201
% -------------------------------------------------------------------------------------------------
% purpose     : plot w and wq distributions for each tissue 
% input       : sim (struct)
% output      : sim (struct) 
% comment     : plot individual figure for each tissue  
% reference   : -
% -------------------------------------------------------------------------------------------------
% date-author : 2013/04 - guillaume.madelin@nyumc.org
%               2018/05 - guillaume.madelin@nyumc.org 
% =================================================================================================
function [sim] = plot_distributions_w_wq_0201(sim)

    % ---- plot w and wq distributions  
    if (sim.display.figure.distribution_w_wq)
        
        % ---- number of tissues to plot (if chemical exchange, don't plot virtual tissue 3)
        if (sim.param.chemical_exchange),  tisn = 2;  
        else,                              tisn = sim.param.n_tissues;  end
        
        % ---- plot figures
        for tis=1:tisn
            figure('name',['w and wq - ' sim.tissue(tis).tissue.name]); 
                hold on;
                    [ax,h1,h2] = plotyy(sim.tissue(tis).w.f.range,sim.tissue(tis).w.dist,sim.tissue(tis).wq.fq.range,sim.tissue(tis).wq.dist);
                hold off; box on;
                set(get(ax(1),'ylabel'),'string','f [%]');  set(h1,'linestyle','-','marker','.','markersize',16);
                set(get(ax(2),'ylabel'),'string','fq [%]');  set(h2,'linestyle','-','marker','o','markersize',7);
                title(['f and fq distributions - ' sim.tissue(tis).tissue.name],'fontweight','bold','fontsize',12,'color',[0.45 0.45 0.45]);
                xlabel('f and fq [hz]');  legend('f','fq');
        end 
        
    end 
    
end
% =================================================================================================
