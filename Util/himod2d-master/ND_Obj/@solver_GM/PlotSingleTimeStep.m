function PlotSingleTimeStep(obj,num)
% function PlotSingleTimeStep(obj,num)
%
%   Plotta secondo i settings grafici scelti (per il default si veda il
%   costruttore)
%   Attenzione se viene cambiato il dt o il obj.t_0 l'istante di tempo a cui si fa
%   riferimento sarà sbagliato.
%
%
if(obj.full2d)
    %plot solo 2D
    figure(obj.currentfig);
    set(gcf,'Color',[1 1 1]);
    if (obj.treD)
        pdeplot(obj.pointsfull2d,obj.segfull2d,obj.trifull2d,'xydata',obj.Solution{num,1},'zdata',obj.Solution{num,1},'mesh','off','colormap','jet');
        caxis([obj.cinf,obj.csup]);
        if(obj.timedipendent)
            title(['t = ',num2str(obj.t_0+(num-1)*obj.dt)]);
        end
        if(obj.frontale)
            view(0,0)
            zlim([obj.cinf,obj.csup])
        else
            view(obj.az,obj.el)
            zlim([obj.cinf,obj.csup])
        end
    else
        pdeplot(obj.pointsfull2d,[],obj.trifull2d,'xydata',obj.Solution{num,1},'colormap','jet','contour','on','mesh','off','levels',obj.levels);
        %caxis([obj.cinf,obj.csup]);
        %WARNING chenge
        caxis([0.0,0.16]);
        if(obj.timedipendent)
            title(['t = ',num2str(obj.t_0+(num-1)*obj.dt)]);
        end
    end
else
    %plot soluzione 1D+2D
    figure(obj.currentfig);
    %set(gcf,'Color',[1 1 1]);
    colormap jet;
    if (obj.treD)
        pdeplot(obj.points2d,obj.seg2d,obj.tri2d,'xydata',obj.Solution{num,1},'zdata',obj.Solution{num,1},'mesh','off','colormap','winter');
        hold on
        if(~obj.krpo)
            pdeplot(obj.points1d,obj.seg1d,obj.tri1d,'xydata',obj.Solution{num,2},'zdata',obj.Solution{num,2},'mesh','off','colormap','winter');
        else
            obj.Plot_1d_krpo(num);
        end
        caxis([obj.cinf,obj.csup]);
        if(obj.frontale)
            view(0,0)
            zlim([obj.cinf,obj.csup])
        else
            view(obj.az,obj.el)
            zlim([obj.cinf,obj.csup])
        end
        if(obj.timedipendent)
            title(['t = ',num2str(obj.t_0+(num-1)*obj.dt)]);
        end
    else
        pdeplot(obj.points2d,[],obj.tri2d,'xydata',obj.Solution{num,1},'colormap','jet','contour','on','mesh','off','levels',obj.levels);
        hold on
 
        if(~obj.krpo)
            pdeplot(obj.points1d,obj.seg1d,obj.tri1d,'xydata',obj.Solution{num,2},'mesh','off','colormap','jet');
        else
            obj.Plot_1d_krpo(num);
        end
        caxis([obj.cinf,obj.csup]);
        if(obj.timedipendent)
            title(['t = ',num2str(obj.t_0+(num-1)*obj.dt)]);
        end
    end
end
%axis equal

set(gca, ...
    'FontUnits','points', ...
    'FontSize',11, ...
    'YTick',[-0.5,0,0.5],...
    'YTickLabel',[0,0.5,1]);
set(gca,'XLim',[0,obj.lun2d+obj.lun1d],'YLim',[-0.5,0.5]);

hold off
end