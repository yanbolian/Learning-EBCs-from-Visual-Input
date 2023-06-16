function plotEBC(root,out,fign, fig_row)
%
% Downloaded from the Github of Hasselmo's lab
% Modified by Yanbo Lian for the purpose of displaying figures

%     %% occupancy circular
%     figure(fign); subplot(1,4,1);
%     % the +pi/2 brings "forwards" to "up"
%     [t2, r2] = meshgrid(wrapTo2Pi(out.params.thetaBins+pi/2), out.params.distanceBins(1:end-1));
%     [x, y] = pol2cart(t2,r2);
%     surface(x,y, out.occ), shading interp
%     hold on
%     set(gca,'XTick',[],'YTick',[])
%     
%     colormap(parula)
%     set(gca, 'YDir','Normal','CLim',[0 prctile(out.occ(:),99)])
%     set(gca,'YDir','Normal')  
%     title('occ')
% %     colorbar
%     axis off; axis square

    
%     %% nspk circular
%     figure(fign); subplot(1,4,2);
%     % the +pi/2 brings "forwards" to "up"
%     [t2, r2] = meshgrid(wrapTo2Pi(out.params.thetaBins+pi/2), out.params.distanceBins(1:end-1));
%     [x, y] = pol2cart(t2,r2);
%     surface(x,y, out.nspk), shading interp
%     hold on
%     set(gca,'XTick',[],'YTick',[])
%     
%     title([])
%     colormap(parula)
%     set(gca, 'YDir','Normal','CLim',[0 prctile(out.nspk(:),99)])
%     set(gca,'YDir','Normal')  
%     title('nspk')
% %     colorbar
%     axis off; axis square
    
    %% scatter of spike directions
    figure(fign);
    ax = subplot(fig_row(1),3,fig_row(2));
    hold on;
    
%     c_size = 5; % default 15;
    c_size = 10; % default 15;
    
    if strcmp(class(root), 'CMBHOME.Session')
        import CMBHOME.Utils.ContinuizeEpochs
        plot(ContinuizeEpochs(root.x),ContinuizeEpochs(root.y),'Color',[.7 .7 .7])
        colormap(ax,hsv)
        xlim([min(ContinuizeEpochs(root.x)) max(ContinuizeEpochs(root.x))]);
        ylim([min(ContinuizeEpochs(root.y)) max(ContinuizeEpochs(root.y))]);
        cx = ContinuizeEpochs(root.cel_x);
        cy = ContinuizeEpochs(root.cel_y);
        ch = ContinuizeEpochs(root.cel_headdir);
        scatter(cx,cy,c_size,ch,'filled')
    else
        plot(root.x,root.y,'Color',[.7 .7 .7])
        colormap(ax,flip(hsv))
        xlim([min(root.x) max(root.x)]); ylim([min(root.y) max(root.y)])
        cx = root.x(root.spike == 1);
        cy = root.y(root.spike == 1);
        ch = root.md(root.spike == 1);
        scatter(cx,cy,c_size,ch,'filled')
    end
%     scatter(out.QP(:,1),out.QP(:,2),30,'k','filled')
    set(ax,'YDir','Normal')
%     title('Trajectory')
    axis off
    axis square
    
    %% ratemap circular
    figure(fign);
    ax = subplot(fig_row(1),3,fig_row(2)+1);
    % the +pi/2 brings "forwards" to "up"
    [t2, r2] = meshgrid(wrapTo2Pi(out.params.thetaBins+pi/2), out.params.distanceBins(1:end-1));
    [x, y] = pol2cart(t2,r2);
    %     h=surface(x,y, out.rm); shading interp
    h = pcolor(x,y, out.rm); shading interp % YL: 20210913 7.22am
    
    hold on
    set(ax,'XTick',[],'YTick',[])
    
    colormap(ax, parula)
    set(ax, 'YDir','Normal','CLim',[0 prctile(out.rm(:), 99)])
    set(ax,'YDir','Normal')
    %     title('rm')
    axis off
    %     title(prctile(out.rm(:), 99))
    %     colorbar
    axis off; axis square
    
    % YL: 20210913 7.14am
    xl = xlim()
    yl = ylim()
    lv = plot([0,0], yl, 'LineWidth', 1, 'Color', 'w');
    lh = plot(xl, [0,0], 'LineWidth', 1, 'Color', 'w');
    lv.LineStyle = '--';
    lh.LineStyle = '--';
    
end
        