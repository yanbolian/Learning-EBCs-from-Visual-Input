function i_plot = plotEBC_tight(root,out,ha,i_plot)
%
% Downloaded from the Github of Hasselmo's lab
% Modified by Yanbo Lian for the purpose of displaying figures

    %% scatter of spike directions
    ax = ha(i_plot); i_plot = i_plot + 1;
    axes(ax);
    hold on;
    
    c_size = 6; % default 15;
%     c_size = 6; % default 15;
    
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
        xlim([min(root.x) max(root.x)]); ylim([min(root.y) max(root.y)]);
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
    ax = ha(i_plot); i_plot = i_plot + 1;
    axes(ax);

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
    xl = xlim();
    yl = ylim();
    lv = plot([0,0], yl, 'LineWidth', 1, 'Color', 'w');
    lh = plot(xl, [0,0], 'LineWidth', 1, 'Color', 'w');
    lv.LineStyle = '--';
    lh.LineStyle = '--';
    
end
        