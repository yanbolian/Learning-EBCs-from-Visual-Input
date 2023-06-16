%% plot EBCs using spike data
% This file is used in main.m
% Author: Yanbo Lian

%% Spatial rate map from spike data
positions = [root(1).x root(1).y] / 100; % unit: cm to m
input_positions = discretize_positions_into_bins(positions, env);

S_test_root = zeros(lca.num_cell, size(positions,1));
for i = 1 : lca.num_cell
    S_test_root(i,:) = root(i).spike;
end
S_test_root = double(S_test_root);

% Spatial rate map normalized by occupancy
rate_sum_map = input_positions * S_test_root'; % L * 100
occupancy_map = sum(input_positions, 2); % L * 1;
occupancy_map(occupancy_map==0) = nan;
rate_maps = rate_sum_map ./ repmat(occupancy_map, 1, lca.num_cell);
rate_maps(isnan(rate_maps)) = 0;
sigma_discrete = 0.03 / env.x_size * env.Nx; % 0.03; Assume x_size:y_size = Nx:Ny
rate_maps_smoothed = reshape(imgaussfilt(reshape(rate_maps,env.Ny,env.Nx,[]), sigma_discrete),env.L,[]);

fig = figure('Name', 'Spatial rate maps', 'numbertitle', 'off');
display_matrix(normalize_matrix(rate_maps_smoothed), 3);
fig.Units = 'normalized'; fig.Position = [0 0.5 0.3 0.3];
colormap(parula); h = colorbar; h.Limits=[0 1];

%% plot figures
n_rows = 10; % number of rows in each figure
i_fig = 10;
for i_cell = 1 : lca.num_cell
    if mod(i_cell, n_rows) == 1
        i_fig = i_fig + 1;
        fig = figure(i_fig); clf
        fig.Units = 'normalized'; fig.Position = [0 0 0.22 0.88]; % figures
        [ha, ~] = tight_subplot(n_rows, 3, [.01 .03],[.1 .01],[.01 .01]);
        i_plot = 1;
    end
    
    % Spatial ratemap
    ax = ha(i_plot); i_plot = i_plot + 1;
    axes(ax);
    rate_map = reshape(rate_maps_smoothed(:,i_cell),env.Ny,env.Nx);
    rate_map = imresize(rate_map, 3, 'lanczos3');
    rate_map = rate_map / max(rate_map(:));
    imagesc(rate_map);
    colormap(ax, 'parula');
    axis image
    axis off
    set(ax, 'YDir', 'normal')
    %     title(['cell #' num2str(i_cell)])
    
    % EBC plot: trajectory plot & egocentric rate map
    if sum(root(i_cell).spike) > 0
        out = EgocentricRatemap(root(i_cell)); % where r is your behavioral/ephys struct
        i_plot = plotEBC_tight(root(i_cell), out, ha, i_plot);
    else
        i_plot = i_plot + 2;
    end
end