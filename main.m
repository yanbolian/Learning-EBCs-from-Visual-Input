%% This script implements V1-RSC model of Lian et al. 2023 J Neurosci paper
% Author: Yanbo Lian

clear
close all
clc

addpath(genpath('functions'))

%% Parameters of the environment and trajectory
env.x_size =  1.25; 1.265; % unit: m
env.y_size =  1.25; 1.265; % unit: m
env.Nx = 100; % number of discrete points that represent the environment
env.Ny = 100; % number of discrete points that represent the environment
env.L = env.Nx * env.Ny;
env.delta_x = env.x_size / (env.Nx-1);
env.delta_y = env.y_size / (env.Ny-1);

% Trajectory
trajectory.data_file_name = '1.25m170degreeSimulatedTrajectory.mat';
trajectory.option = 'virtual rat';
trajectory.dt = 30; % ms

%% load data
NOTE = 'The learning rate is 10% of the original value for the last 25% of the trajectory';
load(trajectory.data_file_name, 'images_data', 'hds_data', 'positions_data', 'images_V1', 'lgn', 'Nx', 'Ny') % wall on the left

images_data = double(images_data) / 255; % scale visual input to 0-1
positions_data = positions_data / 1000; % unit: mm to m

Nt = size(images_V1,1); % number of different positions
num_V1_cell = size(images_V1, 2);

%% Create model RSC cells (parameters of LCA that implements sparse coding)
lca.num_cell = 100; % number of cells in the network
lca.lambda = 0; % firing threshold that controls sparsity
lca.thresh_type = 'soft-non-negative'; % implements non-negative sparse coding
lca.A_eta = 0.1; % learning rate of connection A
lca.s_max = 30; % maximal firing rate
lca.history_flag = 1;
lca.tau = 10; % unit: ms; membrane time constant of model dynamics
lca.dt = 0.5; % unit: ms; each time step
lca.U_eta = lca.dt / lca.tau; % update rate of membrane potential
lca.n_iter = trajectory.dt / lca.dt; % number of iterations required at each position of the trajectory
lca.display_every = 3000; 30; % unit: ms; frequency of generating plot
lca.A = normalize_matrix(rand(num_V1_cell, lca.num_cell)); % Connections between V1 and RSC
lca.n_epoch = Nt-1; % number of training epoch

U_rsc = randn(lca.num_cell, 1); % Membrane potentials of M neurons for batch_size images
S_rsc = max(U_rsc - lca.lambda, 0);% Firing rates (Response) of M neurons for batch_size images

%% main loop of learning V1-RSC connection
display_flag = 0;
A_eta = lca.A_eta;
fig = figure(1);fig.Units = 'normalized'; fig.Position = [0 0 0.8 0.8];

for i_epoch = 2 : lca.n_epoch+1
    
    if i_epoch == floor(0.75 * lca.n_epoch)
        A_eta = 0.1 * A_eta;
    end
    
    % positions of training trajectories
    r = positions_data(i_epoch,:);
    hd = mod(hds_data(i_epoch)+90, 360); % the HD of getting to current position
    disp = -(positions_data(i_epoch-1,:)-r)/trajectory.dt*1000; % ms to s
    disp_vec = disp(1)+disp(2)*1i;
    v = abs(disp_vec); % the speed of getting to the current position
    
    S_v1 = images_V1(i_epoch,:)'; % response of V1 cells
    img = reshape(images_data(i_epoch, :), Nx, Ny)';
    
    [S_rsc, U_rsc, S_rsc_his, ~] = sparse_coding_by_LCA(...
        S_v1, lca.A, lca.lambda, lca.thresh_type, lca.U_eta, lca.s_max, lca.n_iter, lca.history_flag, S_rsc, U_rsc);
    
    % update weights based on the learning rule
    R = S_v1 - lca.A * S_rsc; % calculate residual error
    dA1 = R * S_rsc';
    lca.A = lca.A + A_eta * dA1;
    lca.A = max(lca.A, 0); % A is kept non-negative
    lca.A = normalize_matrix(lca.A, 'L2 norm', 1); % Normalize each column of the connection matrix
    
    
    if (mod(i_epoch-1, lca.display_every/trajectory.dt) == 0)
        fprintf('Learning: %3.1f s\n',i_epoch*trajectory.dt/1000);
        display_flag = 1;
    end
    
    % Plot figures during learning
    if (display_flag == 1)
        display_flag = 0;
        figure(1); clf
        
        % plot the current position, head direction (HD) and movement direction (MD)
        subplot 331
        v = abs(disp_vec); % the speed of getting to current position
        md = mod(rad2deg(angle(disp_vec)),360);
        arrow_length = 0.5 * v;
        q = quiver(r(1),r(2), arrow_length*cos(deg2rad(md)), arrow_length*sin(deg2rad(md)));
        q.MaxHeadSize = 1;
        q.Marker = '.';
        q.MarkerSize = 15;
        q.LineWidth = 1;
        hold on;
        q = quiver(r(1),r(2), arrow_length*cos(deg2rad(hd)), arrow_length*sin(deg2rad(hd)));
        q.MaxHeadSize = 1;
        q.Marker = '.';
        q.MarkerSize = 15;
        q.LineWidth = 1;
        xlim([-0.1 env.x_size+0.1])
        ylim([-0.1 env.y_size+0.1])
        xlabel('x (m)')
        ylabel('y (m)')
        axis square
        set(gca, 'YDir', 'normal')
        set(gca, 'XTick', [0,env.x_size/2,env.x_size], 'XTickLabels', [0,env.x_size/2,env.x_size])
        set(gca, 'YTick', [0,env.y_size/2,env.y_size], 'YTickLabels', [0,env.y_size/2,env.y_size])
        grid on
        legend({'HD','MD'});
        title(['(' num2str(r(1),'%.2f') ',' num2str(r(2),'%.2f') ') ' ...
            num2str(v,'%.2f') 'm/s HD=' num2str(hd,'%.1f') '\circ']);
        
        % plot the current view
        subplot(334)
        imagesc(img,[0 1]);
        colormap(gca, 'gray');
        title('current view')
        axis image
        axis off
        
        % LGN-processed image
        subplot(337)
        img_lgn = imfilter(img,lgn.Ic-lgn.Is)./imfilter(img,lgn.Id);
        img_lgn = img_lgn((1+lgn.size)/2:(Ny-(lgn.size-1)/2),(1+lgn.size)/2:(Nx-(lgn.size-1)/2));
        imagesc(img_lgn);
        colormap(gca, 'gray');
        title('LGN-processed image')
        axis image
        axis off
        
        % Membrane potential and firing rates of RSC cells
        subplot(3,3,5);stem(S_rsc);
        ylabel('S: firing rates')
        subplot(3,3,8);stem(U_rsc);
        ylabel('U: membrane potentials')
        xlabel('RSC cell index')
        
        %
        if lca.history_flag == 1
            subplot(3,3,2);
            plot(S_rsc_his);
            xlabel('Iterations')
            xlim([1 lca.n_iter])
            title('Trajectory of RSC cell firing rates')
        end
        
        % V1-RSC connection
        subplot 133
        display_matrix(lca.A);
        title('A: V1-RSC connection');
        colormap(gray);colorbar
        pause(0.2)
    end
end

%% Spatial rate map normalized by occupancy
epsilon = 1e-16; % A small value to avoid zero division
N_test = Nt - 1; % get rid of the first image
env.Nx = 43; 42; % If the box is 125cm wide, it's divided into 3cm*3cm bins
env.Ny = 43; 42;
env.bin_width = 0.03;
env.L = env.Nx * env.Ny;
input_positions = discretize_positions_into_bins(positions_data(2:Nt,:), env);

U_rsc = randn(lca.num_cell, 1); % membrane potentials of M neurons for batch_size images
S_rsc = max(U_rsc - lca.lambda, 0); % firing rates (Response) of M neurons for batch_size images
S_test = zeros(lca.num_cell,N_test);

for i_test = 1 : N_test
    
    S_v1 = images_V1(i_test+1,:)';
    [S_rsc, U_rsc, S_his, ~] = sparse_coding_by_LCA(...
        S_v1, lca.A, lca.lambda, lca.thresh_type, lca.U_eta, lca.s_max, lca.n_iter, lca.history_flag, S_rsc, U_rsc);
    
    S_test(:,i_test) = S_rsc;
    
    if (mod(i_test, lca.display_every/trajectory.dt) == 0)
        fprintf('Testing: %3.1f s\n', i_test*trajectory.dt/1000);
    end
end

rate_sum_map = input_positions * S_test'; % L * 100
occupancy_map = sum(input_positions, 2); % L * 1;
occupancy_map(occupancy_map==0) = nan;
lca.rate_maps = rate_sum_map ./ repmat(occupancy_map, 1, lca.num_cell);
lca.rate_maps(isnan(lca.rate_maps)) = 0;

% smooth the rate map
sigma_discrete = 0.03 / env.x_size * env.Nx; % smoothing using with 3cm Gaussian filter
lca.rate_maps_smoothed = reshape(imgaussfilt(reshape(lca.rate_maps,env.Ny,env.Nx,[]), sigma_discrete),env.L,[]);

%% Generate spike data of the model that will be used to plot egocentric plots
positions = positions_data(2:end,:) * 100; % unit: cm
hds = hds_data(2:end,1) + 90; % add 90 degrees so that east is 0 degrees
delta_t = trajectory.dt / 1000; % units: second
N = length(hds);

% Set the peak firing rate of the whole population
peak_rate = 30;
S_test = peak_rate / max(S_test(:)) * S_test;

for i_cell = 1 : lca.num_cell
    root(i_cell).x = positions(:,1); % unit: cm
    root(i_cell).y = positions(:,2); % unit: cm
    root(i_cell).ts = delta_t * (1:N)'; %
    root(i_cell).md = deg2rad(mod(hds,360)); % movement (or head) direction: in radians
    root(i_cell).spike = ((S_test(i_cell,:)*delta_t) > rand(1,N))'; % Rate -> spikes
end

%% Save results
save('results.mat', 'env', 'lca', 'trajectory', 'S_test', 'root', 'NOTE');

%% Plot EBCs from spike data: spatial rate maps, spike plots and egocentric rate maps
plot_EBCs;