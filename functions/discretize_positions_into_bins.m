function r_d = discretize_positions_into_bins(r, env)
% 
% r: Nt*2; unit:m
% env: parameters of the environment

delta_x = env.bin_width;
delta_y = env.bin_width;

N = size(r, 1);
r_d = zeros(env.Nx*env.Ny, N);

for i = 1 : N
    r_temp = zeros(env.Ny, env.Nx);
    r_temp(floor(r(i,2)/delta_y)+1, floor(r(i,1)/delta_x)+1) = 1;
    r_d(:,i) = reshape(r_temp, numel(r_temp), 1);
end