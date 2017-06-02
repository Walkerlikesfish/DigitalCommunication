% Digital Communication Lab 7 MPC estimation
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
%% Parameters
clear all
addpath('functions/');
addpath('misc/');
% ------------ basics --------------- %
flags.N_line = 10; % how many points in a line 
flags.N_bins = 501; % the number of bins
flags.N_pts = 128; % the number of points on the figure drawing
flags.BW = 200 * 1e6; % 200M Hz
flags.dt = 1/flags.BW; 
flags.BW2_prop = 0.1; % the propotion for the pre_LPF
% ------- carrier frequency --------- %
flags.f_c = 2.35e9;      % carrier frequency = 2.35Ghz
flags.tap_max = 8;
% selecting the filter: 0-Rectangular; 1-Tuckey(Recommend)
flags.LPFfilter = 1;
% -------- spatial resolution -------- %
flags.theta = [0:pi/20:pi];
flags.phi = [0:pi/40:2*pi];

%% Load datas
% Memo
% the data is stored in the way X-Y(9:-1:0)-Z
data_LOS = load('LOS_2017.mat');
data_NLOS = load('NLOS_2017.mat');

disp('------------------------------------');
disp('Initialise finished!');
disp('------------------------------------');

%% Get the Wideband domain - LOS
[HtMat_LOS] = get_wideBand_reduceBW(data_LOS.Data, flags);
a_n_mat_LOS = calc_a_n(HtMat_LOS, flags);
a_n_toplot = 20*log10(abs(a_n_mat_LOS(:,:,1)));
index = findLocalMaxima(a_n_toplot, -50);
plotImage(flags.phi, flags.theta, a_n_toplot, [-40; -100]);
hold on; plot(flags.phi(index(2)), flags.theta(index(1)), 'r*');
title('LOS angular amplitude'); xlabel('phi'); ylabel('theta');
disp('LOS model');
disp('------------------------------------');

%% Get the Wideband domain - NLOS
[HtMat_NLOS] = get_wideBand_reduceBW(data_NLOS.Data, flags);
a_n_mat_NLOS = calc_a_n(HtMat_NLOS, flags);
a_n_toplot = 20*log10(abs(a_n_mat_NLOS(:,:,1)));
index = findLocalMaxima(a_n_toplot, -50);
plotImage(flags.phi, flags.theta, a_n_toplot, [-40; -100]);
title('NLOS angular amplitude'); xlabel('phi'); ylabel('theta');
disp('NLOS model');
disp('------------------------------------');

%% Generate the h_i(n)
h_mpc_LOS = gen_channel_model(a_n_mat_LOS,flags);
h_mpc_NLOS = gen_channel_model(a_n_mat_NLOS,flags);
mpc_h = {h_mpc_LOS, h_mpc_NLOS};
save('h_mpc.mat','mpc_h');

disp('Generating the h_i');
disp('------------------------------------');

%% Spatial correlation Z direction
cur_tap = 1;
[Rdz_LOS] = calc_spatial_corr_Z(a_n_mat_LOS, flags, cur_tap);
[Rdz_NLOS] = calc_spatial_corr_Z(a_n_mat_NLOS, flags, cur_tap);

%% Spatial correlation X & Y direction
[Rdx_LOS, Rdy_LOS] = calc_spatial_corr_XY(a_n_mat_LOS, flags, cur_tap);
[Rdx_NLOS, Rdy_NLOS] = calc_spatial_corr_XY(a_n_mat_NLOS, flags, cur_tap);

%% Questions
%[1] Beamforming funciton explain

%%
% [mat_beta] = construct_beta(flags);
% [mat_B_LOS] = construct_B(mat_beta, flags);
% [mat_an_LOS] = construct_an(mat_B_LOS, HtMat_LOS, flags);
% plotImage(flags.phi, flags.theta, 20*log10(abs(mat_an_LOS(:,:,1))), [-40; -80])