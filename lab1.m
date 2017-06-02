% Digital Communication Lab 1
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
% Lab 1 
% Channel charateristic estimation 
% -> coherence frequency and PDP
%% Init
clear all
addpath('functions/');
addpath('misc/');

flags.N_line = 10; % how many points in a line 
flags.N_bins = 501; % the number of bins
flags.N_pts = 128; % the number of points on the figure drawing
flags.BW = 200 * 1e6; % 200M Hz
flags.dt = 1/flags.BW; 
flags.BW2_prop = 0.1; % the propotion for the pre_LPF

% selecting the filter: 0-Rectangular; 1-Tuckey(Recommend)
flags.LPFfilter = 1;

% Basic measurement parameters
flags.dist_int = 0.02;
flags.fc = 2.35e9;
flags.vc = 3e8;
flags.lamda = flags.vc/flags.fc;

% Load datas
data_LOS = load('LOS_2017.mat');
data_NLOS = load('NLOS_2017.mat');

disp('------------------------------------');
disp('Initialise finished!');
disp('------------------------------------');

%% Plot the local area - Narrowband
HtMat_NLOS = ifft_3dmat(data_NLOS.Data, flags.N_line, flags.N_bins);
% calculate PDP
[HtMatAmp_NLOS, ~] = calc_PDP(HtMat_NLOS, flags.N_line, flags.N_bins);
HtMatAmp_NLOS = sum(HtMatAmp_NLOS, 3);
HtMatAmp_NLOS = squeeze(sum(HtMatAmp_NLOS, 4));
ht_toplot_NLOS = squeeze(10*log10(HtMatAmp_NLOS));

HtMat_LOS = ifft_3dmat(data_LOS.Data, flags.N_line, flags.N_bins);
% calculate PDP
[HtMatAmp_LOS, ~] = calc_PDP(HtMat_LOS, flags.N_line, flags.N_bins);
HtMatAmp_LOS = sum(HtMatAmp_LOS, 3);
HtMatAmp_LOS = squeeze(sum(HtMatAmp_LOS, 4));
ht_toplot_LOS = squeeze(10*log10(HtMatAmp_LOS));

axis_in_lamda = (1:1:flags.N_line) .* flags.dist_int ./ flags.lamda;
figure(1); surf(axis_in_lamda,axis_in_lamda,ht_toplot_LOS); title('LOS h(t) (dB) in local area');
figure(2); surf(axis_in_lamda,axis_in_lamda,ht_toplot_NLOS); title('NLOS h(t) (dB) in local area');

%% 200Mhz Coherence Frequency estimation
[f_cohr_LOS, PDP_LOS] = calc_coherf_raw(data_LOS.Data, flags);
disp(['LOS 200Mhz BW, f_coherence: ',num2str(f_cohr_LOS), 'Hz']);

[f_cohr_NLOS, PDP_NLOS] = calc_coherf_raw(data_NLOS.Data, flags);
disp(['NLOS 200Mhz BW, f_coherence: ', num2str(f_cohr_NLOS), 'Hz'])

%% 20Mhz Coherence Frenquency estimation
flags.LPFfilter = 0;
[f_cohr_LOS_20Mrec, PDP_LOS_20Mrec] = calc_coherf_filter(data_LOS.Data, flags);
disp(['LOS 20Mhz BW(filtered with rec window), f_coherence: ',num2str(f_cohr_LOS_20Mrec), 'Hz']);
flags.LPFfilter = 1;
[f_cohr_LOS_20Mtuc, PDP_LOS_20Mtuc] = calc_coherf_filter(data_LOS.Data, flags);
disp(['LOS 20Mhz BW(filtered with rec window), f_coherence: ',num2str(f_cohr_LOS_20Mtuc), 'Hz']);

flags.LPFfilter = 0;
[f_cohr_NLOS_20Mrec, PDP_NLOS_20Mrec] = calc_coherf_filter(data_NLOS.Data, flags);
disp(['LOS 20Mhz BW(filtered with Tuckey window), f_coherence: ',num2str(f_cohr_NLOS_20Mrec), 'Hz']);
flags.LPFfilter = 1;
[f_cohr_NLOS_20Mtuc, PDP_NLOS_20Mtuc] = calc_coherf_filter(data_NLOS.Data, flags);
disp(['LOS 20Mhz BW(filtered with Tuckey window), f_coherence: ',num2str(f_cohr_NLOS_20Mtuc), 'Hz']);

%% Questions
% Observe the the limited BW, and fix this problem ? what problem ?
% 
