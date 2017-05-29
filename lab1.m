% Digital Communication Lab 1
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
% Lab 1 
% Channel charateristic estimation 
% -> coherence frequency and PDP
%% Init
clear all

flags.N_line = 10; % how many points in a line 
flags.N_bins = 501; % the number of bins
flags.N_pts = 128; % the number of points on the figure drawing
flags.BW = 200 * 1e6; % 200M Hz
flags.dt = 1/flags.BW; 
flags.BW2_prop = 0.1; % the propotion for the pre_LPF

% selecting the filter: 0-Rectangular; 1-Tuckey(Recommend)
flags.LPFfilter = 1;

% Load datas
data_LOS = load('LOS_2017.mat');
data_NLOS = load('NLOS_2017.mat');

disp('------------------------------------');
disp('Initialise finished!');
disp('------------------------------------');

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
