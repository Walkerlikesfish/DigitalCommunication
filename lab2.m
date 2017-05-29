% Digital Communication Lab 2
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
% Lab 2 
% Channel model statistical estimation
%% Init
clear all

flags.N_line = 10; % how many points in a line 
flags.N_bins = 501; % the number of bins
flags.N_pts = 128; % the number of points on the figure drawing
flags.BW = 200 * 1e6; % 200M Hz
flags.dt = 1/flags.BW; 
flags.BW2_prop = 0.1; % the propotion for the pre_LPF

flags.tap_max = 6;

% selecting the filter: 0-Rectangular; 1-Tuckey(Recommend)
flags.LPFfilter = 1;

% Load datas
data_LOS = load('LOS_2017.mat');
data_NLOS = load('NLOS_2017.mat');

disp('------------------------------------');
disp('Initialise finished!');
disp('------------------------------------');
%% Narrowband Model
% in narrow band, T_s is too large to distinguish different MPC -> 1 tap
% only
% LOS
% [?]LOS should be Rician fading -> ch5.p2 -> how to model estimate
[rayl_fit_nb_los, rice_fit_nb_los] = fit_narrowb(data_LOS.Data, flags);

% NLOS - should be Rayl fading
% [?] why the rice.K do not goes to 0 -> rayl fading
[rayl_fit_nb_nlos, rice_fit_nb_nlos] = fit_narrowb(data_NLOS.Data, flags);

%% Wideband Model
% in wideband model the MPC can be distinguished, thus it is a channel with
% memory. We have to get the statistical parameters for each of the taps

% LOS - Rician fading
[rayl_fit_wb_los, rice_fit_wb_los] = fit_wideb(data_LOS.Data, flags);

%%
% NLOS - rayl fading [?]
% [?] why there is the error when fiting with Rician 
[rayl_fit_wb_nlos] = fit_wideb_rayl(data_NLOS.Data, flags);

%% Generate the channel
% using the stastical to generate the channel function
%

ht.h_LOS_rice = gen_channel(rice_fit_wb_los, flags, 0);
ht.h_NLOS_rayl = gen_channel(rayl_fit_wb_nlos, flags, 1);
save('impulse_response.mat','ht');


%% Questions:
% [1] Narrowband model: how to model LOS ch5.p2
% [2] Why rice.K donot -> 0 for NLOS narrowband
% [3] is it same distribution wideband that in each tap, NLOS->rayl,
% LOS->rice, if so why there is error for fitting with Rician
