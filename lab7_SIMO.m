% Digital Communication Lab 7 SIMO
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% %% Parameters
%% setting of the basic parameters
clear all
close all;
% [Basic Settings]
flags.Nbits = 1024*60;   % number of total bits ready to send
flags.f_c = 2.35e9;      % carrier frequency = 2.35Ghz
flags.BW = 20e6;         % Bandwidth = 20MHz
flags.Nbps = 6;          % 2^6=64 QAM  
flags.N_subcarr = 64;    % number of sub carriers
flags.N_cp = 16;         % length of Cyclic prefix length
% [AWGN Settings]
flags.AWGN = 1;          % shall we turn on the AWGN channel ? 1-yes;0-no
flags.EbN0 = -20:300;     % EbN0 interval
% [MPC settings] - Load the multi-receiver settings
ht_MO = load('h_mpc');
flags.MPCht_LOS = ht_MO.mpc_h{1};
flags.MPCht_NLOS = ht_MO.mpc_h{2};
flags.MPCchoice = 1; % the choice for MPC model: 1:LOS; 2:NLOS; -1:No channel

ht = load('impulse_response.mat');
flags.MPCht = ht.ht;
% [Channel Estimation]
flags.preamble_size = -1; % -1 for initialize meaning no preamble added
flags.preamble_yes = 1;
% [Receiver]
flags.tdEQ = 0; % switch for time domain equalisation : 0-OFF/ 1-ON
% [Time shifting]
flags.STO = 1; % switch for STO: 0-OFF / 1-ON
flags.timeshift = 1; % set the time shifting (unit:[bit])
flags.N_averageWindow = flags.N_cp*2;
% [CFO]
flags.CFO = 1; % switch for CFO: 0-OFF/ 1-ON
flags.f_tx = 10; % transmitter frequency shift unit[Hz]
flags.pilot_channel = [-21,-7,7,21];
% [SIMO]
% the selected_channels indicate which channel(receiver) is activated in
% the simulation. As the receiver matrix is 10x10x10, we use 3 decimal
% digits to indicate the id of the receiver
flags.selected_channels = [111,131,142,153,164,175];

%% [Transmitter] - Generate the source messages

[symbol_cp_s, bits_tx] = SISO_transmitter(flags);
% adding preamble to the tranmitting signals
if flags.preamble_yes == 1
    [preamble_F, flags.preamble_t, trans_bits_ap] = SISO_adding_preamble(symbol_cp_s, flags);
    % update the preamble information into flags
    flags.preamble_size = size(preamble_F);
    flags.preamble_F = preamble_F;
    % update the transmitting signal
    symbol_cp_s = trans_bits_ap; 
end

%% [Channel + receiver] Setup 
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
flags.AWGN = 1; 
flags.EbN0i = 50;
flags.STO = 1; % turn on the STO
flags.timeshift = 1; % set STO 
flags.f_tx = 50; % transmitter frequency shift unit[Hz]

%% [Channel + receiver] iterate through CFO
scan_CFO = [0:1000:50000];
for ii=1:50
    flags.f_tx = scan_CFO(ii);
    
    [ach_bits_mo] = SIMO_channel(symbol_cp_s, flags);
    ach_bits = ach_bits_mo(1,:); 
    % [Receiver] - Time acquisition
    [t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
    % [Receiver] - CFO acquisition
    [df_est, ach_bits] = SISO_estimate_CFO(ach_bits, flags);
    % [Receiver] - Channel estimation independently for each receiver
    n_channels = length(flags.selected_channels);
    hf_est_mo = [];
    ach_bits_mo_np = [];  % the bits stream after clipping the preamables
    for ic=1:n_channels
        cur_cid = flags.selected_channels(ic);
        id_channel = [floor(cur_cid/100), floor(mod(cur_cid,100)/10), mod(cur_cid,10)];
        % estimate the channel function
        cur_ach_bits = ach_bits_mo(ic,:);
        [cur_hf_est] = SISO_ZF_estimator(cur_ach_bits, flags);
        hf_est_mo = [hf_est_mo; cur_hf_est];
        % extract the singal clipping the head
        cur_ach_bits = cur_ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
        ach_bits_mo_np = [ach_bits_mo_np; cur_ach_bits];
    end
    % [Receiver] - Equalizer
    [arec_bits] = SIMO_receiver(ach_bits_mo_np, hf_est_mo, flags, 0);

    % [Result]
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER(ii) = 1-(sum(howcorrect)/flags.Nbits);
    
    disp(['processed  ' num2str(ii) '/' num2str(length(scan_CFO))])
end

figure(2)
semilogy(scan_CFO(1:50), BER,'-gx');
hold on;    
xlabel('CFO (Hz)');
ylabel('Bit Error Rate (BER)');
title('BER vs CFO');
legend('Modulation 64QAM');
grid on   

%% [EbN0]
flags.EbN0 = -20:30;     % EbN0 interval

flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel

flags.AWGN = 1; 
BER_SIMO = zeros(length(flags.EbN0),1);
disp('A little bit patience is required...')
for ii=1:length(flags.EbN0)
    % [Channel] AWGN
    flags.EbN0i = ii;
    flags.EbN0(ii)
    [ach_bits_mo] = SIMO_channel(symbol_cp_s, flags);
    ach_bits = ach_bits_mo(1,:); 
    % [Receiver] - Time acquisition
    [t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
    % [Receiver] - CFO acquisition
    [df_est, ach_bits] = SISO_estimate_CFO(ach_bits, flags);
    % [Receiver] - Channel estimation independently for each receiver
    n_channels = length(flags.selected_channels);
    hf_est_mo = [];
    ach_bits_mo_np = [];  % the bits stream after clipping the preamables
    for ic=1:n_channels
        cur_cid = flags.selected_channels(ic);
        id_channel = [floor(cur_cid/100), floor(mod(cur_cid,100)/10), mod(cur_cid,10)];
        % estimate the channel function
        cur_ach_bits = ach_bits_mo(ic,:);
        [cur_hf_est] = SISO_ZF_estimator(cur_ach_bits, flags);
        hf_est_mo = [hf_est_mo; cur_hf_est];
        % extract the singal clipping the head
        cur_ach_bits = cur_ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
        ach_bits_mo_np = [ach_bits_mo_np; cur_ach_bits];
    end
    % [Receiver] - Equalizer
    [arec_bits] = SIMO_receiver(ach_bits_mo_np, hf_est_mo, flags, 0);
    % [Result]
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER_SIMO(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
end

% decide if you want to overlap the figure or not
figure(3)
semilogy(flags.EbN0, BER_SIMO,'-ro');
hold on;    
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs EbN0');
legend('Modulation 64QAM');
grid on   

%% SISO

BER_SISO = zeros(length(flags.EbN0),1);
disp('A little bit patience is required...')
for ii=1:length(flags.EbN0)
    % [Channel] AWGN
    flags.EbN0i = ii;
    flags.EbN0(ii)
    [ach_bits] = SISO_channel(symbol_cp_s, flags);
    % [Receiver] 
    if flags.preamble_yes == 1
        % [Receiver - Estimate the channel]
        % estimating in frequency domain
        [hf_est] = SISO_ZF_estimator(ach_bits, flags);
        % estimating in time domain
        [ht_est] = SISO_TD_estimator(ach_bits, flags);
        flags.MPCZF = hf_est;
        flags.MPCTD = ht_est;
        % extract the preamble from the signals
        ach_bits = ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
    end
    % [Receiver - Equalisation and estimation]
    [arec_bits] = SISO_receiver(ach_bits, flags, 0);
    % [Result]
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER_SISO(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
end

% decide if you want to overlap the figure or not
figure(3)
semilogy(flags.EbN0, BER_SISO,'-gx');
hold on;    
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs EbN0');
legend('Modulation 64QAM');
grid on   

%% [Test] 

% Channel
flags.EbN0i = 320;
flags.f_tx = 100; % transmitter frequency shift unit[Hz]
[ach_bits_mo] = SIMO_channel(symbol_cp_s, flags);

% Reciever - CTO and CFO acquisition from preamble
% always take the 1st channel as an example channel
% since the local oscillator is the same even for the whole matrix
ach_bits = ach_bits_mo(1,:); 
% [Receiver] - Time acquisition
[t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
% [Receiver] - CFO acquisition
[df_est, ach_bits] = SISO_estimate_CFO(ach_bits, flags);

% [Receiver] - Channel estimation independently for each receiver
n_channels = length(flags.selected_channels);
hf_est_mo = [];
ach_bits_mo_np = [];  % the bits stream after clipping the preamables
for ic=1:n_channels
    cur_cid = flags.selected_channels(ic);
    id_channel = [floor(cur_cid/100), floor(mod(cur_cid,100)/10), mod(cur_cid,10)];
    % estimate the channel function
    cur_ach_bits = ach_bits_mo(ic,:);
    [cur_hf_est] = SISO_ZF_estimator(cur_ach_bits, flags);
    hf_est_mo = [hf_est_mo; cur_hf_est];
    % extract the singal clipping the head
    cur_ach_bits = cur_ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
    ach_bits_mo_np = [ach_bits_mo_np; cur_ach_bits];
end

% [Receiver] - Equalizer
[arec_bits] = SIMO_receiver(ach_bits_mo_np, hf_est_mo, flags, 1);



