% Digital Communication Lab 5
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
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
flags.EbN0 = -20:30;     % EbN0 interval
% [MPC settings]
ht = load('impulse_response.mat');
flags.MPCht = ht.ht;
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
% [Channel Estimation]
flags.preamble_size = -1; % -1 for initialize meaning no preamble added
flags.preamble_yes = 1;
% [Receiver]
flags.tdEQ = 0; % switch for time domain equalisation : 0-OFF/ 1-ON
% [Time shifting]
flags.STO = 1; % switch for STO: 0-OFF / 1-ON
flags.timeshift = 1; % set the time shifting (unit:[bit])
flags.N_averageWindow = flags.N_cp*2;

%% [Transmitter]
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

%% [Channel+Receiver] scan throught the EbN0
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel

flags.AWGN = 1; 
BER = zeros(length(flags.EbN0),1);
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
    BER(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
end

% decide if you want to overlap the figure or not
figure(3)
semilogy(flags.EbN0, BER,'-gx');
hold on;    
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs EbN0');
legend('Modulation 64QAM');
grid on   

%% [Channel + receiver] scan through the STO
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
flags.AWGN = 1; 
flags.EbN0i = 40;
flags.STO = 1; % turn on the STO

disp('Scanning the STO, A little bit patience is required...')
N_max_STO = 10;
flags.MPCchoice = 0;
N_shift = [-50:50];
BER = zeros(length(N_shift),1);

%%
for ii=1:length(N_shift)
    % [STO] - iterate the timeshift
    flags.timeshift = -N_shift(ii);
    % [Channel] AWGN
    [ach_bits] = SISO_channel(symbol_cp_s, flags);
    % [Receiver] - Time acquisition
    [t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
    % [Receiver] - Channel Estimation 
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
    % [Receiver] - Equalisation
    [arec_bits] = SISO_receiver(ach_bits, flags, 0);
    % [Result]
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
    
    disp(['processed  ' num2str(ii) '/' num2str(length(N_shift))])
end

figure(1)
semilogy(N_shift, BER,'-gx');
hold on;    
xlabel('Shift samples (bit)');
ylabel('Bit Error Rate (BER)');
title('BER vs N_shift');
legend('Modulation 64QAM');
grid on   


%% Questions

%% backup
%% Test
% ii = 20;
% flags.timeshift = -50;
% [ach_bits] = SISO_channel(symbol_cp_s, flags);
% % [Receiver] - Time acquisition
% [t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
% % [Receiver] - Channel Estimation 
% if flags.preamble_yes == 1
%     % [Receiver - Estimate the channel]
%     % estimating in frequency domain
%     [hf_est] = SISO_ZF_estimator(ach_bits, flags);
%     % estimating in time domain
%     [ht_est] = SISO_TD_estimator(ach_bits, flags);
%     flags.MPCZF = hf_est;
%     flags.MPCTD = ht_est;
%     % extract the preamble from the signals
%     ach_bits = ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
% end
% % [Receiver] - Equalisation
% [arec_bits] = SISO_receiver(ach_bits, flags, 0);
% % [Result]
% bits_rx = arec_bits;
% howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not