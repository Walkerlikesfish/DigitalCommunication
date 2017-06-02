% Digital Communication Lab 5
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
%% setting of the basic parameters
clear all
addpath('functions/');
addpath('misc/');
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
% [CFO]
flags.CFO = 0; % switch for CFO: 0-OFF/ 1-ON
flags.f_tx = 10; % transmitter frequency shift unit[Hz]
flags.pilot_channel = [-21,-7,7,21];


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

%% [Channel + receiver] scan through the STO
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
flags.AWGN = 1; 
flags.EbN0i = 50;
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
    arec_bits = f_SISO_simulation(symbol_cp_s, flags);
    % [Result]
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
    
    disp(['processed  ' num2str(ii) '/' num2str(length(N_shift))])
end

figure(1)
semilogy(N_shift, BER,'-ro');
hold on;    
xlabel('Shift samples (bit)');
ylabel('Bit Error Rate (BER)');
title('BER vs N_shift');
legend('Modulation 64QAM');
grid on   


%% Questions

%% backup
%% Test
% flags.EbN0i = 300;
% flags.timeshift = 10;
% [ach_bits] = SISO_channel(symbol_cp_s, flags);
% % [Receiver] - Time acquisition
% [t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
% t_est
% figure(1); hold on; plot(abs(ach_bits(1000:1200)), 'b-x');
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
