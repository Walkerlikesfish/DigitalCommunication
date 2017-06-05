% Digital Communication Lab 4
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
% Lab4 - Channel estimation
% NO CFO, NO SCO, only have to estimate the channel
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
flags.EbN0 = -30:50;     % EbN0 interval
% [MPC settings]
ht = load('impulse_response.mat');
flags.MPCht = ht.ht;
flags.MPCchoice = 0; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
% [Channel Estimation]
flags.preamble_size = -1; % -1 for initialize meaning no preamble added
flags.preamble_yes = 1; % [switch] for Preamble: 0-OFF / 1-ON [!]Prequisit for CFO and SCO
% [Receiver]
flags.tdEQ = 0; % [switch] for time domain equalisation : 0-OFF/ 1-ON
% [Time shifting]
flags.STO = 0; % [switch] for STO: 0-OFF / 1-ON
flags.timeshift = 1; % set the time shifting (unit:[bit])
flags.N_averageWindow = flags.N_cp*2;
% [CFO]
flags.CFO = 0; % [switch] for CFO: 0-OFF/ 1-ON
flags.f_tx = 0; % transmitter frequency shift unit[Hz]
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

%% [Channel+Receiver] scan throught the EbN0
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
flags.tdEQ = 1; % [switch] for time domain equalisation : 0-OFF/ 1-ON
flags.AWGN = 1; 
BER = zeros(length(flags.EbN0),1);
disp('A little bit patience is required...')
for ii=1:length(flags.EbN0)
    % [Channel] AWGN
    flags.EbN0i = ii;
    flags.EbN0(ii)
    [arec_bits] = f_SISO_simulation(symbol_cp_s, flags);
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
end

% decide if you want to overlap the figure or not
figure(3)
semilogy(flags.EbN0, BER,'-rx');
hold on;    
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs EbN0');
legend('Modulation 64QAM');
grid on   

%% [Test]
flags.tdEQ = 0;
flags.EbN0i = 50;
[ach_bits] = SISO_channel(symbol_cp_s, flags);
% [Receiver] 
% [Receiver - Estimate the channel]
[hf_est] = SISO_ZF_estimator(ach_bits, flags);
% estimating in time domain
[ht_est] = SISO_TD_estimator(ach_bits, flags);
flags.MPCZF = hf_est; % set the estimated channel to global var
flags.MPCTD = ht_est;

ach_bits = ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
% [Receiver - Equalisation and estimation]
[arec_bits] = SISO_receiver(ach_bits, flags, 1);


%% Questions
% [1] why in practice we do not use the centre frequency <- cuz when
% demodulating, there will be a DC shift corrupting the 0Hz on the baseband
% [2] How to get the Frequency responce from the h[n] which has only 8 taps