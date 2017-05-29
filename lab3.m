% Digital Communication Lab 3
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
% Lab 3 
% Simple system model: NO preamble, NO Synchronization, ONLY include the
% Channel equalization based on the known channel model
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
flags.EbN0 = -50:50;     % EbN0 interval
% [MPC settings]
ht = load('impulse_response.mat');
flags.MPCht = ht.ht;
flags.MPCchoice = -1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel
% [Channel Estimation]
flags.preamble_size = -1; % -1 for initialize meaning no preamble added
flags.preamble_yes = -1;
% [Receiver]
flags.tdEQ = 0; % switch for time domain equalisation : 0-OFF/ 1-ON
% [Time shifting]
flags.STO = 0; % switch for STO: 0-OFF / 1-ON
flags.timeshift = 1; % set the time shifting (unit:[bit])
flags.N_averageWindow = flags.N_cp*2;
% [CFO]
flags.CFO = 0; % switch for CFO: 0-OFF/ 1-ON
flags.f_tx = 10; % transmitter frequency shift unit[Hz]
flags.pilot_channel = [-21,-7,7,21];

%% [Transmitter]
[symbol_cp_s, bits_tx] = SISO_transmitter(flags);
%% [Channel+Receiver] scan throught the EbN0
flags.MPCchoice = 1; % the choice for MPC model: 0:LOS; 1:NLOS; -1:No channel

flags.AWGN = 1; 
BER = zeros(length(flags.EbN0),1);
disp('A little bit patience is required...')
for ii=1:length(flags.EbN0)
    % [Channel] AWGN
    flags.EbN0i = ii;
    flags.EbN0(ii)
    arec_bits = f_SISO_simulation(symbol_cp_s, flags);
    bits_rx = arec_bits;
    howcorrect=(bits_tx==bits_rx);        % check the original signal and the processed signal is equal or not
    BER(ii)=1-(sum(howcorrect)/flags.Nbits);    % Bit Error Rate (BER)
end

% figure
semilogy(flags.EbN0, BER,'-xg');
hold on;    
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs EbN0');
legend('Modulation 64QAM');
grid on   

%% Questions
% [1] why in practice we do not use the centre frequency <- cuz when
% demodulating, there will be a DC shift corrupting the 0Hz on the baseband
% [2] How to get the Frequency responce from the h[n] which has only 8 taps