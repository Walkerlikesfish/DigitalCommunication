function [symbol_cp_s, bit_tx] = SISO_transmitter(flags)
% Simulating the transmitter side for a SISO channel, input the parameter
% group flags, generating a series of bits as the transmitter output
%
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 

bit = randi([0 1],flags.Nbits,1);  % generate a length of N bits random signal
bit_tx= bit;

%% [Transmitter] 64QAM
symb = mapping(bit_tx,flags.Nbps,'qam');
symb_R = real(symb);
symb_I = imag(symb);

%% [Transmitter] S/P
symbol_para=reshape(symb, flags.N_subcarr, []);
% ifft the symbol parallel series
len_p = size(symbol_para,2);

%% [Transmitter] CP
% Cyclic Prefix(CP) is added to mitigate 2 problem:
% 1) The wideband channel has memory, which lead to possibly corruption
% between frames;
% 2) When carrying out the DFT at the receveiver end, the DFT input should
% be periodic, adding the CP making it periodical 
% => the periodic is refering to when the signal is passing through the channel(h(t)), the
% response produced is periodic.
%
symbol_cp_t = zeros(flags.N_subcarr+flags.N_cp, len_p);
for ii=1:len_p
    cur_freq = symbol_para(:,ii);
    % check if pilot need to be injected in certain channels
    if flags.CFO == 1
        ch_mid = flags.N_subcarr/2;
        for ii2=1:4
            cur_freq(flags.pilot_channel(ii2)+ch_mid)=1; % mark the pilot channel
        end
    end
    cur_t = ifft(cur_freq);
    cur_cp = cur_t(end-flags.N_cp+1:end); % extract CP
    cur_t = [cur_cp; cur_t];              % add CP before the symbol
    symbol_cp_t(:,ii) = cur_t;
end

%% [Transmitter]: P/S
symbol_cp_s = reshape(symbol_cp_t, 1, []);

end

