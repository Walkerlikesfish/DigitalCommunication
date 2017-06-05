function [arec_bits] = SISO_receiver(ach_bits, flags, verbal)
% Simulating the receiver side for a SISO channel, input the parameter
% group flags, generating a series of bits as the transmitter output
%
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
if nargin < 3
    verbal = 0;
end

%% Series->Parrallel
r_symbol_cp_t = reshape(ach_bits, flags.N_subcarr+flags.N_cp, []);
% CP^(-1)
r_symbol_t = r_symbol_cp_t(flags.N_cp+1:end,:);
len_p = size(r_symbol_t,2);

%% TD h(t) equalisation
if flags.tdEQ == 1
    h_tf = fft(flags.MPCTD, flags.N_subcarr);
    r_symbol_sp_te = zeros(size(r_symbol_t));
    for ii=1:len_p
        cur_t = r_symbol_t(:,ii);
        cur_tf = fft(cur_t);
        r_symbol_sp_te(:,ii) = cur_tf ./ h_tf;
    end
end

%% FD equalisation
% FFT
flags.tdEQ = 0;

for ii=1:len_p
    cur_t = r_symbol_t(:,ii);
    cur_f = fft(cur_t);
    r_symbol_sp(:,ii) = cur_f;
end


if flags.tdEQ == 0
    % Equalisation in FD
    if flags.MPCchoice == 0
        ht = flags.MPCht.h_LOS_rice;
        hf = fft(ht, flags.N_subcarr);
        % use estimation or use knowledge
        if flags.preamble_yes == 1
            hf = flags.MPCZF;
        end
    elseif flags.MPCchoice == 1
        ht = flags.MPCht.h_NLOS_rayl;
        hf = fft(ht, flags.N_subcarr);
        % use estimation or use knowledge
        if flags.preamble_yes == 1
            hf = flags.MPCZF;
        end
    elseif flags.MPCchoice == -1
        hf = ones(1,flags.N_subcarr);
    end
    hf = hf.';
    
    for ii=1:len_p
        cur_symbol = r_symbol_sp(:,ii);
        r_symbol_aE = cur_symbol ./ hf;
        r_symbol_sp(:,ii) = r_symbol_aE;
    end
end

%% Frequency compensation
ch_mid = flags.N_subcarr/2;
pilot_shift = 0;
if flags.CFO == 1
    for ii=1:len_p
        cur_symbol = r_symbol_sp(:,ii);
        % extract the pilot and estimate the CFO
        for ii2=1:4
            cur_shift = angle(cur_symbol(flags.pilot_channel(ii2)+ch_mid));
            pilot_shift = pilot_shift + cur_shift;
        end
        pilot_shift = pilot_shift/4;
        
        % correct the CFO using the pilot estimation
        cur_symbol = cur_symbol * exp(-1j*cur_shift);
        r_symbol_sp(:,ii) = cur_symbol;
    end
end

%% the rest
% P/S
if flags.tdEQ == 1
    r_symbol_sp = r_symbol_sp_te;
end

r_symbol_f = reshape(r_symbol_sp, 1, []);
r_symbol_R = real(r_symbol_f);
r_symbol_I = imag(r_symbol_f);

if verbal>0
    figure
    scatter(r_symbol_R,r_symbol_I);
end

bits_rx = demapping(r_symbol_f, flags.Nbps);
bits_rx = reshape(bits_rx, 1, []);
bits_rx = bits_rx.';

arec_bits = bits_rx;
end

