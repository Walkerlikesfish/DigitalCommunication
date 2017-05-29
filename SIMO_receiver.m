function [arec_bits] = SIMO_receiver(ach_bits_mo, hf_est_mo, flags, verbal)
% Simulating the receiver side for a SISO channel, input the parameter
% group flags, generating a series of bits as the transmitter output
%
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 
if nargin < 4
    verbal = 0;
end

%% Series->Parrallel

n_recs = size(ach_bits_mo);
n_recs = n_recs(1); % get the number of receivers
r_symbol_cpt_t = cell(1,n_recs);
r_symbol_cpt_f = cell(1,n_recs);

for ic = 1:n_recs
    cur_ach_bits = squeeze(ach_bits_mo(ic, :));
    cur_r_symbol_cp = reshape(cur_ach_bits, flags.N_subcarr+flags.N_cp, []);
    cur_r_symbol_t = cur_r_symbol_cp(flags.N_cp+1:end,:);
    len_p = size(cur_r_symbol_t,2);
    r_symbol_cpt_t{ic} = cur_r_symbol_t;
    cur_r_symbol_f = fft(cur_r_symbol_t);
    r_symbol_cpt_f{ic} = cur_r_symbol_f;
end

%% TD h(t) equalisation
if flags.tdEQ == 1
    disp('Sorry, do not support Time domain equalisation for SIMO system!')
end

%% FD equalisation
if flags.tdEQ == 0
    r_symbol_sp = zeros(flags.N_subcarr, len_p);
    % Equalisation in FD
    hf_est_mo_conj = conj(hf_est_mo);
    hf_est_mo_sum = zeros(1, flags.N_subcarr);
    for ic=1:n_recs
        cur_hf = abs(hf_est_mo(ic,:)).^2;
        hf_est_mo_sum = hf_est_mo_sum + cur_hf;
    end
    
    for ii=1:len_p % iterate through the frames
        cur_symbol_aEs = zeros(1, flags.N_subcarr);
        for ic=1:n_recs % sum through the channels(recvs)
            cur_symbol = r_symbol_cpt_f{ic}(:,ii);
            cur_symbol_aE = cur_symbol .* (hf_est_mo_conj(ic,:)).';
            cur_symbol_aEs = cur_symbol_aEs + cur_symbol_aE.';
        end
        cur_symbol_aEs = cur_symbol_aEs./hf_est_mo_sum;
        r_symbol_sp(:, ii) = cur_symbol_aEs;
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

