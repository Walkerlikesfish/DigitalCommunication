function [hf_est] = SISO_ZF_estimator(r_symbols, flags)
%SISO_ZF_estimator 
%   intputs:
%       -r_symbols: the bit seires coming out from the transmitter
%       -flags:
%   outputs:
%       -hf_est: [N_sub] preamble in frequency domain

N_preamble = 2;
% extract the preamble: always at the begining of the frame
preamble_t = r_symbols(1:(flags.N_subcarr+flags.N_cp)*2);
% CP^(-1) for preamble
preamble_t = preamble_t(flags.N_cp*N_preamble+1:end);
hf_est = zeros(1,flags.N_subcarr);

for ii=1:N_preamble
    cur_one_preamble = preamble_t((ii-1)*flags.N_subcarr+1:ii*flags.N_subcarr);
    cur_one_preamble_f = fft(cur_one_preamble); % r_hat
    h_cur_est = cur_one_preamble_f ./ flags.preamble_F;
    hf_est = hf_est + h_cur_est;
end

% estimate the H(f) by averaging the estimation from N_one_preamble
% preambles
hf_est = hf_est/N_preamble;

end

