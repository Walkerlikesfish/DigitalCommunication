function [ht_est] = SISO_TD_estimator(r_symbols, flags)
%SISO_TD_estimator 
%   intputs:
%       -r_symbols: the bit seires coming out from the transmitter
%       -flags:
%   outputs:
%       -ht_est: [N_sub] preamble in frequency domain

N_preamble = 2;
% extract the preamble: always at the begining of the frame
preamble_t = r_symbols(1:(flags.N_subcarr+flags.N_cp)*2);
% CP^(-1) for preamble
preamble_t = preamble_t(flags.N_cp*N_preamble+1:end);
% ht_est = zeros(1,flags.N_subcarr);

for ii=1:N_preamble
    cur_one_preamble = preamble_t((ii-1)*flags.N_subcarr+1:ii*flags.N_subcarr);
    h_cur_est = cur_one_preamble \ flags.preamble_t;
    if ii == 1
        ht_est = h_cur_est;
    else
        ht_est = ht_est + h_cur_est;
    end
end

% estimate the H(f) by averaging the estimation from N_one_preamble
% preambles
ht_est = ht_est/N_preamble;
size(ht_est);


end

