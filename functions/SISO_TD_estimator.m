function [ht_est] = SISO_TD_estimator(r_symbols, flags)
%SISO_TD_estimator 
%   intputs:
%       -r_symbols: the bit seires coming out from the transmitter
%       -flags:
%   outputs:
%       -ht_est: [N_sub] preamble in frequency domain
N_taps = 6;
N_preamble = 2;
% extract the preamble: always at the begining of the frame
preamble_t = r_symbols(1:(flags.N_subcarr+flags.N_cp)*2);
% CP^(-1) for preamble
preamble_t = preamble_t(flags.N_cp*N_preamble+1:end);

% Conv matrix construction
preamble_mat = zeros(flags.N_subcarr, flags.N_subcarr);
cur_line = flags.preamble_t.';
for ii = 1:flags.N_subcarr
    preamble_mat(ii,:) = cur_line;
    cur_line = circshift(cur_line, 1);
    %cur_line = cur_line.';
end
preamble_mat_r = preamble_mat(:,1:N_taps); % take only the first N_taps 

for ii=1:N_preamble
    cur_one_preamble = preamble_t((ii-1)*flags.N_subcarr+1:ii*flags.N_subcarr);
    %Construct the preamble matrix
    h_cur_est = preamble_mat_r \ cur_one_preamble.';
    if ii == 1
        ht_est = h_cur_est;
    else
        ht_est = ht_est + h_cur_est;
    end
end

% estimate the H(f) by averaging the estimation from N_one_preamble
% preambles
ht_est = ht_est/N_preamble;

end

