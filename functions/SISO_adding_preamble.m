function [one_preamble_F, one_preamble_t, trans_bits_ap] = SISO_adding_preamble(trans_symbols, flags)
%SISO_adding_preamble 
%   intputs:
%       -trans_symbols: the bit seires coming out from the transmitter
%       -flags:
%   outputs:
%       -preamble: [N_sub] preamble in frequency domain

% Estimate amplitude A

% A_est = A_est * 2 * pi;
% [?] is the energy calculated from the TD == FD ? (Pascal theory)->2pi

% Generate -A +A series in Frequency domain
one_preamble_F = rand([1, flags.N_subcarr]);
for ii=1:flags.N_subcarr
    if one_preamble_F(ii)>0.5
        one_preamble_F(ii) = 1;
    else
        one_preamble_F(ii) = -1;
    end
end

% IFFT to time domain -> preamble -> preamble x2 repeat
one_preamble_t = ifft(one_preamble_F);

power_pre = sum(abs(one_preamble_t).^2)/size(one_preamble_t,2);
power_sig = sum(abs(trans_symbols).^2)/size(trans_symbols,2);
prop_pre_sig = power_pre/power_sig;
one_preamble_t = one_preamble_t .* prop_pre_sig;

% calc CP of Preamble
cur_cp = one_preamble_t(end-flags.N_cp*2+1:end); % extract CP
% add preamblex2+CP to signal 
two_preamble_n_cp = [cur_cp one_preamble_t, one_preamble_t];

trans_bits_ap = [two_preamble_n_cp trans_symbols];

end

