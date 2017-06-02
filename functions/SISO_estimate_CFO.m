function [f_est, out_ach_bits] = SISO_estimate_CFO(ach_bits, flags)
%SISO_estimate_CFO 
% try to estimate CFO by preambles
%   intputs:
%       -ach_bits: the bit seires coming out from the transmitter
%       -flags:
%   outputs:
%       -f_est: [N_sub] preamble in frequency domain
%       -out_ach_bits: after the simple compensation on the time domain
Ts = 1/flags.BW;
preamble_1 = ach_bits(flags.N_cp*2+1 : flags.N_cp*2+flags.preamble_size(2));
preamble_2 = ach_bits(flags.N_cp*2+flags.preamble_size(2)+1 : flags.N_cp*2+flags.preamble_size(2)*2);
% figure(10);plot(preamble_1); hold on; plot(preamble_2)

phase_est = angle(sum(preamble_1 .* conj(preamble_2)));
f_est = phase_est / (flags.preamble_size(2)*2*pi*Ts); % 2pi*f_est(Hz)
%f_est = flags.f_tx;
n_bits = length(ach_bits);

comp_phase = exp(1j*f_est*2*pi*[1:n_bits]*Ts);
out_ach_bits = ach_bits .* comp_phase;

end

