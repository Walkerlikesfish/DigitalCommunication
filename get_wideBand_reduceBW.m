function [HtMat] = get_wideBand_reduceBW(data, flags)
% Calculate the coherence frequency given the data
%   Input:
%       data:   the data matrix
%       flags:  all the parameters goes here
%   output:
%       HtMat: time domain H with reduced bandwidth and wideband model

% filter the signal to get the 20Mhz signal
N_windows = ceil(flags.N_bins * flags.BW2_prop);
tukey_alpha = 0.8;
% Rectangular window
lpf_tukey = tukeywin(N_windows, tukey_alpha);
cur_lpf = lpf_tukey.';
% inverse FFT: H(f) -> h(n)
[HtMat, ~] = ifft_3dmat_filter(data, flags.N_line, flags.N_bins, cur_lpf, flags.BW2_prop);

end

