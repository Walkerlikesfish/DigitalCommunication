function [rayl] = fit_wideb_rayl(data, flags)
% Calculate the coherence frequency given the data
%   Input:
%       data:   the data matrix
%       flags:  all the parameters goes here
%   output:
%       rayl_fit: the reyleih fitting object

% filter the signal to get the 20Mhz signal
N_windows = ceil(flags.N_bins * flags.BW2_prop);
tukey_alpha = 0.5;
mid = floor(flags.N_bins/2);
% Rectangular window
lpf_tukey = tukeywin(N_windows, tukey_alpha);
cur_lpf = lpf_tukey.';
% inverse FFT: H(f) -> h(n)
[HtMat, N_bins_f] = ifft_3dmat_filter(data, flags.N_line, flags.N_bins, cur_lpf, flags.BW2_prop);

% only take the first *flags.max_taps* taps
for ii=1:flags.tap_max
    cur_vec = HtMat(:,:,:,ii);
    cur_vec = reshape(abs(cur_vec),1000,1);
    rayl_fit = fitdist(cur_vec,'rayleigh');
    rayl.B(ii) = rayl_fit.B;
end

end

