function [f_coher, PDPMat_cut] = calc_coherf_filter(data, flags)
% Calculate the coherence frequency given the data
%   Input:
%       data:   the data matrix
%       flags:  all the parameters goes here
%   output:
%       f_coher: the estimated f_coherence

N_windows = ceil(flags.N_bins * flags.BW2_prop);
tukey_alpha = 0.5;
mid = floor(flags.N_bins/2);
% Rectangular window
lpf_rec = ones(1,N_windows);
lpf_tukey = tukeywin(N_windows, tukey_alpha);
if flags.LPFfilter == 0 % choose the rectangular window
    cur_lpf = lpf_rec;
elseif flags.LPFfilter == 1 % choose the tuckey window
    cur_lpf = lpf_tukey.';
end

[HtMat, N_bins_f] = ifft_3dmat_filter(data, flags.N_line, flags.N_bins, cur_lpf, flags.BW2_prop);
[HtMatAmp, PDPMat] = calc_PDP(HtMat, flags.N_line, N_bins_f);
PDPMat = PDPMat./(flags.N_line^3);


% Plot figure
[~,locmax] = max(PDPMat);
PDPMat_cut_rec = PDPMat(locmax:end-10);
PDPMat_cut_db_rec = 10*log10(PDPMat_cut_rec);

x_axis = [1:length(PDPMat_cut_rec)] .* flags.dt;

figure
plot(x_axis,PDPMat_cut_db_rec)
xlabel('delay (s)')
ylabel('PDP (dB)')

f_coher = calc_fcohr(PDPMat_cut_rec, flags.dt/flags.BW2_prop);
PDPMat_cut = PDPMat_cut_rec;

end

