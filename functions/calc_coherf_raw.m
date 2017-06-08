function [f_coher, PDPMat_cut] = calc_coherf_raw(data, flags)
% Calculate the coherence frequency given the data
%   Input:
%       data:   the data matrix
%       flags:  all the parameters goes here
%   output:
%       f_coher: the estimated f_coherence

HtMat = ifft_3dmat(data, flags.N_line, flags.N_bins);
% calculate PDP
[HtMatAmp, PDPMat] = calc_PDP(HtMat, flags.N_line, flags.N_bins);
PDPMat = PDPMat./(flags.N_line^3);

[~,locmax] = max(PDPMat); % find the first local maxima
PDPMat_cut = PDPMat(locmax:locmax+flags.N_pts-1); % clip the responce
PDPMat_cut_db = 10*log10(PDPMat_cut); % adapt unit to dB

x_axis = [1:length(PDPMat_cut_db)] .* flags.dt;

figure
plot(x_axis,PDPMat_cut_db)
title('Power Delay Profile - Original')
xlabel('delay (s)')
ylabel('PDP (dB)')

% estimate coherence frequency
f_coher = calc_fcohr(PDPMat_cut, flags.dt);

end

