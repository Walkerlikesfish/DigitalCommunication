function [rayl_fit, rice_fit] = fit_narrowb(data, flags)
% Calculate the coherence frequency given the data
%   Input:
%       data:   the data matrix
%       flags:  all the parameters goes here
%   output:
%       rayl_fit: the reyleih fitting object
%       rice_fit: the reyleih fitting object

HtMat_LOS = ifft_3dmat(data, flags.N_line, flags.N_bins);
HtMat_LOS_sum = Matsum(flags.N_line, HtMat_LOS);
Vec = reshape(abs(HtMat_LOS_sum),1000,1);

% fit the random variables
rayl_fit = fitdist(Vec,'rayleigh');
rice_fit = fitdist(Vec,'rician');

K = rice_fit.s^2 / (2 * rice_fit.sigma^2)

end


