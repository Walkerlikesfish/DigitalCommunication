function [f_cohr] = calc_fcohr(PDPMat_cut, dt)

P_T = sum(PDPMat_cut);
Tau = [0:1:size(PDPMat_cut,2)-1];

Tau_mean = sum(Tau .* PDPMat_cut)/P_T;
Delta_T = sqrt(sum(Tau.^2 .* PDPMat_cut)/P_T - Tau_mean^2) * dt;

f_cohr = 1/(2*pi*Delta_T);

end