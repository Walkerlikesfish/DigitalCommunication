function [nmse] = calc_NMSE(h_est, h_ori)
%CALC_NMSE is the function to calculate the n mean square error
% input:
%   - h_est: estimated transfer function
%   - h_ori: original transfer function
% output:
%   - nmse: NMSE

nmse = 0;
n = length(h_est);
cur_sum = 0;
for ii = 1:n
    cur_sum = cur_sum + (abs(h_est(ii)-h_ori(ii)))^2;
    nmse = nmse + abs(h_ori(ii))^2;
end
nmse = cur_sum/nmse;

end

