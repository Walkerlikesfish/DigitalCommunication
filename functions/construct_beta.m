function [mat_beta] = construct_beta(flags)
%Construct the beta matrix for each position in space
% input:
%   -flags:
% output:
%   -mat_beta:
mat_beta = zeros(length(flags.theta), length(flags.phi), 3);
lamda = 3e8 / flags.f_c;

for i1=1:length(flags.theta)
    for i2=1:length(flags.phi)
        cur_theta = flags.theta(i1);
        cur_phi = flags.phi(i2);
        mat_beta(i1, i2, :) = 2*pi/lamda*[sin(cur_theta)*cos(cur_phi)...
            ,sin(cur_theta)*sin(cur_phi), cos(cur_theta)];
    end
end

end

