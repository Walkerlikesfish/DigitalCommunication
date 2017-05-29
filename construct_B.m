function [mat_B] = construct_B(mat_beta, flags)
% construct the B matric
%

mat_B = zeros(flags.N_line,flags.N_line,flags.N_line, length(flags.theta), length(flags.phi));
for ix = 1:flags.N_line
    for iy = 1:flags.N_line
        for iz = 1:flags.N_line
            for itheta=1:length(flags.theta)
                for iphi=1:length(flags.phi)
                    ir = [ix-1, iy-1, iz-1]; % 
                    cur_beta = squeeze(mat_beta(itheta, iphi, :));            
                    mat_B(ix,iy,iz,itheta,iphi) = exp(-1j*ir*cur_beta);
                end
            end
        end
    end
end

end

