function [mat_a_n] = calc_a_n(hn, flags)

lamda = 3e8 / flags.f_c; %wavelength of the carrier 
mat_a_n = zeros(length(flags.theta), length(flags.phi), flags.tap_max);

for itheta=1:length(flags.theta)
    for iphi=1:length(flags.phi)
        cur_theta = flags.theta(itheta);
        cur_phi = flags.phi(iphi);
        cur_beta = (2*pi/lamda)*[sin(cur_theta)*cos(cur_phi), sin(cur_theta)*sin(cur_phi), cos(cur_theta)];
        cur_nomin = 0;
        cur_denomin = 0;
        for ix=0:9
            for iz=0:9
                for iy=9:-1:0
                    ir = [ix, iy, iz] * 0.02;
                    
                    cur_hi = hn(ix+1,10-iy,iz+1,:);
                    N_h = size(cur_hi, 2);
                    [~, n_start] = max(abs(cur_hi(1:(N_h+1)/2)));
                    cur_h = cur_hi(n_start:n_start+flags.tap_max-1);
                    
                    cur_B = exp(-1*1j*cur_beta*ir');
                    
                    cur_nomin = cur_nomin + cur_h*conj(cur_B);
                    cur_denomin = cur_denomin + (abs(cur_B))^2;
                end
            end
        end
        mat_a_n(itheta, iphi, :) = cur_nomin/cur_denomin;
    end
end

end

