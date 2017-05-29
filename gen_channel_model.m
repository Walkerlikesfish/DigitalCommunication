function [hspace_mat] = gen_channel_model(a_n, flags)
% Generate the channel model (h(n)) based on the a_n(beamforming), the
% channel model also takes into consideration of the MPC

lamda = 3e8 / flags.f_c; %wavelength of the carrier 

hspace_mat = zeros(flags.N_line,flags.N_line,flags.N_line,flags.tap_max);
% iterate through the spacatial space
for itap = 1:flags.tap_max
for ix=0:9
    for iz=0:9
        for iy=9:-1:0
            %sum through the angular space
            ir = [ix, iy, iz] * 0.02;
            phase_rand = rand * 2 * pi;
            for itheta=1:length(flags.theta)
                for iphi=1:length(flags.phi)
                    cur_theta = flags.theta(itheta);
                    cur_phi = flags.phi(iphi);
                    cur_beta = (2*pi/lamda)*[sin(cur_theta)*cos(cur_phi), sin(cur_theta)*sin(cur_phi), cos(cur_theta)];
                    cur_h = a_n(itheta, iphi, itap)*exp(1j*(phase_rand-cur_beta*ir'));
                    hspace_mat(ix+1,10-iy,iz+1,itap) = hspace_mat(ix+1,10-iy,iz+1,itap) + cur_h;
                end
            end
        end
    end
end
end



end

