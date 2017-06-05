function [Rdz] = calc_spatial_corr_Z(a_n, flags, cur_tap)

lamda = 3e8 / flags.f_c; %wavelength of the carrier 
flags.dz_list = [0:lamda/20:lamda*5];
Rdz = zeros(1, length(flags.dz_list));

beta_s = 2*pi/lamda;

for idz=1:length(flags.dz_list)
    cur_dz = flags.dz_list(idz);
    cur_rdz = 0;
    for itheta=1:length(flags.theta)
        rdz_per_theta = 0;
        for iphi=1:length(flags.phi)
            cur_theta = flags.theta(itheta);
            cur_u = beta_s * cos(cur_theta);
            if size(size(a_n),2) == 3
                 tmpv = abs(a_n(itheta, iphi, cur_tap))^2 * exp(1j*cur_u*cur_dz);
            else if size(size(a_n),2) == 2
                     tmpv = abs(a_n(itheta, iphi))^2 * exp(1j*cur_u*cur_dz);
                end
            end
           
            rdz_per_theta = rdz_per_theta + tmpv;
        end
        rdz_per_theta = rdz_per_theta/length(flags.phi); % average through the phi
        cur_rdz = cur_rdz + rdz_per_theta;
    end
    Rdz(idz)=Rdz(idz) + cur_rdz;
end

%Normalize
Rdz_to_plot = real(Rdz);
coef = 1/Rdz_to_plot(1);
Rdz_to_plot = Rdz_to_plot .* coef;

figure(1);plot(flags.dz_list/lamda,real(Rdz_to_plot));title('Spatial Correlation on Z direction');
xlabel('Distance (x lamda)');ylabel('R(dz)');
hold on

ref_line = zeros(1,length(flags.dz_list));
plot(flags.dz_list/lamda,ref_line, '-.')

end

