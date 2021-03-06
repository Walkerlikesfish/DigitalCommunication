function [Rdx, Rdy] = calc_spatial_corr_XY(a_n, flags, cur_tap)

lamda = 3e8 / flags.f_c; %wavelength of the carrier 
flags.dz_list = [0:lamda/20:lamda*5];
beta_s = 2*pi/lamda;

tmp_theta_phi = zeros(length(flags.theta), length(flags.phi));

gama_record = [];
tmpv_record = zeros(1, length(flags.theta)*length(flags.phi));
cnt_tmpv = zeros(1, length(flags.theta)*length(flags.phi));

for idx=1:length(flags.dz_list)
    cur_dx = flags.dz_list(idx);
    cur_rdx = 0;
    for itheta=1:length(flags.theta)
        for iphi=1:length(flags.phi)
            cur_theta = flags.theta(itheta);
            cur_phi = flags.phi(iphi);
            cur_gama = acos(cos(cur_phi)*sin(cur_theta));
            
            if any(gama_record == cur_gama)
                gama_ind = find(gama_record == cur_gama);
            else
                gama_record = [gama_record cur_gama];
                gama_ind = find(gama_record == cur_gama);
            end
            
            cur_u = beta_s * cos(cur_gama);
            if size(size(a_n),2) == 3
                 tmpv = abs(a_n(itheta, iphi, cur_tap))^2 * exp(1j*cur_u*cur_dx);
            else if size(size(a_n),2) == 2
                     tmpv = abs(a_n(itheta, iphi))^2 * exp(1j*cur_u*cur_dx);
                end
            end
            
            tmpv_record(gama_ind) = tmpv_record(gama_ind) + tmpv;
            cnt_tmpv(gama_ind) = cnt_tmpv(gama_ind) + 1; % counting the possible overlaps
        end
    end
    cur_rdx = sum(nonzeros(tmpv_record)./nonzeros(cnt_tmpv));
    Rdx(idx) = cur_rdx;
end

% Normalize
Rdx_to_plot = real(Rdx);
coef = 1/Rdx_to_plot(1);
Rdx_to_plot = Rdx_to_plot .* coef;

figure(1);plot(flags.dz_list/lamda,real(Rdx_to_plot));title('Spatial Correlation on X direction');
xlabel('Distance (lamda)');ylabel('R(dx)');
hold on
% draw ref line
ref_line = zeros(1,length(flags.dz_list));
plot(flags.dz_list/lamda,ref_line, '-.')

gama_record = [];
tmpv_record = zeros(1, length(flags.theta)*length(flags.phi));
cnt_tmpv = zeros(1, length(flags.theta)*length(flags.phi));

for idy=1:length(flags.dz_list)
    cur_dx = flags.dz_list(idy);
    cur_rdy = 0;
    for itheta=1:length(flags.theta)
        for iphi=1:length(flags.phi)
            cur_theta = flags.theta(itheta);
            cur_phi = flags.phi(iphi);
            cur_gama = asin(sin(cur_phi)*sin(cur_theta));
            
            if any(gama_record == cur_gama)
                gama_ind = find(gama_record == cur_gama);
            else
                gama_record = [gama_record cur_gama];
                gama_ind = find(gama_record == cur_gama);
            end
            
            cur_u = beta_s * cos(cur_gama);
            if size(size(a_n),2) == 3
                 tmpv = abs(a_n(itheta, iphi, cur_tap))^2 * exp(1j*cur_u*cur_dx);
            else if size(size(a_n),2) == 2
                     tmpv = abs(a_n(itheta, iphi))^2 * exp(1j*cur_u*cur_dx);
                end
            end
            
            tmpv_record(gama_ind) = tmpv_record(gama_ind) + tmpv;
            cnt_tmpv(gama_ind) = cnt_tmpv(gama_ind) + 1; % counting the possible overlaps
        end
    end
    cur_rdy = sum(nonzeros(tmpv_record(1:900))./nonzeros(cnt_tmpv(1:900)));
    Rdy(idy) = cur_rdy;
end

% Normalize
Rdy_to_plot = real(Rdy);
coef = 1/Rdy_to_plot(1);
Rdy_to_plot = Rdy_to_plot .* coef;

figure(2);plot(flags.dz_list/lamda,real(Rdy_to_plot));title('Spatial Correlation on Y direction');
xlabel('Distance (lamda)');ylabel('R(dx)');
hold on
% draw ref line
ref_line = zeros(1,length(flags.dz_list));
plot(flags.dz_list/lamda,ref_line, '-.')

end

