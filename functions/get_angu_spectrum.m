function [theta_list, phi_list, a_n] = get_angu_spectrum(data,narrow_band,plot_it,save_name, hn)
    N_theta = 25;
    N_phi = 50;
    N_X = 10;
    N_Y = 10;
    N_Z = 10;

    if narrow_band
        N_taps = 1;
    else
        N_taps = 8;
    end

    h_tot = zeros(1, 1, N_taps);

    theta_list = linspace(0, pi, N_theta);
    phi_list = linspace(0, 2*pi, N_phi);
    X_list = 0:9;
    Y_list = 9:-1:0;
    Z_list = 0:9;

    a_n  = zeros(N_theta, N_phi, N_taps);

    c = 3e8;
    fc = 2.35e9;
    lambda = c/fc;
    N_kept = 30;

    for phi_ind = 1:N_phi
        phi_ind
        phi = phi_list(phi_ind);
        for theta_ind = 1:N_theta
            theta = theta_list(theta_ind);
            beta = (2*pi/lambda)*[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
            B_i_tot = 0;

            for X_ind = 1:N_X
                X = X_list(X_ind);
                for Y_ind = 1:N_Y
                    Y = Y_list(Y_ind);
                    for Z_ind = 1:N_Z
                        Z = Z_list(Z_ind);
                        transfer_func = data{X_ind}{Y_ind}{Z_ind};
                        N = size(transfer_func, 2);
                        h = ifft(transfer_func)*sqrt(N);

                        h = h(1:10:end);
                        h = hn(X_ind,Y_ind,Z_ind,:);
                        N_h = size(h, 2);
                        % security: only search in the first half
                        [thing, n_start] = max(abs(h(1:(N_h+1)/2)));
                        h = h(n_start:n_start+N_kept-1);

                        if narrow_band
                            h_tot(1,1,1) = sum(h);         % narrow
                        else
                            h_tot(1,1,:) = h(1:N_taps);
                        end

                        r = [X Y Z]*0.02;

                        B_i = exp(-1*1j*beta*r');

                        a_n(theta_ind,phi_ind,:) = a_n(theta_ind,phi_ind,:)+h_tot(1,1,:)*conj(B_i);

                        B_i_tot = B_i_tot+abs(B_i)^2;
                    end
                end
            end
            a_n(theta_ind,phi_ind,:) = a_n(theta_ind,phi_ind,:)/B_i_tot;
        end
    end

    if plot_it
        if narrow_band
            figure;
            mesh(phi_list, theta_list, abs(a_n(:,:,1)));
        else
            for ii=1:N_taps
                figure;
                mesh(phi_list, theta_list, abs(a_n(:,:,ii)));
            end
        end
    end

    % non-empty filename: save it
    if ~strcmp(save_name, '')
        save(save_name, 'phi_list', 'theta_list', 'a_n')
    end
end
