function [mat_an] = construct_an(mat_B, hn, flags)
% construct_an used to create the beamforming
mat_an = zeros(length(flags.theta),length(flags.phi), flags.tap_max);

for in=1:flags.tap_max
    for itheta=1:length(flags.theta)
        for iphi=1:length(flags.phi)
            % sum through the space 1000 points of |Bi|^2
            cur_Bi = squeeze(mat_B(:,:,:,itheta,iphi));
            cur_Bi = abs(cur_Bi).^2;
            cur_deno = sum(sum(sum(cur_Bi)));
            % sum through the space 1000 points
            cur_Bi = squeeze(mat_B(:,:,:,itheta,iphi));
            cur_Bi = conj(cur_Bi);
            cur_hn = hn(:,:,:,in);
            cur_nomi = sum(sum(sum(cur_hn .* cur_Bi)));
            mat_an(itheta,iphi,in) = cur_nomi/cur_deno;
        end
    end
end

end

