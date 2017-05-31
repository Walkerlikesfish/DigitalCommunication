function [HtMat_f, N_bins_f] = ifft_3dmat_filter_back(Data, N_line, N_bins, cur_lpf, BW2_prop)
N_taps = 51;
N_bins_f = ceil(N_bins * BW2_prop);
HtMat_f = zeros(N_line,N_line,N_line,N_taps+1);
mid = floor(501/2);
rec_f = zeros(1, N_bins);
rec_f(ceil(mid-N_bins_f/2):floor(mid+N_bins_f/2)) = 1;

for ir = 1 : N_line
   for icol = 1:N_line
      for iz = 1:N_line
          curBin = Data{ir}{icol}{iz};
          curBin_f = curBin .* rec_f;
          curBin_f(ceil(mid-N_bins_f/2):floor(mid+N_bins_f/2)) = curBin_f(ceil(mid-N_bins_f/2):floor(mid+N_bins_f/2)) .* cur_lpf; % apply the lpf
          curHt = ifft(curBin_f);
          curHt = downsample(curHt, 10);
          [~,locmax] = max(abs(curHt(1:N_bins_f/2)));
          curHt = curHt(locmax:locmax+N_taps);
          HtMat_f(ir,icol,iz,:) = curHt; % normalisation
      end
   end
end

end