function [HtMat] = ifft_3dmat(Data, N_line, N_bins)
    HtMat = zeros(N_line,N_line,N_line,N_bins);
    for ir = 1 : N_line
        for icol = 1:N_line
            for iz = 1:N_line
                curBin = Data{ir}{icol}{iz};
                % reverse FFT for each space location i
                curHt = ifft(curBin);
                HtMat(ir,icol,iz,:) = curHt;
            end
        end
    end   
end