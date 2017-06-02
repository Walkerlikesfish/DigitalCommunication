function [HtMatAmp, PDPMat] = calc_PDP(HtMat, N_line, N_bins)

HtMatAmp = zeros(N_line,N_line,N_line,N_bins);
for ir = 1 : N_line
    for icol = 1:N_line
        for iz = 1:N_line
            % |h(t)|^2
            HtMatAmp(ir,icol,iz,:) = (abs(squeeze(HtMat(ir,icol,iz,:)))).^2;
        end
    end
end

PDPMat = zeros(1,N_bins);

% sum through the local area 
for ib = 1:N_bins
    curArr = HtMatAmp(:,:,:,ib);
    curAmp = sum(sum(sum(curArr)));
    PDPMat(ib) = curAmp;
end

end