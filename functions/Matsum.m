function [HtMat_sum] = Matsum(N_line,HtMat)
HtMat_sum = zeros(N_line,N_line,N_line,1);
for ir = 1 : N_line
    for icol = 1:N_line
        for iz = 1:N_line
            HtMat_sum(ir,icol,iz) = sum(HtMat(ir,icol,iz,:));
        end
    end
end