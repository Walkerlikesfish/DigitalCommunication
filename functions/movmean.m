function [outarr] = movmean(in_arr, k)
%movmean
% use to calculate moving average on input array in_arr, with window number
% k

len_in = length(in_arr);
outarr = in_arr;
for ii=1:len_in-k
    cur_win = in_arr(ii:ii+k);
    cur_win = mean(cur_win);
    outarr(ii) = cur_win;
end

end

