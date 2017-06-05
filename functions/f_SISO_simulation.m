function [arec_bits] = f_SISO_simulation(symbol_cp_s, flags, verbal)
% f_SIMO_simulation is designed to wrap all simulation together into one
% funciton, simulation from the input_symbol to the bits stream at the
% receiver end.
if nargin < 3
    verbal = 0;
end

[ach_bits] = SISO_channel(symbol_cp_s, flags);
% [Receiver] - Time acquisition
if flags.STO
    %[t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
end
% [Receiver] - CFO acquisition
if flags.CFO
    [df_est, ach_bits] = SISO_estimate_CFO(ach_bits, flags);
end
%df_est_show = df_est/(2*pi) * flags.BW
% [Receiver] - Channel Estimation 
if flags.preamble_yes == 1
    % [Receiver - Estimate the channel]
    % estimating in frequency domain
    [hf_est] = SISO_ZF_estimator(ach_bits, flags);
    % estimating in time domain
    [ht_est] = SISO_TD_estimator(ach_bits, flags);
    flags.MPCZF = hf_est; % set the estimated channel to global var
    flags.MPCTD = ht_est;
    % extract the preamble from the signals
    ach_bits = ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
end
% [Receiver] - Equalisation
[arec_bits] = SISO_receiver(ach_bits, flags, verbal);

end

