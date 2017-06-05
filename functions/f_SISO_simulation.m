function [arec_bits] = f_SISO_simulation(symbol_cp_s, flags, verbal)
% f_SIMO_simulation is designed to wrap all simulation together into one
% funciton, simulation from the input_symbol to the bits stream at the
% receiver end.
if nargin < 3
    verbal = 0;
end

global NMSErec;
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
    % calculate and record the NMSE
    if flags.MPCchoice == 0
        ht = flags.MPCht.h_LOS_rice;
    elseif flags.MPCchoice == 1
        ht = flags.MPCht.h_NLOS_rayl;   
    end
    hf_ori = fft(ht, flags.N_subcarr);
    if flags.tdEQ == 0
        htf_est = hf_est;
    else
        htf_est = fft(ht_est, flags.N_subcarr);
    end
    NMSErec(flags.EbN0i) = calc_NMSE(htf_est, hf_ori); % record the NMSE
    % extract the preamble from the signals
    ach_bits = ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
end
% [Receiver] - Equalisation
[arec_bits] = SISO_receiver(ach_bits, flags, verbal);

end

