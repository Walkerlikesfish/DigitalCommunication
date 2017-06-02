function [arec_bits] = f_SIMO_simulation(symbol_cp_s, flags)
% f_SIMO_simulation is designed to wrap all simulation together into one
% funciton, simulation from the input_symbol to the bits stream at the
% receiver end.

[ach_bits_mo] = SIMO_channel(symbol_cp_s, flags);
ach_bits = ach_bits_mo(1,:); 
% [Receiver] - Time acquisition
[t_est, ach_bits] = SISO_estimate_STO(ach_bits, flags);
% [Receiver] - CFO acquisition
[df_est, ach_bits] = SISO_estimate_CFO(ach_bits, flags);
% [Receiver] - Channel estimation independently for each receiver
n_channels = length(flags.selected_channels);
hf_est_mo = [];
ach_bits_mo_np = [];  % the bits stream after clipping the preamables
for ic=1:n_channels
    cur_cid = flags.selected_channels(ic);
    id_channel = [floor(cur_cid/100), floor(mod(cur_cid,100)/10), mod(cur_cid,10)];
    % estimate the channel function
    cur_ach_bits = ach_bits_mo(ic,:);
    [cur_hf_est] = SISO_ZF_estimator(cur_ach_bits, flags);
    hf_est_mo = [hf_est_mo; cur_hf_est];
    % extract the singal clipping the head
    cur_ach_bits = cur_ach_bits(1,(flags.N_subcarr+flags.N_cp)*2+1: end);
    ach_bits_mo_np = [ach_bits_mo_np; cur_ach_bits];
end
% [Receiver] - Equalizer
[arec_bits] = SIMO_receiver(ach_bits_mo_np, hf_est_mo, flags, 0);

end

