function [ach_bits_mo] = SIMO_channel(trans_bits_in, flags)
% Simulating the transmitter side for a SISO channel, input the parameter
% group flags, generating a series of bits as the transmitter output
%
% VUB BRUFACE
% Yu Liu, Bohan Zhang, Xianjun Mao
% 

n_channels = length(flags.selected_channels);
ach_bits_mo = [];

% iterate through each channel

for ic=1:n_channels
    cur_cid = flags.selected_channels(ic);
    id_channel = [floor(cur_cid/100), floor(mod(cur_cid,100)/10), mod(cur_cid,10)];
    
    % MPC
    if flags.MPCchoice == 1
        ht = squeeze(flags.MPCht_LOS(id_channel(1),id_channel(2),id_channel(3),:));
        trans_bits_h = conv(trans_bits_in, ht);
    elseif flags.MPCchoice == 2
        ht = squeeze(flags.MPCht_NLOS(id_channel(1),id_channel(2),id_channel(3),:));
        trans_bits_h = conv(trans_bits_in, ht);
    elseif flags.MPCchoice == -1
        trans_bits_h = trans_bits_in;
    end
    trans_bits = trans_bits_h;
    
    % AWGN
    if flags.AWGN == 1
        % if AWGN is added to the channel
        % calc signal series energy
        En_signal=(trapz(abs(trans_bits).^2))*(1/flags.BW);
        % calc the per bit Energy
        Eb = En_signal/flags.Nbits/2;
        EbN0_mag = 10^(flags.EbN0(flags.EbN0i)/10);
        N0 = Eb/EbN0_mag;
        P_noise = 2*N0*flags.BW;
        noise_s = sqrt(P_noise/2)*(randn(1,length(trans_bits))+1i*randn(1,length(trans_bits)));
        anoise_bits = trans_bits + noise_s;
    else
        % No AWGN is added
        anoise_bits = trans_bits;
    end

    if flags.MPCchoice ~= -1
        ht_len = length(ht);
        ach_bits = anoise_bits(1:end-ht_len+1);
    else
        ach_bits = anoise_bits;
    end

    %% Comment about STO and CFO 
    % STO and CFO are the same receiver matrix cuz they are sharing the
    % same oscillator
    
    %% Possible Time shifting
    if flags.STO == 1
        ach_bits = ach_bits';
        ach_bits = circshift(ach_bits, flags.timeshift);
        ach_bits = ach_bits';
    end

    %% Possible CFO
    if flags.CFO == 1
        n_bits = length(ach_bits);
        Ts = 1/flags.BW;
        s_shift = exp(1j*2*pi*flags.f_tx*[1:n_bits]*Ts);
        ach_bits = ach_bits .* s_shift;
    end  
    
    ach_bits_mo = [ach_bits_mo; ach_bits];
end

end

