function [h_gen] = gen_channel(distr_model, flags, channel_type)
% Generate the impulse responce of the channle for flags.tap_max taps, with
% the specified distribution
% input:
%   - distr_model: the distribution parameter specified for the target
%   channel_type
%   - flags: flags group
%   - channel_type: 0-rician or 1-reylei

h_gen = zeros(1,flags.tap_max);

if channel_type == 1
    for ii=1:flags.tap_max
        cur_distr = makedist('rayleigh','B',distr_model.B(ii));
        phase_rand = rand * 2 * pi;
        h_gen(ii) = random(cur_distr) * exp(1j * phase_rand);
    end
elseif channel_type == 0
    for ii=1:flags.tap_max
        cur_distr = makedist('rician','s',distr_model.S(ii),'sigma',distr_model.Sigma(ii));
        phase_rand = rand * 2 * pi;
        h_gen(ii) = random(cur_distr) * exp(1j * phase_rand);
    end
else
    disp('Sorry, we do not provide this service !')
end

end

