function [noise] = spectrally_matched_noise(sample_sound)

% [NOISE] = SPECTRALLY_MATCHED_NOISE(SAMPLE_SOUND)
%
% generates noise with same length and spectrum as sample_sound

spec = abs(fft(sample_sound));

N = length(spec);
Nh = ceil(N/2);

phi = rand(size(spec))*2*pi;
x_r = (spec(1:Nh).*[0; cos(phi(2:Nh))])';
x_i = (spec(1:Nh).*[0; sin(phi(2:Nh))])';

if rem(N,2)==0
    x(1:Nh) = x_r + i*x_i;
    x(Nh+1:N) = fliplr([(x_r(2:Nh)-i*x_i(2:Nh)) 1]);
else
    x(1:Nh) = x_r + i*x_i;
    x(Nh+1:N) = fliplr([(x_r(2:Nh)-i*x_i(2:Nh))]);
end

noise = real(ifft(x));

