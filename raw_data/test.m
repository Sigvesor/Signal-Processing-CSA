load('data_healthy.mat')

% park
id = sqrt(2/3) * ia - sqrt(1/6) * ib - sqrt(1/6) * ic;
iq = sqrt(1/2) * ib - sqrt(1/2) * ic;
ip = sqrt(id.^2 + iq.^2);

[freq, val] = mfft(time, ip);

plot(freq, val)


function [f, Y] = mfft(t, y)
N	= length(y); 
Ts	= (t(end) - t(1))/(N - 1);
df	= 1/(N*Ts);
f	= ((0:(N-1))*df).';
Y	=abs(2*fft(y,N)/N); 
end

% note that a FFT amplitude equal to the amplitude signal fed to the FFT, 
%  We need to normalize FFTs by the number of sample points.
% Note that doing this will divide the power between the positive and negative sides, 
%so if you are only going to look at one side of the FFT, you can multiply the FFT by 2, and you'll get the expected magnitude
% study 2 different frequencies, then multiply the FFT by 4

