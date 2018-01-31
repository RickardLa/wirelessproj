function [c] = filter_method(fc, v, Ts, N, Ns, kc)
% This function generates samples of the fading channel 
%using the filter method. 

% Inputs: 
%    fc = Carrier freq.   [Hz]
%    v  = Relative speed  [m/s]
%    Ts = Sample interval [s]
%    N  = Length of filter [samples]
%    Ns = Number of simulated samples (Ns>>2*N+1)
%    kc = LOS-link 

fd = v*fc/3e8;                   % Doppler freq.   [Hz]


% Samples spaces Ts-seconds apart
t = -N*Ts:Ts:(N-1)*Ts;           

% Generate filter based on the Doppler spectrum for Clarke's method
g = besselj(0.25, 2*pi*fd.*abs(t))./(abs(t).^0.25); 
g(t==0) = ((pi*fd)^0.25)/gamma(5/4);

% Truncate filter by multiplying with a window. 
% Hann-window gives most smooth PSD.
w = hann(2*N); 
g_hat = w'.*g; 

% Divide by energy of signal to get unit energy
K = 1/sqrt(sum(abs(g_hat).^2));
g_hat = K.*g_hat;

% Colvolve g_hat with a zero-mean, unit-variance Gaussian r.v.
x = (randn(Ns, 1) + 1i*randn(Ns, 1))*sqrt(1/2);
c = conv(x,g_hat,'same') + kc;


end

