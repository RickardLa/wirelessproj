function [c] = spectrum_method(fc, v, Ts, Ns, kc)
% This function generates samples of the fading channel 
% using the spectrum method. 

% Inputs: 
%    fc = Carrier freq.   [Hz]
%    v  = Relative speed  [m/s]
%    Ts = Sample interval [s]
%    N  = Length of filter [samples]
%    Ns = Number of simulated samples (Ns>>2*N+1)
%    kc = LOS-link 

fd = v*fc/3e8;                   % Doppler freq.   [Hz]
fs = 1/Ts; 

% Set G(f) to be the square root of the theoretical Doppler spectrum for
% Clarke's method
G = @(f) sqrt((1/(pi*fd))./sqrt((1-(f./fd).^2)));     

% Generate Ns samples of Gtil
f = 0:fs/Ns:(Ns-1)*fs/Ns; 
Gtil = G(f) + G(f-fs); 


x = (randn(Ns, 1) + 1i*randn(Ns, 1))*sqrt(1/2);

c = ifft(Gtil*x,Ns) + kc;

% NOTE: We have to scale x so c has unit variance. 
var(c)



end

