%% part1
clc
clf
clear
close


fc = 2e9;               % Carrier freq.   [Hz]
v = 30/3.6;             % Relative speed  [m/s]
fd = v*fc/3e8;          % Doppler freq.   [Hz]
Ts = 1e-4;              % Sample interval [s]
fs = 1/Ts;              % Sample freq.    [Hz]

N = 1e3;                % Length of filter [samples]
Ns = 1e5;               % Number of simulated samples (Ns>>2*N+1) 
kc = 0; 

t = -N*Ts:Ts:(N-1)*Ts;

% Filter
g = besselj(0.25, 2*pi*fd.*abs(t))./(abs(t).^0.25);
g(t==0) = ((pi*fd)^0.25)/gamma(5/4);

    
w = hann(2*N); 
g_hat = w'.*g; 

K = 1/sqrt(sum(abs(g_hat).^2));
g_hat = K.*g_hat;


% Input 
x = (randn(Ns, 1) + 1i*randn(Ns, 1))*sqrt(1/2 );
c = conv(x,g_hat,'same') + kc;


[P, f] = pwelch(c,[] ,[] ,[],fs,'centered'); 
plot(f,10*log10(P),'b')
hold on
Sc = (1/(pi*fd))./sqrt((1-(f./fd).^2));
Sc(f<-fd | f>fd ) = eps;
plot(f,10*log10(Sc),'r')


axis([-300 300 -100 0])
grid on
legend('Simulation', 'Theoretical')
title('PSD')
xlabel('Frequency [Hz]')
ylabel('Power [dB]')
