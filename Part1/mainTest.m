%% filter
clc
clf
clear
close



fc = 2e9;               % Carrier freq.   [Hz]
v = 30/3.6;             % Relative speed  [m/s]
Ts = 1e-4;              % Sample interval [s]
N = 1e3;                % Length of filter [samples]
Ns = 1e5;               % Number of simulated samples (Ns>>2*N+1) 
kc = 0; 



c = filter_method(fc, v, Ts, N, Ns, kc);

fd = v*fc/3e8;          % Doppler freq.   [Hz]
fs = 1/Ts;              % Sample freq.    [Hz]

[P, f] = pwelch(c,[] ,[] ,[],fs,'centered'); 
plot(f,10*log10(P),'b')
hold on
Sc = (1/(pi*fd))./sqrt((1-(f./fd).^2));     % Theoretical
Sc(f<-fd | f>fd ) = eps;
plot(f,10*log10(Sc),'r')


axis([-300 300 -100 0])
grid on
legend('Simulation', 'Theoretical')
title('PSD')
xlabel('Frequency [Hz]')
ylabel('Power [dB]')

%% spectrum
clc
clf
clear
close

fc = 2e9;               % Carrier freq.   [Hz]
v = 30/3.6;             % Relative speed  [m/s]
Ts = 1e-4;              % Sample interval [s]
Ns = 1e5;               % Number of simulated samples (Ns>>2*N+1) 
kc = 0; 


c = spectrum_method(fc, v, Ts, Ns, kc);



fd = v*fc/3e8;          % Doppler freq.   [Hz]
fs = 1/Ts;              % Sample freq.    [Hz]

[P, f] = pwelch(c,[] ,[] ,[],fs,'centered'); 
plot(f,10*log10(P),'b')
hold on
Sc = (1/(pi*fd))./sqrt((1-(f./fd).^2));     % Theoretical
Sc(f<-fd | f>fd ) = eps;
plot(f,10*log10(Sc),'r')


axis([-300 300 -100 0])
grid on
legend('Simulation', 'Theoretical')
title('PSD')
xlabel('Frequency [Hz]')
ylabel('Power [dB]')