clc
clf
clear 
close

% ============ Simulation parameters ========================
N = 64;                 % Number of subcarriers
Nsym = 10;           % Number of OFDM-symbols
M = 4;                  % Modulation order. M-QAM
knorm = 2;              % 4-QAM --> 2,  16-QAM --> 10
Ncp = 4;                 % Length of cyclic prefix. cp > delay spread
Rb = 1.85e6;            % Bit rate
Rs = Rb/(log2(M));      % Symbol rate
fs = 1e6;               % Sample frequency
Ts = 1/fs;              % Sample duration


Pt = 0.1;               % Transmit power 
E = Pt*Ts;              % Energy per symbol
EbN0 = 0:25;            % EbN0 in dB
% ============== Channel parameters =======================

fc = 2e9;               % Carrier frequency
v = 15;                 % Relative speed
c = 3e8;                % Speed of light
fd = fc*(v/c);          % Doppler frequency
tau = [0 4];            % Two taps. Delay 0 and 4.
fdTs = fd*Ts;           % Normalized Doppler frequency
N02 = 2.07e-20;         % Noise spectral density
pathLoss = 10^(-101/10);% Pathloss linear scale
