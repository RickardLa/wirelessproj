clf
clc
clear 
close all

% ============== Specification of system ==================
 
N = 20;                 % Number of subcarriers
N_sym = 2;              % Number of OFDM-symbols
M = 4;                  % Modulation order. M-QAM
knorm = 2;              % 4-QAM --> 2,  16-QAM --> 10
cp = 10;                % Length of cyclic prefix. cp > delay spread
Rb = 1.85e6;            % Bit rate
Rs = Rb/(log2(M));      % Symbol rate
Ts = 1/Rs;              % Sample duration
fc = 2e9;               % Carrier frequency
N0 = 2.07e-20;          % Noise spectral density
Pt = 0.1;               % Transmit power 
E = Pt*Ts;              % Energy per symbol

% ============== Specification of channel =================

v = 15;                 % Relative speed
c = 3e8;                % Speed of light
fd = fc*(v/c);          % Doppler frequency
tau = [0 4];            % Two taps. Delay 0 and 4.
fdTs = fd*Ts;           % Normalized Doppler frequency
% =========================================================
%                   TRANSMITTER
% =========================================================

% =============== Bits to symbols =========================

b = randi([0 1],1,N*N_sym*log2(M));          % Generate bitstream
b_s = reshape(b(:),log2(M),N*N_sym);         % Reshape to matrix
s_b = bi2de(b_s','left-msb');                % Convert to decimal
sp = sqrt(E/knorm)*qammod(s_b,M);            % Modulate and scale symbols


% ============= Transform to time domain to retrieve signal ========
sp = reshape(sp,N,N_sym);                    % Reshape to (N,N_sym)-matrix
z = sqrt(N/Ts)*ifft(sp,[],1);                % IFFT column-wise and scale

% ============ Add cyclic prefix ===================================
z_cp = [z(end-cp+1:end,:); z];               % Add cyclic prefix for each OFDM symbol
z_concat = reshape(z_cp,(N+cp)*N_sym,1);     % Concatenate OFDM symbols to one vector for transmission 

% =========================================================
%                     CHANNEL
% =========================================================

[r,h] = Fading_Channel(z_concat,tau,fdTs);

hs = [h(1,1) 0 0 0 h(1,2)];                  % Add zeros to get correct delay profile
C = diag(fft(hs,N));

% % =============================================
% %             RECEIVER
% % =============================================

% =========== Separate OFDM symbols and remove CP ====================
y_cp = r(1:end-max(tau));                    % Receive transmitted signal. NOTE: No channel added here
y_cp = reshape(y_cp,N+cp,N_sym);
y = y_cp(end-N+1:end,:);                     % Remove cyclic prefix

% ========== Retrieve symbols mapped in freq domain ==================
r = sqrt(Ts/N)*fft(y,[],1); 

% scatterplot(r)
% ========== Demodulate through ML ===================================

rtest = C'/abs(C')^2*r;                      % Equalize r (A*C'*r)
% scatterplot(rtest)
s_hat = qamdemod(sqrt(knorm/E)*rtest,M);     % Demodulate with the factor used in modulation

out = reshape(de2bi(s_hat,'left-msb')',1,N*N_sym*log2(M));
BER = sum(out ~= b)/(N*N_sym*log2(M));


