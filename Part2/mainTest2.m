%%  SSY135 Project part 1

clear all
close all

%-----------
% Parameters
%-----------

Ms = [4 16 64];

BERS = zeros(3,26);

for index = 1:length(Ms)

Fc = 2E9;       % 2 Ghz carrier
Fs = 1E6;       % 1 Mhz sample rate
M = Ms(index);         % size of constellation
Ts = 1/Fs;      % symbol rate
N = 64;         % number of subcarriers, limited by Tcoh 
tBER = 5E-3;    % target BER
Nsym = round((100/tBER)/N*log2(M)); % total number of OFDM symbols
Nsym = 1000;
tau = [0; 4];   % delay taps i n samples
CP = 4;         % cyclic prefix
v = 15;         % speeed in m/s
c = 3E8;        % speed of light
fD = v/(c/Fc);  % Doppler freq
fdTs = fD*Ts;   % normalised Doppler freq
Pavg = 0.1;     % average tx power in watts
N0_2 = 2.07E-20;    % noise spectral density at receiver
EbN0 = [0:1:25];    % Eb/N0 range
%EbN0 = 14;
PL = 101;           % path loss in dB        
PL_linear = 10^(PL/10);    % path loss linear
nbits = N*log2(M);  % bits per OFDM symbol
% not used
rb = 1.85E6;    % target bit rate (1.85 Mbps) 
rs = rb/log2(M); % symbol rate

%---------------
% Begin BER loop
%---------------

BER = zeros(1,length(EbN0));

for i = 1 : length(EbN0)    

    %---------------------
    % Generate random bits
    %---------------------

    b = randi([0, 1], 1, N*Nsym*log2(M));   % bits
    bs = reshape(b(:), log2(M), N*Nsym);    % bits
    sb = bi2de(bs','left-msb');
    
    %-----------------
    % Calculate energy
    %-----------------
    
    Eb = 10.^(EbN0(i)/10)*sqrt(N0_2);
    E = log2(M)*Eb;
    
    %--------------
    % Constellation
    %--------------

    %E = Pavg*Ts; % 0.1 watts linear scale
    knorm = mean(abs(qammod(sb,M).^2))
    sp = sqrt(E/knorm)*qammod(sb,M);
    %scatterplot(sp);
    sp = reshape(sp,[N,Nsym]); % 1 column: 1 OFDM symbol   

    %-------------------------
    % Begin transmission loop 
    %-------------------------

    br = zeros(1,length(b)); % empty container for received bits

    for j = 1 : Nsym

        %--------------------------
        % Grab the next OFDM symbol
        %--------------------------
        
        s = sp(:,j);
        
        %-----------------
        % Take inverse FFT
        %-----------------
        
        z = sqrt(N/1)*ifft(s);
        
        %------------------
        % Add cyclic prefix
        %------------------

        z = [ z(end-CP+1:end).' z.'].';

        %----------------------
        % Transmit over channel
        %---------------------- 

        [y, h] = Fading_Channel(z, tau, fdTs);    
        h = [h(:,1) zeros(N+8,3) h(:,2)];   % t0 0 0 0 t4
        
        %--------------------------
        % Apply noise and path loss
        %--------------------------

        %sigma = sqrt(1/(2*10^(EbN0(1)/10)));

        % n[k] = scaled FFT of w[n] with variance N0
        sigma = sqrt(N0_2);

        n = sqrt(0.5)*(randn(1,length(y))*sigma + 1i*randn(1,length(y))*sigma);
        y = y / sqrt(PL_linear);
        y = y + n';
        
        %---------------------------
        % Remove prefix and take FFT
        %---------------------------

        y = y(CP+1:end);
        r = sqrt(1/N)*fft(y,N); 

        %--------------------------------------------------------------
        % Calculate C matrix (compensates for phase shift) and multiply
        %--------------------------------------------------------------

        C = fft(h(1,:).',N); % create diagonal matrix from one column in FFT
        A = 1./((C.*conj(C)));
        C = diag(C);
        A = diag(A);
        demod = sqrt(PL_linear)*(1/sqrt(E/knorm))*A*C'*r;
        sbrx = qamdemod(demod,M);
        
        sprx = de2bi(sbrx,log2(M),'left-msb');
        bsrx = reshape(sprx', 1, log2(M)*N);
        br((j-1)*nbits+1:j*nbits) = bsrx; 

    end % end Nsym
    
    BER(i) = sum(bitxor(b,br))/length(b);
    
end % end EbN0
%scatterplot(demod);
BERS(index,:) = BER; 
end

figure
semilogy(EbN0,BERS(1,:),'--');
hold on
semilogy(EbN0,BERS(2,:),'--');
semilogy(EbN0,BERS(3,:),'--');
plot([0 25],[tBER tBER],'r--');
xlabel('EbN0 [dB]'); ylabel('BER');
text(1,0.0075,'Target BER','Color','red');
grid on
legend('M = 4','M = 16','M = 64');
title('Uncoded OFDM');

