%%  Part 2 - OFDM Simuation
run('parameters.m'); 


BER = zeros(1,length(EbN0));


for i = 1:length(EbN0)
    TxBits = randi([0 1], 1, N*Nsym*log2(M));        % Random bits 
    TxBits_s = reshape(TxBits(:), log2(M),N*Nsym);   % Reshape to matrix
    TxBits_s = bi2de(TxBits_s', 'left-msb');         % Convert to decimal

    EbN0lin = 10^(EbN0(i)/10);                       % EbN0 in linear scale
    Eb = sqrt(N02)*EbN0lin;                          % Bit energy at receiver
    E = log2(M)*Eb;                                  % Energy per transmitted symbol
    
    TxSymb = sqrt(E/knorm)*qammod(TxBits_s,M);       % Modulate and scale symbols
    TxSymb = reshape(TxSymb,N,Nsym);                 % Reshape to (N,Nsym)-matrix
    
    % For given Eb and E, the OFDM-symbols are ready for transmission. 
    receivedBits = zeros(Nsym,N*log2(M)); 
    for j = 1:Nsym
        z = sqrt(N/Ts)*ifft(TxSymb(:,j));            % Inverse FFT of OFDM-symbols, one by one. 
        z = [z(end-Ncp+1:end); z];                   % Add cyclic prefix
        
        [r,h] = Fading_Channel(z,tau,fdTs);          % Send over channel
        hs = [h(:,1) zeros(N+2*Ncp,3) h(:,2)];
        C = diag(fft(hs,N))';
        
        noise = sqrt(N02)*(randn(1,length(r)) + 1i*randn(1,length(r)));
        r = r*pathLoss+noise;                        % Adding pathloss and Gaussian noise
        r = r(Ncp+1:end);                            % Remove cyclic prefix
        r = sqrt(Ts/N)*fft(r,N);                     % Retrieve symbols
        
        
        requal = C/abs(C).^2*r;                     % Equalize channel
        RxSymb = qamdemod(sqrt(knorm/E)*requal,M); 
        RxBits = de2bi(RxSymb, 'left-msb');
        receivedBits(j,:) = reshape(RxBits,1,N*log2(M));
        

               
    end
    receivedBits = reshape(receivedBits,1,N*Nsym*log2(M));   

    
    BER(i) = sum(receivedBits~=TxBits)/(N*Nsym*log2(M));
    
end
semilogy(EbN0,BER,'or')



