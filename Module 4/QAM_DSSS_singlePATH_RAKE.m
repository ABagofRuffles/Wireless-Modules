    clear all
    clc
    numPkts=1000; %number of packets sent

    N = 200; % number of symbols
    M = 16;   % constellation size
    k = log2(M); % bits per symbol
    
    %%%%%for spreading
    seq=[1;1;1;1;1;0;0;1;1;0;1;0;1]; 
    sequence=seq*2-1; %convert in 1s and -1s
    seqLen=length(sequence);

    %%%%%%%%for RAKE receiver
    numFingers=1; %if using RAKE defines number of fingers
    finger=zeros(seqLen,N,numFingers);
    decodeBitx=zeros(numFingers,N);

    % defining the real and imaginary PAM constellation
    % for 16-QAM
    alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
    alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
    k_16QAM = 1/sqrt(10);

    Eb_N0_dB  = [0:15]; % multiple Es/N0 values
    Es_N0_dB  = Eb_N0_dB + 10*log10(k);
    erTot=zeros(1,length(Eb_N0_dB));

    % Mapping for binary <--> Gray code conversion
    ref = [0:k-1];
    map = bitxor(ref,floor(ref/2));
    [tt ind] = sort(map);                                

    for ii = 1:length(Eb_N0_dB)
    for pktX=1:numPkts    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; % conversion from binary to decimal
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + j*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    %%spreading
    sT=sequence*s;
    s2=reshape(sT,1,seqLen*N);

    % noise
    % -----
    
EsNo=10^(Es_N0_dB(ii)/10);
stanDevNoise=sqrt((seqLen)/(2*EsNo));
    
%     stanDevNoise=sqrt((seqLen));
    noise =stanDevNoise *[randn(1,length(s2)) + j*randn(1,length(s2))]; % white guassian noise, 0dB variance 

h=(1/sqrt(2))*(randn+j*randn);    %equal power
 
%     y1 = s2 + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    y1 = conv(s2,h) + noise; % additive white gaussian noise

    
    
    %%despreading    
[~,I]=sort(abs(h),'descend'); %puts path gains in descending order
for n=1:numFingers
Rxchips=reshape(y1(I(n):end+(I(n)-1)),seqLen,N); %selects stongest path
% finger(:,:,n)=(Rxchips*h(I(n))')/sqrt(abs(h(I(n)))); %selects first 'n' stongest paths
finger(:,:,n)=(Rxchips*h(I(n))')/(abs(h(I(n))))^2; %selects first 'n' stongest paths

decodeBitx(n,:)=(sequence'*finger(:,:,n)); %multiplies by spreading sequence
end
decodeBit=sum(decodeBitx,1); %adds up columns, coherently combining o/p of each correlator

    y=decodeBit/(seqLen*numFingers);
     
% demodulation
    % ------------
    
    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(find(ipHatRe>max(alphaRe))) = max(alphaRe);
    ipHatRe(find(ipHatRe<min(alphaRe))) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(find(ipHatIm>max(alphaIm))) = max(alphaIm);
    ipHatIm(find(ipHatIm<min(alphaIm))) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % counting errors for real and imaginary
    nBitErr(pktX) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;
    erPkt(pktX)=size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1);
    pktInEr=erPkt>0; %returns 1 if current packet has an error in

    end
    erTot(ii)=erTot(ii)+sum(nBitErr); %total errors in all packets

    simBer(ii)=(erTot(ii)/(N*k*numPkts)); %bit error rate
    erTot(ii)=sum(erPkt); %total errors in all packets
    totPktEr(ii)=sum(pktInEr); %total packets in error
    totPktErRate(ii)=totPktEr(ii)/numPkts; %packet error rate
    end

    theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));

    EbNoLin=10.^(Eb_N0_dB/10);
    gamma_c=EbNoLin*k;
    ber = 1/2 * ( 1 - sqrt(gamma_c/k./(1+gamma_c/k)) );
    close all; figure
    semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
    hold on
    semilogy(Eb_N0_dB,ber,'o-')
    hold on
    semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
    axis([0 15 10^-5 1])
    grid on
    legend('theory AWGN','theory Rayleigh', 'simulation');
    xlabel('Eb/No, dB')
    ylabel('Bit Error Rate')

    title('Bit error probability curve for 16-QAM with DSSS, Rayleigh single tap and RAKE rx')

    figure
    semilogy(Eb_N0_dB,totPktErRate,'mx-','LineWidth',2);
    axis([0 15 10^-5 1])
    grid on
    legend('simulation');
    xlabel('Eb/No, dB')
    ylabel('Packet Error Rate')
    title('Packet error probability curve for 16-QAM with DSSS, Rayleigh single tap and RAKE rx')