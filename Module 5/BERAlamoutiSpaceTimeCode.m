% Script for an implementation of a Rayleigh fading channel with Alamouti Space Time Block Coding
% for 4 transmit antenna real modulation.
% BER for a BPSK modulated signal is computed.


clc;
close all;
clear all;

NumberofBits = 10^6; % number of bits or symbols

Eb_N0_dB = [0:25]; % multiple Eb/N0 values


for ii = 1:length(Eb_N0_dB)
    % Transmitter
    ip = rand(2,NumberofBits) > 0.5; % generating 0,1 with equal probability
    currentBit = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0

    % Alamouti STBC
    sCode = zeros(4,NumberofBits);
    sCode(:,1:2:end) = (1/sqrt(2))*reshape(currentBit,4,NumberofBits/2); % [x1 x2  ...]
    sCode(:,2:2:end) = (1/sqrt(2))*(kron(ones(2,NumberofBits/2),[-1;1]).*flipud(reshape(conj(currentBit),4,NumberofBits/2))); % [-x2* x1* ....]

    h = 1/sqrt(2)*[randn(2,NumberofBits) + 1i*randn(2,NumberofBits)]; % Rayleigh channel
    hMod = kron(reshape(h,4,NumberofBits/2),ones(1,2)); % repeating the same channel for two symbols
    n = 1/sqrt(2)*[randn(4,NumberofBits) + 1i*randn(4,NumberofBits)]; % white gaussian noise, 0dB variance

    % Channel and noise Noise addition
    y = sum(hMod.*sCode,2) + 10^(-Eb_N0_dB(ii)/20)*n;

    % Receiver
    yMod = kron(reshape(y,4,NumberofBits/2),ones(1,2)); % [y1 y1 ... ; y2 y2 ...]
    yMod(2,:) = conj(yMod(2,:)); % [y1 y1 ... ; y2* y2*...]

    % forming the equalization matrix
    hEq = zeros(4,NumberofBits);
    hEq(:,[1:2:end]) = reshape(h,4,NumberofBits/2); % [h1 0 ... ; h2 0...]
    hEq(:,[2:2:end]) = kron(ones(2,NumberofBits/2),[1;-1]).*flipud(reshape(h,4,NumberofBits/2)); % [h1 h2 ... ; h2 -h1 ...]
    hEq(1,:) = conj(hEq(1,:)); %  [h1* h2* ... ; h2 -h1 .... ]
    hEqPower = sum(hEq.*conj(hEq),1);
    yHat = sum(hEq.*yMod,1)./hEqPower; % [h1*y1 + h2y2*, h2*y1 -h1y2*, ... ]
    yHat(2:2:end) = conj(yHat(2:2:end));

    % receiver - hard decision decoding
    ipHat = real(yHat) > 0;

    % counting the errors
    nErr(ii) = size(find([ip- ipHat]),2);
end

simBer = nErr/NumberofBits; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5));
p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p));
pAlamouti = 1/2 - 1/2*(1+2./EbN0Lin).^(-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti));

close all
figure(2)
semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2);
semilogy(Eb_N0_dB,theoryBerAlamouti_nTx2_nRx1,'c+-','LineWidth',2);
semilogy(Eb_N0_dB,simBer,'mo-','LineWidth',2);
axis([0 25 10^-5 0.5])
grid on
legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'theory (nTx=2, nRx=1, Alamouti)', 'sim (nTx=2, nRx=1, Alamouti)');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for BPSK modulation with Alamouti STBC (Rayleigh channel)');
