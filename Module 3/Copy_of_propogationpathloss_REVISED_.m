% Free Space Propagation Loss with noise consideration
clc;
close all;
clear all;

fc = input('Enter carrrier frequency(MHz): ');  % carrier frequency
c = 299.792458; % speed of light with implicit 10e6 (Mega) multiplication
d = 1:1:10000;  % distance vector for wave to travel in meters
d0 = 1;
Pt = 50; % typical cell tower transmission power in urban area 50 Watts
Gt = 1; % isotropic antenna transmitter gain
Gr = 1; % isotropic antenna reciever gain
n = 4;  % path loss exponent for obstruction in building from 4-6
mu=0;   % mu
sigma=6;  % standard deviation


Ht = 30:0.007:99.993; % Hieght of Antenna

% Mobile station antenna height gain factor
Ghr = 20*log10(3/3);
Ghr2 = 20*log10(4/3);
Ghr3 = 20*log10(5/3);
Ghr4 = 20*log10(6/3);
Ghr5 = 20*log10(7/3);
Ghr6 = 20*log10(8/3);

Amu = 35;   % Median attenuation relative to free space
Garea = 9;  % Gain due to type of environment, 9 for urban
Ght = 20*log10(Ht/200); % Base station antenna height gain factor

Pr = ((c./(4*pi*d0*fc)).^n).*(Pt*Gt*Gr)*((d0./d).^n);   % Pr formula

noise = lognrnd(mu,sigma,size(d));
PrSh = Pr.*noise;  % Pr formula + Shadowing

rayNoise = abs(sigma*randn(size(d))+1i*sigma*randn(size(d)));
PrShRy = PrSh.*rayNoise; % Pr formula + Shadowing + Rayliegh fading


figure(1);
semilogx(d,10*log10(Pr),'b',d,10*log10(PrSh),'r',d,10*log10(PrShRy),'m');
xlabel('Magnitude of distance in kilometers'); 
ylabel('Recieved Power in dB'); 
title('Free space model'); 
grid on;
legend('Path Loss','Path Loss with Shadowing','Path Loss with Shadowing and Multipath Fading');

figure(2)
subplot(2,1,1)
semilogx(d,10*log10(noise),'.r')
title('Noise scatterplot for lognormal shadowing');
grid on
xlabel('Distance in kilometers')
ylabel('Lognormal Shadowing in dB')

subplot(2,1,2)
semilogx(d,10*log10(rayNoise),'.k')
title('Noise scatterplot for Rayleigh Fading');
grid on
xlabel('Distance in kilometers')
ylabel('Rayleigh distribution random variable in dB')

figure(3)
semilogx(Ht,10*log10(Pr)+Amu-Ght-Ghr-Garea,'k',Ht,10*log10(Pr)+Amu-Ght-Ghr2-Garea,'c',Ht,10*log10(Pr)+Amu-Ght-Ghr3-Garea,'r',Ht,10*log10(Pr)+Amu-Ght-Ghr4-Garea,'b',Ht,10*log10(Pr)+Amu-Ght-Ghr5-Garea,'g',Ht,10*log10(Pr)+Amu-Ght-Ghr6-Garea,'m');
title('Okumura Model for Various Receiver Antenna Heights');
grid on
xlabel('Transmitter antenna Height in meters')
ylabel('Propagation Path loss in dB')
legend('3 meters','4 meters','5 meters','6 meters','7 meters','8 meters');

% QPSK Modulation and Demodulation without consideration of noise

% Using Matlab, generate the QPSK waveform with the following parameters: 
% Symbol rate: 2000 symbols/second; 
% carrier frequency: 10000 Hz, 
% bit sequence to be transmitted: 1 0 1 1 0 0 1 1. 
% The in-phase and quadrature waveforms should be separately plotted. 
% Use sufficiently high sampling rate to ensure smooth waveforms.

data = [1 0 1 1 0 0 1 1]; % information to be transmitted

figure(4)
stem(data, 'linewidth', 3), grid on;
title('Information before Transmiting');
axis([ 0 (length(data) + 1) 0 1.5]);

data_NZR = 2 * data-1; % Data Represented at NZR form for QPSK modulation
s_p_data = reshape(data_NZR, 2, length(data)/2)  % S/P convertion of data


br = 1e3; %Let the transmission bit rate be 1000 bps for 2000 symbols per second
T = 1/br; % bit duration
t = T/99:T/99:T; % Time vector for one bit information



% XXXXXXXXXXXXXXXXXXXXXXX QPSK modulation  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
y = [];
yInPhase = [];
yQuadrature = [];

for(i = 1: length(data)/2)
    y1 = s_p_data(1, i) * cos(2*pi*fc*t); % inphase component
    y2 = s_p_data(2, i) * sin(2*pi*fc*t) ;% Quadrature component
    
    yInPhase = [yInPhase y1]; % inphase signal vector
    yQuadrature = [yQuadrature y2]; %quadrature signal vector
    
    y = [y y1+y2]; % modulated signal vector
end

Tx_sig = y; % transmitting signal after modulation
tt = T/99:T/99:(T*length(data))/2;

figure(5)
subplot(3,1,1);
plot(tt, yInPhase, 'linewidth', 3), grid on;
title('Waveform for inphase component in QPSK modulation');
xlabel('Time[sec]');
ylabel('Amplitude[Volts]');

subplot(3,1,2);
plot(tt, yQuadrature, 'linewidth', 3), grid on;
title('Waveform for Quadrature component in QPSK modulation');
xlabel('Time [sec]');
ylabel('Amplitude [Volts]');


subplot(3,1,3);
plot(tt, Tx_sig, 'r', 'linewidth', 3), grid on;
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('Time[sec]');
ylabel('Amplitude[Volts]');



% QPSK demodulation
Rx_data = [];
Rx_sig = Tx_sig; % Received signal

for(i=1: 1:length(data)/2)

    %%XXXXXX inphase coherent vector XXXXXXX
    Z_in = Rx_sig((i-1)*length(t)+1:i*length(t)).*cos(2*pi*fc*t); 
    % above line indicates multiplication of received & inphase carred signal
    
    Z_in_intg = (trapz(t,Z_in)) * (2/T);% integration using trapizodial rull
    if(Z_in_intg > 0) % Decession Maker
        Rx_in_data = 1;
    else
        Rx_in_data = 0; 
    end
    
    %%XXXXXX Quadrature coherent vector XXXXXX
    Z_qd = Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*fc*t);
    %above line indicates multiplication of received & Quadphase carred signal
    
    Z_qd_intg = (trapz(t,Z_qd)) * (2/T);%integration using trapizodial rull
        if (Z_qd_intg > 0)% Decession Maker
        Rx_qd_data = 1;
        else
        Rx_qd_data = 0; 
        end
        
        
        Rx_data = [Rx_data  Rx_in_data  Rx_qd_data]; % Received Data vector
end

figure(6)
subplot(3,1,1);
plot(tt, Z_in, 'linewidth', 3), grid on;
title('Waveform for inphase component in QPSK modulation');
xlabel('Time[sec]');
ylabel('Amplitude[Volts]');

subplot(3,1,2);
plot(tt, Z_qd, 'linewidth', 3), grid on;
title('Waveform for Quadrature component in QPSK modulation');
xlabel('Time [sec]');
ylabel('Amplitude [Volts]');


subplot(3,1,3);
plot(tt, Rx_sig, 'r', 'linewidth', 3), grid on;
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('Time[sec]');
ylabel('Amplitude[Volts]');

figure(7)
stem(Rx_data, 'linewidth', 3) 
title('Information after Receiveing');
axis([0 (length(data) + 1) 0 1.5]), grid on;

% End of program
    



