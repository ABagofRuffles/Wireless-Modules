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