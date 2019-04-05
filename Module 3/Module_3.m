clc;
close all;
clear all;

% Propogation Path Loss
Dstart = 1; %Start Distance = 1m
Dresolution = 1; %1 meter resolution
Values = 100000; %total values
D =(Dstart:Dresolution:(Dstart+Dresolution*(Values-1))); %Distance Range
F = 2.4 * 10^9; %Frequency (2.4 GHz)
c = 299792458; %Speed of Light (m/s)
SIG = 100; %Tx Signal (100dB)
G = 10; %Gain = (10dB)
mu = 0; %mean = 0 (normalized about 0)
sigma = 6; %Shadowing Variance (dB) (6-12dB normal)
B = 0.5; %Rayleigh b
i = 1;
while i < (Values + 1)
Dt = Dstart + Dresolution * (i-1); %Distance (1m to 1000m)
PL(i) = log((c./(4*pi*Dresolution*F)).^2) + log((Dresolution./i).^2); %Path Loss (dB)
S(i) = normrnd(mu,sigma); %Shadowing (Random Value)
PLS(i) = PL(i) + log(S(i)); %Path Loss + Shadowing (dB)
R(i) = raylrnd(i); %Rayleigh Distribution Fading (dB)
PLSR(i) = log(R(i)) + PLS(i); %Multipath Fading + Path Loss + Shadowing (dB)
i = i + 1;
end

% Figure Generation
X_Values = log(D);
Y1_Values = PL;
Y2_Values = PLS;
Y3_Values = PLSR;
figure(1);
plot(X_Values, Y1_Values,'DisplayName', 'Path Loss');
hold on;

xlabel('Log(Distance)'); 
grid on; 
legend('Location', 'eastoutside');
hold off;

% Individual Figures
figure(2);
plot1 = subplot(3,1,1);
X_Values = D;
Y_Values = PL;
plot(X_Values, Y_Values);
ylabel('SNR (dB)'); 
xlabel('Distance (m)');
title('Propogation Path Loss vs. Distance');

plot2 = subplot(3,1,2);
X_Values = D;
Y_Values = S;
plot(X_Values, Y_Values);
ylabel('Shadowing (dB)'); 
xlabel('Distance (m)');
title('Shadowing Path Loss vs. Distance');

plot3 = subplot(3,1,3);
X_Values = D;
Y_Values = R;
plot(X_Values, Y_Values);
ylabel('Rayleigh Fading (dB)'); 
xlabel('Distance (m)');
title('Rayleigh Fading Loss vs. Distance');
