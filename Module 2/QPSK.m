% QPSK Modulation and Demodulation without consideration of noise

% Using Matlab, generate the QPSK waveform with the following parameters: 
% Symbol rate: 2000 symbols/second; 
% carrier frequency: 10000 Hz, 
% bit sequence to be transmitted: 1 0 1 1 0 0 1 1. 
% The in-phase and quadrature waveforms should be separately plotted. 
% Use sufficiently high sampling rate to ensure smooth waveforms.
% Submit your waveforms with Matlab code.


clc;
clear all;
close all;

data = [1 0 1 1 0 0 1 1]; % information to be transmitted

figure(1)
stem(data, 'linewidth', 3), grid on;
title('Information before Transmiting');
axis([ 0 (length(data) + 1) 0 1.5]);

data_NZR = 2 * data-1; % Data Represented at NZR form for QPSK modulation
s_p_data = reshape(data_NZR, 2, length(data)/2)  % S/P convertion of data


br = 1e3; %Let the transmission bit rate be 1000 bps for 2000 symbols per second
fc = 10.^3; % minimum carrier frequency
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

figure(2)
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
    % above line indicat multiplication of received & inphase carred signal
    
    Z_in_intg = (trapz(t,Z_in)) * (2/T);% integration using trapizodial rull
    if(Z_in_intg > 0) % Decession Maker
        Rx_in_data = 1;
    else
        Rx_in_data = 0; 
    end
    
    %%XXXXXX Quadrature coherent vector XXXXXX
    Z_qd = Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*fc*t);
    %above line indicat multiplication ofreceived & Quadphase carred signal
    
    Z_qd_intg = (trapz(t,Z_qd)) * (2/T);%integration using trapizodial rull
        if (Z_qd_intg > 0)% Decession Maker
        Rx_qd_data = 1;
        else
        Rx_qd_data = 0; 
        end
        
        
        Rx_data = [Rx_data  Rx_in_data  Rx_qd_data]; % Received Data vector
end


figure(3)
stem(Rx_data, 'linewidth', 3) 
title('Information after Receiveing');
axis([0 (length(data) + 1) 0 1.5]), grid on;

% End of program
    



