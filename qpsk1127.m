clc;
clear;
close all;

% Given binary data for transmission
data = [1 0 0 0 0 1 1 1];

% Extend the data length for testing
for k = 1:1:10
    data = [data data]; % doubles the length of data
end

% System parameters
% p = 10;               % Transmitter power
% z = 73;                 % Typical 1/2 wave dipole impedance
n = length(data);       % Number of bits to transmit
A = 3;          % Amplitude of transmitted signal
freq = 4e6;             % Frequency of carrier (4 MHz)
fs = 1e8;               % Sampling frequency (100 MHz)
B = 10;              % Background noise amplitude

% Time vector for one symbol (2 bits per symbol)
t_symbol = 0:1/fs:2e-6; % Duration of one symbol (2 microseconds)

% Carrier signals for QPSK
carrier_I = cos(2 * pi * freq * t_symbol); % In-phase carrier
carrier_Q = sin(2 * pi * freq * t_symbol); % Quadrature carrier

% Modulation (QPSK)
modulated_signal = [];
for i = 1:2:n
    % Map data bits to QPSK symbol
    if data(i) == 1 && data(i+1) == 0
        I = A / sqrt(2); % Symbol 45° (binary 10)
        Q = A / sqrt(2);
    elseif data(i) == 1 && data(i+1) == 1
        I = -A / sqrt(2); % Symbol 135° (binary 11)
        Q = A / sqrt(2);
    elseif data(i) == 0 && data(i+1) == 1
        I = -A / sqrt(2); % Symbol 225° (binary 01)
        Q = -A / sqrt(2);
    else
        I = A / sqrt(2); % Symbol 315° (binary 00)
        Q = -A / sqrt(2);
    end
    % Modulate the carrier with I and Q components
    modulated_symbol = I * carrier_I - Q * carrier_Q;
    modulated_signal = [modulated_signal modulated_symbol];
end

% Adding Gaussian noise
SNR_dB = .3;                           % Desired Signal-to-Noise Ratio in dB
sigma = B / sqrt(2 * 10^(SNR_dB / 10)); % Noise standard deviation
noise = sigma * randn(size(modulated_signal)); % Generate noise

noisy_signal = modulated_signal + noise; % Add noise to signal

% Demodulation
received_data = [];
constellation_points = []; % Store constellation points for plotting
for i = 1:n/2
    % Extract the ith symbol segment (two bits per symbol)
    symbol_segment = noisy_signal((i-1)*length(t_symbol)+1:i*length(t_symbol));
    
    % Demodulate the I and Q components
    I_component = sum(symbol_segment .* carrier_I) / length(t_symbol);
    Q_component = -sum(symbol_segment .* carrier_Q) / length(t_symbol);
    
    % Normalize I and Q components
    I_component = I_component / A;
    Q_component = Q_component / A;
    
    % Save constellation points for plotting
    constellation_points = [constellation_points; I_component + 1j*Q_component];
    
    % Decision thresholds for QPSK
    if I_component > 0 && Q_component > 0
        received_data = [received_data 1 0]; % Symbol 45° (binary 10)
    elseif I_component < 0 && Q_component > 0
        received_data = [received_data 1 1]; % Symbol 135° (binary 11)
    elseif I_component < 0 && Q_component < 0
        received_data = [received_data 0 1]; % Symbol 225° (binary 01)
    else
        received_data = [received_data 0 0]; % Symbol 315° (binary 00)
    end
end

% Calculate Bit Error Rate (BER)
errors = sum(data ~= received_data);
BER = errors / n;

% Display results
disp("Transmitted Data:");
disp(data(1:16)); % Show the first 16 bits for brevity
disp("Received Data:");
disp(received_data(1:16)); % Show the first 16 bits for brevity
disp("Bit Error Rate (BER):");
disp(BER);

% FFT of the received signal
NFFT = 2^nextpow2(length(noisy_signal)); % FFT length (next power of 2)
received_fft = fft(noisy_signal, NFFT).^2;
f = fs * (0:(NFFT/2)) / NFFT; % Frequency vector for one-sided spectrum

% Plot the signals, constellation diagram, and FFT magnitude
figure;

subplot(6,1,1);
stem(data(1:16), 'filled'); title('Original Data');
grid on;

subplot(6,1,2);
plot(modulated_signal(1:1600)); title('Modulated Signal (QPSK)');
ylim([-1.5*A 1.5*A]); grid on;

subplot(6,1,3);
plot(noisy_signal(1:1600)); title('Noisy Signal');
ylim([-1.5*A 1.5*A]); grid on;

subplot(6,1,4);
scatter(real(constellation_points), imag(constellation_points), 'filled');
title('Constellation Diagram');
xlabel('In-Phase'); ylabel('Quadrature');
grid on;

subplot(6,1,5);
plot(f, abs(received_fft(1:NFFT/2+1))); 
title('Relative Power Spectrum of Transmitted Signal');
xlabel('Frequency (Hz)'); ylabel('|FFT(Received Signal)|');
grid on;

subplot(6,1,6);
bar(1, BER); title('Bit Error Rate (BER)');
ylim([0 0.5]);
grid on;
