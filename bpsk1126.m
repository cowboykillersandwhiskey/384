% BPSK Modulation and Demodulation with Enhanced Plotting
clc;
close all;
clear;

% Given binary data for transmission
data = [1 0 0 0 0 1 1 1]; 

for k = 1:1:10
    data = [data data]; % doubles the length of data
end                     % this gives us more bits of data to play with

n = length(data);       % Number of bits to transmit
A = 3;                  % Amplitude of transmitted signal 
freq = 4e6;             % Frequency of carrier (4 MHz)
fs = 1e8;               % Sampling frequency (100 MHz)
B = 10;                 % Background noise amplitude

% Time vector for one bit duration
t_bit = 0:1/fs:1e-6;    % Duration per bit, 1 microsecond per bit

% Carrier signal for one bit
carrier = A * cos(2 * pi * freq * t_bit);

% Modulate the signal (BPSK Modulation)
modulated_signal = [];
for i = 1:n
    if data(i) == 1
        modulated_bit = carrier;
    else
        modulated_bit = -carrier;
    end
    modulated_signal = [modulated_signal modulated_bit];
end

% Adding Gaussian Noise
SNR_dB = 2; % Desired Signal-to-Noise Ratio in dB
sigma = B / sqrt(2 * 10^(SNR_dB / 10)); % Noise standard deviation
noise = sigma * randn(size(modulated_signal)); % Generate noise

noisy_signal = modulated_signal + noise; % Add noise to signal

% Demodulation
received_data = [];
constellation_points = []; % Store constellation points for plotting
for i = 1:n
    % Extract the ith bit segment
    bit_segment = noisy_signal((i-1)*length(t_bit)+1:i*length(t_bit));
    
    % Multiply by carrier to demodulate 
    demodulated_bit = bit_segment .* carrier;
    
    % Integrate over the bit duration and make decision
    decision_metric = sum(demodulated_bit);
    constellation_points = [constellation_points decision_metric]; % Save for constellation plot
    if decision_metric > 0
        received_data(i) = 1;
    else
        received_data(i) = 0;
    end
end

% Calculate Bit Error Rate (BER)
errors = sum(data ~= received_data);
BER = errors / n;

% Display results
disp("Transmitted Data:");
disp(data);
disp("Received Data:");
disp(received_data);
disp("Bit Error Rate (BER):");
disp(BER);

% FFT of the received signal
NFFT = 2^nextpow2(length(noisy_signal)); % FFT length (next power of 2)
received_fft = fft(noisy_signal, NFFT).^2;
f = fs * (0:(NFFT/2)) / NFFT; % Frequency vector for one-sided spectrum

% Plot the signals, power, constellation, and FFT magnitude
figure;

% Adjust global figure properties
set(gcf, 'Color', 'w', 'Position', [100, 100, 1200, 800]); % Set background to white and resize

% Original Data
subplot(6, 1, 1);
stem(data(1:16), 'filled', 'LineWidth', 1.5, 'MarkerSize', 8);
title('Original Data', 'FontSize', 14, 'FontWeight', 'Bold');
xlabel('Bit Index', 'FontSize', 12); ylabel('Value', 'FontSize', 12);
ylim([-0.5 1.5]); xticks(1:16);
grid on;

% Modulated Signal
subplot(6, 1, 2);
plot(modulated_signal(1:1600), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]); % Blue color
title('Modulated Signal (BPSK)', 'FontSize', 14, 'FontWeight', 'Bold');
xlabel('Sample Index', 'FontSize', 12); ylabel('Amplitude', 'FontSize', 12);
ylim([-1.5*A 1.5*A]);
grid on;

% Noisy Signal
subplot(6, 1, 3);
plot(noisy_signal(1:1600), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]); % Orange color
title('Noisy Signal', 'FontSize', 14, 'FontWeight', 'Bold');
xlabel('Sample Index', 'FontSize', 12); ylabel('Amplitude', 'FontSize', 12);
ylim([-1.5*A 1.5*A]);
grid on;

% % Constellation Diagram
% subplot(6, 1, 4);
% scatter(real(constellation_points), zeros(size(constellation_points)), 50, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]); % Green color
% title('Constellation Diagram', 'FontSize', 14, 'FontWeight', 'Bold');
% xlabel('In-Phase Component', 'FontSize', 12); ylabel('Quadrature Component', 'FontSize', 12);
% xlim([-A*1.5 A*1.5]);
% grid on;

% FFT Magnitude Plot
subplot(6, 1, 5);
plot(f, abs(received_fft(1:NFFT/2+1)), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]); % Blue color
title('Relative Power Spectrum of Transmitted Signal', 'FontSize', 14, 'FontWeight', 'Bold');
xlabel('Frequency (Hz)', 'FontSize', 12); ylabel('|FFT(Received Signal)|', 'FontSize', 12);
xlim([0 max(f)]);
grid on;

% % Bit Error Rate (BER)
% subplot(6, 1, 6);
% bar(1, BER, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
% title('Bit Error Rate (BER)', 'FontSize', 14, 'FontWeight', 'Bold');
% xlabel('Simulation', 'FontSize', 12); ylabel('BER', 'FontSize', 12);
% ylim([0 1]);
% grid on;

% Adjust subplot spacing
sgtitle('BPSK Modulation and Analysis', 'FontSize', 19, 'FontWeight', 'Bold'); % Add global title
