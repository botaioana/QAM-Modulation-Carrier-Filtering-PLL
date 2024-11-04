%% Tema 9: Generare constelatie QAM pe purtatoare (alegem patratica),
% extragere purtatoare prin ridicare la puterea 4, filtrare, limitare si
% aplicare semnal extras la un circuit PLL
clear all;
n = 4; % numarul de biti/simbol
M = 2^n; % numarul de fazori (puncte in constelatie)
L = sqrt(M); % nivele per constelatie patratica
A0 = 10; % unitatea elementara a amplitudinii
N = 1000; % Numar de simboluri de transmis

% Generarea nivelelor de amplitudine pe axele I si Q
Ik = zeros(1, L);
Qk = zeros(1, L);
for i = 1:L
    Ik(i) = (2*(i-1) + 1 - L) * A0;
    Qk(i) = (2*(i-1) + 1 - L) * A0;
end

% Generarea codului Gray pentru axele I si Q
gray_I = generate_gray(log2(L)); % sau n/2
gray_Q = generate_gray(log2(L));

% Crearea constelației QAM și afișarea simbolurilor și bitilor Gray
arr = zeros(1, M);
labels = cell(1, M); % Etichete pentru codurile Gray
k = 1;
for i = 1:L
    for j = 1:L
        arr(k) = complex(Ik(i), Qk(j)); % Punctul complex din constelatie
        labels{k} = strcat(gray_I{i}, gray_Q{j}); % Bitii Gray asociati
        k = k + 1;
    end
end

% Plotarea constelatiei QAM
figure;
plot(real(arr), imag(arr), 'o'); % Doar punctele
grid on;
xlabel('Ik = Ak * cos(\Phi_k)'), ylabel('Qk = Ak * sin(\Phi_k)');
title(sprintf('Constelatia %d-QAM patratica', M));

% Adaugarea bitilor Gray langa fiecare punct din constelatie
for i = 1:length(arr)
    text(real(arr(i)), imag(arr(i)), labels{i}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

%% Generarea semnalului QAM

% Generare date aleatorii
data_bits = randi([0 1], N, n); % Generarea de biti aleatori (N simboluri, fiecare simbol n biti)
symbol_indices = bi2de(data_bits, 'left-msb') + 1; % Maparea bitilor la indicii simbolurilor
tx_symbols = arr(symbol_indices); % Asocierea simbolurilor din constelatie

% Parametri pentru semnalul modulat
fc = 1 * 1e3; % Frecventa purtatoare
fs = 15e3; % Frecventa de esantionare (increased to avoid aliasing issues)
t = (0:N-1) / fs; % Vector de timp

% Modulare QAM (componentele In-phase si Quadrature)
I_signal = real(tx_symbols) .* cos(2 * pi * fc * t);
Q_signal = imag(tx_symbols) .* sin(2 * pi * fc * t);
qam_signal = I_signal - Q_signal; % Semnalul QAM modulat

% Plotare semnal modulat QAM
figure;
plot(t, qam_signal);
title('Semnal QAM Modulat');
xlabel('Timp (s)');
ylabel('Amplitudine');

%% Raise the Signal to the Fourth Power
carrier = qam_signal.^4;
carrier = carrier - mean(carrier); % Remove DC offset

% Manually assign f_p4 (the frequency after raising to the fourth power)
f_p4 = 4 * fc;

% Chebyshev Bandpass Filter around the assigned f_p4
[b, a] = cheby1(4, 0.5, [f_p4 * 0.98, f_p4 * 1.02] / (fs / 2), 'bandpass'); % Design Chebyshev filter
carrier_iso = filter(b, a, carrier);

% Plot the filtered carrier signal and its Fourier Transform
f_carrier_iso = fft(carrier_iso);
f = (0:length(f_carrier_iso)-1) * (fs / length(f_carrier_iso)); % Frequency vector

figure;
subplot(221);
plot(t, carrier);
title('Extracted Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(222);
plot(f, abs(fft(carrier)));
title('Fourier Transform of Extracted Carrier Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);

subplot(223);
plot(t, carrier_iso);
title('Chebyshev Filtered Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(224);
plot(f, abs(f_carrier_iso));
title('Fourier Transform of Chebyshev Filtered Carrier Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);

%% Limiting the Filtered Signal
limited_carrier = sign(carrier_iso); % Binary limiter to isolate phase

% Plot the limited carrier signal
figure;
subplot(2, 1, 1);
plot(t, limited_carrier);
title('Limited Carrier Signal (Input to PLL)');
xlabel('Time (s)');
ylabel('Amplitude');

%% PLL (Phase-Locked Loop) for Carrier Recovery
pll = comm.CarrierSynchronizer( ...
    'Modulation', 'QAM', ...
    'SamplesPerSymbol', 1, ...
    'DampingFactor', 0.707, ... % Critical damping
    'NormalizedLoopBandwidth', 0.01); % Bandwidth to lock onto carrier

% Apply the limited signal to the PLL
recovered_carrier = pll(limited_carrier.'); % Ensure column vector for PLL

% Plot the recovered carrier signal from PLL
subplot(2, 1, 2);
plot(t, real(recovered_carrier));
title('PLL Output (Recovered Carrier)');
xlabel('Time (s)');
ylabel('Amplitude');

% Optional: Fourier Transform of PLL output for analysis
f_recovered_carrier = fft(recovered_carrier);
figure;
plot(f, abs(f_recovered_carrier));
title('Fourier Transform of Recovered Carrier Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);

%% Performance Analysis

% Carrier Frequency Offset (CFO) Error
estimated_frequency = f_p4 / 4;
recovered_frequency = abs(mean(diff(angle(recovered_carrier)))) * fs / (2 * pi);
CFO_error = abs(estimated_frequency - recovered_frequency);
fprintf('Carrier Frequency Offset (CFO) Error: %.2f Hz\n', CFO_error);

% Phase Error Calculation
ideal_carrier = cos(2 * pi * estimated_frequency * t);
phase_error = mean(abs(angle(ideal_carrier) - angle(recovered_carrier)));
fprintf('Mean Phase Error: %.4f radians\n', phase_error);

% Signal-to-Noise Ratio (SNR)
% Calculate signal and noise power
signal_power = mean(abs(qam_signal).^2);  % Scalar power of the transmitted QAM signal
noise_power = mean(abs(recovered_carrier(:) - ideal_carrier(:)).^2);  % Scalar power of the noise

% Calculate SNR in dB
SNR = 10 * log10(signal_power / noise_power);
fprintf('SNR of Recovered Signal: %.2f dB\n', SNR);

% Plot Spectral Purity of Recovered Carrier
figure;
plot(f, abs(f_recovered_carrier));
title('Spectral Purity of Recovered Carrier');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 fs/2]);
