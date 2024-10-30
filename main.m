%% Tema 9: Generare constelatie QAM pe purtatoare (alegem patratica),
%extragere purtatoare prin ridicare la puterea 4, filtrare, limitare si
%aplicare semnal extras la un circuit PLL
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
gray_I = generate_gray(log2(L))%sau n/2
gray_Q = generate_gray(log2(L))

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
figure,
plot(real(arr), imag(arr), 'o'), % Doar punctele
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
fs = 10e3; % Frecventa de esantionare
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

%% Ridicare la puterea 4
carrier=qam_signal.^4;
carrier=carrier - mean(carrier);
% Fourier Transform of the extracted signal
f_carrier = fft(carrier);
f = (0:length(f_carrier)-1)*(fs/length(f_carrier));
[f_max, f_indice]=max(f_carrier,[],"all");
f_p4=f(f_indice)
%f_p4=4 * 1e3;
% Bandpass filter around the carrier frequency after raising to the fourth power
carrier_iso = bandpass(carrier, [f_p4 * 0.98, f_p4 * 1.02], fs);
% Fourier Transform of the filtered signal
f_carrier_iso = fft(carrier_iso);

% Define frequency axis
f = (0:length(f_carrier)-1)*(fs/length(f_carrier)); % Frequency vector

% Plot the extracted carrier signal and its Fourier Transform
figure;
subplot(221);
plot(t, carrier);
title('Extracted Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(222);
plot(f, abs(f_carrier)); % Use abs() to get the magnitude
title('Fourier Transform of Extracted Carrier Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Plot the bandpass filtered carrier signal and its Fourier Transform
subplot(223);
plot(t, carrier_iso); % Plot the filtered signal
title('Bandpass Filtered Carrier Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(224);
plot(f, abs(f_carrier_iso)); % Use abs() to get the magnitude
title('Fourier Transform of Filtered Carrier Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
