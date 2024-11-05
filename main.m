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
    Ik(i) = (2*(i-1) + 1 - L) * A0; % Calculam nivelele de amplitudine pe axa I
    Qk(i) = (2*(i-1) + 1 - L) * A0; % Calculam nivelele de amplitudine pe axa Q
end

% Generarea codului Gray pentru axele I si Q
gray_I = generate_gray(log2(L)); % sau n/2 - generare cod Gray pentru axa I
gray_Q = generate_gray(log2(L)); % generare cod Gray pentru axa Q

% Crearea constelației QAM și afișarea simbolurilor și bitilor Gray
arr = zeros(1, M); % Vector pentru stocarea punctelor din constelatie
labels = cell(1, M); % Etichete pentru codurile Gray
k = 1;
for i = 1:L
    for j = 1:L
        arr(k) = complex(Ik(i), Qk(j)); % Punctul complex din constelatie
        labels{k} = strcat(gray_I{i}, gray_Q{j}); % Bitii Gray asociati fiecarui punct
        k = k + 1;
    end
end

% Plotarea constelatiei QAM
figure;
plot(real(arr), imag(arr), 'o'); % Doar punctele din constelatie
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
fs = 15e3; % Frecventa de esantionare (pentru a evita aliasing-ul)
t = (0:N-1) / fs; % Vector de timp

% Modulare QAM (componentele In-phase si Quadrature)
I_signal = real(tx_symbols) .* cos(2 * pi * fc * t); % Componenta In-phase
Q_signal = imag(tx_symbols) .* sin(2 * pi * fc * t); % Componenta Quadrature
qam_signal = I_signal - Q_signal; % Semnalul QAM modulat

% Plotare semnal modulat QAM
figure;
plot(t, qam_signal);
title('Semnal QAM Modulat');
xlabel('Timp (s)');
ylabel('Amplitudine');

%% Ridicarea semnalului la puterea a patra
carrier = qam_signal.^4; % Ridicare la puterea a patra pentru a extrage purtatoarea
carrier = carrier - mean(carrier); % Eliminarea offset-ului DC

% Atribuirea manuala a frecventei dupa ridicare la puterea a patra
f_p4 = 4 * fc;

% Filtrare cu filtru Chebyshev in jurul frecventei f_p4
[b, a] = cheby1(4, 0.5, [f_p4 * 0.98, f_p4 * 1.02] / (fs / 2), 'bandpass'); % Design filtru Chebyshev
carrier_iso = filter(b, a, carrier); % Aplicarea filtrului asupra semnalului

% Plotarea semnalului filtrat si transformata Fourier a acestuia
f_carrier_iso = fft(carrier_iso);
f = (0:length(f_carrier_iso)-1) * (fs / length(f_carrier_iso)); % Vector de frecventa

figure;
subplot(221);
plot(t, carrier);
title('Semnal Purtatoare Extras');
xlabel('Timp (s)');
ylabel('Amplitudine');

subplot(222);
plot(f, abs(fft(carrier)));
title('Transformata Fourier a Semnalului Purtatoare Extras');
xlabel('Frecventa (Hz)');
ylabel('Magnitudine');
xlim([0 fs/2]);

subplot(223);
plot(t, carrier_iso);
title('Semnal Purtatoare Filtrat Chebyshev');
xlabel('Timp (s)');
ylabel('Amplitudine');

subplot(224);
plot(f, abs(f_carrier_iso));
title('Transformata Fourier a Semnalului Purtatoare Filtrat');
xlabel('Frecventa (Hz)');
ylabel('Magnitudine');
xlim([0 fs/2]);

%% Limitarea semnalului filtrat
limited_carrier = sign(carrier_iso); % Limitare binara pentru izolarea fazei

% Plotarea semnalului limitat
figure;
subplot(2, 1, 1);
plot(t, limited_carrier);
title('Semnal Purtatoare Limitat (Input pentru PLL)');
xlabel('Timp (s)');
ylabel('Amplitudine');

%% PLL (Phase-Locked Loop) pentru recuperarea purtatoarei
pll = comm.CarrierSynchronizer( ...
    'Modulation', 'QAM', ...
    'SamplesPerSymbol', 1, ...
    'DampingFactor', 0.707, ... % Amortizare critica
    'NormalizedLoopBandwidth', 0.01); % Largimea de banda pentru a se bloca pe purtatoare

% Aplicarea semnalului limitat la PLL
recovered_carrier = pll(limited_carrier.'); % Asigurarea ca vectorul este coloana pentru PLL

% Plotarea semnalului recuperat de la PLL
subplot(2, 1, 2);
plot(t, real(recovered_carrier));
title('Iesire PLL (Purtatoare Recuperata)');
xlabel('Timp (s)');
ylabel('Amplitudine');

% Optional: Transformata Fourier a iesirii PLL pentru analiza
f_recovered_carrier = fft(recovered_carrier);
figure;
plot(f, abs(f_recovered_carrier));
title('Transformata Fourier a Semnalului Purtatoare Recuperat');
xlabel('Frecventa (Hz)');
ylabel('Magnitudine');
xlim([0 fs/2]);

%% Analiza Performantei

% Eroarea de Offset a Frecventei Purtatoarei (CFO)
estimated_frequency = f_p4 / 4;
recovered_frequency = abs(mean(diff(angle(recovered_carrier)))) * fs / (2 * pi);
CFO_error = abs(estimated_frequency - recovered_frequency);
fprintf('Eroarea de Offset a Frecventei Purtatoarei (CFO): %.2f Hz\n', CFO_error);

% Calculul Erorii de Faza
ideal_carrier = cos(2 * pi * estimated_frequency * t);
phase_error = mean(abs(angle(ideal_carrier) - angle(recovered_carrier)));
fprintf('Eroarea Medie de Faza: %.4f radiani\n', phase_error);

% Rata Semnal-Zgomot (SNR)
% Calculul puterii semnalului si a zgomotului
signal_power = mean(abs(qam_signal).^2);  % Puterea semnalului QAM transmis
noise_power = mean(abs(recovered_carrier(:) - ideal_carrier(:)).^2);  % Puterea zgomotului

% Calculul SNR in dB
SNR = 10 * log10(signal_power / noise_power);
fprintf('SNR-ul Semnalului Recuperat: %.2f dB\n', SNR);

% Plotarea puritatii spectrale a purtatoarei recuperate
figure;
plot(f, abs(f_recovered_carrier));
title('Puritatea Spectrala a Purtatoarei Recuperate');
xlabel('Frecventa (Hz)');
ylabel('Magnitudine');
xlim([0 fs/2]);
