%Tema 9: Generare constelatie QAM pe purtatoare (alegem patratica),
%extragere purtatoare prin ridicare la puterea 4, filtrare, limitare si
%aplicare semnal extrasla un circuit PLL
clear all;
n = 4; % numarul de biti/simbol
M = 2^n; % numarul de fazori (puncte in constelatie)
L = sqrt(M); % nivele per constelatie patratica
A0 = 1; % unitatea elementara a amplitudinii

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
