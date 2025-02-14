eps_0 = 8.85e-12;        
mu_0 = 4 * pi * 1e-7;    
eps_r = 4.4;               
mu_r = 1;                
f = linspace(1e8, 10e9, 1000); 
d = 0.002;                
R = 1;                  
L = 0.94e-7;            
C0 = 4.31e-14;              
x_values = linspace(0, 0.7, 20); 
cases = {'linear', 'cubic'};
colors = {'r', 'b'};
SNR = 30; % Livello di rumore in dB

%% Gamma a diverse umidità
figure('Name','Gamma vs Frequenza');
hold on;

% Valori di umidità
selected_x = [ 0,0.5, 1];
for idx = 1:3
    x = selected_x(idx);
    C = C0 * (1 + x^3); % Caso cubico, alla fine sono simili la differenza viene apprezzata solo quantitativamente. 
    
    % Calcolo Gamma
   [gamma, ~] = calculate_gamma(f, C, R, L, d, eps_r, mu_r);
    
    plot(f/1e9, gamma, 'LineWidth', 1.5, 'DisplayName', sprintf('H = %.1f', x));
end

xlabel('Frequenza (GHz)');
ylabel('|\Gamma|');
title('Coefficiente di Riflessione a Diverse Umidità');
legend('Location', 'best');
grid on;
xlim([1 5]);

%% Rette di Calibrazione con Rumore
figure('Name','Calibrazione con Rumore');
hold on;

for c = 1:length(cases)
    f_res = zeros(size(x_values));
    f_res_noisy = zeros(size(x_values));
    
    for idx = 1:length(x_values)
        x = x_values(idx);
        if strcmp(cases{c}, 'linear')
            C = C0 * (1 +  x);
        else
            C = C0 * (1 +   x^3);
        end
        
        % Calcolo Gamma 
        [gamma, ~] = calculate_gamma(f, C, R, L, d, eps_r, mu_r);
        [gamma_min, idx_min] = max(gamma);
        f_res(idx) = f(idx_min);
        
        % Aggiunta rumore
        gamma_noisy = awgn(gamma, SNR, 'measured');
        [~, idx_min_noisy] = max(gamma_noisy);
        f_res_noisy(idx) = f(idx_min_noisy);
    end
    
    % Plot dati e dati corrotti
    plot(x_values, f_res/1e9, [colors{c} '-'], 'LineWidth', 2, 'DisplayName', [cases{c} ' pulito']);
    plot(x_values, f_res_noisy/1e9, [colors{c} 'o'], 'MarkerSize', 6, 'DisplayName', [cases{c} ' rumoroso']);
    errorbar(x_values,f_res_noisy/1e9,f_res/1e9/10^2) % con SNR di 40dB le barre di errore non si vedono
end

xlabel('x (normalizzato)');
ylabel('Frequenza di Ris. (GHz)');
title('Curve di Calibrazione con Rumore Gaussiano');
legend('Location', 'best');
grid on;

%% Grafico 3: Combinazione
figure('Position', [100 100 1200 500])

%Gamma corrotto
subplot(1,2,1);
hold on;
selected_x = [ 0,0.5, 1];
for x = selected_x
    C = C0 * (1 + x);
    
    [gamma, ~] = calculate_gamma(f, C, R, L, d, eps_r, mu_r);
    gamma_noisy = awgn(gamma, SNR, 'measured');
    plot(f/1e9, gamma_noisy, 'LineWidth', 1.2, 'DisplayName', sprintf('H = %.1f', x));
end
xlabel('Frequenza (GHz)');
ylabel('|\Gamma| Corrotto');
title('Gamma con Rumore');
legend('Location', 'best');
grid on;

%Rette di calibrazione
subplot(1,2,2);
hold on;

for c = 1:length(cases)
    f_res = zeros(size(x_values));
    f_res_noisy = zeros(size(x_values));
    
    for idx = 1:length(x_values)
        x = x_values(idx);
        if strcmp(cases{c}, 'linear')
            C = C0 * (1 + x);
        else
            C = C0 * (1 + x^3);
        end
        
        % Calcolo Gamma 
        [gamma, ~] = calculate_gamma(f, C, R, L, d, eps_r, mu_r);
        [gamma_min, idx_min] = max(gamma);
        f_res(idx) = f(idx_min);
        
        % Aggiunta rumore
        gamma_noisy = awgn(gamma, SNR, 'measured');
        [~, idx_min_noisy] = max(gamma_noisy);
        f_res_noisy(idx) = f(idx_min_noisy);
    end
    
    % Plot dati e dati corrotti
    plot(x_values, f_res/1e9, [colors{c} '-'], 'LineWidth', 2, 'DisplayName', [cases{c} ' pulito']);
    plot(x_values, f_res_noisy/1e9, [colors{c} 'o'], 'MarkerSize', 6, 'DisplayName', [cases{c} ' rumoroso']);
end

xlabel('x (normalizzato)');
ylabel('Frequenza di Ris. (GHz)');
title('Curve di Calibrazione con Rumore Gaussiano');
legend('Location', 'best');
grid on;

title('Confronto Calibrazioni');
grid on;

%% Funzione ausiliaria
function [gamma, f_res] = calculate_gamma(f, C, R, L, d, eps_r, mu_r)
    omega = 2 * pi * f;
    beta = omega .* sqrt(4*pi*1e-7*mu_r * 8.85e-12*eps_r);
    Zin = sqrt(4*pi*1e-7*mu_r/(8.85e-12*eps_r)) .* (sqrt(4*pi*1e-7*mu_r/(8.85e-12*eps_r)) + 1i*377*tan(beta*d)) ./ (377 + 1i*sqrt(4*pi*1e-7*mu_r/(8.85e-12*eps_r)).*tan(beta*d));
    Zrlc = R + 1i*omega*L - 1i./(omega*C);
    Ztot = 1 ./ (1./Zin + 1./Zrlc);
    gamma = abs((Ztot - 377) ./ (Ztot + 377));
    [~, idx_min] = min(gamma);
    f_res = f(idx_min);
end




