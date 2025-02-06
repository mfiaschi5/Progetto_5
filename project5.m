clear all; close all; clc;

eps_0 = 8.85e-12;        
mu_0 = 4 * pi * 1e-7;    

eps_r = 4;               
mu_r = 1;                
f = linspace(1e8, 10e9, 1000); 
d = 0.01;                
R = 50;                  
L = 0.94e-10;
x=linspace(0.1,10,100);
C = 3e-11*x^3;               

f0 = 1 / (2 * pi * sqrt(L * C));
fprintf('Frequenza di risonanza f0 = %.2f GHz\n', f0 / 1e9);

omega = 2 * pi * f;                         
beta = omega .* sqrt(mu_0 * mu_r * eps_0 * eps_r); 
Z0 = 377 * ones(size(f));                   
Zd = sqrt((mu_0 * mu_r) / (eps_0 * eps_r)) * ones(size(f));

j = 1i;
Zin = Zd .* (Zd + j .* Z0 .* tan(beta .* d)) ./ (Zd + j .* Z0 ./ tan(beta .* d));
Zrlc = R + j .* omega .* L - j ./ (omega .* C);
Ztot = 1 ./ (1 ./ Zin + 1 ./ Zrlc);
gamma = abs((Zd - Ztot) ./ (Zd + Ztot));
gamma = gamma + 0.01*randn(size(f));
figure;
plot(f / 1e9, gamma, 'LineWidth', 1.5);
xlabel('Frequenza (GHz)');
ylabel('Coefficiente di riflessione |\Gamma|');
title('Coefficiente di riflessione della struttura');
grid on;


