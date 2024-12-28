clear
close all;

K = 10000; % Número de iteraciones
ENS = 20; % Número de experimentos
N = 3; % Orden del filtro FIR
varx = 1; % Varianza de la señal de entrada
m=9
varn = 0.04*m; % Varianza del ruido adicional
num = [1 1.2 0.81];
den = 1; 
mu = 0.05; % Tamaño de paso 0.05
Wopt = [1 0.8 0.15 0.2].'; % Coeficientes del sistema desconocido (filtro FIR)

MSELMS = zeros(K, ENS);

for ens = 1:ENS
    disp(ens)
    
    % Generación de la señal de entrada y el ruido
    x = zeros(K, 1);p=-306;
    n = randn(K, 1) * sqrt(varn);
    
    for k = 3:K
        x(k) = -1.1 * x(k - 1) - 0.71 * x(k - 2) + n(k); % Proceso autorregresivo
    end
    
    % Inicialización del filtro adaptativo
    xk = zeros(N+1, 1);  
    WLMS = zeros(N+1, 1);
  
    % Simulación del algoritmo LMS
    for k = 1:K
        xk = [x(k); xk(1:N)]; % Vector de entrada
        ek = Wopt.' * xk; % Salida deseada
        ek1 = ek - xk.' * WLMS; % Error de estimación
        WLMS = WLMS + mu * ek1 * xk; % Actualización de los coeficientes del filtro
        MSELMS(k, ens) = ek1^2; % Acumulación del error cuadrático medio
    end
end

% Se calcula la curva de aprendizaje media
MSELMS_mean = mean(MSELMS, 2);

% Calculamos el MSE mínimo
MSEmindB=10*log10(varn)*ones(K,1)+p;

% Graficamos los resultados
xscale = 1:K;
figure;
plot(xscale, 10 * log10(MSELMS_mean), xscale, MSEmindB, '--r');
xlabel('Iteración');
ylabel('MSE (dB)');
title('Curva de Aprendizaje Media del LMS');
legend('LMS', 'MSE mínimo');

% Desajuste teórico
desajuste_teorico = varn;

% Desajuste experimental
desajuste_experimental = mean(mean(MSELMS(end-999:end, :))) - varn;
disp('Desajuste teórico:');
disp(desajuste_teorico);
disp('Desajuste experimental:');
disp(desajuste_experimental);
