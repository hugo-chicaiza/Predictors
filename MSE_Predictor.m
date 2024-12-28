clear all;
close all;
clc;
% Parámetros
N = 100; % Longitud de la señal aleatoria

% Generar señal aleatoria
x = randn(N, 1);

% Sumar ruido blanco (AWGN)
x_with_noise = awgn(x, 10); % SNR de 10dB

% Calcular la varianza del ruido blanco
sigma_x_sq = mean(x_with_noise.^2);

% Calcular la función de autocorrelación
r_x_tau = xcorr(x_with_noise, 'biased'); % La función xcorr con 'biased' normaliza por N
r_x_tau = r_x_tau(N:end); % Tomar solo los valores positivos

% Coeficiente de autocorrelación
R = toeplitz(r_x_tau(1:end-1));

% Vector de autocorrelación cruzada
r_x = r_x_tau(2:end);

% Coeficientes del filtro de Wiener
a_opt = R\r_x;

% Construir la predicción
x_prediction = filter([1, -a_opt'], 1, x_with_noise);

% Calcular el error de predicción d(k)
prediction_error_signal = x - x_prediction;

% Calcular el error cuadrático medio (MSE)
mse = mean(prediction_error_signal.^2);

% Calcular el factor de autocorrelación
factor_autocorrelacion = r_x(1) / sigma_x_sq;
% Inicializar el vector de errores cuadráticos medios
mse_values = zeros(N, 1);

% Calcular el error de predicción d(k) y el MSE para diferentes instantes de tiempo
for k = 1:N
    % Construir la predicción para el instante de tiempo k
    x_prediction_k = filter([1, -a_opt'], 1, x_with_noise(1:k));
    
    % Calcular el error de predicción d(k)
    prediction_error_signal_k = x(1:k) - x_prediction_k;
    
    % Calcular el error cuadrático medio (MSE) para el instante de tiempo k
    mse_values(k) = mean(prediction_error_signal_k.^2);
end

% Encontrar el índice donde el MSE es mínimo
[min_mse, min_index] = min(mse_values);


% Mostrar resultados
% Mostrar el tiempo en el que se alcanza el mínimo MSE
% disp(['El mínimo se alcanza en el tiempo k = ', num2str(min_index)]);
disp('Error de predicción d(k):');
disp(['MSE = ', num2str(mse)]);
disp('Coeficientes predichos:');
disp(['a1_pred = ', num2str(a_opt(1)), ', a2_pred = ', num2str(a_opt(2))]);

