clear all;
close all;
clc;

% Parámetros del proceso estacionario
sigma_x_squared = 1; % Varianza del ruido blanco
syms a1 a2 lambda;
% Definir la función objetivo J(a1, a2) = a1^2 + a2^2
J = a1^2 + a2^2;
% Definir la función de restricción g(a1, a2) = a1 + a2 - 1
g = a1 + a2 - 1;
% Definir la función Lagrangiana L(a1, a2, lambda)
L = J - lambda * g;
% Tomar las derivadas parciales de L con respecto a a1, a2 y lambda
dL_da1 = diff(L, a1);
dL_da2 = diff(L, a2);
dL_dlambda = diff(L, lambda);
% Resolver el sistema de ecuaciones resultante
[a1_sol, a2_sol, lambda_sol] = solve(dL_da1 == 0, dL_da2 == 0, dL_dlambda == 0, a1, a2, lambda);
% Mostrar los resultados
disp(['a1 = ', char(a1_sol)]);
disp(['a2 = ', char(a2_sol)]);
% Coeficientes óptimos del predictor
a1_opt = double(a1_sol);
a2_opt = double(a2_sol);

% Tiempo de simulación
N = 100;

% Generación de la señal x(k) como ruido blanco
x = sqrt(sigma_x_squared) * randn(1, N);

% Inicialización del error de predicción
d = zeros(1, N);

% Cálculo del error de predicción
for k = 3:N
    d(k) = x(k) - a1_opt*x(k-1) - a2_opt*x(k-2);
end

% Calcular el error de predicción
error_prediccion = mean(d(3:end).^2);

% Encontrar el índice donde el error de predicción alcanza su mínimo
[min_error, min_index] = min(d);

% Frecuencia de muestreo (ejemplo)
%fs = 1000; % 1000 muestras por segundo (ejemplo)


% Mostrar resultados
%fprintf('El valor mínimo del error de predicción es %f en k %d y en el tiempo %f segundos.\n', min_error, min_index, tiempo_minimo);
fprintf('El error de predicción es %f.\n', error_prediccion);
fprintf('Los coeficientes óptimos del predictor son: a1 = %f, a2 = %f\n', a1_opt, a2_opt);
% % Calcular el tiempo correspondiente al índice del valor mínimo
% tiempo_minimo = min_index / fs; % en segundos

% % Graficar el error de predicción
% figure;
% stem(d);
% hold on;
% stem(min_index, min_error, 'r', 'LineWidth', 2);
% hold off;
% title('Error de predicción');
% xlabel('k');
% ylabel('d(k)');
% legend('Error de predicción', 'Valor mínimo');
