clear
close all;
K = 10000; % interacciones
ENS = 1000;
N = 3;
varx = 1; % POTENCIA DE LA SEÑAL
varn = 1e-8; % POTENCIA DEL RUIDO 
num = 1;
den = [1 1.2 .81];
mu = 0.05;
mun = 0.5;
lambda = 0.989;
Wopt = [1 .9 .1 .2].'; % coefientes a identificar

MSELMS = zeros(K,1);
for ens = 1:ENS
    disp(ens);
    % Generating input and reference signals
    x = randn(K,1);
    x = filter(num, den, x);
    x = x / std(x);
    x = x * sqrt(varx);
    n = randn(K,1);
    n = filter(Wopt, 1, x) + n;
    n = n / std(n);
    n = n * sqrt(varn);
    d = filter(Wopt, 1, x) + n;
    
    % INICIALIZACIONES
    xk = zeros(N+1,1);  
    WLMS = zeros(N+1,1);
  
    for k = 1:K
        xk = [x(k)
              xk(1:N)];
        % ALGORITMO LMS
        ek1 = d(k) - xk.' * WLMS;
        WLMS = WLMS + mu * ek1 * xk;
        MSELMS(k) = MSELMS(k) + ek1^2;
    end % for ending k
end % for ending ens

MSELMS = MSELMS / ENS;
MSEmindB = 10*log10(varn) * ones(K,1);

% SIMULATION RESULTS
xscale = 1:K;
xlabel('k'); ylabel('dB');
title('MSE of conventional LMS');
plot(xscale, 10*log10(MSELMS), xscale, MSEmindB, '--r');
legend('LMS', 'MSE_{min}')

