% supressão de ruído utilizando ondaletas

load noisdopp

dopp_den = wdenoise(noisdopp,DenoisingMethod='UniversalThreshold',Wavelet='db6');
dopp2_den = wdenoise(noisdopp,DenoisingMethod='SURE');
dopp3_den = wdenoise(noisdopp,DenoisingMethod='Minimax');
dopp4_den = wdenoise(noisdopp,DenoisingMethod='Bayes')

subplot(2,2,1); plot(dopp_den); title('Limiar Universal')
subplot(2,2,2); plot(dopp2_den); title('SureShrink')
subplot(2,2,3); plot(dopp3_den); title('Minimax')
subplot(2,2,4); plot(dopp4_den); title('BayesShrink')

% Parâmetros
N = 1024;
explosao_prob = 0.01;
ruido_std = 0.3;
bump_width = 5;

sinal = randn(1, N);

for i = 1:N
    if rand < explosao_prob
        start_idx = max(1, i - bump_width);
        end_idx = min(N, i + bump_width);
        sinal(start_idx:end_idx) = sinal(start_idx:end_idx) + 10;
    end
end

% aqui é adicionado o ruído gaussiano
ruido = ruido_std * randn(1, N);
sinal_exemplo_1 = sinal + ruido;

dopp_den = wdenoise(sinal_exemplo_1,DenoisingMethod='UniversalThreshold',Wavelet='db6');
dopp2_den = wdenoise(sinal_exemplo_1,DenoisingMethod='SURE');
dopp3_den = wdenoise(sinal_exemplo_1,DenoisingMethod='Minimax');
dopp4_den = wdenoise(sinal_exemplo_1,DenoisingMethod='Bayes');

subplot(3,2,1); plot(sinal_exemplo_1); title('Original')
subplot(3,2,2); plot(dopp_den); title('Limiar Universal')
subplot(3,2,3); plot(dopp2_den); title('SureShrink')
subplot(3,2,4); plot(dopp3_den); title('Minimax');
subplot(3,2,5); plot(dopp4_den); title('BayesShrink');

N = 1024;
t = linspace(0, 10, N);
trend = sin(t) + 0.5 * cos(0.5 * t); % Tendências fortes de subida e queda
explosao_prob = 0.05;
ruido_std = 0.2;
bump_width = 5;

sinal = trend + randn(1, N);

for i = 1:N
    if rand < explosao_prob
        start_idx = max(1, i - bump_width);
        end_idx = min(N, i + bump_width);
        sinal(start_idx:end_idx) = sinal(start_idx:end_idx) + 10;
    end
end

ruido = ruido_std * randn(1, N);
sinal_exemplo_1 = sinal + ruido;

plot(sinal_exemplo_1);
title('Sinal com Tendências Acentuadas, Bumps Quadrados e Ruído Gaussiano');
xlabel('Amostras');
ylabel('Amplitude');

