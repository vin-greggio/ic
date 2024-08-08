#Aqui, é implementado o método de Ogden & Parzen, teste de hipóteses, para seleção de thresholds, utilizando a linguagem R.

library(wavethresh)

set.seed(2024)
N <- 1024
t <- seq(0, 10, length.out = N)
signal <- sin(t) + 0.76 * cos(2 * t)
noise <- rnorm(N, sd = 0.45)
noisy_signal <- signal + noise

wavelet_decomp <- wd(noisy_signal)

wavelet_decomp <- TOthreshda1(wavelet_decomp, alpha = 0.05)

denoised_signal <- wr(wavelet_decomp)
par(mfrow = c(3, 1))

plot(t, signal, type = "l", col = "blue", main = "Sinal Original", xlab = "Tempo", ylab = "Amplitude")
plot(t, noisy_signal, type = "l", col = "red", main = "Sinal com Ruído", xlab = "Tempo", ylab = "Amplitude")
plot(t, denoised_signal, type = "l", col = "green", main = "Sinal com Ruído retirado, usando alpha = 0.05", xlab = "Tempo", ylab = "Amplitude")
