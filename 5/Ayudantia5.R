## @knitr previa
set.seed(300823)

## @knitr setBootstrap
n <- 300
mu <- .2
sigma <- .2
x <- rnorm(300, mean = mu, sd = sigma)

## @knitr estadisticas
est1 <- function(x) {
  # First Statistic, lambda estimator by the sample mean
  quantile(x, 0.975)
}

est2 <- function(x) {
  # Second Statistic, lambda estimator by the sample standard deviation
  quantile(x, 0.025)
}

muEst <- mean(x)
sigmaEst <- sd(x)
## @knitr bootstrap

k <- 100
m <- 30

meanEst1 <- numeric(k)

for (i in 1:k){
  #z <- sample(x, size = m, replace = FALSE)
  z <- rnorm(m, mean = muEst, sd = sigmaEst)
  meanEst1[i] <- est1(z)
}

## @knitr comparativo

S1 <- quantile(meanEst1, 0.02)
S2 <- quantile(meanEst1, 0.98)

T1 <- mean(meanEst1) + sd(meanEst1) / sqrt(k) * qt(0.02, k - 1)
T2 <- mean(meanEst1) + sd(meanEst1) / sqrt(k) * qt(0.98, k - 1)

print(sprintf("El cuantil empírico 0.02 es %2.3f", S1))
print(sprintf("El cuantil empírico 0.98 es %2.3f", S2))
print(sprintf("El cuantil T 0.02 es %2.3f", T1))
print(sprintf("El cuantil T 0.98 es %2.3f", T2))
