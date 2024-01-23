## @knitr inic
#Jose Miguel Saavedra Aguilar
library(ordinal)
set.seed(157017)

## @knitr t1inicializacion
# Creamos un vector de 100 puntos en el soporte.
# En nuestro caso, tomamos de 0 a 10.
x <- seq(0, 10, length = 100)

## @knitr t1ej1
# y es la función de densidad elegida evaluada en los puntos de x
# Para obtener la función de densidad, se toma "d"+nombre de la distribución,
# por ejemplo dnorm o dexp
y1 <- dexp(x, rate = 5)
y2 <- dexp(x, rate = 0.5)
y3 <- dexp(x, rate = 2)


## @knitr t1ej2
# f es la función de distribución asociada evaluada en los puntos de x
# Para obtener la función de densidad, se toma "p"+nombre de la distribución,
# por ejemplo pnorm o pexp.
f1 <- pexp(x, rate = 5)
f2 <- pexp(x, rate = 0.5)
f3 <- pexp(x, rate = 2)


## @knitr t1ej3
# Tomamos n=100 y simulamos n datos pseudoaleatorios de la distribución elegida
# Para obtener los datos aleatorios, se toma "r"+nombre de la distribución,
# por ejemplo rnorm o rexp
n <- 100
z <- rexp(n, rate = 0.5)

## @knitr t1ej4
#Ahora vamos a ordenar los datos como indica la tarea.
z2 <- sort(z)
# Inicializamos k =i/n+1
k <- numeric()

for (i in 1:n) {
  k[i] <- i / (n + 1)
}

## @knitr t1gr1
plot(x, y1, type = "l", col = "red", xlab = "x", ylab = "f(x)")
lines(x, y2, col = "blue")
lines(x, y3, col = "green")

## @knitr t1gr2
plot(x, f1, type = "l", col = "red", xlab = "x", ylab = "f(x)")
lines(x, f2, col = "blue")
lines(x, f3, col = "green")

## @knitr t1gr3
hist(z, breaks = "FD", freq = FALSE)
lines(x, y2)

## @knitr t1gr4
#Graficamos los puntos y comparamos con la función de densidad teórica
plot(z2, k, col = "red", type = "l")
lines(x, f2, col = "blue")

## @knitr t2ej1
#Ejemplo 2
set.seed(157017)
prob <- c(.1, 0.3, 0.5, 0.7, 1)
frec <- findInterval(runif(1000), prob)
table(frec) / 1000
set.seed(157017)
frecSample <- sample(x = 0:4, size = 1000, replace = TRUE,
                     prob = c(0.1, 0.2, 0.2, 0.2, 0.3))
table(frecSample) / 1000

## @knitr t2ej2
#Ejemplo 3
p <- function(n) {
  values <- sample(1:100, n, replace = TRUE)
  prob <- (2 * values) / (100 * 101)
  return(prob)
}
head(p(10000), n = 30)

## @knitr t2ej3
#Ejemplo 4
mayor <- 0
for (i in 1:10000){
  reclamaciones <- sum(rbinom(1000, 1, 0.09245))
  montos <- sum(rgamma(reclamaciones, 7000, 1))
  if (montos > 500000) {
    mayor <- mayor + 1
  }
}
mayor / 10000


## @knitr t3estimadores
est1 <- function(x) {
  mean(x)
}

est2 <- function(x) {
  n <- length(x)
  sum((x - mean(x))^2) / n
}

est3 <- function(x) {
  n <- length(x)
  -0.5 + sqrt(.25 + sum(x^2) / n)
}

## @knitr t3ej2
lambda <- 3
n1 <- 10
n2 <- 20
n3 <- 40
n4 <- 80
n5 <- 200

x1 <- rpois(n1, lambda)
x2 <- rpois(n2, lambda)
x3 <- rpois(n3, lambda)
x4 <- rpois(n4, lambda)
x5 <- rpois(n5, lambda)

e11 <- est1(x1)
e12 <- est2(x1)
e13 <- est3(x1)
e21 <- est1(x2)
e22 <- est2(x2)
e23 <- est3(x2)
e31 <- est1(x3)
e32 <- est2(x3)
e33 <- est3(x3)
e41 <- est1(x4)
e42 <- est2(x4)
e43 <- est3(x4)
e51 <- est1(x5)
e52 <- est2(x5)
e53 <- est3(x5)

print("  n |  est1,  est2,  est3")
print(sprintf("%3i | %1.3f, %1.3f, %1.3f", n1, e11, e12, e13))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f", n2, e21, e22, e23))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f", n3, e31, e32, e33))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f", n4, e41, e42, e43))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f", n5, e51, e52, e53))


## @knitr t3ej3dfn
n1 <- 15
x1 <- rnorm(n1, mean = 60, sd = 5)
k1 <- 1:n1 / (n1)


## @knitr t3ej3a

estMu <- function(x) {
  mean(x)
}

estSigma <- function(x) {
  n <- length(x)
  sqrt(sum((x - mean(x))^2) / n)
}

mu1 <- estMu(x1)
sigma1 <- estSigma(x1)
k1 <- 1 : n1 / (n1)

## @knitr t3ej3b
x <- seq(40, 80, length = 500)
plot(x, pnorm(x, mean = 60, sd = 5), type = "l", col = "red",
     xlab = "x", ylab = "F(x)")
lines(x, pnorm(x, mean = mu1, sd = sigma1), lty = 2, col = "blue")
lines(c(0, sort(x1), 100), c(0,k1,1), col = "green", type = "s")
legend("bottomright", legend = c("Densidad teórica", "Densidad estimada",
                                 "Densidad Empírica"),
       col = c("red", "blue", "green"), lty = c(1, 2, 1))

## @knitr t3ej3dfn2
n2 <- 30
x2 <- rnorm(n2, mean = 60, sd = 5)
n3 <- 100
x3 <- rnorm(n3, mean = 60, sd = 5)
k2 <- 1:n2 / (n2)
k3 <- 1:n3 / (n3)

## @knitr t3ej3c
mu2 <- estMu(x2)
sigma2 <- estSigma(x2)
mu3 <- estMu(x3)
sigma3 <- estSigma(x3)

## @knitr t3ej3c1
plot(x, pnorm(x, mean = 60, sd = 5), type = "l", col = "red",
     xlab = "x", ylab = "F(x)")
lines(x, pnorm(x, mean = mu2, sd = sigma2), lty = 2, col = "blue")
lines(c(0, sort(x2), 100), c(0, k2, 1), col = "green", type = "s")
legend("bottomright", legend = c("Densidad teórica", "Densidad estimada",
                                 "Densidad Empírica"),
       col = c("red", "blue", "green"), lty = c(1, 2, 1))

## @knitr t3ej3c2
plot(x, pnorm(x, mean = 60, sd = 5), type = "l", col = "red", xlab = "x",
     ylab = "F(x)")
lines(x, pnorm(x, mean = mu3, sd = sigma3), lty = 2, col = "blue")
lines(c(0, sort(x3), 100), c(0, k3, 1), col = "green", type = "s")
legend("bottomright",
       legend = c("Densidad teórica", "Densidad estimada", "Densidad Empírica"),
       col = c("red", "blue", "green"), lty = c(1, 2, 1))


## @knitr t4setBootstrap
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)

## @knitr t4estadisticas
est1 <- function(x) {
  # First Statistic, lambda estimator by the sample mean
  1 / mean(x)
}

est2 <- function(x) {
  # Second Statistic, lambda estimator by the sample standard deviation
  n <- length(x)
  mX <- mean(x)
  sqrt((n - 1) / sum((x - mX)^2))
}


## @knitr t4bootstrap

k <- 100
m <- 30

lambdaEst1 <- numeric(k)
lambdaEst2 <- numeric(k)

for(i in 1:k){
    z <- sample(x, size = m, replace = FALSE)
    lambdaEst1[i] <- est1(z)
    lambdaEst2[i] <- est2(z)
}

## @knitr t4comparativo

S1 <- sd(lambdaEst1)
S2 <- sd(lambdaEst2)

print(sprintf("La desviación estándar del estimador 1 es %2.3f", S1))
print(sprintf("La desviación estándar del estimador 2 es %2.3f", S2))


## @knitr t5setBootstrap
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)

## @knitr t5estadisticas
est1 <- function(x) {
  # First Statistic, lambda estimator by the sample mean
  1 / mean(x)
}

est2 <- function(x) {
  # Second Statistic, lambda estimator by the sample standard deviation
  n <- length(x)
  mX <- mean(x)
  sqrt((n - 1) / sum((x - mX)^2))
}


## @knitr t5bootstrap
r <- mean(x)
k <- 100
m <- 30

lambdaEst1 <- numeric(k)
lambdaEst2 <- numeric(k)

for (i in 1:k){
  z <- rexp(m, rate = 1 / r)
  lambdaEst1[i] <- est1(z)
  lambdaEst2[i] <- est2(z)
}

## @knitr t5comparativo

S1 <- sd(lambdaEst1)
S2 <- sd(lambdaEst2)

print(sprintf("La desviación estándar del estimador 1 es %2.3f", S1))
print(sprintf("La desviación estándar del estimador 2 es %2.3f", S2))


## @knitr t6setBootstrap
n <- 300
mu <- .2
sigma <- .2
x <- rnorm(300, mean = mu, sd = sigma)

## @knitr t6estadisticas
est1 <- function(x) {
  # First Statistic, lambda estimator by the sample mean
  quantile(x, 0.975)
}


muEst <- mean(x)
sigmaEst <- sd(x)

## @knitr t6bootstrap

k <- 100
m <- 30

meanEst1 <- numeric(k)

for (i in 1:k){
  #z <- sample(x, size = m, replace = FALSE)
  z <- rnorm(m, mean = muEst, sd = sigmaEst)
  meanEst1[i] <- est1(z)
}

## @knitr t6comparativo

S1 <- quantile(meanEst1, 0.02)
S2 <- quantile(meanEst1, 0.98)

T1 <- mean(meanEst1) + sd(meanEst1) / sqrt(k) * qt(0.02, k - 1)
T2 <- mean(meanEst1) + sd(meanEst1) / sqrt(k) * qt(0.98, k - 1)

print(sprintf("El cuantil empírico 0.02 es %2.3f", S1))
print(sprintf("El cuantil empírico 0.98 es %2.3f", S2))
print(sprintf("El cuantil T 0.02 es %2.3f", T1))
print(sprintf("El cuantil T 0.98 es %2.3f", T2))