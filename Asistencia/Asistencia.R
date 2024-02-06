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


## @knitr t4bootstrap
r <- mean(x)
k <- 100
m <- 30
q1 <- 0.025
q2 <- 0.975

q1Est <- numeric(k)
q2Est <- numeric(k)

for (i in 1:k){
  z <- sample(x,size = m, replace = FALSE)
  q1Est[i] <- quantile(z,q1)
  q2Est[i] <- quantile(z,q2)
}

## @knitr t4comparativo

Q1 <- mean(q1Est)
Q2 <- mean(q2Est)

S1 <- sd(q1Est)
S2 <- sd(q2Est)

print(sprintf("La media del cuantil %1.3f es %2.3f", q1, Q1))
print(sprintf("La media del cuantil %1.3f es %2.3f", q2, Q2))
print(sprintf("La desviación estándar del cuantil %1.3f es %2.3f", q1, S1))
print(sprintf("La desviación estándar del cuantil %1.3f es %2.3f", q2, S2))

## @knitr t5setBootstrap
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)


## @knitr t5bootstrap
r <- mean(x)
k <- 100
m <- 30
q1 <- 0.025
q2 <- 0.975

q1Est <- numeric(k)
q2Est <- numeric(k)

for (i in 1:k){
  z <- rexp(m, rate = 1 / r)
  rBoot <- mean(z)
  q1Est[i] <- qexp(q1, rate = 1 / rBoot)
  q2Est[i] <- qexp(q2, rate = 1 / rBoot)
}

## @knitr t5comparativo

Q1 <- mean(q1Est)
Q2 <- mean(q2Est)

S1 <- sd(q1Est)
S2 <- sd(q2Est)

print(sprintf("La media del cuantil %1.3f es %2.3f", q1, Q1))
print(sprintf("La media del cuantil %1.3f es %2.3f", q2, Q2))
print(sprintf("La desviación estándar del cuantil %1.3f es %2.3f", q1, S1))
print(sprintf("La desviación estándar del cuantil %1.3f es %2.3f", q2, S2))

## @knitr t6setBootstrap
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)


## @knitr t6bootstrap
r <- mean(x)
k <- 100
m <- 30
q1 <- 0.025
q2 <- 0.975

q1Est <- numeric(k)
q2Est <- numeric(k)

for (i in 1:k){
  z <- rexp(m, rate = 1 / r)
  rBoot <- mean(z)
  q1Est[i] <- qexp(q1, rate = 1 / rBoot)
  q2Est[i] <- qexp(q2, rate = 1 / rBoot)
}

## @knitr t6normal

Q1 <- mean(q1Est)
Q2 <- mean(q2Est)

S1 <- sd(q1Est)  # R uses ddof=1 by default for unbiased standard deviation
S2 <- sd(q2Est)

T1L <- Q1 - S1 * qnorm(0.975)
T1R <- Q1 + S1 * qnorm(0.975)

T2L <- Q2 - S2 * qnorm(0.975)
T2R <- Q2 + S2 * qnorm(0.975)

q1_real <- qexp(q1, rate = lambda)
q2_real <- qexp(q2, rate = lambda)

print("Por el método normal:")
print(sprintf("El intervalo de confianza para el cuantil %.3f es (%.5f, %.5f)", q1, T1L, T1R))
if (T1L < q1_real && q1_real < T1R) {
  print(sprintf("El cuantil %.3f está en el intervalo de confianza", q1))
} else {
  print(sprintf("El cuantil %.3f no está en el intervalo de confianza", q1))
}

print(sprintf("El intervalo de confianza para el cuantil %.3f es (%.5f, %.5f)", q2, T2L, T2R))
if (T2L < q2_real && q2_real < T2R) {
  print(sprintf("El cuantil %.3f está en el intervalo de confianza", q2))
} else {
  print(sprintf("El cuantil %.3f no está en el intervalo de confianza", q2))
}

## @knitr t6pivotal

q1_N <- qexp(q1, rate = 1 / r)
q2_N <- qexp(q2, rate = 1 / r)

q1_R <- quantile(q1Est, 0.025)
q1_L <- quantile(q1Est, 0.975)
q2_R <- quantile(q2Est, 0.025)
q2_L <- quantile(q2Est, 0.975)

P1L <- 2 * q1_N - q1_L
P1R <- 2 * q1_N - q1_R
P2L <- 2 * q2_N - q2_L
P2R <- 2 * q2_N - q2_R

print("Por el método pivotal:")
print(sprintf("El intervalo de confianza para el cuantil %.3f es (%.5f, %.5f)", q1, P1L, P1R))
if (P1L < q1_real && q1_real < P1R) {
  print(sprintf("El cuantil %.3f está en el intervalo de confianza", q1))
} else {
  print(sprintf("El cuantil %.3f no está en el intervalo de confianza", q1))
}

print(sprintf("El intervalo de confianza para el cuantil %.3f es (%.5f, %.5f)", q2, P2L, P2R))
if (P2L < q2_real && q2_real < P2R) {
  print(sprintf("El cuantil %.3f está en el intervalo de confianza", q2))
} else {
  print(sprintf("El cuantil %.3f no está en el intervalo de confianza", q2))
}
