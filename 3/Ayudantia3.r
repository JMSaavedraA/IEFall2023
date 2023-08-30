## @knitr estimadores
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

estMu <- function(x) {
    mean(x)
}

estSigma <- function(x) {
    n <- length(x)
    sqrt(sum((x - mean(x))^2) / n)
}

## @knitr ej2
set.seed(25082023)
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
print(sprintf("%3i | %1.3f, %1.3f, %1.3f",n1,e11,e12,e13))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f",n2,e21,e22,e23))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f",n3,e31,e32,e33))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f",n4,e41,e42,e43))
print(sprintf("%3i | %1.3f, %1.3f, %1.3f",n5,e51,e52,e53))


## @knitr ej3dfn
n1 <- 15
x1 <- rnorm(n1, mean = 60, sd = 5)
n2 <- 30
x2 <- rnorm(n2, mean = 60, sd = 5)
n3 <- 100
x3 <- rnorm(n3, mean = 60, sd = 5)
x <- seq(40, 80, length = 500)
k1 <- 1:n1 / (n1)
k2 <- 1:n2 / (n2)
k3 <- 1:n3 / (n3)

## @knitr ej3a
mu1 <- estMu(x1)
sigma1 <- estSigma(x1)

## @knitr ej3b
plot(x,pnorm(x,mean = 60, sd = 5),type="l",col="red",xlab = "x", ylab = "F(x)")
lines(x,pnorm(x,mean = mu1, sd = sigma1),lty = 2, col = "blue")
lines(c(0,sort(x1),100),c(0,k1,1),col = "green",type = "s")
legend("bottomright",legend = c("Densidad teórica", "Densidad estimada", "Densidad Empírica"),col = c("red","blue", "green"),lty = c(1,2,1))

## @knitr ej3c
mu2 <- estMu(x2)
sigma2 <- estSigma(x2)
mu3 <- estMu(x3)
sigma3 <- estSigma(x3)

## @knitr ej3c1
plot(x,pnorm(x,mean = 60, sd = 5),type="l",col="red",xlab = "x", ylab = "F(x)")
lines(x,pnorm(x,mean = mu2, sd = sigma2),lty = 2, col = "blue")
lines(c(0,sort(x2),100),c(0,k2,1),col = "green",type = "s")
legend("bottomright",legend = c("Densidad teórica", "Densidad estimada", "Densidad Empírica"),col = c("red","blue", "green"),lty = c(1,2,1))

## @knitr ej3c2
plot(x,pnorm(x,mean = 60, sd = 5),type="l",col="red",xlab = "x", ylab = "F(x)")
lines(x,pnorm(x,mean = mu3, sd = sigma3),lty = 2, col = "blue")
lines(c(0,sort(x3),100),c(0,k3,1),col = "green",type = "s")
legend("bottomright",legend = c("Densidad teórica", "Densidad estimada", "Densidad Empírica"),col = c("red","blue", "green"),lty = c(1,2,1))