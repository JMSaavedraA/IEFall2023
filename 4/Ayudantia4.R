## @knitr previa
set.seed(300823)

## @knitr setBootstrap
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)

## @knitr estadisticas
est1 <- function(x) {
    # First Statistic, lambda estimator by the sample mean
    1/mean(x)
}

est2 <- function(x) {
    # Second Statistic, lambda estimator by the sample standard deviation
    n <- length(x)
    mX <- mean(x)
    sqrt((n-1) / sum((x - mX)^2))
}


## @knitr bootstrap

k <- 100
m <- 30

lambdaEst1 <- numeric(k)
lambdaEst2 <- numeric(k)

for(i in 1:k){
    z <- sample(x, size = m, replace = FALSE)
    lambdaEst1[i] <- est1(z)
    lambdaEst2[i] <- est2(z)
}

## @knitr comparativo

S1 <- sd(lambdaEst1)
S2 <- sd(lambdaEst2)

print(sprintf("La desviación estándar del estimador 1 es %2.3f",S1))
print(sprintf("La desviación estándar del estimador 2 es %2.3f",S2))