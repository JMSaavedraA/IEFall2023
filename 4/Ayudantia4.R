## @knitr previa
set.seed(300823)
n <- 300
lambda <- 5
x <- rexp(300, rate = lambda)

est1 <- function(x) {
    # First Statistic, lambda estimator by the sample mean
    1/mean(x)
}

est2 <- function(x) {
    # Second Statistic, lambda estimator by the sample standard deviation
    n <- length(x)
    sqrt((n-1) / sum((x - mean(x))^2))
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

