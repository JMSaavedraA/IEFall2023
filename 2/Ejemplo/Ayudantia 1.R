## @knitr inicializacion
# Creamos un vector de 100 puntos en el soporte.
# En nuestro caso, tomamos de 0 a 10.
x <- seq(0, 10, length = 100)

## @knitr ej1
# y es la función de densidad elegida evaluada en los puntos de x
# Para obtener la función de densidad, se toma "d"+nombre de la distribución, por ejemplo dnorm o dexp
y1 <- dexp(x, rate = 5)
y2 <- dexp(x, rate = 0.5)
y3 <- dexp(x, rate = 2)


## @knitr ej2
# f es la función de distribución asociada evaluada en los puntos de x
# Para obtener la función de densidad, se toma "p"+nombre de la distribución, por ejemplo pnorm o pexp
f1 <- pexp(x, rate = 5)
f2 <- pexp(x, rate = 0.5)
f3 <- pexp(x, rate = 2)


## @knitr ej3
# Tomamos n=100 y simulamos n datos pseudoaleatorios de la distribución elegida
# Para obtener los datos aleatorios, se toma "r"+nombre de la distribución, por ejemplo rnorm o rexp
n <- 100
z <- rexp(n, rate = 0.5)

## @knitr ej4
#Ahora vamos a ordenar los datos como indica la tarea.
z2 <- sort(z)
# Inicializamos k =i/n+1
k <- numeric()

for (i in 1:n) {
  k[i] <- i / (n + 1)
}