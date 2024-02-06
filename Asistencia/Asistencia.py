# Jose Miguel Saavedra Aguilar
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon, norm

random.seed(157017)

# Creamos un vector de 100 puntos en el soporte.
# En nuestro caso, tomamos de 0 a 10.
x = np.linspace(0, 10, num=100)

# y es la función de densidad elegida evaluada en los puntos de x
# Para obtener la función de densidad, se toma "pdf" + nombre de la distribución,
# por ejemplo expon.pdf
y1 = expon.pdf(x, scale=1/5)
y2 = expon.pdf(x, scale=1/0.5)
y3 = expon.pdf(x, scale=1/2)

# f es la función de distribución asociada evaluada en los puntos de x
# Para obtener la función de distribución, se toma "cdf" + nombre de la distribución,
# por ejemplo expon.cdf.
f1 = expon.cdf(x, scale=1/5)
f2 = expon.cdf(x, scale=1/0.5)
f3 = expon.cdf(x, scale=1/2)

# Tomamos n=100 y simulamos n datos pseudoaleatorios de la distribución elegida
# Para obtener los datos aleatorios, se toma "rvs" + nombre de la distribución,
# por ejemplo expon.rvs
n = 100
z = expon.rvs(scale=1/0.5, size=n)

# Ordenamos los datos como indica la tarea
z2 = np.sort(z)

# Inicializamos k = i/n+1
k = np.zeros(n)
for i in range(n):
    k[i] = (i+1) / (n+1)

# Gráficas
# Gráfico 1
plt.figure()
plt.plot(x, y1, 'red', label='rate = 5')
plt.plot(x, y2, 'blue', label='rate = 0.5')
plt.plot(x, y3, 'green', label='rate = 2')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.show()

# Gráfico 2
plt.figure()
plt.plot(x, f1, 'red', label='rate = 5')
plt.plot(x, f2, 'blue', label='rate = 0.5')
plt.plot(x, f3, 'green', label='rate = 2')
plt.xlabel('x')
plt.ylabel('F(x)')
plt.legend()
plt.show()

# Gráfico 3
plt.figure()
plt.hist(z, bins='auto', density=True)
plt.plot(x, y2)
plt.xlabel('x')
plt.ylabel('Frequency')
plt.title('Histogram of z')
plt.show()

# Gráfico 4
plt.figure()
plt.plot(z2, k, 'red')
plt.plot(x, f2, 'blue')
plt.xlabel('z')
plt.ylabel('k')
plt.title('Empirical Quantile Function')
plt.tight_layout()
plt.show()


# Ejemplo 2

prob = [0.1, 0.3, 0.5, 0.7, 1]
frec = np.histogram(np.random.uniform(size=1000), bins=np.cumsum(prob))[0] / 1000

np.random.seed(157017)
frecSample = np.histogram(np.random.choice(np.arange(5), size=1000, replace=True, p=[0.1, 0.2, 0.2, 0.2, 0.3]))[0] / 1000

print(f"Ejemplo 2:\n{frec}")
print(f"Ejemplo 2 (using sample):\n{frecSample}")

# Ejemplo 3
def p(n):
    values = np.random.choice(np.arange(1, 101), size=n, replace=True)
    prob = (2 * values) / (100 * 101)
    return prob

print("\nEjemplo 3:")
print(p(10000)[:30])

# Ejemplo 4
mayor = sum(np.sum(np.random.binomial(1, 0.09245, 1000)) for _ in range(10000))
montos = sum(np.sum(np.random.gamma(7000, 1, reclamaciones)) for reclamaciones in np.random.binomial(1000, 0.09245, 10000))
result = mayor / 10000
print(f"\nEjemplo 4:\n{result}")

# Tarea 3
# Ejercicio 1. Estimadores
def est1(x):
    return np.mean(x)

def est2(x):
    n = len(x)
    return np.sum((x - np.mean(x))**2) / n

def est3(x):
    n = len(x)
    return -0.5 + np.sqrt(0.25 + np.sum(x**2) / n)

# Ejercicio 1. Simulación
lambda_val = 3
n_values = [10, 20, 40, 80, 200]

for n in n_values:
    x = np.random.poisson(lambda_val, n)
    e1 = est1(x)
    e2 = est2(x)
    e3 = est3(x)
    print(f"n = {n:3} | est1: {e1:.3f}, est2: {e2:.3f}, est3: {e3:.3f}")

# Ejercicio 2. Definiciones
n1 = 15
x1 = np.random.normal(loc=60, scale=5, size=n1)
k1 = np.arange(1, n1 + 1) / n1

# Ejercicio 2. Estimadores
mu1 = est1(x1)
sigma1 = est2(x1)

# Ejercicio 2. Gráficas
x = np.linspace(40, 80, 500)
plt.figure()
plt.plot(x, norm.cdf(x, loc=60, scale=5), label="Densidad teórica", color="red")
plt.plot(x, norm.cdf(x, loc=mu1, scale=sigma1), label="Densidad estimada", linestyle="--", color="blue")
plt.step(np.concatenate(([0], np.sort(x1), [100])), np.concatenate(([0], k1, [1])), label="Densidad Empírica", color="green", linestyle="-.")
plt.xlabel("x")
plt.ylabel("F(x)")
plt.legend()
plt.show()

# Ejercicio 2b. Definiciones
n2 = 30
x2 = np.random.normal(loc=60, scale=5, size=n2)
n3 = 100
x3 = np.random.normal(loc=60, scale=5, size=n3)
k2 = np.arange(1, n2 + 1) / n2
k3 = np.arange(1, n3 + 1) / n3

# Ejercicio 2b. Estimadores
mu2, sigma2 = est1(x2), est2(x2)
mu3, sigma3 = est1(x3), est2(x3)

# Ejercicio 2c. Gráfica 1
plt.figure()
plt.plot(x, norm.cdf(x, loc=60, scale=5), label="Densidad teórica", color="red")
plt.plot(x, norm.cdf(x, loc=mu2, scale=sigma2), label="Densidad estimada", linestyle="--", color="blue")
plt.step(np.concatenate(([0], np.sort(x2), [100])), np.concatenate(([0], k2, [1])), label="Densidad Empírica", color="green", linestyle="-.")
plt.xlabel("x")
plt.ylabel("F(x)")
plt.legend()
plt.show()

# Ejercicio 2c. Gráfica 2
plt.figure()
plt.plot(x, norm.cdf(x, loc=60, scale=5), label="Densidad teórica", color="red")
plt.plot(x, norm.cdf(x, loc=mu3, scale=sigma3), label="Densidad estimada", linestyle="--", color="blue")
plt.step(np.concatenate(([0], np.sort(x3), [100])), np.concatenate(([0], k3, [1])), label="Densidad Empírica", color="green", linestyle="-.")
plt.xlabel("x")
plt.ylabel("F(x)")
plt.legend()
plt.show()


# Bootstrap no paramétrico. Parámetros de la distribución real
n = 300
lambda_val = 5
x = np.random.exponential(scale=1/lambda_val, size=n)

# Bootstrap no paramétrico. Parámetros de bootstrap
r = np.mean(x)
k = 100
m = 30
q1 = 0.025
q2 = 0.975

q1_est = np.zeros(k)
q2_est = np.zeros(k)

# Bootstrap no paramétrico. Bootstrap

for i in range(k):
    z = np.random.choice(x, size=m, replace=False)
    q1_est[i] = np.percentile(z, q=q1)
    q2_est[i] = np.percentile(z, q=q2)

# Bootstrap no paramétrico. Resultados
Q1 = np.mean(q1_est)
Q2 = np.mean(q2_est)

S1 = np.std(q1_est, ddof=1)  # ddof=1 for unbiased standard deviation
S2 = np.std(q2_est, ddof=1)

print(f"La media del cuantil {q1:.3f} es {Q1:.3f}")
print(f"La media del cuantil {q2:.3f} es {Q2:.3f}")
print(f"La desviación estándar del cuantil {q1:.3f} es {S1:.3f}")
print(f"La desviación estándar del cuantil {q2:.3f} es {S2:.3f}")

# Bootstrap paramétrico. Parámetros de la distribución real

n = 300
lambda_val = 5
x = np.random.exponential(scale=1/lambda_val, size=n)

# Bootstrap paramétrico. Parámetros de bootstrap
r = np.mean(x)
k = 100
m = 30
q1 = 0.025
q2 = 0.975

q1_est = np.zeros(k)
q2_est = np.zeros(k)

# Bootstrap paramétrico. Bootstrap
for i in range(k):
    z = np.random.exponential(scale=r, size=m)
    r_boot = np.mean(z)
    q1_est[i] = expon.ppf(q1, scale=r_boot)
    q2_est[i] = expon.ppf(q2, scale=r_boot)

# Bootstrap paramétrico. Resultados

Q1 = np.mean(q1_est)
Q2 = np.mean(q2_est)

S1 = np.std(q1_est, ddof=1)  # ddof=1 for unbiased standard deviation
S2 = np.std(q2_est, ddof=1)

print(f"La media del cuantil {q1:.3f} es {Q1:.3f}")
print(f"La media del cuantil {q2:.3f} es {Q2:.3f}")
print(f"La desviación estándar del cuantil {q1:.3f} es {S1:.3f}")
print(f"La desviación estándar del cuantil {q2:.3f} es {S2:.3f}")


# Bootstrap paramétrico. Parámetros de la distribución real
n = 300
lambda_val = 5
x = np.random.exponential(scale=1/lambda_val, size=n)

# Bootstrap paramétrico. Parámetros de bootstrap
r = np.mean(x)
k = 1000
m = 100
q1 = 0.025
q2 = 0.975

q1_est = np.zeros(k)
q2_est = np.zeros(k)

# Bootstrap paramétrico. Bootstrap

for i in range(1, k+1):
    z = np.random.exponential(scale=r, size=m)
    r_boot = np.mean(z)
    q1_est[i-1] = expon.ppf(q1, scale=r_boot)
    q2_est[i-1] = expon.ppf(q2, scale=r_boot)

# Bootstrap paramétrico. Resultados incluyendo intervalos de confianza
Q1 = np.mean(q1_est)
Q2 = np.mean(q2_est)

S1 = np.std(q1_est, ddof=1)  # ddof=1 for unbiased standard deviation
S2 = np.std(q2_est, ddof=1)

T1L = Q1 - S1 * norm.ppf(0.975)
T1R = Q1 + S1 * norm.ppf(0.975)

T2L = Q2 - S2 * norm.ppf(0.975)
T2R = Q2 + S2 * norm.ppf(0.975)

q1_real = expon.ppf(q1, scale=1/lambda_val)
q2_real = expon.ppf(q2, scale=1/lambda_val)

print("Por el método normal:")
print(f"El intervalo de confianza para el cuantil {q1:.3f} es ({T1L:.5f}, {T1R:.5f})")
if T1L < q1_real < T1R:
    print(f"El cuantil {q1:.3f} está en el intervalo de confianza")
else:
    print(f"El cuantil {q1:.3f} no está en el intervalo de confianza")

print(f"El intervalo de confianza para el cuantil {q2:.3f} es ({T2L:.5f}, {T2R:.5f})")
if T2L < q2_real < T2R:
    print(f"El cuantil {q2:.3f} está en el intervalo de confianza")
else:
    print(f"El cuantil {q2:.3f} no está en el intervalo de confianza")

q1_N = expon.ppf(q1, scale=r)
q2_N = expon.ppf(q2, scale=r)

q1_R = np.percentile(q1_est, q=2.5)
q1_L = np.percentile(q1_est, q=97.5)
q2_R = np.percentile(q2_est, q=2.5)
q2_L = np.percentile(q2_est, q=97.5)

P1L = 2 * q1_N - q1_L
P1R = 2 * q1_N - q1_R
P2L = 2 * q2_N - q2_L
P2R = 2 * q2_N - q2_R

print("Por el método pivotal:")
print(f"El intervalo de confianza para el cuantil {q1:.3f} es ({P1L:.5f}, {P1R:.5f})")
if P1L < q1_real < P1R:
    print(f"El cuantil {q1:.3f} está en el intervalo de confianza")
else:
    print(f"El cuantil {q1:.3f} no está en el intervalo de confianza")

print(f"El intervalo de confianza para el cuantil {q2:.3f} es ({P2L:.5f}, {P2R:.5f})")
if P2L < q2_real < P2R:
    print(f"El cuantil {q2:.3f} está en el intervalo de confianza")
else:
    print(f"El cuantil {q2:.3f} no está en el intervalo de confianza")