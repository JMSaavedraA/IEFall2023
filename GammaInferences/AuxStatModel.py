import numpy as np
from scipy.special import digamma, gamma
from scipy.optimize import minimize_scalar
##  Funciones que se usarán:                                               ####
#   -------------------------------- --
# Estadística Anderson Darling para comparación de modelos. 
# Se le alimenta el vector ZZ de datos ordenados transformados con la
# distribución estimada: ZZ=Fhat(xsort).
def AndersonDarling(ZZ):
    nn = len(ZZ)
    ns = np.arange(1, nn + 1)
    aux1 = (2 * ns - 1)
    aux2 = (2 * nn + 1 - (2 * ns))
    AD = -nn - (1 / nn) * np.sum((aux1 * np.log(ZZ)) + aux2 * np.log(1 - ZZ))
    return AD
#--------------- --
# Estadística propuesta por Eloísa del promedio ponderado de distancias
# verticales en la gráfica de probabilidad.
def WAVD(ZZ):
    nn = len(ZZ)
    ns = np.arange(1, nn + 1)
    aens = nn + 1 - ns
    Ebetas = ns / (nn + 1)
    varian = (ns * aens) / (((nn + 1) ** 2) * (nn + 2))
    sali = np.sum(np.abs(ZZ - Ebetas) / np.sqrt(varian)) / nn
    return sali
#--------------- --
# Función para hallar la raíz de primera derivada de lp(alfa)=0 para el
# parámetro de forma alfa de la distribución Gamma.
def Lgampsi(ALPHA, T1, T3, N):
    if ALPHA > 0:
        sali = np.log(ALPHA * N / T1) - digamma(ALPHA) + (T3 / N)
    else:
        sali = 55555
    return sali
#------------------------- --
# Log verosimilitud perfil del parámetro forma Gamma:
def lpgamalfa(alfas, T1, T3, N):
    A = np.sum(alfas <= 0)
    if A > 1:
        print('Warning alpha must be positive')
        lp = -14
    else:
        lp = -N * gamma(alfas) + (N * alfas * np.log(alfas)) - N * alfas * np.log(T1 / N) + alfas * (T3 - N)
    return lp
#------------------------- --
# Para encontrar numericamente intervalos verosim perfil de nivel c 
# para parámetro alfa de la Gamma se calcula rp(alpha)-lnc 
def rpgamalfalnc(alfas, T1, T3, N, ALPHAMLE, CC):
    rsali = lpgamalfa(alfas, T1, T3, N) - lpgamalfa(ALPHAMLE, T1, T3, N) - np.log(CC)
    return rsali
#------------------------- --
# Logverosimilitud (alpha,mu) para la distribución Gamma 
def lvgamalfamu(vec2, T1, T3, N):
    ALFA = vec2[0]
    MU = vec2[1]
    if ALFA > 0 and MU > 0:
        lv = -N * gamma(ALFA) + N * ALFA * np.log(ALFA) + (ALFA * T3) - (ALFA * T1 / MU) - (N * ALFA * np.log(MU))
    else:
        lv = -50
    return lv
#------------------------- --
# Para calcular numéricamente la logperfil de la media Gamma:
def lpgamamu(vec1, T1, T3, N, MU):
    ALFA = vec1
    if ALFA > 0 and MU > 0:
        lv = -N * gamma(ALFA) + N * ALFA * np.log(ALFA) + (ALFA * T3) - (ALFA * T1 / MU) - (N * ALFA * np.log(MU))
    else:
        lv = -50
    return lv
#------------------------- --
# Para obtener intervalos perfiles de la media Gamma, se calcula
#  rp(mu)-lnc:
def rpgamamulnc(MU, T1, T3, N, ALPHAMLE, MUMLE, CC):
    result = minimize_scalar(lambda ALFA: lpgamamu(ALFA, T1, T3, N, MU), bounds=(0.001, ALPHAMLE * 15), method='bounded')
    lpgmu = result.fun
    sali = lpgmu - lpgamamu(ALPHAMLE, T1, T3, N, MUMLE) - np.log(CC)
    return sali
#----------------------------------
# Se calculan niveles de verosimilitud de intervalos y regiones para 
# inferencias con la distribución Gamma (uno y dos parámetros)
def cgamlevels(n):
    # Obtained jointly with Elías Mercado for his BA Thesis 
    #for a single Gamma parameter:
    cg90 = 0.2585 - (1 / (-0.5171 + (1.6667 * n)))
    cg95 = 0.1465 - (1 / (2.029 + (2.084 * n)))
    cg99 = 0.0362 - (1 / (14.609 + (4.886 * n)))
    cgamint = {"90%": cg90, "95%": cg95, "99%": cg99}
    #for both Gamma parameters:
    cgr90 = 0.1 - (1 / (1.060 + (4.617 * n)))
    cgr95 = 0.05 - (1 / (4.932 + (7.030 * n)))
    cgr99 = 0.01 - (1 / (51.82 + (21.66 * n)))
    cgamr = {"90%": cgr90, "95%": cgr95, "99%": cgr99}

    cesgam = {"cgamint": cgamint, "cgamr": cgamr}
    return cesgam
#---------------------------
# Se calculan niveles de verosimilitud de intervalos y regiones para 
# inferencias con la distribución Normal (uno y dos parámetros)
def cnorlevels(n):
    #para mu:
    cc90 = 0.2585 - (1 / (-0.4743 + (1.9029 * n)))
    cc95 = 0.1465 - (1 / (0.8469 + (2.3544 * n)))
    cc99 = 0.0362 - (1 / (8.753 + (5.551 * n)))
    ccesmu = {"90%": cc90, "95%": cc95, "99%": cc99}
    #para sigma:
    cs90 = 0.2585 - (1 / (0.775 + (1.641 * n)))
    cs95 = 0.1465 - (1 / (2.757 + (2.082 * n)))
    cs99 = 0.0362 - (1 / (18.294 + (5.229 * n)))
    ccessig = {"90%": cs90, "95%": cs95, "99%": cs99}
    #para (mu,sigma)
    ctous90 = 0.10 - (1 / (1.11 + (4.655 * n)))
    ctous95 = 0.05 - (1 / (1.61 + (7.27 * n)))
    ctous99 = 0.01 - (1 / (19.51 + (23.59 * n)))
    ctous = {"90%": ctous90, "95%": ctous95, "99%": ctous99}

    cesn = {"cmu": ccesmu, "csig": ccessig, "cmusig": ctous}
    return cesn
#---------------------------
# Para inferencias con distribuciones Weibull y Gumbel de mínimos: 
def K1bet(BETA, XX):
    if BETA > 0:
        k1 = np.sum(np.exp(BETA * XX))
    else:
        print("Error in K1: Beta must be positive")
        k1 = 0
    return k1
#---------------------------
def K2bet(BETA, XX):
    if BETA > 0:
        k2 = np.sum(XX * np.exp(BETA * XX))
    else:
        print("Error in K2: Beta must be positive")
        k2 = 0
    return k2
#---------------------------
def K3bet(BETA, XX):
    if BETA > 0:
        k3 = np.sum((XX ** 2) * np.exp(BETA * XX))
    else:
        print("Error in K3: Beta must be positive")
        k3 = 0
    return k3
#---------------------------
#  Para encontrar el emv del parámetro de forma beta Weibull dando con la 
# raíz de la primera derivada de la perfil de beta igual a cero:
# Los argumentos son beta, los logdatos XX, s suma TT y el tamaño muestral NN.
def Weibetmle(BETA, XX, TT, NN):
    if BETA > 0:
        k1 = K1bet(BETA, XX)
        k2 = K2bet(BETA, XX)
        sali = (NN / BETA) + TT - (NN * k2 / k1)
    else:
        sali = 5555
    return sali
#  ----------------------- -- 
#  Se calcula logverosimilitud Weibull (a,beta), beta es forma y aa logescala.
#  Los argumentos son VEC2=(a,beta), logdata XX, sumaTT, y tamaño mtra NN.
def lvabetaWei(VEC2, XX, TT, NN):
    AA = VEC2[0]
    BETA = VEC2[1]
    if BETA > 0:
        k1 = K1bet(BETA, XX)
        sali = NN * np.log(BETA) + (BETA * TT) - (NN * AA * BETA) - np.exp(-AA * BETA) * k1
    else:
        sali = -999999999999
    return sali
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Normal  rp(mu)
def rpmu(mu, t1, t2, n):
    HH = t2 - (t1 ** 2) / n
    Hmu = t2 - (2 * t1 * mu) + (n * (mu ** 2))
    rp = -(n / 2) * (np.log(Hmu) - np.log(HH))
    return rp

