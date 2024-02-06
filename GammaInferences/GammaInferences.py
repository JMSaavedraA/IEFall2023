####    INFERENCIAS GAMMA ____ ####
#   SE ILUSTRAN INFERENCIAS POSIBLES PARA PARAMETROS DE LA GAMMA:
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.special import polygamma
import AuxStatModel as ASM

# Ejemplo que queremos cargar
file_path = 'Ej1.pickle'

# Cargar los datos
with open(file_path, 'rb') as f:
    X = pickle.load(f)

xsort = np.sort(X) # datos ordenados
n = len(X)

#    Estadísticas suficientes:                                   ####
t1 = np.sum(X)  # Para distribuciones Exp, Normal, Lognormal, Gamma e Inv Gaussiana
t2 = np.sum(X ** 2) # Para la Normal y Lognormal
t3 = np.sum(np.log(X))  # Para la Gamma
# EMVs Gamma                                                           ####
# Estimador momentos forma Gamma:
alphamom = (t1**2) / ((n * t2) - (t1**2))
# EMV de media Gamma, Exponencial y Normal:
mumle = t1 / n
# EMV de parámetros de forma Gamma, alphamle,  y escala, betagamle:
root_func = lambda alph: ASM.Lgampsi(alph, t1, t3, n)
result = root_scalar(root_func, bracket=[alphamom / 10, alphamom * 10], xtol=0.000001)
alphamle = result.root
betagamle = mumle / alphamle

# Definir la función trigamma
trigamma = lambda x: polygamma(1, x)

#   Gamma Rp(alpha)
Psi1 = trigamma(alphamle)  # trigamma function evaluated at alphamle
infoalfa = n * Psi1 - (n / alphamle)
a1 = alphamle / 5
a2 = alphamle * 2
alfs = np.linspace(a1, a2, 200)

# Good approximating function proposed by Eloisa:
rpstalfa = -(9 * infoalfa * (alphamle**2) / 2) * ((1 - ((alfs / alphamle)**(1 / 3)))**2)
Rpstalfa = np.exp(rpstalfa)

# The exact Rp(alpha) is calculated here
rpalfa = ASM.lpgamalfa(alfs, t1, t3, n) - ASM.lpgamalfa(alphamle, t1, t3, n)
Rpalfa = np.exp(rpalfa)

# Exact Intervals alpha Gamma
# Computing exact profile likelihood intervals for alpha of levels c
cesgam = ASM.cgamlevels(n)["cgamint"]

alfizq90 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle / 7, alphamle], args=(t1, t3, n, alphamle, cesgam["90%"]), method='brentq'
).root

alfder90 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle, alphamle * 5], args=(t1, t3, n, alphamle, cesgam["90%"]), method='brentq'
).root

alfizq95 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle / 7, alphamle], args=(t1, t3, n, alphamle, cesgam["95%"]), method='brentq'
).root
alfder95 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle, alphamle * 5], args=(t1, t3, n, alphamle, cesgam["95%"]), method='brentq'
).root

alfizq99 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle / 7, alphamle], args=(t1, t3, n, alphamle, cesgam["99%"]), method='brentq'
).root

alfder99 = root_scalar(
    ASM.rpgamalfalnc, bracket=[alphamle, alphamle * 5], args=(t1, t3, n, alphamle, cesgam["99%"]), method='brentq'
).root

alfizq = np.array([alfizq90, alfizq95, alfizq99])
alfder = np.array([alfder90, alfder95, alfder99])
IntgamAlpha = np.column_stack((alfizq, alfder))

plt.figure()
plt.plot(alfs, Rpalfa, label=r"$R_p(\alpha)$", color='black', linewidth=1)
plt.plot(alfs, Rpstalfa, label=r"$R_{pst}(\alpha)$", linestyle='dashed', color='red')
#cvec = np.array(list(cesgam.values())) * np.array([[1, 1, 1], [1, 1, 1]])
CesMuS = np.array(list(cesgam.values()))
cvec = np.array([x * np.ones((2)) for x in CesMuS])

plt.plot(IntgamAlpha[0, :], cvec[0, :], linestyle='solid', color='blue')
plt.plot(IntgamAlpha[1, :], cvec[1, :], linestyle='solid', color='green')
plt.plot(IntgamAlpha[2, :], cvec[2, :], linestyle='solid', color='purple')
plt.title(r"$R_p(\alpha)$")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$R_p(\alpha)$")
plt.ylim(0, 1)
plt.legend()
plt.show()

##  ............................................................................
# Proposed Alpha Intervals vs exact:                        ####
# Proposed by Eloisa and Elías:

###### CHECK FROM HERE ###################
import numpy as np
from scipy.optimize import minimize_scalar

# Assuming mumle, t1, t3, n are already defined

# Function for log-likelihood of alpha in Gamma distribution
def lpgamalpha(alpha, t1, t3, n):
    return -n * np.log(alpha) - n * alpha * np.log(mumle) + (alpha - 1) * t3 - t1 / alpha

# Function for log-likelihood of mu in Gamma distribution
def lpgamamu(alpha, t1, t3, n, mu):
    return -n * np.log(alpha) - n * alpha * np.log(mu) + (alpha - 1) * t3 - t1 / alpha

# Search interval for alpha
alpha_interval = [0.001, alphamle * 15]

# Gamma profile likelihood intervals for alpha
ALF1s = alphamle * ((1 - np.sqrt(-2 * np.log(0.9) / infoalfa) / (3 * alphamle)) ** 3)
ALF2s = alphamle * ((1 + np.sqrt(-2 * np.log(0.9) / infoalfa) / (3 * alphamle)) ** 3)

# Display exact and approximate intervals for alpha
result_intervals = np.column_stack((alfizq, alfder, ALF1s, ALF2s))
print("Exact and Approximate Intervals for Alpha:")
print(result_intervals)

# Gamma profile likelihoods for mu
MU1 = mumle / 1.2
MU2 = mumle * 1.2
MM = 250
musis = np.linspace(MU1, MU2, MM)
lvecpgmu = np.zeros(MM)
alpharmle = np.zeros(MM)

for j in range(MM):
    interv = [0.001, alphamle * 15]  # Search interval for alpha
    MUU = musis[j]
    result = minimize_scalar(lambda alpha: -lpgamamu(alpha, t1, t3, n, MUU),
                             bounds=(min(interv), max(interv)),
                             method='bounded')
    lvecpgmu[j] = -result.fun
    alpharmle[j] = result.x

# Calculate relative likelihood for mu in Gamma distribution
Rpgmunum = np.exp(lvecpgmu - lpgamamu(alphamle, t1, t3, n, mumle))
