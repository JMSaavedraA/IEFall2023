#######     Loglikelihood for interval censored data    ###########

import numpy as np
from scipy.stats import weibull_min, norm, lognorm, gamma, expon, gumbel_l, invgauss

# Normal loglikelihood for interval censored data in matrix DATMAT (nx2)
def lvnorCEN(vec2, DATMAT):
    epsi = 0.000001
    mu = vec2[0]
    sigma = vec2[1]
    if sigma > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = norm.cdf(DATMAT[:, 1], loc=mu, scale=sigma) - norm.cdf(DATMAT[:, 0], loc=mu, scale=sigma)
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv
#------------------------
# Lognormal loglikelihood for interval censored data in matrix DATMAT (nx2)
def lvLognorCEN(vec2, DATMAT):
    epsi = 0.000001
    muLN = vec2[0]
    sigmaLN = vec2[1]
    if sigmaLN > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = lognorm.cdf(DATMAT[:, 1], s=sigmaLN, scale=np.exp(muLN)) - lognorm.cdf(
            DATMAT[:, 0], s=sigmaLN, scale=np.exp(muLN)
        )
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv
#-------------------------
# Gamma loglikelihood for interval censored data in matrix DATMAT (nx2)
def lvGamCEN(vec2, DATMAT):
    epsi = 0.000001
    alfa = vec2[0]
    mugam = vec2[1]
    betagam = mugam / alfa
    if alfa > 0 and mugam > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = gamma.cdf(DATMAT[:, 1], a=alfa, scale=betagam) - gamma.cdf(
            DATMAT[:, 0], a=alfa, scale=betagam
        )
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv

#------------------------
# Exponential loglikelihood for interval censored data in matrix DATMAT (nx2)
# This distribution is parametrized in terms of its mean theta

def lvExpCEN(vec1, DATMAT):
    epsi = 0.000001
    theta = vec1[0]
    if theta > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = expon.cdf(DATMAT[:, 1], scale=1/theta) - expon.cdf(
            DATMAT[:, 0], scale=1/theta
        )
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv

#------------------------
# Weibull of minima loglikelihood for interval censored data in matrix DATMAT (nx2)

def lvWeiCEN(vec2, DATMAT):
    epsi = 0.000001
    betaWei = vec2[0]  # shape parameter
    sigWei = vec2[1]  # scale parameter
    if betaWei > 0 and sigWei > 0:
        difi = weibull_min.cdf(DATMAT[:, 1], c=betaWei, scale=sigWei) - weibull_min.cdf(
            DATMAT[:, 0], c=betaWei, scale=sigWei
        )
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv

#------------------------
# Gumbel of minima loglikelihood for interval censored data in matrix DATMAT (nx2)

def lvGuminCEN(vec2, DATMAT):
    epsi = 0.000001
    agum = vec2[0]  # location parameter
    bgum = vec2[1]  # scale parameter
    if bgum > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = gumbel_l.cdf(DATMAT[:, 1], loc=agum, scale=bgum) - gumbel_l.cdf(
            DATMAT[:, 0], loc=agum, scale=bgum
        )
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv

#---------------------- --
#--------Inverse Gaussian log likelihood with interval censoring: --- --

def lvGauInvCEN(vec2, DATMAT):
    epsi = 0.000001
    muGI = vec2[0]  # mean parameter
    lamGI = vec2[1]  # shape parameter
    if lamGI > 0 and muGI > 0:
        # Prob. of each obs. of falling in its corresponding interval:
        difi = invgauss.cdf(DATMAT[:, 1], muGI, scale=lamGI) - invgauss.cdf(
            DATMAT[:, 0], muGI, scale=lamGI
        )
        # very small values are replaced before taking logarithm:
        difnew = np.where(difi < epsi, epsi, difi)
        lv = np.sum(np.log(difnew))
    else:
        lv = -99999999999999999  # discouraging impossible values for parameters
    return lv
