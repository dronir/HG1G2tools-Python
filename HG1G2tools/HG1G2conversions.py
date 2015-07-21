
from numpy import zeros, log10, exp

def a1a2a3_to_HG1G2(params):
    res = zeros(3)
    x = params.sum()
    res[0] = -2.5 * log10(x)
    res[1] = params[0] / x
    res[2] = params[1] / x
    return res

def HG1G2_to_a1a2a3(params):
    res = zeros(3)
    res += 10**(-0.4 * params[0])
    res[0] *= params[1]
    res[1] *= params[2]
    res[2] *= (1 - params[1] - params[2])
    return res

def a1a2_to_HG12(params):
    res = zeros(2)
    res[0] = -2.5 * log10(params[0])
    res[1] = params[1] / params[0]
    return res

def HG12_to_a1a2(params):
    res = zeros(2)
    res[0] = exp(-0.9210340371976184 * params[0])
    res[1] = res[0] * params[1]
    return res
