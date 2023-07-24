from numpy import exp


def dvdt(v, m, h, n, a, b, r, jinj, jsin):
    gNa, gK, gV = 30, 23, 0.05
    gA, gh = 16, 12
    ENa, EK, EV = 90, -100, -70
    EA, Eh = -90, -32.9
    eq = gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) + gA * a * b * (EA-v) + gh * r * (Eh-v) + gV * (EV-v) + jsin + jinj
    return eq


def dmdt(v, m):
    minf = x_inf('m', v)
    taum = tau_x('m', v)
    eq = (minf - m) / taum
    return eq


def dhdt(v, h):
    hinf = x_inf('h', v)
    tauh= tau_x('h', v)
    eq = (hinf - h) / tauh
    return eq


def dndt(v, n):
    ninf = x_inf('n', v)
    taun = tau_x('n', v)
    eq = (ninf - n) / taun
    return eq


def dadt(v, a):
    ainf = x_inf('a', v)
    taua = tau_x('a', v)
    eq = (ainf - a) / taua
    return eq


def dbdt(v, b):
    binf = x_inf('b', v)
    taub = tau_x('b', v)
    eq = (binf - b) / taub
    return eq


def drdt(v, r):
    rinf = x_inf('r', v)
    taur = tau_x('r', v)
    eq = (rinf - r) / taur
    return eq


def x_inf(x, v):
    alpham = -0.1 * (v + 38) / (exp(-(v + 38) / 10) - 1)
    alphah = 0.07 * exp(-(v + 63) / 20)
    alphan = 0.018 * (v - 25) / (1 - exp(-(v - 25) / 25))
    betam = 4 * exp(-(v + 65) / 18)
    betah = 1 / (1 + exp(-(v + 33) / 10))
    betan = 0.0036 * (v - 35) / (exp((v - 35) / 12) - 1)
    eq = {
        'm': alpham / (alpham+betam),
        'h': alphah / (alphah+betah),
        'n': alphan / (alphan+betan),
        'a': 1 / (1+exp(-(v+14)/16.6)),
        'b': 1 / (1+exp((v+71)/7.3)),
        'r': 1 / (1+exp((v+84)/10.2))
    }
    return eq[x]


def tau_x(x, v):
    alpham = -0.1 * (v + 38) / (exp(-(v + 38) / 10) - 1)
    alphah = 0.07 * exp(-(v + 63) / 20)
    alphan = 0.018 * (v - 25) / (1 - exp(-(v - 25) / 25))
    betam = 4 * exp(-(v + 65) / 18)
    betah = 1 / (1 + exp(-(v + 33) / 10))
    betan = 0.0036 * (v - 35) / (exp((v - 35) / 12) - 1)
    eq = {
        'm': 1 / (alpham+betam),
        'h': 1 / (alphah+betah),
        'n': 1 / (alphan+betan),
        'a': 5,
        'b': 1 / (0.000009/exp((v-26)/18.5) + 0.014/(0.2+exp(-(v+70)/11))),
        'r': 1 / (exp(-14.59-0.086*v) + exp(-1.87+0.0701*v))
    }
    return eq[x]
