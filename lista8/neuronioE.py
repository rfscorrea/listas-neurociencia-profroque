from numpy import exp


def dvdt(v, h, n, jinj, jsin):
    Cm = 1
    gNa, gK, gV = 100, 80, 0.1
    ENa, EK, EV = 50, -100, -67
    minf = x_inf('m', v)
    eq = (gNa * minf**3 * h * (ENa-v) + gK * n**4 * (EK-v) + gV*(EV-v) + jsin + jinj) / Cm
    return eq


def dhdt(v, h):
    hinf = x_inf('h', v)
    tauh = tau_x('h', v)
    eq = (hinf - h) / tauh
    return eq


def dndt(v, n):
    ninf = x_inf('n', v)
    taun = tau_x('n', v)
    eq = (ninf - n) / taun
    return eq


def x_inf(x, v):
    a = alpha(x, v)
    b = beta(x, v)
    eq = a / (a + b)
    return eq


def tau_x(x, v):
    a = alpha(x, v)
    b = beta(x, v)
    eq = 1 / (a + b)
    return eq


def alpha(x, v):
    eq = {
        'm': 0.32 * (v+54) / (1-exp(-(v+54)/4)),
        'h': 0.128 * exp(-(v+50)/18),
        'n': 0.032 * (v+52) / (1-exp(-(v+52)/5))
    }
    return eq[x]


def beta(x, v):
    eq = {
        'm': 0.28 * (v+27) / (exp((v+27)/5)-1),
        'h': 4 / (1+exp(-(v+27)/5)),
        'n': 0.5 * exp(-(v+57)/40)
    }
    return eq[x]
