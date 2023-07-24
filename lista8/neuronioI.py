from numpy import exp


def dvdt(v, h, n, jinj, jsin):
    Cm = 1
    gNa, gK, gV = 35, 9, 0.1
    ENa, EK, EV = 55, -90, -65
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
    eq = 0.2 / (a + b)
    return eq


def alpha(x, v):
    eq = {
        'm': 0.1 * (v+35) / (1-exp(-(v+35)/10)),
        'h': 0.07 * exp(-(v+58)/10),
        'n': 0.01 * (v+34) / (1-exp(-0.1*(v+34)))
    }
    return eq[x]


def beta(x, v):
    eq = {
        'm': 4 * exp(-(v+60)/18),
        'h': 1 / exp(-0.1*(v+28)+1),
        'n': 0.125 * exp(-(v+44)/80)
    }
    return eq[x]

