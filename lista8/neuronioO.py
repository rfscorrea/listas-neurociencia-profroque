from numpy import tanh, exp


def dvdt(v, m, h, n, a, b, r, j):
    Cm = 1.3
    gNa = 30
    gK = 23
    gV = 0.05
    gA = 16
    gh = 12
    ENa = 90
    EK = -100
    EV = -70
    EA = -90
    Eh = -32.9
    eq = (gNa * m**3 * h*(ENa-v) + gK * n**4 * (EK-v) + gA*a*b*(EA-v) + gh*r*(Eh-v) + gV*(EV-v) + j)/Cm
    return eq


def dmdt(v, m):
    eq = (x_inf('m', v) - m)/tau_x('m', v)
    return eq


def dhdt(v ,h):
    eq = (x_inf('h', v) - h) / tau_x('h', v)
    return eq


def dndt(v, n):
    eq = (x_inf('n', v) - n) / tau_x('n', v)
    return eq


def dadt(v, a):
    eq = (x_inf('a', v) - a) / tau_x('a', v)
    return eq


def dbdt(v, b):
    eq = (x_inf('b', v) - b) / tau_x('b', v)
    return eq


def drdt(v, r):
    eq = (x_inf('r', v) - r) / tau_x('r', v)
    return eq


def dsdt(v, s):
    taus = 0.2
    taud = 20
    eq = ((1+tanh(v/4))/2) * ((1-s)/taus) - (s/taud)
    return eq


def x_inf(x, v):
    eq = {
        'm': alpha('m', v)/(alpha('m',v)+beta('m',v)),
        'h': alpha('h', v)/(alpha('h',v)+beta('h',v)),
        'n': alpha('n', v)/(alpha('n',v)+beta('n',v)),
        'a': 1/(1+exp(-(v+14)/16.6)),
        'b': 1/(1+exp((v+71)/7.3)),
        'r': 1/(1+exp((v+84)/10.2))
    }
    return eq[x]


def tau_x(x, v):
    eq = {
        'm': 1/(alpha('m',v)+beta('m',v)),
        'h': 1/(alpha('h',v)+beta('h',v)),
        'n': 1/(alpha('n',v)+beta('n',v)),
        'a': 5,
        'b': 1/( ( 0.000009/exp((v-26)/18.5) ) + ( 0.014/(0.2+exp(-(v+70)/11)) ) ),
        'r': 1/( ( exp(-14.59-0.086*v) ) + ( exp(-1.87+0.0701*v) ) )
    }
    return eq[x]


def alpha(x, v):
    eq = {
        'm': -0.1*(v+38)/(exp(-(v+38)/10)-1),
        'h': 0.07*exp(-(v+63)/20),
        'n': 0.018*(v-25)/(1-exp(-(v-25)/25)),
    }
    return eq[x]


def beta(x, v):
    eq = {
        'm': 4*exp(-(v+65)/18),
        'h': 1/(1+exp(-(v+33)/10)),
        'n': 0.0036*(v-35)/(exp((v-35)/12)-1),
    }
    return eq[x]
