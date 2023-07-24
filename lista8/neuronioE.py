from numpy import exp, tanh


def dvdt(v,h,n,j):
    C_m = 1
    g_Na = 100
    g_K = 80
    g_V = 0.1
    E_Na = 50
    E_K = -100
    E_V = -67
    eq = (g_Na *h* x_inf('m',v)**3 * (E_Na-v) + g_K * n**4 * (E_K-v) + g_V * (E_V-v) + j)/C_m
    return eq


def dhdt(v, h):
    eq = (x_inf('h', v) - h) / tau_x('h', v)
    return eq


def dndt(v, n):
    eq = (x_inf('n', v) - n) / tau_x('n', v)
    return eq


def dsdt(v, s):
    taus = 0.1
    taud = 3
    eq = (1+tanh(v/4)/2) * (1-s)/taus - s/taud
    return eq


def x_inf(x, v):
    return alpha(x, v) / (alpha(x, v)+beta(x, v))


def tau_x(x, v):
    return 1 / (alpha(x, v)+beta(x, v))


def alpha(x, v):
    eq = {
        'm': 0.32*(v+54)/(1-exp(-(v+54)/4)),
        'h': 0.128*exp(-(v+50)/18),
        'n': 0.032*(v+52)/(1-exp(-(v+52)/5))
    }
    return eq[x]


def beta(x, v):
    eq = {
        'm': 0.28*(v+27)/(exp((v+27)/5)-1),
        'h': 4/(1+exp(-(v+27)/5)),
        'n': 0.5*exp(-(v+57)/40)
    }
    return eq[x]
