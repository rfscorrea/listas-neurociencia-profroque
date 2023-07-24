from numpy import exp, tanh


def dvdt(v,h,n,j):
    C_m = 1
    g_Na = 35
    g_K = 9
    g_V = 0.1
    E_Na = 55
    E_K = -90
    E_V = -65
    eq = (g_Na *h* x_inf('m',v)**3 * (E_Na-v) + g_K * n**4 * (E_K-v) + g_V * (E_V-v) + j)/C_m
    return eq


def dhdt(v, h):
    eq = (x_inf('h', v) - h) / tau_x('h', v)
    return eq


def dndt(v, n):
    eq = (x_inf('n', v) - n) / tau_x('n', v)
    return eq


def dsdt(v, s):
    taus = 0.3
    taud = 9
    eq = (1+tanh(v/4)/2) * (1-s)/taus - s/taud
    return eq


def x_inf(x, v):
    return alpha(x, v) / (alpha(x, v)+beta(x, v))


def tau_x(x, v):
    return 0.2 / (alpha(x, v)+beta(x, v))


def alpha(x, v):
    eq = {
        'm': 0.1*(v+35)/(1-exp(-(v+35)/10)),
        'h': 0.07*exp(-(v+58)/20),
        'n': 0.01*(v+34)/(1-exp(-0.1*(v+34)))
    }
    return eq[x]


def beta(x, v):
    eq = {
        'm': 4*exp(-(v+60)/18),
        'h': 1/(exp(-0.1*(v+28))+1),
        'n': 0.125*exp(-(v+44)/80)
    }
    return eq[x]
