from numpy import exp, tanh


def neuronioO(v, h, n, m, a, b, r, s, v_pre, j_sin, j_inj):
    C_m = 1.3 #uF/cm2
    tau_s = 0.2 #ms
    tau_d = 20   #ms
    g_Na = 30  #mS/cm2
    g_V = 0.05   #mS/cm2
    g_K = 23   #mS/cm2
    g_A = 16    #mS/cm2
    g_h = 12    #mS/cm2
    E_K = -100  #mV
    E_Na = 90   #mV
    E_V = -70   #mV
    E_A = -90   #mV
    E_h = -32.9 #mV
    E_rev = -80   #mV

    dvdt = (g_Na * m**3 * h * (E_Na - v) + g_K * n**4 * (E_K - v) + g_A*a*b*(E_A - v)
            + g_h*r*(E_h - v) + g_V*(E_V - v) + j_sin + j_inj)/C_m
    dmdt = (x_inf('m', v) - m) / tau_x('m', v)
    dhdt = (x_inf('h', v) - m) / tau_x('h', v)
    dndt = (x_inf('n', v) - m) / tau_x('n', v)
    dadt = (1/(1+exp(-(v+14)/16.6)) - a) / tau_x('a', v)
    dbdt = (1/(1+exp((v+71)/7.3)) - b) / tau_x('b', v)
    drdt = (1/(1+exp((v+84)/10.2)) - r) / tau_x('r', v)
    dsdt = ((1 + tanh(v_pre / 4)) / 2) * ((1 - s) / tau_s) - (s / tau_d)

    return dvdt, dmdt, dhdt, dndt, dadt, dbdt, drdt, dsdt


def tau_x(x: str, v: float):
    taus = {
        'm': 1 / (alpha('m', v) + beta('m', v)),
        'h': 1 / (alpha('h', v) + beta('h', v)),
        'n': 1 / (alpha('n', v) + beta('n', v)),
        'a': 5,
        'b': 1 / ((0.000009/exp((v-26)/18.5))+(0.014/(0.2+exp(-(v+70)/11)))),
        'r': 1 / (exp(-14.59-0.086*v)+exp(-1.87+0.0701*v))

    }
    return taus[x]


def x_inf(x: str, v: float):
    return alpha(x, v) / alpha(x, v) + beta(x, v)


def alpha(x: str, v: float):
    alphas = {
        'm': -0.1*(v+38)/(exp(-(v+38)/10)-1) if v != -38 else 1/exp(-(v+38)/10),
        'h': 0.07*exp(-(v+63)/20),
        'n': 0.018*(v-25)/(1-exp(-(v-25)/25)) if v != 25 else 0.018*25/exp(-(v-25)/25),
    }
    return alphas[x]


def beta(x: str, v: float):
    betas = {
        'm': 4*exp(-(v+65)/18),
        'h': 1/(exp(-(v+33)/10)+1),
        'n': 0.0036*(v-35)/(exp((v-35)/12) - 1) if v != 35 else 0.0036*12/exp((v-35)/12)
    }
    return betas[x]
