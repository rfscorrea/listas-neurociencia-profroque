from numpy import exp, tanh


def neuronioI(v, h, n, s, v_pre, j_sin, j_inj):
    C_m = 1 #uF/cm2
    tau_s = 0.3 #ms
    tau_d = 9   #ms
    g_Na = 35  #mS/cm2
    g_V = 0.1   #mS/cm2
    g_K = 9   #mS/cm2
    E_K = -90  #mV
    E_Na = 55   #mV
    E_V = -65   #mV
    E_rev = -80   #mV

    dvdt = (g_Na * x_inf('m', v)**3 * h * (E_Na - v) + g_K * n**4 * (E_K - v) + g_V*(E_V - v) + j_sin + j_inj)/C_m
    dhdt = (x_inf('h', v) - h) * (alpha('h', v) + beta('h', v)) / 0.2
    dndt = (x_inf('n', v) - n) * (alpha('n', v) + beta('n', v)) / 0.2
    dsdt = ((1 + tanh(v_pre / 4)) / 2) * ((1 - s) / tau_s) - (s / tau_d)

    return dvdt, dhdt, dndt, dsdt


def x_inf(x: str, v: float):
    return alpha(x, v) / alpha(x, v) + beta(x, v)


def alpha(x: str, v: float):
    alphas = {
        'm': 0.1*(v+35)/(1-exp(-(v+35)/10)) if v != -54 else 0.1*10/exp(-(v+35)/10),
        'h': 0.07*exp(-(v+58)/20),
        'n': 0.01*(v+34)/(1-exp(-0.1*(v+34))) if v != -52 else 0.01/0.1*exp(-0.1*(v+34)),
    }
    return alphas[x]


def beta(x: str, v: float):
    betas = {
        'm': 4*exp(-(v+60)/18),
        'h': 1/(exp(-0.1*(v+28))+1),
        'n': 0.125*exp(-(v+44)/80)
    }
    return betas[x]
