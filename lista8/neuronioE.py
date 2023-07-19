from numpy import exp, tanh


def neuronioE(v, h, n, s, v_pre, j_sin, j_inj):
    C_m = 1 #uF/cm2
    tau_s = 0.1 #ms
    tau_d = 3   #ms
    g_Na = 100  #mS/cm2
    g_V = 0.1   #mS/cm2
    g_K = 80    #mS/cm2
    E_K = -100  #mV
    E_Na = 50   #mV
    E_V = -67   #mV
    E_rev = 0   #mV

    dvdt = (g_Na * x_inf('m', v)**3 * h * (E_Na - v) + g_K * n**4 * (E_K - v) + g_V*(E_V - v) + j_sin + j_inj)/C_m
    dhdt = (x_inf('h', v) - h) * (alpha('alpha_h', v) + beta('beta_h', v))
    dndt = (x_inf('n', v) - n) * (alpha('alpha_n', v) + beta('beta_n', v))
    dsdt = ((1 + tanh(v_pre / 4)) / 2) * ((1 - s) / tau_s) - (s / tau_d)

    return dvdt, dhdt, dndt, dsdt



def x_inf(var: str, v: float):
    return alpha(f'alpha_{var}', v) / alpha(f'alpha_{var}', v) + beta(f'beta_{var}', v)


def alpha(x: str, v: float):
    alphas = {
        'm': 0.32*(v+54)/(1-exp(-(v+54)/4)) if v != -54 else 0.32*4/exp(-(v+54)/4),
        'h': 0.128*exp(-(v+50)/18),
        'n': 0.032*(v+52)/(1-exp(-(v+52)/5)) if v != -52 else 0.32*5/exp(-(v+52)/5),
    }
    return alphas[x]


def beta(x: str, v: float):
    betas = {
        'm': 0.28 * (v + 27) / (exp((v + 27)/5) - 1) if v != -27 else 0.28*5/exp(-(v+27)/5),
        'h': 4/(1+exp(-(v+27)/5)),
        'n': 0.5*exp(-(v+57)/40)
    }
    return betas[x]
