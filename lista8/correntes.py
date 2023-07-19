def J_inj(j, t, ti, tf):
    return j*(t >= ti) - j*(t >= tf)


def J_sin(g_sin, s, E_rev, v):
    return g_sin * s * (E_rev - v)
