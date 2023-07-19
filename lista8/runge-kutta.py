"""
t = np.arange(0, tmax, dt)
J_inj = j*(t >= ti) - j*(t >= tf)

j_inj = J_inj[k]
j_sin = j_sin(*args)
v_E, h_E, n_E = v_E[k], h_E[k], n_E[k]
v_I, h_I, n_I = v_I[k], h_I[k], n_I[k]

equações do neurônio E = neuronioE(v_E, h_E n_E, j_sin, j_inj)
equações do neurônio I = neuronioI(v_I, h_I, n_I, j_sin, j_inj)

dsdt = dsdt(*args)
"""
