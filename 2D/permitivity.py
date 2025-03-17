import numpy as np
def Drude(lam, omega_p, gamma):
    w=2.0*np.pi*3e8/lam
    eps = 1 - omega_p**2 / (w**2 + 1j*gamma*w)
    return eps