import numpy as np
import scipy
import scipy.integrate as scint
import matplotlib.pyplot as plt

# Constants
Mue = 2
R0 = 7.72*10**8 / Mue
M0 = 5.67*10**33 / Mue**2
Mch = 5.836 / Mue**2
Rho0 = 9.74*10**5 * Mue

def solver():