import numpy as np
import scipy
import scipy.integrate as scint
import matplotlib.pyplot as plt

# Constants
Mue = 2
R0 = 7.72e8 / Mue #cm
M0 = 5.67e33 / Mue**2 #g
Mch = 5.836 / Mue**2 #g
Rho0 = 9.74e5 * Mue #g/cm**3
MSun = 1.989e33 #g
RSun = 6.963e10 #cm

RhoC = np.logspace(-1, 6.4, 10)

def Equations8n9(r, y): 

    D0M, D0Rho = y # D0 before the variables refer to dimensionaless (would have done 0D if that was possible)
    x = D0Rho**(1/3)
    gammaX = x**2/(3*np.sqrt(1+x**2))
    dRhodR = - D0M*D0Rho/(gammaX*r**2)
    dMdR = r**2*D0Rho
    return dRhodR, dMdR

def stateSolver(rho):
    intState = [0, rho]
 
    # WDrLimit = 9.74*10**8 #radius of white dwarf in m
    rRange = [1e-10, 1e10]
    def denisityCheck(r, y):
        return y[1]
    result = scint.solve_ivp(Equations8n9, rRange, intState, rtol=1e-3, atol=1e-8, events=denisityCheck)
    return result


for rho in RhoC:
    solution = stateSolver(rho)

    print(solution)
    