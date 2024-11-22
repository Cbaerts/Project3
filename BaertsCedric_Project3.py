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
    return  dMdR, dRhodR

def stateSolver(rho):
    intState = [0, rho]
 
    # WDrLimit = 9.74*10**8 #radius of white dwarf in m
    rMin = 1e-4
    result = scint.solve_ivp(Equations8n9, [rMin, 1000], intState, events=lambda r, y: y[1])
    return result

fig, ax = plt.subplots()
for rho in RhoC:
    solution = stateSolver(rho)
    mass = solution.y[0] * M0
    density = solution.y[1] * Rho0
    radius = solution.t * R0
    ax.plot(mass, radius, label=f"RhoC = {rho:.1e}")

plt.ylabel("Radius (cm)")
plt.xlabel("Mass (g)")
plt.title("Mass vs Radius")
plt.legend()
plt.show()

