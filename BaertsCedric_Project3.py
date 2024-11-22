import numpy as np
import scipy
import scipy.integrate as scint
import matplotlib.pyplot as plt
from astropy import units as u

## Question 1 ##

# Constants
Mue = 2
R0 = 7.72e8 / Mue #cm
M0 = 5.67e33 / Mue**2 #g
Mch = 5.836 / Mue**2 #g
Rho0 = 9.74e5 * Mue #g/cm**3
MSun = 1.989e33 #g
RSun = 6.963e10 #cm

RhoC = np.logspace(-1, 6.4, 10)
RhoCP3 = [RhoC[3], RhoC[6], RhoC[9]]
MassP3 = []
RadiusP3 = []

def Equations8n9(r, y): 

    D0M, D0Rho = y # D0 before the variables refer to dimensionaless (would have done 0D if that was possible)
    x = D0Rho**(1/3)
    gammaX = x**2/(3*np.sqrt(1+x**2))
    dRhodR = - D0M*D0Rho/(gammaX*r**2)
    dMdR = r**2*D0Rho
    return  dMdR, dRhodR

def stateSolver(rho, methos):
    intState = [0, rho]
    rMin = 1e-4

    result = scint.solve_ivp(Equations8n9, [rMin, 1000], intState, events=lambda r, y: y[1], method=methos)
    return result

## Question 2 ##

fig, ax = plt.subplots()
for rho in RhoC:
    solution = stateSolver(rho, 'RK45')
    mass = (solution.y[0] * M0) / u.solMass.to(u.g) 
    density = solution.y[1] * Rho0
    radius = (solution.t * R0) / u.solRad.to(u.cm) 
    ax.plot(mass, radius, label=f"RhoC = {rho:.1e}")
    if rho in RhoCP3:
        RadiusP3.append(radius)
        MassP3.append(mass)
plt.ylabel("Radius (Solar Radii))")
plt.xlabel("Mass (Solar Masses)")
plt.title("Mass vs Radius")
plt.axvline(Mch, label="Chandrasekhar Mass Limit")
plt.legend()
plt.show()

## Question 3 ##
i = 0
fig, ax = plt.subplots()
for rho in RhoCP3:
    solution1 = stateSolver(rho, 'RK23')
    mass = (solution.y[0]*M0)/u.solMass.to(u.g) 
    density = solution.y[1]*Rho0
    radius = (solution.t*R0)/u.solRad.to(u.cm) 
    ax.plot(mass, radius, label=f"RhoC = {rho:.1e}, Method = RK23")
    ax.plot(MassP3[i], RadiusP3[i], label=f"RhoC = {rho:.1e}, Method = RK45", linestyle='dotted')
    i += 1

plt.ylabel("Radius (Solar Radii))")
plt.xlabel("Mass (Solar Masses)")
plt.title("Mass vs Radius")
plt.axvline(Mch, label="Chandrasekhar Mass Limit")
plt.legend()
plt.show()