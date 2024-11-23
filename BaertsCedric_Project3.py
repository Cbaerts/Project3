import numpy as np
import scipy
import scipy.integrate as scint
import matplotlib.pyplot as plt
from astropy import units as u
import csv

## Question 1 ##

# Constants
Mue = 2
R0 = 7.72e8 / Mue #cm
M0 = 5.67e33 / Mue**2 #g
Mch = 5.836 / Mue**2 #g
Rho0 = 9.74e5 * Mue #g/cm**3
MSun = 1.989e33 #g
RSun = 6.963e10 #cm

# RhoC is in logspace since we want to see how the mass of various magnitudes affect this process (linspace would be around 5/6 orders of magnitude for all values.)
RhoC = np.logspace(-1, 6.4, 10)
# P3 suffix is for arrays used for Question 3
RhoCP3 = [RhoC[3], RhoC[6], RhoC[9]]
MassP3 = []
RadiusP3 = []


def Equations8n9(r, y): 
    '''
The function takes in the radius, mass and density at those radii then computes Equations 8 & 9 from the project outline. 
These equations represent the derivatives of the density and the mass as a function of the radii.

Parameters
----------
r: float
    The current radius of the star (cm)
y: ndarray
    States the dimensionless values for the mass and density at the radii (g), and (g/cm^3)

Return
------
dMdR: float
    The Change in the mass over the change in the radius. Equation 8
dRhodR: float
    The Change in the density over the change in the radius. Equation 9
    '''
    # D0 before the variables refer to dimensionaless (would have done 0D if that was possible)
    # All equations are given in the project 3 details. 
    # y is an array which contains the dimensionless mass and density
    D0M, D0Rho = y 
    x = D0Rho**(1/3)
    gammaX = x**2/(3*np.sqrt(1+x**2))
    dRhodR = - D0M*D0Rho/(gammaX*r**2)
    dMdR = r**2*D0Rho
    return  dMdR, dRhodR

def stateSolver(rho, methos):
    '''
The function takes in the density(RhoC) and the integration method. 
Solves the initial value problem based on the initial state and the given density.

Parameters
----------
rho: float
    The current density of the star (g)
methos: str
    The integration method for the solve_ivp function

Return
------
result: float
    The solution to the IVP when the density becomes zero.
    '''
    # Initial Conditions
    intState = [0, rho]
    # Minimum radius to avoid dividing by zero
    rMin = 1e-4
    # result of the IVP. Will stop iterating when the density becomes zero (controlled by the "events=lambda r, y: y[1]")
    result = scint.solve_ivp(Equations8n9, [rMin, 1000], intState, events=lambda r, y: y[1], method=methos)
    return result

## Question 2 ##
# Creates arrays for Question 4
radP4 = []
massP4 = []
# Allow for multiple lines on the plot
fig, ax = plt.subplots()
# Iterate for every density in the RhoC array.
for rho in RhoC:
    # Calls the function to run the IVP for each of the rho values. RK45 is the default integration method.
    solution = stateSolver(rho, 'RK45')
    # Sets the 0-D solution arrays of the IVP to their actual value in solar masses, g/cm^3, and solar radii respectively. 
    mass = (solution.y[0] * M0) / u.solMass.to(u.g) 
    density = solution.y[1] * Rho0
    radius = (solution.t * R0) / u.solRad.to(u.cm) 
    # Appends whole array to be used in question 4
    radP4.append(radius)
    massP4.append(mass)
    ax.plot(mass, radius, label=f"RhoC = {rho:.1e}")
    # Appends the array solutions to be used in question 3
    if rho in RhoCP3:
        RadiusP3.append(radius)
        MassP3.append(mass)
# Plot answers
plt.ylabel("Radius (Solar Radii))")
plt.xlabel("Mass (Solar Masses)")
plt.title("Mass vs Radius")
plt.axvline(Mch, label="Chandrasekhar Mass Limit")
plt.legend()
plt.show()
plt.close()

## Question 3 ##
i = 0
fig, ax = plt.subplots()
# Performs same function as for statement at line 88 but for RhoCP3
for rho in RhoCP3:
    solution = stateSolver(rho, 'RK23')
    mass = (solution.y[0]*M0)/u.solMass.to(u.g) 
    density = solution.y[1]*Rho0
    radius = (solution.t*R0)/u.solRad.to(u.cm) 
    # Plots solution to ivp
    ax.plot(mass, radius, label=f"RhoC = {rho:.1e}, Method = RK23")
    # Dotted previous solution for easy comparison
    ax.plot(MassP3[i], RadiusP3[i], label=f"RhoC = {rho:.1e}, Method = RK45", linestyle='dotted')
    i += 1
# Plots solutions of the new integration method compared to the default method
plt.ylabel("Radius (Solar Radii))")
plt.xlabel("Mass (Solar Masses)")
plt.title("Mass vs Radius")
plt.axvline(Mch, label="Chandrasekhar Mass Limit")
plt.legend()
plt.show()
plt.close()

## Question 4 ##
# Name data file to be used later
datafile = 'wd_mass_radius.csv'
# Create arrays
WDMass, WDRadius, WDMassUnc, WDRadUnc = [], [], [], []
# Open the datafile and read it
with open(datafile, 'r') as data:
        # Skip line 1 (Very boring line)
        for line in data.readlines()[1:]:
            # Strips usless line parts (commas) then seperates them there.
            cols = line.strip().split(',')
            # Appends the values to their respective arrays
            WDMass.append(float(cols[0]))
            WDMassUnc.append(float(cols[1]))
            WDRadius.append(float(cols[2]))
            WDRadUnc.append(float(cols[3]))

fig, ax = plt.subplots()
# Plots solutions from Question 2
for i in range(len(massP4)):
    ax.plot(massP4[i], radP4[i], label=f"RhoC = {RhoC[i]:.1e}")
# Plot errorbars and points for the given observed white dwarf data
plt.errorbar(WDMass, WDRadius, xerr=(WDMassUnc), yerr=(WDRadUnc), linestyle='', color='black', zorder=-1)
plt.scatter(WDMass, WDRadius, s=10, label="Observational Data", color='dodgerblue', zorder=1)
plt.title('Plot of the Mass vs Radii of White Dwarfs in Binary Star Systems')
plt.xlabel('Mass (Solar Masses)')
plt.ylabel('Radius (Solar Radii)')
plt.legend()
plt.axvline(Mch, label="Chandrasekhar Mass Limit")
plt.title('Plot of the Mass vs Radii of White Dwarfs in Binary Star Systems')
plt.show()
plt.close()

