"""
Shock and Detonation Toolbox Demo Program

Calculates the constant pressure explosion state.

################################################################################
Theory, numerical methods and applications are described in the following report:

    Numerical Solution Methods for Shock and Detonation Jump Conditions, S.
    Browne, J. Ziegler, and J. E. Shepherd, GALCIT Report FM2006.006 - R3,
    California Institute of Technology Revised August, 2017

Please cite this report and the website if you use these routines. 

Please refer to LICENCE.txt or the above report for copyright and disclaimers.

http://shepherd.caltech.edu/EDL/public/sdt/SD_Toolbox/

################################################################################ 
Updated January 2018
Tested with: 
    Python 3.5 and 3.6, Cantera 2.3
Under these operating systems:
    Windows 10, Linux (Debian 9)
"""

import cantera as ct
from sdtoolbox.thermo import soundspeed_eq,soundspeed_fr

print('demo_HPstate')
## set initial state, composition, and gas object
transport = False # not all mechanisms have transport data defined
P1 = ct.one_atm; T1 = 300 
q = 'H2:2 O2:1  N2:3.76'    
mech = 'Mevel-2009-PCI32.cti'
# gas = ct.Solution(mech,'keromnes2013_mix')
gas = ct.Solution(mech)
gas.TPX = T1,P1,q
rho1 = gas.density
if transport:
    mu1 = gas.viscosity
    nu1 = mu1/rho1
    kcond1 = gas.thermal_conductivity
    kdiff1 = kcond1/(rho1*gas.cp_mass)

## Find HP state
gas.equilibrate('HP')
## Evaluate properties of gas object 
T2 = gas.T
P2 = gas.P
rho2 = gas.density
V2 = 1/rho2
S2 = gas.entropy_mass
x2 = gas.X
c2_eq = soundspeed_eq(gas)
c2_fr = soundspeed_fr(gas)
gamma2_fr =  c2_fr**2*rho2/P2
gamma2_eq =  c2_eq**2*rho2/P2
if transport:
    mu = gas.viscosity
    nu = mu/rho2
    kcond = gas.thermal_conductivity
    kdiff = kcond/(rho2*gas.cp_mass)
    Pr = mu*gas.cp_mass/kcond

## Print out
print('HP computation for '+mech+' with composition '+q)
print('Initial State')
print('   Pressure '+str(P1)+' (Pa)')
print('   Temperature '+str(T1)+' (K)')
print('   Density '+str(rho1)+' (kg/m3)')
print('   Specific heat at constant pressure '+str(gas.cp_mass)+' (J/kg K)')
if transport:
    print('   viscosity '+str(mu1)+' (kg/m s)')
    print('   viscosity (kinematic)'+str(nu1)+' (m2/s)')
    print('   thermal conductivity '+str(kcond1)+' (W/m K)')
    print('   thermal diffusivity '+str(kdiff1)+' (m2/s)')

print('HP State');
print('   Pressure '+str(P2)+' (Pa)')
print('   Temperature '+str(T2)+' (K)')
print('   Density '+str(rho2)+' (kg/m3)')
print('   Specific heat at constant pressure '+str(gas.cp_mass)+' (J/kg K)')
print('   Expansion  ratio '+str(rho1/rho2)+' (kg/m3)')
print('   Entropy '+str(S2)+' (J/kg-K)')
print('   Sound speed (frozen) '+str(c2_fr)+' (m/s)')
print('   Sound speed (equilibrium) '+str(c2_eq)+' (m/s)')
print('   gamma2 (based on frozen sound speed) '+str(gamma2_fr)+' (m/s)')
print('   gamma2 (based on equilibrium sound speed) '+str(gamma2_eq)+' (m/s)')
if transport:
    print('   viscosity '+str(mu)+' (kg/m s)')
    print('   viscosity (kinematic)'+str(nu)+' (m2/s)')
    print('   thermal conductivity '+str(kcond)+' (W/m K)')
    print('   thermal diffusivity '+str(kdiff)+' (m2/s)')
    print('   specific heat Cp '+str(gas.cp_mass)+' (J/kg K)')
    print('   Prandtl number '+str(Pr))

print('Reactants at HP temperature')
gas.TPX = T2,P1,q
rho1 = gas.density
print('   Temperature '+str(T2)+' (K)')
print('   Density '+str(rho1)+' (kg/m3)')
print('   Specific heat at constant pressure '+str(gas.cp_mass)+' (J/kg K)')
if transport:
    mu1 = gas.viscosity
    nu1 = mu1/rho1;
    kcond1 = gas.thermal_conductivity
    kdiff1 = kcond1/(rho1*gas.cp_mass)
    print('   viscosity '+str(mu1)+' (kg/m s)')
    print('   viscosity (kinematic)'+str(nu1)+' (m2/s)')
    print('   thermal conductivity '+str(kcond1)+' (W/m K)')
    print('   thermal diffusivity '+str(kdiff1)+' (m2/s)')

