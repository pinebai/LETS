#testing Cantera 2.3 versions Thermo routines of SDT with Python 3.5.2
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import math

import sdtoolbox as sdt

mech = 'gri30.cti'
#mech_xml = 'GRI30.xml'
gas = ct.Solution(mech)
##q = 'C2H2:2.0,  O2:5.0, AR:21.0'
##gas.X ='C2H2:2.0,  O2:5.0, AR:21.0'
##gas.X = 'C2H2:0.333,O2:0.167,AR:0.5'
##T0=298.0
##P0=5000
##gas.TP = T0, P0
##gas()
##rho = gas.density
##gas.equilibrate('UV')
##gas()
##ae = soundspeed_eq(gas)
##print('equilibrium sound speed ', ae)
##af = soundspeed_fr(gas)
##print('frozen sound speed ', af)
##G_eq = gruneisen_eq(gas)
##print('equilibrium gruneisen coefficient ', G_eq)
##G_fr = gruneisen_fr(gas)
####print('frozen gruneisen coefficient ', G_fr)
##[Ucj, r2] = CJspeed(P0, T0, q, mech, 0)
##print(Ucj, r2)
##gas = PostShock_eq(Ucj, P0, T0, q, mech)
##gas()
##comp1 = 'C2H2:0.2,O2:0.2,AR:0.6'
q = 'C2H2:0.333,O2:0.167,AR:0.5'
T0 = 300.
P0 = 10000.
gas.TPX = T0, P0, 'C2H2:0.333,O2:0.167,AR:0.5'
rho0 = gas.density
h0 = gas.enthalpy_mass
gas()
[Ucj, r2] = sdt.CJspeed(P0, T0, q, mech, 0)
print('CJ speed', Ucj,' fit correlation', r2)
gas = sdt.PostShock_eq(Ucj, P0, T0, q, mech)
gas()
rho2 = gas.density
h2 = gas.enthalpy_mass
p2 = gas.P
# equilibrium sound speed
ae = sdt.soundspeed_eq(gas)
af = sdt.soundspeed_fr(gas)
print('equilibrium sound speed', ae,' frozen sound speed', af)
# downstream velocity
w2 = Ucj*rho0/rho2
err = np.abs(w2-ae)/w2
print('downstream velocity',w2,' error ', err*100, ' %')
P2 = P0 + rho0*Ucj*Ucj*(1 - rho0/rho2)
err = np.abs(P2-p2)/p2
print('downstream pressure',P2,' error ',err*100,' %')
H2 = h0 + 0.5*Ucj*Ucj*(1 - (rho0/rho2)**2)
err = np.abs(H2-h2)/h2
print('downstream enthalpy',H2,' error ',100*err,' %')
