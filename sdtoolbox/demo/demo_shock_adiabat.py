"""
Shock and Detonation Toolbox Demo Program

Generates the points on a frozen shock adiabat and creates an output file.

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
from sdtoolbox.thermo import soundspeed_fr
from sdtoolbox.postshock import PostShock_fr
from numpy import linspace
import datetime

print('demo_shock_adiabat')
# set initial state of gas, avoid using precise stoichiometric ratios. 
P1 = 100000; T1 = 300; 
q = 'O2:1 N2:3.76'  
mech = 'gri30_highT.cti' 
gas1 = ct.Solution(mech)
gas1.TPX = T1,P1,q

print('Computing frozen postshock state for '+q+' using '+mech)
Umin = soundspeed_fr(gas1) + 50

# Evaluate the frozen state of the gas behind a shock wave traveling at speeds up to a maximum value. 

i = 1;
Umax = 2000.; Ustep = 100.; # Ustep will be approximated
nsteps = int((Umax-Umin)/Ustep)

speed = []; P = []; R = []; T = []; a = []; gamma = []; w2 = []; u2 = [];
for speeds in  linspace(Umin,Umax,num=nsteps):
    speed.append(speeds)
    gas = PostShock_fr(speed[-1],P1, T1, q, mech)
    P.append(gas.P)
    R.append(gas.density)
    T.append(gas.T)
    a.append(soundspeed_fr(gas))
    gamma.append(a[-1]**2*R[-1]/P[-1])
    w2.append(gas1.density*speed[-1]/R[-1])
    u2.append(speed[-1] - w2[-1])



##################################################################################################
# Create output file in text format with Tecplot-style variable labels for
# columns
##################################################################################################
fn = 'shock_adiabat.txt'
d = datetime.date.today()
fid = open(fn, 'w')
fid.write('# Shock adiabat\n')
fid.write('# Calculation run on %s\n' % d )
fid.write('# Initial conditions\n')
fid.write('# Temperature (K) %4.1f\n' % T1)
fid.write('# Pressure (Pa) %2.1f\n' % P1)
fid.write('# Density (kg/m^3) %1.4e\n' % gas1.density)
fid.write('# Initial species mole fractions: '+q+'\n')
fid.write('# Reaction mechanism: '+mech+'\n\n')
fid.write('Variables = "shock speed (m/s)", "density (kg/m3)", "Temperature (K)", "Pressure (atm)", "sound speed - fr (m/s)", "particle speed (m/s)", "shock-fixed speed (m/s)", "gamma (fr)"\n')

for i in range(nsteps):
    fid.write('%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e\n' % (speed[i], R[i], T[i], P[i], a[i], u2[i], w2[i], gamma[i]))

fid.close()
