"""
Shock and Detonation Toolbox Demo Program

Calculates post-relected-shock state for a specified shock speed and a specified 
initial mixture.  In this demo, both shocks are non-reactive, i.e. the computed states 
behind both the incident and reflected shocks are FROZEN states.  
The reflected state is computed by the SDToolbox function reflected_fr.

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
from sdtoolbox.postshock import PostShock_fr
from sdtoolbox.reflections import reflected_fr
from sdtoolbox.thermo import soundspeed_fr

# Initial state specification:
# P1 = Initial Pressure  
# T1 = Initial Temperature 
# U = Shock Speed 
# q = Initial Composition 
# mech = Cantera mechanism File name
P1 = 100000; T1 = 300; P1atm = P1/ct.one_atm;
q = 'H2:2 O2:1 N2:3.76';    
mech = 'h2air_highT.cti';               
gas1 = ct.Solution(mech)
gas1.TPX = T1, P1, q

print ('Initial state: ' + q + ', P1 = %.2f atm,  T1 = %.2f K' % (P1atm,T1) )
print ('Mechanism: ' + mech)

# create gas objects for other states
gas2 = ct.Solution(mech)
gas3 = ct.Solution(mech)

# compute minimum incident wave speed
a_fr = soundspeed_fr(gas1)
# incident wave must be greater than or equal to frozen sound speed for
# frozen shock wave computations
UI = 3.*a_fr

print('Incident shock speed UI = %.2f m/s' % (UI))

# compute postshock gas state object gas2
gas2 = PostShock_fr(UI, P1, T1, q, mech);
P2 = gas2.P/ct.one_atm;

print ('Frozen Post-Incident-Shock State')
print ('T2 = %.2f K, P2 = %.2f atm' % (gas2.T,P2))

# compute reflected shock post-shock state gas3
[p3,UR,gas3]= reflected_fr(gas1,gas2,gas3,UI);
# Outputs:
# p3 - pressure behind reflected wave
# UR = Reflected shock speed relative to reflecting surface
# gas3 = gas object with properties of postshock state

P3 = gas3.P/ct.one_atm
print ('Frozen Post-Reflected-Shock State')
print ('T3 = %.2f K,  P3 = %.2f atm' % (gas3.T,P3))
print ("Reflected Wave Speed = %.2f m/s" % (UR))

# gas states
print('Incident gas state')
gas1()
print('Post-incident-shock gas state')
gas2()
print('Post-reflected-schock gas state')
gas3()