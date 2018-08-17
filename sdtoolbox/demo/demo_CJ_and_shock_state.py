"""
Shock and Detonation Toolbox Demo Program

Calculates 2 reflection conditions
  1) equilibrium post-initial-shock state behind a shock traveling at CJ speed (CJ state) 
     followed by equilibrium post-reflected-shock state
  2) frozen post-initial-shock state behind a CJ wave 
     followed by frozen post-reflected-shock state
 
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
from sdtoolbox.postshock import CJspeed, PostShock_eq, PostShock_fr
from sdtoolbox.reflections import reflected_eq, reflected_fr
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr

# Initial state specification:
# P1 = Initial Pressure  
# T1 = Initial Temperature 
# U = Shock Speed 
# q = Initial Composition 
# mech = Cantera mechanism file name

# stoichiometric H2-air detonation at nominal atmospheric conditions
P1 = 100000.
T1 = 300
q ='H2:2 O2:1 N2:3.76'
mech = 'gri30_highT.cti'
gas_initial = ct.Solution(mech)
gas_initial.TPX = T1, P1, q
rho_1 = gas_initial.density

# compute CJ speed
cj_speed = CJspeed(P1, T1, q, mech)  

# compute equilibrium CJ state parameters
gas_cj = PostShock_eq(cj_speed, P1, T1, q, mech)
ae = soundspeed_eq(gas_cj)
af = soundspeed_fr(gas_cj)
rho_2 = gas_cj.density
gammae = ae**2*rho_2/gas_cj.P
gammaf = af**2*rho_2/gas_cj.P
w2 = cj_speed*rho_1/rho_2
u2 = cj_speed-w2
print ('CJ computation for ' + mech + ' with composition ' + q )
print ('Initial conditions: P1 = %.3e Pa & T1 = %.2f K'  % (P1,T1)  )
print ('CJ Speed   %.1f m/s' % cj_speed)
print ('CJ State')
print ('   Pressure   %.3e Pa' % gas_cj.P)
print ('   Temperature  %.1f K' % gas_cj.T)
print ('   Density  %.3f kg/m3' % gas_cj.density)
print ('   Entropy  %.3e J/K' % gas_cj.entropy_mass)
print ('   w2 (wave frame) %.1f m/s' % w2)
print ('   u2 (lab frame) %.1f m/s' % u2)
print ('   c2 frozen %.1f m/s' % af)
print ('   c2 equilbrium %.1f m/s' % ae)
print ('   gamma2 frozen %.3f ' % gammaf)
print ('   gamma2 equilbrium %.3f ' % gammae)

# compute equilibrium reflected CJ parameters
gas_reflected_eq = ct.Solution(mech)
[p3,UR,gas_reflected_eq] = reflected_eq(gas_initial,gas_cj,gas_reflected_eq,cj_speed)
print ('Reflected CJ shock (equilibrium) computation ')
print ('   Reflected wave speed %.1f (m/s)' % UR)
print ('   Pressure %.3e (Pa)'% gas_reflected_eq.P)
print ('   Temperature %.1f (K)' % gas_reflected_eq.T)
print ('   Density  %.3f (kg/m3)' % gas_reflected_eq.density)

# compute frozen CJ shock state parameters
gas_vn = PostShock_fr(cj_speed, P1, T1, q, mech)
ae = soundspeed_eq(gas_vn)
af = soundspeed_fr(gas_vn)
rho_2 = gas_vn.density
gammae = ae**2*rho_2/gas_vn.P
gammaf = af**2*rho_2/gas_vn.P
w2 = cj_speed*rho_1/rho_2
u2 = cj_speed-w2
print ('VN (frozen postshock) State')
print ('   Pressure   %.3e Pa' % gas_vn.P)
print ('   Temperature  %.1f K' % gas_vn.T)
print ('   Density  %.3f kg/m3' % gas_vn.density)
print ('   Entropy  %.3e J/K' % gas_vn.entropy_mass)
print ('   w2 (wave frame) %.1f m/s' % w2)
print ('   u2 (lab frame) %.1f m/s' % u2)
print ('   c2 frozen %.1f m/s' % af)
print ('   c2 equilbrium %.1f m/s' % ae)
print ('   gamma2 frozen %.3f ' % gammaf)
print ('   gamma2 equilbrium %.3f ' % gammae)

# compute frozen reflected CJ parameters
gas_reflected_fr = ct.Solution(mech)
[p3,UR,gas_reflected_fr] = reflected_fr(gas_initial,gas_vn,gas_reflected_fr,cj_speed)
print ('Reflected CJ shock (frozen) computation ')
print ('   Reflected wave speed %.1f (m/s)' % UR)
print ('   Pressure %.3e (Pa)'% gas_reflected_fr.P)
print ('   Temperature %.1f (K)' % gas_reflected_fr.T)
print ('   Density  %.3f (kg/m3)' % gas_reflected_fr.density)
