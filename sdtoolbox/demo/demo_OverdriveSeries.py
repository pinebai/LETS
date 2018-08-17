"""
Shock and Detonation Toolbox Demo Program

This is a demostration of how to vary the Overdrive (U/UCJ)
in a loop for constant volume explosions and the ZND detonation simulations.

NOTE: The convention in this demo is overdrive = shock_speed/CJspeed,
       Another convention is (U/UCJ)^2, so be careful when using this demo

Here the explosion functions are called for varying conditions, new properties are calculated, 
plots are displayed, and output files are written.
 
WARNING: due to the inefficiency of the Scipy ODE solvers as implemented here,
         cvsolve and zndsolve are much slower in the Python version of SDToolbox
         as compared to MATLAB. As a result, this demo takes a LONG TIME to complete
         due to looping these functions. Some of the outputs, particularly the 
         activation energy estimates, may also differ from MATLAB. See the docstrings
         for the cv and znd modules for further discussion.

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

from sdtoolbox.postshock import CJspeed, PostShock_fr, PostShock_eq
from sdtoolbox.cv import cvsolve
from sdtoolbox.znd import zndsolve
import cantera as ct
import numpy as np
import datetime
import matplotlib.pyplot as plt
# Use warning filters on cvsolve/zndsolve to suppress output of 'excess work' warnings
# if taking too long to solve
import warnings

print('Overdrive (U/UCJ) Series')

mech = 'h2air_highT.cti';   
gas  = ct.Solution(mech); 
gas1 = ct.Solution(mech);
nsp  = gas.n_species;

npoints = 10

P1          = np.zeros(npoints,float)
overdrive   = np.zeros(npoints,float)
exo_time_CV = np.zeros(npoints,float)
ind_time_CV = np.zeros(npoints,float)
theta_effective_CV = np.zeros(npoints,float)
T2   = np.zeros(npoints,float)
P2   = np.zeros(npoints,float)
rho2 = np.zeros(npoints,float)
Ts   = np.zeros(npoints,float)
Ps   = np.zeros(npoints,float)
exo_time_ZND    = np.zeros(npoints,float)
exo_len_ZND     = np.zeros(npoints,float)
ind_time_ZND    = np.zeros(npoints,float)
ind_len_ZND     = np.zeros(npoints,float)
Tf_ZND          = np.zeros(npoints,float)
theta_effective_ZND = np.zeros(npoints,float)

# find Hydrogen nitrogen, and oxygen indices
ih2  = gas.species_index('H2')
io2  = gas.species_index('O2')
in2  = gas.species_index('N2')
x = 'H2:2.0,O2:1.0,N2:3.76'  
T1 = 300
P1 = ct.one_atm; P1atm = P1/ct.one_atm

print('Initial Conditions')
print(x + ', Pressure = %.2f atm, Temperature = %.2f K' % (P1atm,T1))
print('For %s values of overdrive' % (npoints))
cj_speed = CJspeed(P1, T1, x, mech)


for i in range(npoints):
    overdrive[i] = (1.0 +0.6/npoints*(i)) #Overdrive = U/Ucj
    ratio = overdrive[i]
    print('%i : Overdrive = %.2f ' % (i+1,ratio))

    # Find post shock state for given speed
    gas.TPX = T1,P1,x
    gas = PostShock_fr(cj_speed*overdrive[i], P1, T1, x, mech)
    Ts[i] = gas.T # frozen shock temperature   
    Ps[i] = gas.P # frozen shock pressure
    
    ### Constant Volume Explosion Data ###
    # Solve constant volume explosion ODEs
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        CVout = cvsolve(gas,dt=1e-8)
    exo_time_CV[i] = CVout['exo_time']
    ind_time_CV[i] = CVout['ind_time']    

    
    ### ZND Detonation Data ###
    gas.TPX = T1,P1,x
    gas = PostShock_fr(cj_speed*overdrive[i], P1, T1, x, mech)
    # Solve znd detonation ODEs
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        ZNDout = zndsolve(gas,gas1,cj_speed*overdrive[i])
    exo_time_ZND[i] = ZNDout['exo_time_ZND']
    exo_len_ZND[i] = ZNDout['exo_len_ZND']
    ind_time_ZND[i] = ZNDout['ind_time_ZND']
    ind_len_ZND[i] = ZNDout['ind_len_ZND']
    Tf_ZND[i] = ZNDout['T'][-1]
    
    ### Calculate CJ state properties ###
    gas = PostShock_eq(cj_speed*overdrive[i], P1, T1, x, mech);
    T2[i] = gas.T
    P2[i] = gas.P
    rho2[i] = gas.density
   
    # Approximate the effective activation energy using finite differences
    factor = 0.02 
    Ta = Ts[i]*(1.0+factor)
    gas.TPX = Ta,Ps[i],x
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        CVout1 = cvsolve(gas,dt=1e-8)
    taua = CVout1['ind_time']
    Tb = Ts[i]*(1.0-factor)
    gas.TPX = Tb,Ps[i],x
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        CVout2 = cvsolve(gas,dt=1e-8)
    taub = CVout2['ind_time']
    # Approximate effective activation energy for CV explosion
    if taua == 0 and taub == 0:
        theta_effective_CV[i] = 0
    else:
        theta_effective_CV[i] = 1.0/Ts[i]*((np.log(taua)-np.log(taub))/((1.0/Ta)-(1.0/Tb)))
    
    gas = PostShock_fr(cj_speed*overdrive[i]*(1.0+factor), P1, T1, x, mech)
    Ta= gas.T
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        ZNDout1 = zndsolve(gas,gas1,cj_speed*overdrive[i]*(1.0+factor))
    ind_len_a = ZNDout1['ind_len_ZND']
    gas = PostShock_fr(cj_speed*overdrive[i]*(1.0-factor),P1, T1, x, mech)
    Tb = gas.T
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=UserWarning)
        ZNDout2 = zndsolve(gas,gas1,cj_speed*overdrive[i]*(1.0-factor))   
    ind_len_b = ZNDout2['ind_len_ZND']
    # Approximate effective activation energy for ZND Detonation
    if ind_len_a == 0 and ind_len_b ==0:
        theta_effective_ZND[i] = 0        
    else:
        theta_effective_ZND[i] = 1/Ts[i]*((np.log(ind_len_a)-np.log(ind_len_b))/((1/Ta)-(1/Tb)))
    
## Make plots
fontsize=12
plt.figure(1)

plt.subplot(221)
plt.plot(overdrive,T2,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel('Temperature (K)',fontsize=fontsize)
plt.title('Post CJ State Temperature',fontsize=fontsize)
plt.xlim(xmin=1)

plt.subplot(222)
plt.plot(overdrive,P2,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel('Pressure (Pa)',fontsize=fontsize)
plt.title('Post CJ State Pressure',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.xlim(xmin=1)

plt.subplot(223);
plt.plot(overdrive,rho2,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'Density (kg/m$^3$)',fontsize=fontsize)
plt.title('Post CJ State Density',fontsize=fontsize)
plt.xlim(xmin=1)

plt.tight_layout()

# Plots for the Induction Zone (times and lengths)

plt.figure(2)
plt.subplot(221)
plt.plot(overdrive,ind_time_CV,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\tau_{CV_i}$ (s)',fontsize=fontsize)
plt.title('Induction time for CV explosion',fontsize=fontsize)
plt.xlim(xmin=1)

plt.subplot(222)
plt.plot(overdrive,ind_time_ZND,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\tau_{ZND_i}$ (s)',fontsize=fontsize)
plt.title('Induction time for ZND detonation',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.xlim(xmin=1)

plt.subplot(224)
plt.plot(overdrive,ind_len_ZND,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\Delta_{ZND_i}$ (m)',fontsize=fontsize)
plt.title('Induction length for ZND detonation',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.xlim(xmin=1)

plt.tight_layout()

# Plots for the Exothermic Zone (times and lengths)
plt.figure(3)
plt.subplot(221)
plt.plot(overdrive,exo_time_CV,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\tau_{CV_e}$ (s)',fontsize=fontsize)
plt.title('Exothermic time for CV explosion',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

plt.subplot(222)
plt.plot(overdrive,exo_time_ZND,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\tau_{ZND_e}$ (s)',fontsize=fontsize)
plt.title('Exothermic time for ZND detonation',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

plt.subplot(224)
plt.plot(overdrive,exo_len_ZND,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\Delta_{ZND_e}$ (m)',fontsize=fontsize)
plt.title('Exothermic length for ZND detonation',fontsize=fontsize)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

plt.tight_layout()

# Ts and Tf for ZND detonation
plt.figure(4)
plt.subplot(121)
plt.plot(overdrive,Ts,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$T_s$ (K)',fontsize=fontsize)
plt.title('Frozen shock Temperature for ZND detonation',fontsize=fontsize)

plt.subplot(122)
plt.plot(overdrive,Tf_ZND,'k')
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$T_f$ (K)',fontsize=fontsize)
plt.title('Final Temperature for ZND detonation',fontsize=fontsize)

plt.tight_layout()

# Approximation of the effective activation energy for CV explosion
plt.figure(5)
plt.plot(overdrive,theta_effective_CV)
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\Theta_{CV}$ (J)',fontsize=fontsize);
plt.title('Effective activation energy for CV explosion',fontsize=fontsize);
plt.tight_layout()

# Approximation of the effective activation energy for ZND detonation
plt.figure(6)
plt.plot(overdrive, theta_effective_ZND)
plt.xlabel(r'Overdrive $= U/U_{CJ}$',fontsize=fontsize)
plt.ylabel(r'$\Theta_{ZND}$ (J)',fontsize=fontsize)
plt.title('Effective activation energy for ZND detonation',fontsize=fontsize)
plt.tight_layout()

plt.show()

##################################################################################################
# CREATE OUTPUT TEXT FILE
##################################################################################################
d = datetime.date.today()
fn = 'Overdrive_Series.txt'
outfile = open(fn, 'w')
outfile.write('# CONSTANT VOLUME EXPLOSION\n')
outfile.write('# CALCULATION RUN ON %s\n\n' % d)
outfile.write('# INITIAL CONDITIONS\n')
outfile.write('# TEMPERATURE (K) %0.2f\n' % T1)
outfile.write('# PRESSURE (Pa) %0.2f\n' % P1)
q = 'H2:2 O2:1 N2:3.76'
outfile.write( '# SPECIES MOLE FRACTIONS: %s \n' % q )
outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
outfile.write('Variables = "Overdrive(U/UCJ)", "Temperature state 2", "Pressure state 2", "density state 2", "Temperature Post Shock", "Pressure Post Shock", "Induction time CV", "Exothermic time CV",  "Effective Activation Energy CV "  "Effective Activation Energy ZND ", "Induction  time ZND", "Induction length ZND", "Exothermic time ZND", "Exothermic length ZND", "Final Temperature ZND"\n')
for i in range(npoints):
     outfile.write('%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e \t %14.5e  \t %14.5e \n' 
                   % (overdrive[i], T2[i], P2[i], rho2[i], Ts[i], Ps[i], ind_time_CV[i], exo_time_CV[i], theta_effective_CV[i], theta_effective_ZND[i], 
                      ind_time_ZND[i], ind_len_ZND[i], exo_time_ZND[i], exo_len_ZND[i], Tf_ZND[i]))
outfile.close()

