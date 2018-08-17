"""
Shock and Detonation Toolbox Demo Program

Generate plots and output files for a ZND detonation with the shock front
traveling at the CJ speed.

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

from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve
from sdtoolbox.utilities import CJspeed_plot, znd_plot, znd_fileout
import cantera as ct
import warnings

P1 = 100000 
T1 = 300
q = 'H2:2 O2:1 N2:3.76'
mech = 'h2air_highT.cti'
file_name = 'h2air'

# Find CJ speed and related data, make CJ diagnostic plots
cj_speed,R2,plot_data = CJspeed(P1,T1,q,mech,fullOutput=True)
CJspeed_plot(plot_data,cj_speed)

# Set up gas object
gas1 = ct.Solution(mech)
gas1.TPX = T1,P1,q

# Find post shock state for given speed
gas = PostShock_fr(cj_speed, P1, T1, q, mech)

# Solve ZND ODEs, make ZND plots
with warnings.catch_warnings():
    warnings.simplefilter('ignore',category=UserWarning)
    znd_out = zndsolve(gas,gas1,cj_speed,dt=1e-9,absTol=1e-7,t_end=1e-5)
znd_plot(znd_out)
znd_fileout(file_name,znd_out)

print('Reaction zone pulse width (exothermic length) = %.4g m' % znd_out['exo_len_ZND'])
print('Reaction zone induction length = %.4g m' % znd_out['ind_len_ZND'])
print('Reaction zone pulse time (exothermic time) = %.4g s' % znd_out['exo_time_ZND'])
print('Reaction zone induction time = %.4g s' % znd_out['ind_time_ZND'])