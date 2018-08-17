"""
Shock and Detonation Toolbox
"utilities" module

Utility tools for creating pre-defined plots and output files from the outputs
of the CJspeed, cvsolve, and zndsolve functions.
 
This module defines the following functions:

    CJspeed_plot
    cv_plot
    znd_plot
    znd_fileout
    
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

import numpy as np
import matplotlib.pyplot as plt
from cantera import one_atm

def CJspeed_plot(plot_data,cj_speed):
    """
    Creates two plots of the CJspeed fitting routine: both display density ratio
    vs. speed. The first is very 'zoomed in' around the minimum, and shows the
    quadratic fit plotted through the calculated points. The second shows the 
    same fit on a wider scale, with the minimum and its corresponding speed
    indicated.
    
    FUNCTION SYNTAX:
        CJspeed_plot(plot_data,cj_speed)
        
    INPUT:
        plot_data = tuple (rr,w1,dnew,a,b,c) produced by sdt.postshock.CJspeed
                    rr = density ratio
                    w1 = speed
                    dnew = minimum density
                    a,b,c = quadratic fit coefficients
                    
        cj_speed = CJ speed from same calculation as plot_data
        
    OUTPUT:
        (none, but displays plots)
        
    """
    # Unpack tuple of plot data
    rr,w1,dnew,a,b,c = plot_data
    
    # Generate plots
    xmin = np.min(rr)
    xmax = np.max(rr)
    x = np.linspace(xmin,xmax)
    y = a*x*x+b*x+c
    fig, ax = plt.subplots()
    line1, = ax.plot(rr,w1,marker='s',linestyle='none')
    line2, = ax.plot(x,y)
    ax.ticklabel_format(style='plain',useOffset=False)
    ax.set_title('CJspeed fitting routine output, CJ speed ='+str(cj_speed))
    ax.set_xlabel('density ratio')
    ax.set_ylabel('speed (m/s)')
    
    x = np.linspace(1.5,2.0)
    y = a*x*x+b*x+c
    fig, ay = plt.subplots()
    line1, = ay.plot(x,y)
    line1, = ay.plot(dnew,cj_speed,marker='s',linestyle='none')
    ay.set_title('CJspeed fitting routine output, CJ speed ='+str(cj_speed))
    ay.set_xlabel('density ratio')
    ay.set_ylabel('speed (m/s)')
    plt.show()
    
def cv_plot(cv_output):
    """
    Creates two subplots from the solution to a CV explosion (in sdt.cv.cvsolve):
    Temperature vs. time, and pressure vs. time.
    
    FUNCTION SYNTAX:
        cv_plot(cv_output)
        
    INPUT:
        cv_output: dictionary of outputs produced by sdt.cv.cvsolve.
        
    OUTPUT:
        (none, but displays plots)
        
    """
    ###########################################################
	 # PLOT TEMPERATURE PROFILE - MAX TIME = MAX TEMP + 10%
    ###########################################################    
    k = cv_output['T'].argmax()
    
    if cv_output['time'][k] == 0:
        maxt = cv_output['ind_time']*5
    elif cv_output['time'][k] >= cv_output['ind_time']*50:
        maxt = cv_output['ind_time']*5
    else:
        maxt = cv_output['time'][k] + 0.1*cv_output['time'][k]

    mint = 0
    	
    maxT = max(cv_output['T'])+0.1*min(cv_output['T'])
    minT = min(cv_output['T'])-0.1*min(cv_output['T'])
    maxP = max(cv_output['P'])+0.1*min(cv_output['P'])
    minP = min(cv_output['P'])-0.1*min(cv_output['P'])
    	
    fig = plt.figure()
    fig.suptitle('Constant volume explosion',fontsize=12)
	
    # Temperature as a function of time
    ax1 = fig.add_subplot(121)
    ax1.plot(cv_output['time'],cv_output['T'])
    ax1.set_xlabel('Time (s)',fontsize=12)
    ax1.set_ylabel('Temperature (K)',fontsize=12)
    ax1.set_xlim((mint,maxt))
    ax1.set_ylim((minT,maxT))
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    	
    # Pressure as a function of time
    ax2 = fig.add_subplot(122)
    ax2.plot(cv_output['time'],cv_output['P'])
    ax2.set_xlabel('Time (s)',fontsize=12)
    ax2.set_ylabel('Pressure (Pa)',fontsize=12)
    ax2.set_xlim((mint,maxt))
    ax2.set_ylim((minP,maxP))
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    plt.tight_layout()
    plt.show()
    
    print('Max T = '+str(max(cv_output['T']))+' K, Max P = '+str(max(cv_output['P'])/one_atm)
            +' atm, Induction time = '+str(cv_output['ind_time'])+' s')


    
def znd_plot(znd_output):
    """
    Creates two subplots from the solution to a ZND detonation (in sdt.znd.zndsolve):
    Temperature vs. distance, and pressure vs. distance.
    
    FUNCTION SYNTAX:
        znd_plot(znd_output)
        
    INPUT:
        znd_output: dictionary of outputs produced by sdt.znd.zndsolve.
        
    OUTPUT:
        (none, but displays plots)
        
    """
    k = znd_output['T'].argmax()
    
    if znd_output['time'][k] == 0:
        maxx = znd_output['ind_len_ZND']*5
    else:
        maxx = znd_output['distance'][k]

    minx = 0
    	
    maxT = max(znd_output['T'])+0.1*min(znd_output['T'])
    minT = min(znd_output['T'])-0.1*min(znd_output['T'])
    maxP = max(znd_output['P'])+0.1*min(znd_output['P'])
    minP = min(znd_output['P'])-0.1*min(znd_output['P'])
    	
    fig = plt.figure()
    fig.suptitle('ZND structure',fontsize=12)
	
    # Temperature as a function of position
    print('Delta = %.3g m' % znd_output['ind_len_ZND'])
    print('Final T = %.5g K; Max T = %.5g K' % (znd_output['T'][-1],max(znd_output['T'])))
    ax1 = fig.add_subplot(121) #subplot(1,2,1)
    ax1.plot(znd_output['distance'],znd_output['T'])
    ax1.set_xlabel('Distance (m)',fontsize=12)
    ax1.set_ylabel('Temperature (K)',fontsize=12)
    ax1.set_xlim((minx,maxx))
    ax1.set_ylim((minT,maxT))
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    	
    # Pressure as a function of position
    print('Final T = %.3g K; Max T = %.3g K' % (znd_output['P'][-1],max(znd_output['P'])))
    ax2 = fig.add_subplot(122) #subplot(1,2,1)
    ax2.plot(znd_output['distance'],znd_output['P'])
    ax2.set_xlabel('Distance (s)',fontsize=12)
    ax2.set_ylabel('Pressure (Pa)',fontsize=12)
    ax2.set_xlim((minx,maxx))
    ax2.set_ylim((minP,maxP))
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    
    plt.tight_layout()
    plt.show()
    
def znd_fileout(fname,znd_output):
    """
    Creates 2 formatted text files to store the output of the solution to a ZND
    detonation (from sdt.znd.zndsolve). Includes a timestamp of when the file was created,
    input conditions, and tab-separated columns of output data.
    
    FUNCTION SYNTAX:
        znd_fileout(fname,znd_output)
        
    INPUT:
        fname: filename (appended by detonation velocity and '_znd' or '_znd2') (str)
        znd_output: dictionary of outputs produced by sdt.znd.zndsolve.
        
    OUTPUT:
        (none, but generates files)
    """
    import datetime
    from cantera import one_atm

    P1 = znd_output['gas1'].P
    T1 = znd_output['gas1'].T
    r1 = znd_output['gas1'].density
    
    U1 = znd_output['U1']

    fid = open(fname+'_'+str(U1)+'_znd.txt','w')
    d = datetime.date.today().strftime("%B %d, %Y")

    P = P1/one_atm

    fid.write('# ZND: DETONATION STRUCTURE CALCULATION\n')
    fid.write('# CALCULATION RUN ON %s\n\n' % d)
	
    fid.write('# INITIAL CONDITIONS\n')
    fid.write('# TEMPERATURE (K) %4.1f\n' % T1)
    fid.write('# PRESSURE (atm) %2.1f\n' % P)
    fid.write('# DENSITY (kg/m^3) %1.4e\n' % r1)
	
    fid.write('# SHOCK SPEED (m/s) %5.2f\n\n' % U1)
	
    fid.write('# Induction Times\n')
    fid.write('# Time to Peak Thermicity =   %1.4e s\n' % znd_output['ind_time_ZND'])
    fid.write('# Distance to Peak Thermicity =   %1.4e m\n' % znd_output['ind_len_ZND'])

    fid.write('\n# Exothermic/Reaction Times\n')
    fid.write('# Exothermic Pulse Time =   %1.4e s\n' % znd_output['exo_time_ZND'])
    fid.write('# Exothermic Pulse Distance =   %1.4e m\n' % znd_output['exo_len_ZND'])
	
    fid.write('# REACTION ZONE STRUCTURE\n\n')
	
    fid.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    fid.write('Variables = "Distance (m)", "Mach Number", "Time (s)", "Pressure (Pa)", "Temperature (K)", "Density (kg/m^3)", "Thermicity (1/s)"\n')
	
   
    for val in zip(znd_output['distance'],znd_output['M'],znd_output['time'],znd_output['P'],
                   znd_output['T'],znd_output['rho'],znd_output['thermicity']):
        fid.write('%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e\n' % val)
	
    fid.close()
	
    fid = open(fname+'_'+str(U1)+'_znd2.txt','w')

    fid.write('# ZND: DETONATION STRUCTURE CALCULATION\n')
    fid.write('# CALCULATION RUN ON %s\n\n' % d)

    fid.write('# INITIAL CONDITIONS\n')
    fid.write('# TEMPERATURE (K) %4.1f\n' % T1)
    fid.write('# PRESSURE (atm) %2.1f\n' % P)
    fid.write('# DENSITY (kg/m^3) %1.4e\n' % r1)
	
    fid.write('# SHOCK SPEED (M/S) %5.2f\n\n' % U1)
    fid.write('# REACTION ZONE STRUCTURE\n\n')

    fid.write('# Induction Times\n')
    fid.write('# Time to Peak Thermicity =   %1.4e s\n' % znd_output['ind_time_ZND'])
    fid.write('# Distance to Peak Thermicity =   %1.4e m\n' % znd_output['ind_len_ZND'])

    fid.write('\n# Exothermic/Reaction Times\n')
    fid.write('# Time of Reaction based on Thermicity Pulse Width =   %1.4e s\n' % znd_output['exo_time_ZND'])
    fid.write('# Length of Reaction based on Thermicity Pulse Width =   %1.4e m\n' % znd_output['exo_len_ZND'])
	
    fid.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    fid.write('Variables = "Distance (m)", "Velocity (m/s)", "Sound Speed (m/s)", "Gamma", "Weight (kg/mol)","c^2-U^2 (m/s)"\n')
	
    
    for val in zip(znd_output['distance'],znd_output['U'],znd_output['af'],znd_output['g'],
                   znd_output['wt'],znd_output['sonic']):
        fid.write('%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e\n'  % val)
	
    fid.close()
    