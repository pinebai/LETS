"""
Shock and Detonation Toolbox
"cv" module

Calculates constant-volume explosions.
 
This module defines the following functions:

    cvsys
    cvsolve
    
#### WARNING ####
MATLAB's proprietary solvers (e.g. ode15s) are superior to the solvers in Python's
scipy.integrate module. The latter have difficulty resolving the rapidly-changing
regions in the CV ODE system. As a result, the version of cvsolve in the MATLAB
toolbox is much faster than the version here. This version will produce the correct
output, including the stiff regions of the solution, but it requires about 2 orders
of magnitude more timesteps than MATLAB. It is also prone to failing to converge
if the timestep, end-time, or tolerances are not set correctly, so some manual
adjustment can be required, whereas MATLAB usually converges rapidly with the
default settings.  


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
import numpy as np
from scipy.integrate import ode

def cvsys(t,y,gas):
    """
    Evaluates the system of ordinary differential equations for an adiabatic, 
    constant-volume, zero-dimensional reactor. 
    It assumes that the 'gas' object represents a reacting ideal gas mixture.

    INPUT:
        t = time
        y = solution array [temperature, species mass 1, 2, ...]
        gas = working gas object
    
    OUTPUT:
        An array containing time derivatives of:
            temperature and species mass fractions, 
        formatted in a way that the integrator in cvsolve can recognize.
        
    """
    # Set the state of the gas, based on the current solution vector.
    gas.TDY = y[0],gas.density,y[1:]
    
    # energy/temperature equation (terms separated for clarity)  
    a = gas.standard_enthalpies_RT - np.ones(gas.n_species)
    b = gas.net_production_rates /(gas.density * gas.cv_mass)
    dTdt = -gas.T * ct.gas_constant * np.dot(a,b)
    
    # species equations
    dYdt = gas.net_production_rates*gas.molecular_weights/gas.density
    
    return np.hstack((dTdt, dYdt))
    
   
def cvsolve(gas,t_end=1e-6,dt=1e-9,backend='lsoda',relTol=1e-5,absTol=1e-12):
    """
    Solves the ODE system defined in cvsys, taking the gas object input as the
    initial state. As discussed in the warning for this module, the integrator
    has difficulty with the stiff region of this system. In order to resolve it
    without using an extremely small timestep for the whole domain, the system
    is first solved with the dt specified by the user. 
    
    The region of rapid temperature rise is then detected, and the system restarted
    with finer resolution and solved just across this region. This process is repeated
    once more to achieve adequate resolution for correct calculation of quantities
    such as the exothermic pulse width  and induction time. This is quite a crude
    and inefficient process compared to the equivalent function in MATLAB, and
    can probably be improved on considerably.
    
    FUNCTION SYNTAX:
        output = cvsolve(gas,**kwargs)
    
    INPUT:
        gas = working gas object
        
    OPTIONAL INPUT:
        t_end = end time for integration, in sec
        dt = time step for integration, in sec
        backend = solver backend used for the integrator function
        relTol = relative tolerance
        absTol = absolute tolerance
    
    OUTPUT:
        output = a dictionary containing the following results:
            exo_time = pulse width (in secs) of temperature gradient (using 1/2 max)
            ind_time = time to maximum temperature gradient
            ind_len = distance to maximum temperature gradient
            ind_time_10 = time to 10% of maximum temperature gradient
            ind_time_90 = time to 90% of maximum temperature gradient
            time = array of time
            T = temperature profile array
            P = pressure profile array
            species = species profiles (list of lists)
            
    """
    y0 = np.hstack((gas.T,gas.Y))
    r0 = gas.density
    # Set the integration parameters for the ODE solver
    tel = [0.,t_end] # Timespan
    
    # Initialize output arrays. Use a dictionary to emulate the data structure
    # used in the corresponding MATLAB code.
    output = {}
    
    time = [tel[0]]
    T = [gas.T]
    species = [gas.Y] # a list of numpy arrays
 
    # Set up ODE system for integration
    solver = ode(cvsys)
    solver.set_integrator(backend, method='bdf',rtol=relTol, atol=absTol)   
    solver.set_initial_value(y0,tel[0]).set_f_params(gas)
    
    # Call the integrator to march in time - the equation set is defined in cvsys() 
    while solver.successful() and solver.t < tel[1]:
        solver.integrate(solver.t + dt)
        
        #############################################################################################
        # Extract TIME, TEMPERATURE and MASS FRACTION arrays from integrator output
        #############################################################################################
        time.append(solver.t)
        T.append(solver.y[0])
        species.append(solver.y[1:])
    del solver
    
    #### RESOLUTION OF STIFF REGION ####  
    riseidx = np.ediff1d(T).argmax() # beginning of 'gap' in temperature rise
    riset = time[riseidx+1]-time[riseidx]
    steps = 1e4
    stiffdt = riset/steps
    
    y0 = np.hstack((T[riseidx],species[riseidx]))
    
    solver = ode(cvsys)
    solver.set_integrator(backend, method='bdf',rtol=relTol, atol=absTol)   
    solver.set_initial_value(y0,time[riseidx]).set_f_params(gas)
    
    stifftime = [time[riseidx]]; stiffT = [T[riseidx]]; stiffspecies = [species[riseidx]]
    
    while solver.successful() and solver.t < time[riseidx+1]:
        solver.integrate(solver.t + stiffdt)
        stifftime.append(solver.t)
        stiffT.append(solver.y[0])
        stiffspecies.append(solver.y[1:])
    del solver
    
    
    # Repeat process (resolve stiffer region within stiff region)
    riseidx2 = np.ediff1d(stiffT).argmax()
    riset2 = stifftime[riseidx2+1]-stifftime[riseidx2]
    steps2 = 1e4
    stiffdt2 = riset2/steps2
    
    y0 = np.hstack((stiffT[riseidx2],stiffspecies[riseidx2]))
    
    solver = ode(cvsys)
    solver.set_integrator(backend, method='bdf',rtol=relTol, atol=absTol)   
    solver.set_initial_value(y0,stifftime[riseidx2]).set_f_params(gas)
    
    stifftime2 = [stifftime[riseidx2]]; stiffT2 = [stiffT[riseidx2]]; stiffspecies2 = [stiffspecies[riseidx2]]
    
    while solver.successful() and solver.t < stifftime[riseidx2+1]:
        solver.integrate(solver.t + stiffdt2)
        stifftime2.append(solver.t)
        stiffT2.append(solver.y[0])
        stiffspecies2.append(solver.y[1:])
    del solver
        
     # Splice stiff regions back into main output arrays
    stifftime[riseidx2+1:riseidx2+1] = stifftime2
    stiffT[riseidx2+1:riseidx2+1] = stiffT2
    stiffspecies[riseidx2+1:riseidx2+1] = stiffspecies2 
        
    time[riseidx+1:riseidx+1] = stifftime
    T[riseidx+1:riseidx+1] = stiffT
    species[riseidx+1:riseidx+1] = stiffspecies
    ####################################
    
    
    output['time'] = np.asarray(time)
    output['T'] = np.asarray(T)
    output['species'] = species
    
    # Initialize additional output matrices where needed
    b = len(time)
    
    output['ind_time'] = 0
    output['ind_time_90'] = 0
    output['ind_time_10'] = 0
    output['exo_time'] = 0
    
    output['P'] = np.zeros(b)
    
    temp_grad = np.zeros(b)
    #############################################################################
    # Extract PRESSSURE and TEMPERATURE GRADIENT
    #############################################################################
    
    # Have to loop for operations involving the working gas object
    for i,T in enumerate(output['T']):
        gas.TDY = T,r0,output['species'][i]
        
        wt = gas.mean_molecular_weight        
        s = 0
        for z in range(gas.n_species):
            w = gas.molecular_weights[z]
            e = ct.gas_constant*T*(gas.standard_enthalpies_RT[z]/w - 1/wt)
            s = s + e*w*gas.net_production_rates[z]
            
        temp_grad[i] = -s/(r0*gas.cv_mass)
        output['P'][i] = gas.P

    n = temp_grad.argmax()
    
    if n == b:
        print('Error: Maximum temperature gradient occurs at the end of the reaction zone')
        print('       Your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong')
        print(' ')
        output['ind_time'] = output.time[b] 
        output['ind_time_10'] = output.time[b]
        output['ind_time_90'] = output.time[b] 
        output['exo_time'] = 0
        print('Induction Time: '+str(output['ind_time']))
        print('Exothermic Pulse Time: '+str(output['exo_time']))
        return output
    elif n == 0:
        print('Error: Maximum temperature gradient occurs at the beginning of the reaction zone')
        print('       Your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong')
        print(' ')
        output['ind_time'] = output.time[0] 
        output['ind_time_10'] = output.time[0]
        output['ind_time_90'] = output.time[0] 
        output['exo_time'] = 0 
        print('Induction Time: '+str(output['ind_time']))
        print('Exothermic Pulse Time: '+str(output['exo_time']))
        return output
    else:
        output['ind_time'] = output['time'][n]
        
        k = 0
        MAX10 = 0.1*max(temp_grad)
        d = temp_grad[0]        
        while d < MAX10:
            k = k + 1
            d = temp_grad[k]
        output['ind_time_10'] = output['time'][k]
        
        k = 0
        MAX90 = 0.9*max(temp_grad)
        d = temp_grad[0]
        while d < MAX90:
            k = k + 1
            d = temp_grad[k]
        output['ind_time_90'] = output['time'][k]

        # find exothermic time
        half_T_flag1 = 0
        half_T_flag2 = 0
        # Go into a loop to find two times when temperature is half its maximum
        for j,tgrad in enumerate(list(temp_grad)):
            if half_T_flag1 == 0:
                if tgrad > 0.5*max(temp_grad):
                    half_T_flag1 = 1
                    tstep1 = j
                    
            elif half_T_flag2 == 0:
                if tgrad < 0.5*max(temp_grad):
                    half_T_flag2 = 1
                    tstep2 = j
                else:
                    tstep2 = 0


    # Exothermic time for CV explosion
    if tstep2 == 0:
        print('Error: No pulse in the temperature gradient')
        print('       Your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong') 
        output['exo_time'] = 0
    else:
        output['exo_time'] = output['time'][tstep2] - output['time'][tstep1]

     
    return output