"""
Shock and Detonation Toolbox
"znd" module

Calculates ZND explosions.
 
This module defines the following functions:

    zndsys
    zndsolve
    
#### WARNING ####
As noted in the warning in the cv module, MATLAB's proprietary solvers (e.g. ode15s) 
are superior to the solvers in Python's scipy.integrate module. 
The ZND ODE system tends to solve more easily than the corresponding system for CV explosions, 
but the MATLAB version still outperforms this one.


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
import sdtoolbox as sdt
from scipy.integrate import ode

def zndsys(t, y, gas, U1, r1, PSC):
    """
    Set of ODEs to solve ZND Detonation Problem. Note that the pressure inputs and
    outputs in the ODE system are scaled to bring their magnitude closer to the
    other variables. Scaling is with respect to the initial pressure value, PSC.

    INPUT:
        t = time
        y = solution array [scaled pressure, density, position, species mass 1, 2, ..]
        gas = working gas object
        U1 = shock velocity (m/s)
        r1 = initial density (kg/m^3)
        PSC = postshock pressure value (Pa) (unscaled, i.e. original magnitude)
    
    OUTPUT:
        An array containing time derivatives of:
            (scaled) pressure, density, distance and species mass fractions, 
        formatted in a way that the integrator in cvsolve can recognize.
        
    """
    # Can't just do gas.DY like in MATLAB. Can do gas.TDY, use a dummy T the first pass
    dummyT = 300.
    gas.TDY = dummyT,y[1],y[3:] 
    T = (y[0]*PSC/y[1])*(gas.mean_molecular_weight/ct.gas_constant)
    gas.TDY = T,y[1],y[3:] 
    c = sdt.thermo.soundspeed_fr(gas)
    U = U1*r1/gas.density
    M = U/c
    eta = 1-M**2 
    
    sigmadot = getThermicity(gas)
    Pdot = -gas.density*U**2*sigmadot/eta/PSC # This is a scaled pressure, i.e. divided by PSC
    rdot = -gas.density*sigmadot/eta
    
    dYdt = gas.net_production_rates*gas.molecular_weights/gas.density    

    return np.hstack((Pdot, rdot, U, dYdt))


def getThermicity(gas):
    """
    Returns the thermicity = sum ( (w/wi-hsi/(cp*T))*dyidt ). Used by zndsys.
    
    FUNCTION SYNTAX:
        thermicity = getThermicity(gas)
        
    INPUT:
        gas = Cantera gas object (not modified by this function)
        
    OUTPUT:
        thermicity (1/s)
    """
    thermicity=0.0
    for i in range(gas.n_species):
        w=gas.mean_molecular_weight  # mean molecular weight of the mixture
        wi=gas.molecular_weights[i]  # molecular weight of the ith species
        R = ct.gas_constant          # universal gas constant
        hsi = gas.standard_enthalpies_RT[i]*R*gas.T/wi #s pecific enthalpy of the ith species
        cp = gas.cp_mass               # mixture averaged cp, frozen 
        dyidt = gas.net_production_rates[i]*wi/gas.density # dYi/dt: rate of change of the mass fraction of the ith species
        thermicity =  thermicity + (w/wi-hsi/(cp*gas.T))*dyidt
        
    return thermicity


def zndsolve(gas,gas1,U1,t_end=1e-3,dt=1e-8,backend='lsoda',relTol=1e-5,absTol=1e-8):
    """
    ZND Model Detonation Struction Computation
    Solves the set of ODEs defined in zndsys.
    
    FUNCTION SYNTAX:
    output = zndsolve(gas,gas1,U1,**kwargs)
    
    INPUT
        gas = Cantera gas object - postshock state
        gas1 = Cantera gas object - initial state
        U1 = shock velocity (m/s)
        
    OPTIONAL INPUT:
        t_end = end time for integration, in sec
        dt = time step for integration, in sec
        backend = solver backend used for the integrator function
        relTol = relative tolerance
        absTol = absolute tolerance
    
    
    OUTPUT:
        
        output = a dictionary containing the following results:            
            ind_time_ZND = time to maximum thermicity gradient
            ind_len_ZND = distance to maximum thermicity gradient
            exo_time_ZND = pulse width (in secs) of thermicity  (using 1/2 max)
            ind_time_ZND = pulse width (in meters) of thermicity (using 1/2 max)
            
            time = time array
            distance = distance array
            
            T = temperature array
            P = pressure array
            rho = density array
            U = velocity array
            thermicity = thermicity array
            species = species profiles (list of lists)
            
            M = Mach number array
            af = frozen sound speed array
            g = gamma (cp/cv) array
            wt = mean molecular weight array
            sonic = sonic parameter (c^2-U^2) array
    """
    ###########################################################
    # Define initial information
    ###########################################################
    r1 = gas1.density
    PSC = gas.P # PSC = postshock pressure, scaled. Actually stores the original unscaled pressure.
    # --> used so pressure is of similar magnitude to other variables in ODE system,
    #     then PSC is retained to convert scaled pressure solution back to actual values.
    
    x_start = 0.
    y0 = np.hstack((1,gas.density,x_start,gas.Y)) # scaled pressure starts at 1, i.e. PSC/PSC

    # Set the integration parameters for the ODE solver
    tel = [0.,t_end] # Timespan
    
    # Initialize output matrices.
    output = {}
    
    time = [tel[0]]
    distance = [x_start]
    P = [PSC]
    rho = [gas.density]
    species = [gas.Y]
 
    # Set up ODE system for integration
    extraParams = (gas, U1, r1, PSC)
    solver = ode(zndsys)
    solver.set_integrator(backend, method='bdf', rtol=relTol, atol=absTol)
    solver.set_initial_value(y0,tel[0]).set_f_params(*extraParams)
    
    # Call the integrator to march in time - the equation set is defined in zndsys() 
    while solver.successful() and solver.t < tel[1]:
        solver.integrate(solver.t + dt)
        
        #############################################################################################
        # Extract TIME, PRESSURE, DENSITY, DISTANCE and MASS FRACTION arrays from integrator output
        #############################################################################################        
        time.append(solver.t)
        P.append(solver.y[0]*PSC) # rescale pressure to true value
        rho.append(solver.y[1])
        distance.append(solver.y[2])
        species.append(solver.y[3:])
    del solver
        
    output['time'] = time
    output['P'] = P
    output['rho'] = rho
    output['distance'] = distance
    output['species'] = species
    
    b = len(time)
    
    # Initialize additional output matrices where needed
    output['ind_len_ZND'] = 0
    output['ind_time_ZND'] = 0
    output['exo_len_ZND'] = 0
    output['exo_time_ZND'] = 0
    
    output['T'] = np.zeros(b)
    output['U'] = np.zeros(b)
    output['thermicity'] = np.zeros(b)
    output['af'] = np.zeros(b)
    output['g'] = np.zeros(b)
    output['wt'] = np.zeros(b)
    
    
    #############################################################################
    # Extract TEMPERATURE, WEIGHT, GAMMA, SOUND SPEED, VELOCITY, MACH NUMBER, 
    # c^2-U^2, THERMICITY, and TEMPERATURE GRADIENT
    #############################################################################
    
    # Have to loop for operations involving the working gas object
    for i,_ in enumerate(output['P']):
        gas.DPY = output['rho'][i],output['P'][i],output['species'][i]
        af = sdt.thermo.soundspeed_fr(gas) # Use sdt version instead of ct
        U = U1*r1/gas.density
       
        output['T'][i] = gas.T
        output['U'][i] = U
        output['thermicity'][i] = getThermicity(gas)
        output['af'][i] = af
        output['g'][i] = gas.cp/gas.cv
        output['wt'][i] = gas.mean_molecular_weight
        
    
    # Vectorize operations where possible    
    output['M'] = output['U']/output['af']
    eta = 1- output['M']**2
    output['sonic'] = eta*output['af']**2
    
    
    ################################################################################################
    # Find INDUCTION TIME and LENGTH based on MAXIMUM THERMICITY
    ################################################################################################
    n = output['thermicity'].argmax()
    
    output['ind_time_ZND'] = output['time'][n]
    output['ind_len_ZND'] = output['distance'][n]
    
    #######################################################
    # Check for eigenvalue detonation
    #######################################################
    
    if n == b:
        print('Error: Maximum thermicity occurs at the end of the reaction zone')
        print('       You may have an eigenvalue detonation, your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong')
        print(' ')
        print('Mach Number (end of reaction): '+str(output['M'][b])+' - if close to 1, check for eigenvalue detonation')
        output['ind_time_ZND'] = output['time'][b]
        output['ind_len_ZND'] = output['distance'][b] 
        output['exo_time_ZND'] = 0 
        output['exo_len_ZND'] = 0 
        print('Induction Time: '+str(output['ind_time_ZND']))
        print('Exothermic Pulse Time: '+str(output['exo_time_ZND']))
        return output
    
    elif n == 0:
        print('Error: Maximum thermicity occurs at the beginning of the reaction zone')
        print('       You may have an eigenvalue detonation, your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong')
        print(' ')
        print('Mach Number (end of reaction): '+str(output['M'][b])+' - if close to 1, check for eigenvalue detonation')
        output['ind_time_ZND'] = output['time'][0]
        output['ind_len_ZND'] = output['distance'][0] 
        output['exo_time_ZND'] = 0 
        output['exo_len_ZND'] = 0 
        print('Induction Time: '+str(output['ind_time_ZND']))
        print('Exothermic Pulse Time: '+str(output['exo_time_ZND']))
        return output
    
    else:
        max_sigmadot = max(output['thermicity'])
        half_sigmadot_flag1 = 0
        half_sigmadot_flag2 = 0
        # Go into a loop to find two times when sigma_dot is half its maximum
        for j,thermicity in enumerate(list(output['thermicity'])):
            if half_sigmadot_flag1 == 0:
                if thermicity > 0.5*max_sigmadot:
                    half_sigmadot_flag1 = 1
                    tstep1 = j
                    
            elif half_sigmadot_flag2 == 0:
                if thermicity < 0.5*max_sigmadot:
                    half_sigmadot_flag2 = 1
                    tstep2 = j
                else:
                    tstep2 = 0
                    
    
    if tstep2 == 0:
        print('Error: No pulse in the thermicity')
        print('       You may have an eigenvalue detonation, your final integration length may be too short,')
        print('       your mixture may be too rich/lean, or something else may be wrong') 
        output['exo_time_ZND'] = 0
        output['exo_len_ZND'] = 0  
    else:
        output['exo_time_ZND'] = output['time'][tstep2] - output['time'][tstep1]; 
        output['exo_len_ZND'] = output['distance'][tstep2] - output['distance'][tstep1]
        
    
    # Append extra data used to make output file (via znd_fileout)
    output['gas1'] = gas1
    output['U1'] = U1
                    
    
    return output
    