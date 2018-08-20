# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:50:04 2017

@author: jmlawson
"""
##############################################################################
# Import external modules
import sdtoolbox as sdt
from scipy.integrate import quad
from numpy import logspace,log10,sin,cos,arctan,deg2rad,rad2deg,linspace,insert,delete,zeros,sqrt

# Import Cantera with a filter applied to prevent warnings being printed to the console
# (most commonly, temperature out-of-bounds during iteration)
import warnings
# Suppress the 'divide by zero' warning encountered in the trig calls for wedge
# Numpy actually handles this case and returns the correct value for the limit as the arg tends to inf
warnings.filterwarnings("ignore", category=RuntimeWarning)
with warnings.catch_warnings():
    warnings.simplefilter(action="ignore")
    import cantera as ct


from contextlib import contextmanager
import sys, os
        
@contextmanager
def suppress_stdout():
    """
    A function that can be used to suppress another function from writing to stdout (the console).
    Does so by redirecting output to null device (i.e. simply discards output).
    Useful for suppressing messages from functions that clutter the console during a calculation.
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout


##############################################################################
# Classes
class RealGas(ct.Solution):
    """
    This class represents a real gas state. It inherits from the Cantera Solution class,
    so has all the properties and methods of that class (e.g. P,T,h; equilibrate).
    
    Additionally, has extra properties vel,a,M (velocity, sound speed and Mach number)
    with the option to also have Rem (Reynolds number per unit length) if the Cantera input file
    includes transport data.
    
    Create a member of this class using the keyword arguments P,T,X,equil,vel,gas alongside the
    usual 'mech' argument of the parent class, e.g.:
        
        exampleGas = RealGas('air.cti',P=1e5,T=300,X='default',equil=True,vel=2e3,gas='Air')
        
        The keyword arguments represent:
            P: pressure (Pa) [float]
            T: temperature (K) [float]
            X: mole fraction [can give either list format accepted by Cantera, e.g. 'N2:0.7,O2:0.3' or [0.7,0,0,0.3],
                              or can give 'default', which takes the default X as defined in the cti file]
            equil: whether to equilibrate the gas when calculating its sound speed [boolean]
            vel: velocity (m/s) [float]
            gas: the name that will appear in output files to refer to this gas [str]
    """
    def __init__(self,mech,**kwargs):
        self.equil = kwargs['equil']
        self.vel = kwargs['vel']
        self.mech = mech
        ct.Solution.__init__(self,mech)
        
        with warnings.catch_warnings():
            # Suppress FutureWarning about elementwise comparison
            warnings.simplefilter(action="ignore",category=FutureWarning) 
            if (kwargs['X'] != 'default') and (len(self.X) != 1):
                self.TPX = kwargs['T'],kwargs['P'],kwargs['X']
            else:
                self.TP = kwargs['T'],kwargs['P']
        
        self.gas = kwargs['gas'] # 'nice' name for output file
        self.extra()
                
    def extra(self):
        """
        Calculates extra properties beyond those provided by the parent Cantera class,
        namely sound speed, Mach number and total enthalpy (using the provided velocity and static enthalpy).
        
        The sound speed is calculated with or without equilibrating the gas, depending on
        the 'equil' property of the object.
        
        This method is called upon initialization, and should be called again whenever the object's
        properties are altered, to make sure the values remain correct.
        """
        #self.a = equilSoundSpeeds(self)
        if self.equil:
            self.a = sdt.thermo.soundspeed_eq(self)
        else:
            self.a = sdt.thermo.soundspeed_fr(self)
        
        self.M = self.vel/self.a
        self.calc_href()
        self.h0 = self.h - self.href + 0.5*self.vel**2

    def reynolds(self):
        """
        Defines a Reynolds number per unit length, i.e. rho*U/mu
        
        Because this method requires a viscosity, only call this if the Cantera input file
        includes transport properties; otherwise, errors will ensue.
        """
        self.Rem = self.vel*self.density/self.viscosity
        
    def calc_href(self):
        Tref = 298.15
        gas = ct.Solution(self.mech)
        gas.TPX = Tref,ct.one_atm,self.X
        gas.equilibrate('TP')
        self.href = gas.h


##############################################################################
# Functions  
def riemann_integrand(P,s,gas):
    """
    Holding entropy constant, returns 1/(rho*a) as a function of P.
    Note: the gas object (could be RealGas or Cantera) gets modified every time 
    this function is called, because those classes do not create a separate copy
    when passed to a function; rather, a handle to the global object is given.
    
    Note that this assumes a left-facing expansion (i.e. calculated along a C+ characteristic),
    hence the negative sign on the integrand.
    """
    # Suppress Cantera from printing to command line while iterating
    with suppress_stdout(): 
        gas.SP = s,P    
        gas.equilibrate('SP') 
        gas.extra()
    return -1/(gas.a*gas.density) # negative on C+ (lfw)

def PostExp_eq(P2,P1,u1,S1,gas_down):
    """
    This is a custom function that parallels SDToolbox's PostShock_eq() function, for expansion waves
    rather than shock waves.
    
    For a known upstream state and velocity, returns downstream state and velocity
    for an unsteady isentropic expansion to a specified pressure. Assumes that the expansion
    is in quasi-equilibrium, i.e. equilibrates the gas at every iteration point through the expansion wave.
    
    gas_down is a RealGas object representing the downstream gas state, 
    When passed to this function it should be initially set to the known conditions 
    of the upstream state, but will be modified by the subroutine riemann_integrand,
    eventually obtaining downstream properties when this function converges.
    """
    numbreaks = 1 # initial number of breakpoints to try
    maxbreaks = 100 # maximum number of breakpoints
   
    # The integration method (scipy.integrate.quad) can fail to converge if the function isn't well-behaved
    # or changes too quickly. Can manually set breakpoints to divide the integration range into smaller pieces
    # where the changes are less. This takes longer, however. So the approach here is to start with a single breakpoint
    # (i.e. just split the range in 2) and increase the number of breakpoints if the integration fails.
    # So far, have never found a condition that still won't converge with 30 breakpoints.
    try:
        # Turn full_output on to suppress warning messages - now they get written to (discarded) log dictionary
        
        # Integrate along characteristic from P1 to P2, holding entropy constant at S1
        all_out = quad(riemann_integrand,P1,P2,args=(S1,gas_down),full_output=1)
    except RuntimeError:
        while True: # do indefinitely
            breakpoints = logspace(log10(P1),log10(P2),num=numbreaks+2)
            try:
                all_out = quad(riemann_integrand,P1,P2,args=(S1,gas_down),points=breakpoints[1:-2],full_output=1)
                break
            except:
                numbreaks += 1
                print('Numbreaks set to '+str(numbreaks))
                if numbreaks > maxbreaks:
                    # Should be a RuntimeError but want to distinguish it from other Cantera errors that might
                    # be thrown elsewhere
                    raise ValueError('Riemann integration still failed after max. number of breakpoints reached')
                    break
    
    # The first argument of the output is the integral, which is the velocity change            
    delta_u = all_out[0]
    u2 = u1 + delta_u
    # Make sure gas_down is equilibrated to the correct conditions in case integrate.quad doesn't end on P2 for its last iteration
    gas_down.SP = S1,P2
    gas_down.equilibrate('SP')
    gas_down.vel = u2
    gas_down.extra()
    # gas_down now fully represents the final downstream (post-expansion) state
    return


def real_polar_eqn(Ms,*parameters):
    """
    This is the equation that represents the shocktube problem for real gases. It is
    passed to the numerical solver in real_polar_solve(). Since the solver is a minimization
    algorithm, this function returns the absolute difference between the post-shock and post-expansion
    velocities (which should be zero, across the contact surface).
    
    The first input argument is the (unknown) shock Mach number which the solver will try to optimize.
    The additional input parameters are 4 gas objects representing the 4 gas states in the shocktube problem.
    The pre-shock and pre-expansion gases (e.g. states 1 and 4) are known and fixed. The post-shock and post-expansion
    gases are initially unknown, but as mentioned in other docstrings in this library, the gas objects are essentially
    global in scope, so as this function is repeatedly called while the solver iterates, their properties will be altered
    until they converge to their correct states.
    
    This function works by calculating the post-shock state that would result from the provided shock Mach number
    (using PostShock_eq). It then forces the pressures to match across the contact surface (i.e. sets P2=P3) 
    then finds the velocity that would result from an unsteady expansion to that value (using PostExp_eq).
    
    Note that it is assumed that both post-wave states are in equilibrium. Future additions could easily add
    a toggle here that optionally calls PostShock_fr (and a new PostExp_fr) if frozen conditions need exploring.
    """
    print('Iterating on pressure-velocity polar equation...')
    # Note we don't want to modify gas 1 or gas 4 in any way. gas3 will be modified by
    # PostExp_eq, and should initially be set with equal conditions to gas4.
    preShockGas,postShockGas,postExpGas,preExpGas = parameters # unpack tuple
    
    us = Ms*preShockGas.a
    # PostShock_eq generates a Cantera Solution object, but we want our new
    # gas states to be represented as RealGas objects (children of Solution objects)
    # To avoid modifying PostShock_eq to do this, create a separate RealGas object
    # that clones all the properties calculated post-shock
    postShockGas_dummy = sdt.postshock.PostShock_eq(us,preShockGas.P,preShockGas.T,preShockGas.X,preShockGas.ID+'.cti')
    u2 = us - preShockGas.density/postShockGas_dummy.density*(us-preShockGas.vel)
    postShockGas.vel = u2
    
    # Copy values from dummy to actual postShockGas object
    if len(postShockGas_dummy.X)>1: # Cantera throws error if X has length 1 
    # (fixed in future version of Cantera after I submitted a bug report)
    # Can remove this "if" block and just use .TPX once this is confirmed fixed
        postShockGas.TPX = postShockGas_dummy.T,postShockGas_dummy.P,postShockGas_dummy.X
    else:
        postShockGas.TP = postShockGas_dummy.T,postShockGas_dummy.P
    
    postShockGas.extra()
    
    PostExp_eq(postShockGas.P,preExpGas.P,preExpGas.vel,preExpGas.s,postExpGas)

    return abs(postShockGas.vel-postExpGas.vel)
  
    
def real_polar_solve(preShockGas,postShockGas,postExpGas,preExpGas,Ms_perf):
    """
    Top-level function for solving the real-gas shocktube equation.
    The name comes from the fact that we are essentially finding states 2 and 3
    at the intersection of the shock and expansion polars in the pressure-velocity
    plane (done more explicity for the perfect gas case, where the equations for the polars
    are known explicitly).
    
    For arguments, takes 4 gas objects representing the 4 gas states in the shocktube problem
    (could be gases 1-4, but could also be gases 2, 5, 6, and 7 if it is the second 'shocktube'
    in the overall expansion tube problem). The final argument is the initial guess for Ms, 
    given to the solver. We should usually use the solution from the corresponding perfect-gas
    problem.
    """
    # Prepare tuple of parameters to give to polar_eqn
    parameters = (preShockGas,postShockGas,postExpGas,preExpGas)
    
    from scipy.optimize import fminbound
    window = [3.,3.] # below and above perfect Mach number
    lower_bound = max([1.,Ms_perf-window[0]])
    upper_bound = Ms_perf+window[1]    
    
    # After a lot of testing, this was found to be the most robust of the numerical solvers
    # trialled. Doesn't take an initial guess directly, but instead finds the minimum value of 
    # the objective function within 2 bounds, which we construct somewhat arbitrarily about Ms_perf
    Ms = fminbound(real_polar_eqn,lower_bound,upper_bound,args=parameters)
    
    # Ensure sound speeds, Mach numbers and other extras updated
    postShockGas.extra()
    postExpGas.extra()
    
    # Calculate slopes of head and tail characteristics of expansion wave
    dxdt = exp_slopes(preExpGas,postExpGas)
    
    return Ms, dxdt

def real_polar_eqn_inv(drivP,*parameters):
    """
    An alternative version of real_polar_eqn where Ms is known and the pre-expansion
    pressure is unknown. Used for the shock-specified mode of calculation.
    """

    preShockGas,postShockGas,postExpGas,preExpGas,Ms = parameters # unpack tuple
    
    us = Ms*preShockGas.a
    
    postShockGas_dummy = sdt.postshock.PostShock_eq(us,preShockGas.P,preShockGas.T,preShockGas.X,preShockGas.ID+'.cti')
    u2 = us - preShockGas.density/postShockGas_dummy.density*(us-preShockGas.vel)
    postShockGas.vel = u2    
    
    # Copy values from dummy to actual postShockGas object
    if len(postShockGas_dummy.X)>1: # Cantera throws error if X has length 1 
    # (fixed in future version of Cantera after I submitted a bug report)
        postShockGas.TPX = postShockGas_dummy.T,postShockGas_dummy.P,postShockGas_dummy.X
    else:
        postShockGas.TP = postShockGas_dummy.T,postShockGas_dummy.P
        
    # Match across contact surface
    postExpGas.TP = postExpGas.T,postShockGas.P
    postExpGas.vel = postShockGas.vel
    # Set (unknown) driver pressure to current value for iteration. Fix T.
    preExpGas.TP = preExpGas.T,drivP
    
    # Equilibrate, update additional properties
    preExpGas.extra()
    postExpGas.extra()
    postShockGas.extra()
    
    
    PostExp_eq(postExpGas.P,drivP,preExpGas.vel,preExpGas.s,postExpGas)
    
    return abs(postShockGas.vel-postExpGas.vel)

def real_polar_solve_inv(preShockGas,postShockGas,postExpGas,preExpGas,Ms,drivP_perf):
    """
    Alternative version of real_polar_solve used for the 'inverse' problem, i.e. when the
    shock is specified but the driver pressure is unknown. Calls real_polar_eqn_inv instead
    of real_polar_eqn. Similar to the original, uses the solution from the corresponding
    perfect-gas problem as an initial guess. Has an extra argument since Ms is still required.
    """
    # Prepare tuple of parameters to give to polar_eqn
    parameters = (preShockGas,postShockGas,postExpGas,preExpGas,Ms)
    
    from scipy.optimize import fminbound
    window = 2. # scaling on drivP_perf
    lower_bound = drivP_perf/window
    upper_bound = drivP_perf*window
    
    drivP = fminbound(real_polar_eqn_inv,lower_bound,upper_bound,args=parameters)
    # Ensure driver pressure set
    preExpGas.TP = preExpGas.T,drivP
    
    # Ensure sound speeds, Mach numbers and other extras updated
    preExpGas.extra()
    postShockGas.extra()
    postExpGas.extra()
    
    dxdt = exp_slopes(preExpGas,postExpGas)
    
    return drivP, dxdt

def exp_slopes(gas_up,gas_down):
    """
    Calculates the slopes of the head and tail of a left-facing expansion,
    given the known gas states upstream (the head, to the left) 
    and downstream (the tail,to the right).
    """
    dxdt = (gas_up.vel-gas_up.a,gas_down.vel-gas_down.a)
    return dxdt

def detonation_setup(preDetState):
    """
    This function takes a known pre-detonation state (as represented by a tuple of values, detailed in comments), 
    and calculates the post-detonation state along with the CJ detonation wave speed, using functions from SDToolbox.
    
    It returns RealGas objects representing the pre- and post-detonation states, as well as PerfectGas approximations
    of these states (since to solve the full shocktube problem, a perfect-gas solution is required for an initial guess).
    The last output is the CJ speed (in m/s).
    
    Note that for brevity, state 4 and 4a refer to the pre- and post-detonation states, respectively.
    """
    from exp_perf_lib import PerfectGas
    mech,P,T,comp,name = preDetState
    # Don't use equilibrate in rgas4, so that we get the correct pre-detonation density and don't alter the mole fractions
    # (i.e. call frozen sound speed instead)
    rgas4 = RealGas(mech,P=P,T=T,vel=0,X=comp,gas=name,equil=False)
    rgas4a = RealGas(mech,P=P,T=T,vel=0,X=comp,gas=name,equil=True) # initially, a copy of the pre-det state
    
    cj_speed = sdt.postshock.CJspeed(P,T,comp,mech)
    # Use the wave speed to get post-detonation state (Cantera Solution class)
    cjGas = sdt.postshock.PostShock_eq(cj_speed, P, T, comp, mech) 
    
    # Clone post-detonation state into RealGas object
    rgas4a.TPX = cjGas.T,cjGas.P,cjGas.X
    rgas4a.vel = (rgas4.density/rgas4a.density-1)*cj_speed
    rgas4a.extra()
    
    # Clone RealGas objects into PerfectGas objects
    # Convert enthalpy J/kg -> kJ/kg, pressure Pa -> kPa
    # Maintain temperature in K, Mw in g/mol (amu)
    toPerfect4 = (rgas4.cv,rgas4.cp,rgas4.mean_molecular_weight,rgas4.h/1e3)
    toPerfect4a = (rgas4a.cv,rgas4a.cp,rgas4a.mean_molecular_weight,rgas4a.h/1e3)
    pgas4 = PerfectGas(rgas4a.gas,rgas4.P/1e3,rgas4.T,0,None,toPerfect4)
    pgas4a = PerfectGas(rgas4.gas,rgas4a.P/1e3,rgas4a.T,rgas4a.vel,None,toPerfect4a)
    
    return pgas4a,rgas4a,pgas4,rgas4,cj_speed

def oblique_eqn(beta_range,freeGas,gasReturn=False):
    """
    Given a list of oblique shock angles (beta_range) and a RealGas object representing the known
    freestream state, return the deflection angle the flow will experience passing through the stationary
    oblique shock. 
    
    By default, only returns the list of deflection angles (theta) but if the optional keyword
    argument gasReturn=True, will also return the post-oblique shock thermodynamic state (as a Cantera Solution class)
    as well as the velocity components (u1t=u2t, u2n -- representing tangential and normal components wrt shock)
    
    Simply takes the normal component of the freestream and calculates the post-shock state for a normal (1D) shock.
    Uses trigonometry and the fact that the post-shock tangential velocity is unaffected to find the deflection angle.
    """
    theta = zeros(len(beta_range))
    for i,beta in enumerate(beta_range):
        u1n = freeGas.vel*sin(beta)
        u1t = freeGas.vel*cos(beta)
        postOblGas_dummy = sdt.postshock.PostShock_eq(u1n,freeGas.P,freeGas.T,freeGas.X,freeGas.ID+'.cti')
        u2n = freeGas.density/postOblGas_dummy.density*u1n    
        theta[i] = beta - arctan(u2n/u1t)
    
    if gasReturn:    
        return theta,postOblGas_dummy,u1t,u2n
    else:
        return theta

def wedge(freeGas,th_wedge,postOblGas):
    """
    Solves for the flow over a wedge in a RealGas. Inputs are:
        - freeGas: RealGas object representing known freestream state
        - th_wedge: theta (deflection angle) of wedge, in deg
        - postOblGas: RealGas object representing post-oblique shock state. Initially unknown,
                      can just provide as a clone of freeGas. Will converge to correct state.
                      
    Initially, attempted to write this function using numerical solvers and the corresponding
    perfect-gas solution as an initial gas. However, due to several factors, this wasn't very robust:
        - there are two values of beta for a given theta
        - above some maximum theta the shock detaches and there is no solution
        - the periodicity of the trigonometric functions used
        
    Instead, used brute-force/mesh-refinement technique. Considerably slower, but much more robust.
    Essentially exactly the same as the same function in the exp_perf_lib (see description there for
    further details).
    """
    n_coarse = 100
    n_fine = 500

    print('Coarse search for detachment angle...')
    coarse_beta = linspace(0,90,num=n_coarse) # deg    
    coarse_th = rad2deg(oblique_eqn(deg2rad(coarse_beta),freeGas)) # deg
    coarse_max_idx = coarse_th.argmax()
    
    print('Fine search for detachment angle...')
    fine_max_beta_ins = linspace(coarse_beta[coarse_max_idx-1],coarse_beta[coarse_max_idx+1],num=n_fine)
    a = delete(coarse_beta,range(coarse_max_idx-1,coarse_max_idx+2))
    fine_max_beta = insert(a,coarse_max_idx-1,fine_max_beta_ins)
    fine_max_th = rad2deg(oblique_eqn(deg2rad(fine_max_beta),freeGas))
    
    fine_max_idx = fine_max_th.argmax()
    th_max = fine_max_th[fine_max_idx]
    
    if th_wedge < th_max: # attached
        print('Coarse search for shock angle...')
        coarse_target_idx = (abs(fine_max_th[:fine_max_idx]-th_wedge)).argmin()
        fine_target_beta_ins = linspace(fine_max_beta[coarse_target_idx-1],fine_max_beta[coarse_target_idx+1],num=n_fine)
        b = delete(fine_max_beta,range(coarse_target_idx-1,coarse_target_idx+2))
        fine_target_beta = insert(b,coarse_target_idx-1,fine_target_beta_ins)
        
        print('Fine search for shock angle...')
        fine_target_th = rad2deg(oblique_eqn(deg2rad(fine_target_beta),freeGas))
        fine_target_idx = (abs(fine_target_th-th_wedge)).argmin()
        beta = fine_target_beta[fine_target_idx]
    else: # detached
        beta = None
    
    if beta is not None:
        print('Obtaining post-oblique shock thermodynamic state...')
        # Treat oblique shock as normal shock by taking only normal component
        # of freestream Mach number across stationary shock
        _,postOblGas_dummy,u1t,u2n = oblique_eqn(deg2rad([beta]),freeGas,gasReturn=True)
        u2 = sqrt(u1t**2+u2n**2)
        
        if len(postOblGas_dummy.X)>1: # Cantera throws error if X has length 1 
            postOblGas.TPX = postOblGas_dummy.T,postOblGas_dummy.P,postOblGas_dummy.X
        else:
            postOblGas.TP = postOblGas_dummy.T,postOblGas_dummy.P
        
        postOblGas.vel = u2
        postOblGas.extra()
        postOblGas.reynolds()
        
    return beta,th_max

def boundary_layer(preShockGas,postShockGas,ushock,length):
    
    # Define shorter names
#    rho5 = preShockGas.density
    rho6 = postShockGas.density
#    p5 = preShockGas.P
#    p6 = postShockGas.P
#    T6 = postShockGas.T
#    mu5 = preShockGas.viscosity
    mu6 = postShockGas.viscosity
    u67 = postShockGas.vel
    g5 = preShockGas.cp/preShockGas.cv
    
    uw = ushock
    ue0 = uw - u67
    
    W = uw/ue0
#    Z = (g5+1)/(g5-1)
    
    Re_lm = ue0*(W-1)**2*rho6*length/mu6
    
    Pr=0.72
    r0=Pr**(1/3)
    Z=(g5+1)/(g5-1)
    hr_he0=1+(((W-1)**2)/(Z*W-1))*r0
    he0_hw=(Z*W-1)/(W*(Z-W))
    hr_hw=hr_he0*he0_hw
    b=hr_hw-1
    ch=hr_hw-he0_hw
    #print(['b',b])
    #print(['ch',ch])
    #print('***Edit Tracker***')
    #print([Pr, r0, Z])
    
    
    return Re_lm