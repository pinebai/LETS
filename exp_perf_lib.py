# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 08:40:09 2017

@author: LawsonJoelM
"""
from numpy import linspace,sin,cos,tan,deg2rad,rad2deg,sqrt,insert,delete
import warnings
# Suppress the 'divide by zero' warning encountered in the trig calls for wedge
# Numpy actually handles this case and returns the correct value for the limit as the arg tends to inf
warnings.filterwarnings("ignore", category=RuntimeWarning)

##############################################################################
# Classes
class PerfectGas:
    """
    Represents the state of a perfect gas.
    Can initialize from scratch by getting data from an Excel file, or can copy data
    from an existing gas object.
    
    Notes on some of the non-obvious initialization arguments:
        - gas: name of the gas (str) 
                -- if not copying, must match the name of a gas in the gas properties database
        - df: dataframe containing gas property data (pandas dataframe) 
                -- if copying, not required, so can pass anything here (neatest is just to put None)
        -copy: if not provided, defaults to None. Otherwise, must be a tuple containing (cv,cp,Mw,h)
               from the gas object you are copying from
    """    
    def __init__(self,gas,pressure,temperature,velocity,df,copy=None):
        self.gas = gas # string
        self.P = pressure # kPa, float
        self.T = temperature # K, float
        self.vel = velocity # m/s, float
        self.copy = copy # For passing to clones
        
        if copy is None:
            # Get perfect gas data from file
            self.pgData(df)
        else:
            # Get perfect gas data from another perfect or real gas
            cv,cp,Mw,h = copy
            self.Cv = cv
            self.Cp = cp
            self.Mw = Mw
            self.gamma = cp/cv
            self.h = h
#            self.Hf = Hf
         
        self.calc_properties()     
        
    def pgData(self,df):
        """
        Appends perfect gas thermodynamic data from the dataframe (originating from an Excel file) to itself.
        """
        self.data = df.loc[self.gas]
        # Note the Excel file currently uses J/kg.K for R but kJ/kg.K for Cp and Cv
        
        # for syntactical convenience and compability when an gas object tries to copy
        # itself from another already-copied gas object
        self.gamma = self.data.gamma
        self.Cv = self.data.Cv
        self.Cp = self.data.Cp
        self.Mw = self.data.Mw
#        self.Hf = self.data.Hf
       
        
    def calc_properties(self):
        """
        Uses known conditions T,P,u and the imported thermodynamic constants to
        calculate sound speed, Mach number, density, enthalpy and corresponding
        stagnation properties.
        """
        from cantera import gas_constant
        Ru = gas_constant/1000 # universal gas constant, J/mol.K
        Tref = 298.15 # reference temperature for enthalpies. Must match reference temperature for formation enthalpies.
        from numpy import sqrt
        
#        if self.copy is None: # already provided if copied from another gas
#        self.href = 1e3*self.Hf/self.Mw # specific formation enthalpy, kJ/kg
        self.h = self.Cp*(self.T-Tref) # specific enthalpy, kJ/kg
            
        self.R = Ru/(self.Mw/1e3)
        self.a = sqrt(self.gamma*self.R*self.T) # calculate speed of sound
        self.M = self.vel/self.a # Mach number
        self.density = self.P*1e3/(self.R*self.T) # density, kg/m^3
        
        self.h0 = self.h + (0.5*self.vel**2)/1e3 # kJ/kg, from Tref
        

##############################################################################
# Functions    
def polar_solve(preShockGas,preExpGas,pg_filename):
    """
    Inputs:
        - PerfectGas objects for the pre-shock and pre-expansion states
    Outputs:
        - PerfectGas objects for the post-shock and post-expansion states
        - Shock Mach number
        - Slopes of the head and tail of the expansion fan
        
    First calls fsolve on polar_eqn to find the velocity ordinate of the
    intersection of the p-u polars.
    """
    # Prepare tuple of parameters to give to polar_eqn
    polar_parameters = (preShockGas.vel,preExpGas.vel,
                        preShockGas.P,preExpGas.P,
                        preShockGas.gamma,preExpGas.gamma,
                        preShockGas.a,preExpGas.a)
    
    from scipy.optimize import fsolve
    
    for u0 in range(100,10000,500): # attempt initial guesses until soln found
        soln = fsolve(polar_eqn,u0,args=polar_parameters,full_output=True)
        if soln[2] == 1:
            postu = soln[0][0] # velocity post-shock/expansion
            break
        else:
            postu = -1 # solution not found: so far this output is not handled anywhere!
    
    postP = polar_P(preExpGas.P,preExpGas.vel,postu,preExpGas.gamma,preExpGas.a)
    Ms = Ms_from_Pratio(postP/preShockGas.P,preShockGas.gamma)
    
    postShockGas = PerfectGas(preShockGas.gas,postP,
                             shock_jump_T(preShockGas.T,Ms,preShockGas.gamma),
                             postu,pg_filename)
    
    postExpGas, dxdt = unsteady_exp(preExpGas,postu,postP,pg_filename)      

    return Ms, postShockGas, postExpGas, dxdt     

    
def polar_eqn(u,*parameters):
    """
    Defines implicitly the intersection between a left-facing expansion wave
    and a right-facing shock in the p-u plane. Given to a solver function
    in order to obtain the value of u at the intersection.
    
    Here state 1 refers to the pre-shock state, and state 4 to the pre-expansion
    state.
    """
    from numpy import sqrt
    u1,u4,P1,P4,g1,g4,a1,a4 = parameters # unpacks tuple
    
    x = (g1+1)/4*(u-u1)/a1 # dummy variable for a term in the expression for P2
    P2 = P1*(1+g1*(u-u1)/a1*(x+sqrt(x**2+1))) # post-shock pressure
    
    P3 = P4*(1-(g4-1)/2*(u-u4)/a4)**(2*g4/(g4-1)) # post-expansion pressure
    
    return P2-P3 # the solver will try to find P2-P3=0 as desired

    
def polar_P(P1,u1,u2,g,a1):
    """
    Inputs:
        - P1: pressure in state 1
        - u1: velocity in state 1
        - u2: velocity in state 2
        - g: gamma (same in both states)
        - a1: sound speed in state 1 
    Outputs:
        - P2: pressure in state 2
    This is a simple function to back out the pressure ordinate at the intersection
    of the shock and expansion polars once the velocity ordinate has been
    solved for. Uses the expansion since it is a simpler polar than the shock.
    """
    P2 = P1*(1-(g-1)/2*(u2-u1)/a1)**(2*g/(g-1))
    return P2

    
def shock_jump_T(T1,Ms,g):
    """
    Inputs:
        - T1: pre-shock temperature
        - Ms: shock Mach number
        - g: gamma
    Outputs:
        - T2: post-shock temperature
    """
    T2 = T1*(g*Ms**2-(g-1)/2)*((g-1)/2*Ms**2+1)/((g+1)*Ms/2)**2
    
    return T2

def shock_jump(preShockGas,Ms):
    """
    Calculates post-shock state for a shock of known Mach number propagating into
    a stationary upstream gas of known thermodynamic state.
    
    Inputs:
        - preShockGas
        - Ms: shock Mach number
        - pg_filename
    Outputs:
        - postShockGas
    """
    # Use 'u' for velocities in lab frame (moving shock)
    # Use 'w' for velocities in shock frame (fixed shock)
    us = Ms*preShockGas.a # shock speed
    w1 = us
    g = preShockGas.gamma
    
    T2 = preShockGas.T*(g*Ms**2-(g-1)/2)*((g-1)/2*Ms**2+1)/((g+1)*Ms/2)**2
    P2 = preShockGas.P*(1+2*g/(g+1)*(Ms**2-1))
    w2 = w1*(2+(g-1)*Ms**2)/((g+1)*Ms**2)
    
    u2 = us-w2
    
    thermo_copy = (preShockGas.Cv,preShockGas.Cp,preShockGas.Mw,preShockGas.h)
    postShockGas = PerfectGas(preShockGas.gas,P2,T2,u2,None,copy=thermo_copy)
    
    return postShockGas
    
    
def Ms_from_Pratio(Pratio,g):
    """
    Inputs:
        - Pratio: the pressure ratio across the shock: P2/P1
        - g: gamma
    Outputs:
        - Ms: shock Mach number
    """
    from numpy import sqrt
    Ms = sqrt(1/(2*g)*((g+1)*Pratio+g-1))
    return Ms
      

def unsteady_exp(knownGas,u_target,P_target,pg_filename):
    """
    Takes a known gas state upstream of (i.e. unprocessed by) the unsteady expansion.
    The expansion wave is assumed to be left-facing, so a C+ characteristic is
    used to connect the states. It is assumed that the downstream velocity and
    pressure are already known by matching across the contact surface.
    Thus we calculate the downstream sound speed and thereby temperature.
    
    A new PerfectGas object is generated to represent the downstream state,
    as well as returning the slopes dx/dt of the head and tail of the wave.    
    """ 
    a_target = knownGas.a + (knownGas.gamma-1)/2*(knownGas.vel-u_target)
    T_target = a_target**2/(knownGas.gamma*knownGas.R)
    
    # same gas species
    if knownGas.copy is None:
        targetGas = PerfectGas(knownGas.gas,P_target,T_target,u_target,pg_filename)
    else:
        targetGas = PerfectGas(knownGas.gas,P_target,T_target,u_target,None,knownGas.copy)
    
    if knownGas.vel < u_target:
        dxdt = (knownGas.vel-knownGas.a,u_target-a_target) # head, tail
    else:
        dxdt = (u_target-a_target,knownGas.vel-knownGas.a)
    
    return targetGas,dxdt


def wedge(freeGas,th_wedge):
    """
    Calculates flow over a wedge, based on the current thermodynamic state of
    the freestream gas object. Takes the wedge angle theta as an input.
    
    Returns the shock angle beta, or else indicates flow is detached.
    If attached, also returns the new thermodynamic state, including the
    new velocity (parallel to the wedge angle).
    
    Based on the perfect gas theta-beta-M eqn., i.e.
    
    tan(th) = 2*cot(B) * (M^2 sin^2(B) - 1) / (M^2 (g + cos(2B)) + 2)
    
    All angles in degrees.
    
    Works by calculating entire beta-theta curve (for 0 < beta < 90 deg) for a known M.
    The full curve contains n_coarse data points.
    Searches for the peak of this curve. Refines the curve using n_fine data points in the
    vicinity of the coarse maximum. Searches again and finds the fine maximum.
    
    If the wedge angle exceeds this maximum, flow must be detached, so function returns.
    Otherwise, repeats the coarse search followed by refinement about the wedge angle target.
    Searches only to the left of the maximum (to ensure the physical solution is the one returned),
    and finds the oblique shock angle corresponding to the wedge angle.
    
    Once oblique shock angle determined, calculates the post-shock flow.
    """ 
    n_coarse = 500
    n_fine = 5000

 
    coarse_beta = linspace(0,90,num=n_coarse) # deg    
    coarse_th = rad2deg(th_B_M(deg2rad(coarse_beta),freeGas.M,freeGas.gamma)) # deg
    coarse_max_idx = coarse_th.argmax()
    
    fine_max_beta_ins = linspace(coarse_beta[coarse_max_idx-1],coarse_beta[coarse_max_idx+1],num=n_fine)
    a = delete(coarse_beta,range(coarse_max_idx-1,coarse_max_idx+2))
    fine_max_beta = insert(a,coarse_max_idx-1,fine_max_beta_ins)
    fine_max_th = rad2deg(th_B_M(deg2rad(fine_max_beta),freeGas.M,freeGas.gamma))
    
    fine_max_idx = fine_max_th.argmax()
    th_max = fine_max_th[fine_max_idx]
    
    if th_wedge < th_max:
        coarse_target_idx = (abs(fine_max_th[:fine_max_idx]-th_wedge)).argmin()
        fine_target_beta_ins = linspace(fine_max_beta[coarse_target_idx-1],fine_max_beta[coarse_target_idx+1],num=n_fine)
        b = delete(fine_max_beta,range(coarse_target_idx-1,coarse_target_idx+2))
        fine_target_beta = insert(b,coarse_target_idx-1,fine_target_beta_ins)
        
        fine_target_th = rad2deg(th_B_M(deg2rad(fine_target_beta),freeGas.M,freeGas.gamma))
        fine_target_idx = (abs(fine_target_th-th_wedge)).argmin()
        beta = fine_target_beta[fine_target_idx]
    else:
        beta = None
        postOblGas = None
    
    if beta is not None:
        # Treat oblique shock as normal shock by taking only normal component
        # of freestream Mach number across stationary shock
        postOblGas = shock_jump(freeGas,freeGas.M*sin(deg2rad(beta)))
        # Fix normal velocity ref. frame.
        postOblGas.vel = freeGas.vel*sin(deg2rad(beta))-postOblGas.vel
        # Get total velocity by including unchanged tangential component.
        postOblGas.vel = sqrt(postOblGas.vel**2+(freeGas.vel*cos(deg2rad(beta)))**2)
        postOblGas.calc_properties()
   
    return beta,th_max,postOblGas
        
def th_B_M(B,M,g):
    """
    The analytical function relating theta, beta and Mach number for perfect-gas
    flow over a wedge.
    """
    th = 2/tan(B)*(M**2*(sin(B))**2-1)/(M**2*(g+cos(2*B))+2)
    return th


        
        