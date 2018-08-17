# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 17:44:53 2016

@author: LawsonJoelM
"""
import numpy as np
##############################################################################
# Construction of x-t diagram
def exp_tube_xt(us,dxdt,gas2,gas4,u7,lengths,axhandle):
    """
    Generate an x-t diagram for the expansion tube.
    Use the diagram to calculate test time and termination mode.
    
    This diagram is plotted on an existing subplot whose handle is passed
    to this function.
    
    Inputs:
        - us: shock speeds ([1,2]-tuple or [1,3]-tuple -- passive/det mode)
        - dxdt: expansion head/tail slopes ([2,2]-tuple)
        - gas2: gas state 2 (post-primary shock)
        - gas4: gas state 4 (driver gas)
        - u7: velocity of gas 7 (test section freestream)
        - lengths: lengths of facility sections ([1,3]-tuple)
        - axhandle: handle of matplotlib axis to plot on
        
    Outputs:
        - t_test: test time (ms) (float)
        - termination: mode of termination (by head or tail) (str)
        - t_contact2: arrival time of 2nd c.s. at test section wrt time when
                        secondary diaphragm bursts (ms) (float)
        - optimal: [optimal test time (ms), optimal test length (m)]
                    highest possible test time for given conditions, occurring
                    when secondary expansion head and tail arrive simultaneously,
                    and acceleration section length that would cause this to occur.
        
    """
    from numpy import abs as npabs
    # Calculate times in milliseconds: use conversion factor
    conv = 1000.
    
    # First pair of waves
    t_head1 = -lengths[0]/dxdt[0][0]*conv # driver length/dxdt1_head
    t_shock1 = lengths[1]/us[0]*conv # driven length/u_shock1
    t_contact1 = (lengths[1]+lengths[2])/gas2.vel*conv # driven + exp length/u_contact1
    
    # Second pair of waves. Times all wrt t_shock1
    t_tail2 = lengths[2]/dxdt[1][1]*conv
    t_shock2 = lengths[2]/us[1]*conv
    t_contact2 = lengths[2]/u7*conv
    
    if dxdt[0][1] < 0: # primary tail propagates left
        t_tail1 = -lengths[0]/dxdt[0][1]*conv # driver length/dxdt_tail
        axhandle.plot([lengths[0],0],[0,t_tail1],'c',label='Expansion tail') 
    else: # primary tail propagates right
        t_tail1 = (lengths[1]+lengths[2])/dxdt[0][1]*conv
        axhandle.plot([lengths[0],sum(lengths)],[0,t_tail1],
        'c',label='Expansion tail')
    
    termination = 'Tail'
    
    if dxdt[1][0] < 0: # secondary head propagates left
        t_head2 = -(lengths[0]+lengths[1])/dxdt[1][0]*conv # driver length/dxdt_tail
        axhandle.plot([lengths[0]+lengths[1],0],[t_shock1,t_shock1+t_head2],'b',label='Expansion head')
        # Calculate test time
        t_test=t_tail2-t_contact2
        cross2 = None
    else: # secondary head propagates right:
        # Find coordinates of intersection with primary contact surface
        delta_t,delta_x = cs_head_int(dxdt[1][0],us[0],gas2.vel,lengths[1])
        # Calculate similarity solution
        t,x,cross2 = similarity(gas2,delta_t,5*t_contact1/conv,lengths[2],dxdt[1][1])
        # Prepare trajectory data from raw similarity solution
        t_traj2 = t*conv + t_shock1
        x_traj2 = x + lengths[0] + lengths[1]
        # Find index of closest point to intersection
        idx2 = (npabs(t_traj2-(t_shock1+delta_t*conv))).argmin()
        # Plot truncated trajectory (after intersection only)
        axhandle.plot(x_traj2[idx2:],t_traj2[idx2:],'b',label='Expansion head')
        # Plot original trajectory (before intersection only)
        axhandle.plot([lengths[0]+lengths[1],lengths[0]+lengths[1]+delta_x],
                       [t_shock1,t_shock1+delta_t*conv],'b')
        # Find index of closest point to test section
        idx3 = (npabs(x-lengths[2])).argmin()
        t_head2 = t[idx3]*conv # in ms, referenced from primary shock arrival
        # Indicate arrival of reflected head at test section
        #axhandle.plot(sum(lengths),t_head2,'ro')
        # Calculate test time
        t_test=min([t_tail2-t_contact2,t_head2-t_contact2])
        if t_head2 < t_tail2:
            termination = 'Head'
       
    # Plot by connecting endpoints
    # First pair of waves
    if len(us) == 2:
        # Don't plot first expansion head if detonation
        axhandle.plot([lengths[0],0],[0,t_head1],'b')
        # Primary head reflection from end wall
        t,x,cross1 = similarity(gas4,t_head1/conv,10*t_head1/conv,lengths[1]+lengths[2],dxdt[0][1])
        t_traj1 = t*conv
        x_traj1 = x + lengths[0]
        # Find index of closest point to endwall
        idx1 = (npabs(x_traj1)).argmin()
        axhandle.plot(x_traj1[idx1:],t_traj1[idx1:],'b')
        
        
    axhandle.plot([lengths[0],lengths[0]+lengths[1]],[0,t_shock1],'k',label='Shock') 
    axhandle.plot([lengths[0],sum(lengths)],[0,t_contact1],'g--',label='Contact surface')
    
    # Second pair of waves (don't need to label again for legend)
    axhandle.plot([lengths[0]+lengths[1],sum(lengths)],[t_shock1,t_shock1+t_shock2],'k')
    axhandle.plot([lengths[0]+lengths[1],sum(lengths)],[t_shock1,t_shock1+t_contact2],'g--')
    axhandle.plot([lengths[0]+lengths[1],sum(lengths)],[t_shock1,t_shock1+t_tail2],'c')
    
    # Indicate test time    
    axhandle.plot([sum(lengths),sum(lengths)],
                   [t_shock1+t_contact2,t_shock1+t_contact2+t_test],'r',
                   linewidth=3,solid_capstyle='round',clip_on=False)
    
    # Indicate reflection intersection
    #axhandle.plot(lengths[0]+lengths[1]+delta_x,t_shock1+delta_t*conv,'ro')
    
    # Indicate intersection of secondary expansion head and tail
    #axhandle.plot(x_traj2[cross2],t_traj2[cross2],'go')
    
    # Optimal test time: given when the secondary expansion head and tail intersect
    # Find the difference between this time and the arrival time of the contact surface
    # at the same position.
    if cross2 is not None:
        x_opt = x_traj2[cross2]; t_opt = t_traj2[cross2]
        L_opt = x_opt-(lengths[0]+lengths[1]) # optimal length of acceleration section
        t_contact2_opt = L_opt/u7*conv # wrt t_shock1
        t_test_opt = t_opt-(t_contact2_opt+t_shock1)
        optimal = (t_test_opt,L_opt)
    else:
        optimal = None
    
    # 'Contact surface' to represent secondary diaphragm
    axhandle.plot([lengths[0]+lengths[1],lengths[0]+lengths[1]],
                   [0,t_shock1],'g--')
    
    # Detonation wave (use length of 'us' argument as a proxy flag for detonation mode)
    if len(us) == 3:
        t_det = lengths[0]/us[2]*conv # driver length/dxdt1_head
        axhandle.plot([lengths[0],0],[0,t_det],'r',label='Detonation')
    
    
    # Aesthetics
    axhandle.set_xlabel('Length (m)')
    axhandle.set_ylabel('Time (ms)')
    axhandle.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    axhandle.set_xlim([0,sum(lengths)])
    
    
    oldTicks = axhandle.get_xticks()
    newTicks = np.sort(np.append(oldTicks,np.cumsum(lengths)))
    axhandle.set_xticks(newTicks)
    
    axhandle.set_xlim([0,sum(lengths)])
    axhandle.set_ylim([0,1.2*(t_tail2+t_shock1)])
     
    return t_test, termination, t_contact2, optimal


def cs_head_int(u_head,u_shock,u_cs,l_test):
    """
    Uses trigonometry to find the intersection of the primary contact surface and
    the downstream-propogating secondary expansion head.
    
    Inputs:
        u_head: Speed of secondary exp. head
        u_shock: Speed of primary shock
        u_cs: Speed of primary contact surface
        l_test: Length of test section
        
        Note we require u_shock > u_cs > u_head for there to be an intersection
        
    Outputs:
        delta_t,delta_x: t,x coordinates of intersection, referenced from arrival
                         of primary shock at secondary diaphragm.
                         
    For meanings of intermediate variables for angles/lengths, refer to written notes.
    """
    from math import sin,cos,atan
    
    alpha = atan(1/u_shock)
    beta = atan(1/u_cs)
    gamma = atan(1/u_head)
    B = l_test/cos(alpha)
    C = B*sin(beta-alpha)/sin(gamma-beta)
    delta_t = C*sin(gamma)
    delta_x = C*cos(gamma)
    
    return delta_t,delta_x

def similarity(gas2,delta_t1,tbound,xbound,dxdt_tail):
    """
    This function uses a similarity solution to calculate the trajectory of the
    reflected secondary expansion head as it interacts with the rest of the expansion.
    Can also be used for reflection of primary expansion from driver end wall.
    
    Inputs:
        gas2: The gas object representing state 2
        delta_t1: The time delay calculated by the function cs_head_int, i.e.
                  from shock arrival until intersection of exp. head and c.s.
        tbound: An upper bound on the timespan since exp. head arrival 
                  time not known a priori. Also want to capture crossing point,
                  so needs to be sufficiently long.
        xbound: A lower bound on the distance to continue extending the solution linearly
                once the characteristics have crossed. Should be at least longer than the
                facility.
        dxdt_tail: Slope of the expansion tail.
                  
    All times in seconds, all lengths in metres.
        
    Outputs:
        t,x: Matching vectors containing coordinates of trajectory, referenced 
             from intersection of secondary exp. head and primary c.s.
        idx: Index of the above vectors where the characteristics cross each other.
    """      
    from numpy import sign,roll,linspace, append as npapp

    
    points = 5000
    # Don't start timespan from zero (otherwise get a NaN in first element of x)
    tfull = linspace(1e-7,tbound,num=points)
    
    eta = gas2.vel/gas2.a + 2/(gas2.gamma-1)*(1-(gas2.gamma+1)/2*(tfull/delta_t1)**
                               (2*(1-gas2.gamma)/(gas2.gamma+1)))
    
    xfull = eta*tfull*gas2.a

    tail_t = xfull/dxdt_tail # straight-line trajectory of tail, sharing x-coords of sim. soln
    
    # Look for sign change in tail_t-tfull to indicate a crossing of the secondary head and tail
    tsign = sign(tail_t-tfull)
    signchange = ((roll(tsign, 1) - tsign) != 0).astype(int) # array of 0 and 1, 1 indicates a sign change
    signchange[0] = 0 # avoid detecting an intersection at origin
    idx = next((i for i, x in enumerate(signchange) if x), None) # gets index of first sign change (should be only one)

    if idx is not None:
        exitslope = (tfull[idx]-tfull[idx-1])/(xfull[idx]-xfull[idx-1])
        exittraj_x = linspace(xfull[idx],max(xbound,5*xfull[idx]),num=points)
        exittraj_t = tfull[idx]+exitslope*(exittraj_x-xfull[idx])
        
        x = npapp(xfull[:idx],exittraj_x)                  
        t = npapp(tfull[:idx],exittraj_t)
    else: # leave unchanged
        t = tfull
        x = xfull
    
    return t,x,idx
##############################################################################
# Formatting the output file
def stateVar_formatted(stateVar,indepDict,P_units):
        """
        For a given state variable out of the 9 that define the facility,
        i.e. P,T,L for driver, test, acceleration, returns 3 formatted strings:
            
            1) nameUnitStr: '<state name> [<state unit>]
            2) valueUnitStr: '<state value> <state unit>'
            3) fullStateStr = '<state name> = <state value> <state unit>'
            
        For P, unit and value displayed matches the unit selected in the UI.
        """
        drivbar,testtorr,acceltorr = P_units
        if (    (stateVar == 'drivP' and drivbar) 
             or (stateVar == 'testP' and testtorr)
             or (stateVar == 'accelP' and acceltorr)):
            # Select the non-SI unit (bar or torr)
            stateUnit = indepDict[stateVar][2][1]
        elif stateVar[-1] == 'P':
            # Select the SI unit (kPa)
            stateUnit = indepDict[stateVar][2][0]
        else:
            # Don't need nested index here since only one option for T or L unit             
            stateUnit = indepDict[stateVar][2]
            
        if stateVar[-1] == 'P' and indepDict[stateVar][3] is not None:
            dispValue = str(indepDict[stateVar][3])
        else:
            dispValue = str(indepDict[stateVar][0])
        
        nameUnitStr = indepDict[stateVar][1]+' ['+stateUnit+']'
        valueUnitStr = dispValue+' '+stateUnit
        fullStateStr = indepDict[stateVar][1]+' = '+dispValue+' '+stateUnit
            
        return nameUnitStr,valueUnitStr,fullStateStr

def write_shock(shocks,outfile):
    """
    Defines how the shock speeds are formatted in output file
    """
    outfile.write('## Shock speeds:\n\n')
    outfile.write('Primary shock Mach: '+str(shocks[0])+'\n')
    outfile.write('Secondary shock Mach: '+str(shocks[1])+'\n')
    outfile.write('Primary shock velocity: '+str(shocks[2])+' m/s\n')
    outfile.write('Secondary shock velocity: '+str(shocks[3])+' m/s\n')
    outfile.write('\n*************************\n*************************\n')
    return

def write_modes(mode,term,ssr,outfile):
    """
    Defines how the operating modes are formatted in output file
    """
    outfile.write('## Operating modes:\n\n')
    outfile.write('Enthalpy mode: '+mode+'\n')
    outfile.write('Sound speed ratio a3/a2: '+str(ssr)+'\n')
    outfile.write('Test time terminated by: '+term+'\n')
    outfile.write('\n*************************\n*************************\n')
    return

def write_key(gas7,test_time,outfile,real):
    """
    Defines the key test parameters and formats them for output file
    """
    outfile.write('## Test section freestream conditions:\n\n')
    outfile.write('Test time: '+str(test_time*1e3)+' us\n')
    
    # attribute names for PerfectGas object
    attribute = ['P','T','vel','M','h0'] 
    # corresponding descriptions and units for each attribute
    attrDesc = ['Pressure','Temperature','Velocity','Mach number',
                'Total enthalpy']
    attrUnit = ['kPa','K','m/s','','MJ/kg']
    
    if real:
        attribute.extend(['Rem','thermal_conductivity'])
        attrDesc.extend(['Re/m','Thermal conductivity'])
        attrUnit.extend(['1/m','W/m.K'])
        conv = [1e3,1.,1.,1.,1e6,1.,1.]
    else:
        conv = [1.,1.,1.,1.,1e3]
    
    for i,attr in enumerate(attribute):
            outfile.write(attrDesc[i]+': '+str(getattr(gas7,attr)/conv[i])+
            ' '+attrUnit[i]+'\n')
    outfile.write('\n*************************\n*************************\n')
    return

def write_facility(indepDict,gases,lengths,outfile,modes,P_units,detMix):
    """
    Defines how the facility initial setup is formatted in output file
    """
    
    # modes[0] = True/False = Real/Perfect
    # modes[1] = True/False = Detonation/Passive
    if modes[0]:
        conv = 1000.
    else:
        conv = 1.
        
    outfile.write('## Facility initial conditions:\n\n')
    
    if modes[1]:
        outfile.write('Driver mode: Detonation\n\n')
        accelGas = gases[5]
    else:
        outfile.write('Driver mode: Passive\n\n')
        accelGas = gases[4]
        
    
    outfile.write('Driver section:\n')
    outfile.write('Gas: '+str(gases[3].gas)+'\n')
    if modes[1]:
        outfile.write('Initial mole ratio: '+detMix[0]+'\n')
        outfile.write('Initial equivalence ratio: '+str(detMix[1])+'\n')
    outfile.write('Pressure: '+str(gases[3].P/conv)+' kPa')
    if P_units[0]:
        _,rawDrivP,__ = stateVar_formatted('drivP',indepDict,P_units)
        outfile.write('\t('+rawDrivP+')')
    outfile.write('\n')
    outfile.write('Temperature: '+str(gases[3].T)+' K\n')
    outfile.write('Length: '+str(lengths[0])+' m\n\n')
    
    outfile.write('Test section:\n')
    outfile.write('Gas: '+str(gases[0].gas)+'\n')
    outfile.write('Pressure: '+str(gases[0].P/conv)+' kPa')
    if P_units[1]:
        _,rawTestP,__ = stateVar_formatted('testP',indepDict,P_units)
        outfile.write('\t('+rawTestP+')')
    outfile.write('\n')
    outfile.write('Temperature: '+str(gases[0].T)+' K\n')
    outfile.write('Length: '+str(lengths[1])+' m\n\n')
    
    
    outfile.write('Acceleration section:\n')
    outfile.write('Gas: '+str(accelGas.gas)+'\n')
    outfile.write('Pressure: '+str(accelGas.P/conv)+' kPa')
    if P_units[2]:
        _,rawAccelP,__ = stateVar_formatted('accelP',indepDict,P_units)
        outfile.write('\t('+rawAccelP+')')
    outfile.write('\n')
    outfile.write('Temperature: '+str(accelGas.T)+' K\n')
    outfile.write('Length: '+str(lengths[2])+' m')
    
    outfile.write('\n*************************\n*************************\n')
    return

def write_perf_states(pGases,outfile,gasDesc):
    """
    Writes out all perfect gas states (7 or 8 states for passive/det mode)
    """
    # attribute names for PerfectGas object
    attribute = ['gas','P','T','density','vel','M','a','h','gamma'] 
    # corresponding descriptions and units for each attribute
    attrDesc = ['Species','Pressure','Temperature','Density','Velocity','Mach number',
                'Sound speed','Enthalpy','Gamma']
    attrUnit = ['','kPa','K','kg/m^3','m/s','','m/s','kJ/kg','']
    
    if len(pGases) == 7: # Passive
        gasNumber = ['1','2','3','4','5','6','7']
    elif len(pGases) == 8: # Detonation
        gasNumber = ['1','2','3','4','4a','5','6','7']
    
    outfile.write('## Perfect gas states:\n\n')
    for i,gas in enumerate(pGases):
        outfile.write('Gas '+gasNumber[i]+' ('+gasDesc[i]+'):\n\n')
        for j,attr in enumerate(attribute):
            outfile.write(attrDesc[j]+': '+str(getattr(gas,attr))+
            ' '+attrUnit[j]+'\n')
        outfile.write('\n*************************\n')
    
    outfile.write('*************************\n')
        
    return

def write_real_states(rGases,outfile,gasDesc): 
    """
    Writes out all real gas states (7 or 8 states for passive/det mode).
    
    Some of this is done using the Cantera method "report", which prints the thermodynamic
    state in the familiar tabulated format.
    """
    # attribute names for RealGas object - most included in Cantera report
    attribute = ['gas','vel','M','a'] 
    # corresponding descriptions and units for each attribute
    attrDesc = ['Species','Velocity','Mach number','Sound speed']
    attrUnit = ['','m/s','','m/s']
    
    if len(rGases) == 7: # Passive
        gasNumber = ['1','2','3','4','5','6','7']
    elif len(rGases) == 8: # Detonation
        gasNumber = ['1','2','3','4','4a','5','6','7']
    
    outfile.write('## Real gas states:\n\n')
    for i,gas in enumerate(rGases):
        outfile.write('Gas '+gasNumber[i]+' ('+gasDesc[i]+'):\n\n')
        for j,attr in enumerate(attribute):
            outfile.write(attrDesc[j]+': '+str(getattr(gas,attr))+
            ' '+attrUnit[j]+'\n')
        outfile.write(gas.report())
        outfile.write('\n*************************\n')
    
    outfile.write('*************************\n')
        
    return

def write_perf_wedge(pGeom,outfile):
    """
    Defines how the results for a perfect-gas flow over a wedge are formatted in output file
    """
    # attribute names for PerfectGas object
    attribute = ['gas','P','T','density','vel','M','a','h','gamma'] 
    # corresponding descriptions and units for each attribute
    attrDesc = ['Species','Pressure','Temperature','Density','Velocity','Mach number',
                'Sound speed','Enthalpy','Gamma']
    attrUnit = ['','kPa','K','kg/m^3','m/s','','m/s','kJ/kg','']
    
    outfile.write('## Flow over wedge in freestream:\n\n')
    outfile.write('Wedge angle: '+str(pGeom[4])+' deg\n')
    outfile.write('Angle of wedge for detachment: '+str(pGeom[1])+' deg\n')
    if pGeom[0] is None:
        outfile.write('Shock detached.\n\n')  
    else:
        outfile.write('Shock attached.\n')
        outfile.write('Oblique shock angle: '+str(pGeom[0])+' deg\n\n')
        outfile.write('Gas state behind oblique shock:\n\n')
        for j,attr in enumerate(attribute):
                outfile.write(attrDesc[j]+': '+str(getattr(pGeom[2],attr))+
                ' '+attrUnit[j]+'\n')
    outfile.write('\n*************************\n')
    
def write_real_wedge(rGeom,outfile):
    """
    Defines how the results for a real-gas flow over a wedge are formatted in output file
    """
    attribute = ['gas','vel','M','a'] 
    attrDesc = ['Species','Velocity','Mach number','Sound speed']
    attrUnit = ['','m/s','','m/s']
    
    outfile.write('## Flow over wedge in freestream:\n\n')
    outfile.write('Wedge angle: '+str(rGeom[4])+' deg\n')
    outfile.write('Angle of wedge for detachment: '+str(rGeom[1])+' deg\n')
    if rGeom[0] is None:
        outfile.write('Shock detached.\n\n')  
    else:
        outfile.write('Shock attached.\n')
        outfile.write('Oblique shock angle: '+str(rGeom[0])+' deg\n\n')
        outfile.write('Gas state behind oblique shock:\n\n')
        for j,attr in enumerate(attribute):
            outfile.write(attrDesc[j]+': '+str(getattr(rGeom[2],attr))+
            ' '+attrUnit[j]+'\n')
        outfile.write(rGeom[2].report())
    outfile.write('\n*************************\n')
      
##############################################################################
# Simple utility functions  
def torr2kpa(torr):
    """
    Converts torr to kPa
    """
    from cantera import one_atm
    kpa = one_atm/1000/760*torr
    return kpa

    
def bar2kpa(bar):
    """
    Converts bar to kPa
    """
    kpa = 100*bar
    return kpa

##############################################################################
# Testing functions
    
# Define temporary functions here to be called from UI when Calculate is pressed. In this way,
# can use all the functionality of the UI and test the output of some new function on a defined
# run condition.
    
def BL_Re_lm(all_data,allL):
    from exp_real_lib import boundary_layer
    import numpy as np
    gas5 = all_data['Real gases'][4]
    gas6 = all_data['Real gases'][5]
    us2 = all_data['Real shock'][3]
    t_contact2 = all_data['Contact arrival'][1]/1e3
    accelL = allL[2]
    
    # Use ideal separation of secondary shock and secondary contact surface
    # at the time when the C.S. reaches the end of the acceleration tube
    # as an upper bound on the separation distance
    lm_ideal = us2*t_contact2-accelL
    print('Ideal separation length: '+str(round(lm_ideal,3))+' m')

    # Use Mirels' expression for Re_lm
    Re_lm = boundary_layer(gas5,gas6,us2,lm_ideal)
    print('Re_lm: '+str(round(Re_lm/1e6,3))+' x10^6')
    
    # Assume transition Re varies linearly on a loglog plot of Ms vs. Re_t between
    # 1 < Ms < 9, as per Mirels. Significant increase seen experimentally at higher Ms,
    # to values > 10^7
    Ms2 = us2/gas5.a
    Ms_bounds = np.log10(np.array([1,9]))
    Re_t_bounds = np.log10(np.array([0.5,4])*1e6)
    
    if (np.log10(Ms2)>Ms_bounds[0]) and (np.log10(Ms2)<Ms_bounds[1]):
        logRe_t = (np.log10(Ms2)-Ms_bounds[0])/(Ms_bounds[1]-Ms_bounds[0])*(Re_t_bounds[1]-Re_t_bounds[0])+Re_t_bounds[0]
        Re_t = 10**(logRe_t)
        print('Re_t: '+str(round(Re_t/1e6,3))+' x10^6')
    else:
        print('Re_t ~ 10^7')
    
    
    return
    

