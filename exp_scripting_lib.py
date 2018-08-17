#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:46:41 2018

@author: jmlawson
"""
import exp_perf_lib as epl
import exp_real_lib as erl
import exp_common_lib as ecl
import numpy as np
from matplotlib.figure import Figure
import os, sys
      
def single_solve(inputDict,silence=False):  
    
    # Optionally silence output
    old_target = sys.stdout
    new_target = open(os.devnull, "w")
    if silence:
        sys.stdout = new_target
    
    
    gasDesc_passive = ['Initial test gas','Post-shock test gas','Post-exp. driver gas',
                            'Initial driver gas','Initial accel. gas',
                            'Post-shock accel. gas','Final test gas']
    
    gasDesc_det = ['Initial test gas','Post-shock test gas','Post-exp. driver gas',
                        'Pre-det. driver','Post-det. driver','Initial accel. gas',
                        'Post-shock accel. gas','Final test gas']
    
    # Unpack:
    drivGasStr,testGasStr,accelGasStr = inputDict['gas']
    detCompStr = inputDict['detComp']
    drivP,testP,accelP = inputDict['P']
    drivT,testT,accelT = inputDict['T']
    drivL,testL,accelL = inputDict['L']
    drivu,testu,accelu = inputDict['u']
    passive_df,det_df = inputDict['df']
    primShock = inputDict['primShock']
    geom = inputDict['geom']
    
    if inputDict['gasModel'] == 'Perfect':
        perfect = True
    elif inputDict['gasModel'] == 'Real':
        perfect = False
    else:
        print('Invalid gas model: must be ''Perfect'' or ''Real''')
        return
    
    
    if inputDict['driverMode'] == 'Passive':
        passive = True
    elif inputDict['driverMode'] == 'Detonation':
        passive = False
    else:
        print('Invalid driver mode: must be ''Passive'' or ''Detonation''')
        return
        
    if inputDict['driverSpec'] == 'Pressure':
        pressureSpecified = True
    elif inputDict['driverSpec'] == 'ShockMach':
        pressureSpecified = False
        mach = True
    elif inputDict['driverSpec'] == 'ShockSpeed':
        pressureSpecified = False
        mach = False
    else:
        print('Invalid driver specified input: must be ''Pressure'',''ShockMach'', or ''ShockSpeed''')
    
    
    ####################################################################
    # Create objects for perfect driver,test and acceleration gases
    ####################################################################
    
    # If operating in detonation mode, pDriverGas and rDriverGas represent the
    # post-detonation gas, i.e. state 4a.
    # In this case, we need to create a RealGas object representing the pre-detonation
    # mixture, then use SDT to solve for the post-detonation state. From this, we
    # can create an approximate PerfectGas object with fixed gamma etc.
    try:
        if pressureSpecified and passive:
            pDriverGas = epl.PerfectGas(drivGasStr,drivP,drivT,drivu,passive_df)
        elif pressureSpecified and not passive:
            print('Solving for post-detonation state...')
            preDetState = (det_df.loc[drivGasStr,'Cantera'],drivP*1e3,drivT,detCompStr,drivGasStr) # mech,P,T,comp,name
            pDriverGas,rDriverGas,pPreDetGas,rPreDetGas,cj_speed = erl.detonation_setup(preDetState)              
        else:
            # If in shock-specified mode, set an arbitrary driver pressure,
            # in case the pressure entry on other tab was invalid/blank
            pDriverGas = epl.PerfectGas(drivGasStr,1.,drivT,drivu,passive_df)
        
        pTestGas = epl.PerfectGas(testGasStr,testP,testT,testu,passive_df)
        pAccelGas = epl.PerfectGas(accelGasStr,accelP,accelT,accelu,passive_df)
    except AttributeError:
        dataErrormsg = 'Gas data file format incorrect.'
        print(dataErrormsg)
        return
    
    ####################################################################
    # Begin the solving process for perfect gas
    ####################################################################
    if pressureSpecified: # Driver-specified mode
        # First polar intersection: use states 1 & 4 to find 2 & 3
        Ms1_perf,pGas2,pGas3,pdxdt1 = epl.polar_solve(pTestGas,pDriverGas,passive_df)
        us1_perf = Ms1_perf*pTestGas.a
        solvedDrivP_perf = pDriverGas.P
    else: # Shock-specified mode
        if mach:
            Ms1_perf = primShock
            us1_perf = Ms1_perf*pTestGas.a
        else:
            Ms1_perf = primShock/pTestGas.a
            us1_perf = primShock
        pGas2 = epl.shock_jump(pTestGas,Ms1_perf) # fully defined State 2
        pGas3,pdxdt1 = epl.unsteady_exp(pDriverGas,pGas2.vel,pGas2.P,passive_df) # fully defined State 3
        solvedDrivP_perf = epl.polar_P(pGas3.P,pGas3.vel,pDriverGas.vel,pGas3.gamma,pGas3.a) # kPa
        pDriverGas.P = solvedDrivP_perf
    
    
    # Second polar intersection: use states 2 & 5 to find 6 & 7
    Ms2_perf,pGas6,pGas7,pdxdt2 = epl.polar_solve(pAccelGas,pGas2,passive_df)
    us2_perf = Ms2_perf*pAccelGas.a
    
    ####################################################################
    # If required, create objects for corresponding real gases and solve real system
    ####################################################################
    if not perfect:
        conv = 1e3
        
        # Creation. Note Cantera requires Pa, not kPa
        if pressureSpecified and passive:
            # Already created above for detonation
            rDriverGas = erl.RealGas(passive_df.loc[drivGasStr,'Cantera'],vel=drivu,T=drivT,
                                     P=drivP*conv,X=passive_df.loc[drivGasStr,'X'],gas=drivGasStr,equil=True)
        elif not pressureSpecified:
            # For real case, don't use arbitrary pressure. Instead, use solved driver pressure from perfect case above
            rDriverGas = erl.RealGas(passive_df.loc[drivGasStr,'Cantera'],vel=drivu,T=drivT,
                                     P=solvedDrivP_perf*conv,X=passive_df.loc[drivGasStr,'X'],gas=drivGasStr,equil=True)    
        
        rTestGas = erl.RealGas(passive_df.loc[testGasStr,'Cantera'],vel=testu,T=testT,
                                 P=testP*conv,X=passive_df.loc[testGasStr,'X'],gas=testGasStr,equil=True)
        rAccelGas = erl.RealGas(passive_df.loc[accelGasStr,'Cantera'],vel=accelu,T=accelT,
                                 P=accelP*conv,X=passive_df.loc[accelGasStr,'X'],gas=accelGasStr,equil=True)
       
        # Initialize states 2, 3 and 6 as copies of 1, 4, and 5 respectively
        if passive:
            rGas3 = erl.RealGas(passive_df.loc[drivGasStr,'Cantera'],vel=rDriverGas.vel,T=rDriverGas.T,
                                 P=rDriverGas.P,X=rDriverGas.X,gas=rDriverGas.gas,equil=True)
        else:
            rGas3 = erl.RealGas(det_df.loc[drivGasStr,'Cantera'],vel=rDriverGas.vel,T=rDriverGas.T,
                                 P=rDriverGas.P,X=rDriverGas.X,gas=rDriverGas.gas,equil=True)   
                    
        rGas2 = erl.RealGas(passive_df.loc[testGasStr,'Cantera'],vel=testu,T=testT,
                                 P=testP*conv,X=passive_df.loc[testGasStr,'X'],gas=testGasStr,equil=True)
        rGas6 = erl.RealGas(passive_df.loc[accelGasStr,'Cantera'],vel=accelu,T=accelT,
                                 P=accelP*conv,X=passive_df.loc[accelGasStr,'X'],gas=accelGasStr,equil=True)
        
        print('Gases 1-6 initialized!')
        # Solving            
        if pressureSpecified:
            try:
                print('Beginning first shocktube problem...')
                Ms1_real, rdxdt1 = erl.real_polar_solve(rTestGas,rGas2,rGas3,rDriverGas,Ms1_perf)  
                print('First shocktube problem completed!')
            except RuntimeError as ctError:
                ctErrormsg = '{0}'.format(ctError)
                # Format error message by removing * and \n so it fits in UI
                ctErrormsg = ctErrormsg.translate(None,'*\n') 
                print(ctErrormsg+' in first shocktube problem')
                return
            except ValueError as valError:
                print(valError.args[0]+' in first shocktube problem')

                return
            us1_real = Ms1_real*rTestGas.a
            
        else:
            if mach:
                Ms1_real = primShock
                us1_real = Ms1_real*rTestGas.a
            else:
                Ms1_real = primShock/rTestGas.a
                us1_real = primShock
            # Use the perfect driver gas solution as an initial guess for the real solution
            solvedDrivP_real,rdxdt1 = erl.real_polar_solve_inv(rTestGas,rGas2,rGas3,rDriverGas,
                                                               Ms1_real,solvedDrivP_perf*conv)
            
        # Initialize state 7 as a copy of state 2
        rGas7 = erl.RealGas(passive_df.loc[testGasStr,'Cantera'],vel=rGas2.vel,T=rGas2.T,
                                 P=rGas2.P,X=rGas2.X,gas=testGasStr,equil=True)
        print('Gas 7 initialized!')
        
        try:
            print('Beginning second shocktube problem...')
            Ms2_real, rdxdt2 = erl.real_polar_solve(rAccelGas,rGas6,rGas7,rGas2,Ms2_perf)
            print('Second shocktube problem completed!')
        except RuntimeError as ctError:
            ctErrormsg = '{0}'.format(ctError)
            # Format Cantera error message by removing * and \n so it fits in UI
            ctErrormsg = ctErrormsg.translate(None,'*\n') 
            print(ctErrormsg+' in second shocktube problem')
            return
        except ValueError as valError:
            print(valError.args[0]+' in second shocktube problem')
            return
        us2_real = Ms2_real*rAccelGas.a
    
    
    # Postprocessing (much of this could be spun out into separate methods)
    ####################################################################
    
    # Calculate sound speed ratios and enthalpy mode
    ssr_perf = pGas3.a/pGas2.a 
    if ssr_perf < 1:
        mode_perf = 'High'
    else:
        mode_perf = 'Low'           
    ssr = [ssr_perf]
    mode = [mode_perf]
    
    if not perfect:
        rGas7.reynolds()
        ssr_real = rGas3.a/rGas2.a # sound speed ratio
        ssr.append(ssr_real)
        if ssr_real < 1:
            mode_real = 'High'
        else:
            mode_real = 'Low'
        mode.append(mode_real)
    
    # Identify regions of min/max temperature in facility
    if passive:
        pTarray = np.asarray([pTestGas.T,pGas2.T,pGas3.T,pDriverGas.T,pAccelGas.T,pGas6.T,pGas7.T])
        if not perfect:
            rTarray = np.asarray([rTestGas.T,rGas2.T,rGas3.T,rDriverGas.T,rAccelGas.T,rGas6.T,rGas7.T])
        gasDesc = gasDesc_passive
    elif not passive:
        pTarray = np.asarray([pTestGas.T,pGas2.T,pGas3.T,rPreDetGas.T,pDriverGas.T,pAccelGas.T,pGas6.T,pGas7.T])
        if not perfect:
            rTarray = np.asarray([rTestGas.T,rGas2.T,rGas3.T,rPreDetGas.T,rDriverGas.T,rAccelGas.T,rGas6.T,rGas7.T])
        gasDesc = gasDesc_det
   
    pMinTidx = pTarray.argmin()
    pMaxTidx = pTarray.argmax()        
    pExtremeT = (pTarray[pMinTidx],gasDesc[pMinTidx],pTarray[pMaxTidx],gasDesc[pMaxTidx])
    rExtremeT = ()
    if not perfect:          
        rMinTidx = rTarray.argmin()
        rMaxTidx = rTarray.argmax()
        rExtremeT = (rTarray[rMinTidx],gasDesc[rMinTidx],rTarray[rMaxTidx],gasDesc[rMaxTidx])
    
    # Use freestream solution to calculate flow over test body geometry, if required
    # If additional geometries are added, should probably spin this out into own method
    pGeomOut = None
    rGeomOut = None
    if geom[0] == 'Wedge':
        pGeomOut = list(epl.wedge(pGas7,geom[1]))
        pGeomOut.extend(geom)
        if not perfect:
            print('Solving for flow over wedge...')
            rPostOblGas = erl.RealGas(passive_df.loc[testGasStr,'Cantera'],vel=rGas7.vel,T=rGas7.T,
                                 P=rGas7.P,X=rGas7.X,gas=testGasStr,equil=True)
            print('Post-oblique shock gas initialized!')
            rGeomOut = list(erl.wedge(rGas7,geom[1],rPostOblGas))
            rGeomOut.extend([rPostOblGas])
            rGeomOut.extend(geom)
            print('Wedge flow problem complete!')
   
    # Compile outputs for use by UI display
    if passive:
        us_perf = (us1_perf,us2_perf)
#        Ms_perf = (Ms1_perf,Ms2_perf)
    else:
        # Include detonation wave speed and Mach number in these tuples
        # Calculate Mach number based on sound speed in pre-detonation gas mix
        # Use perfect or real versions of this gas as appropriate
        us_perf = (us1_perf,us2_perf,cj_speed)
    pdxdt = (pdxdt1,pdxdt2)
    
    
    if not perfect:
        if passive:
            us_real = (us1_real,us2_real)
        else:
            us_real = (us1_real,us2_real,cj_speed)
        rdxdt = (rdxdt1,rdxdt2)
    
        

    # Generate x-t diagram and display output to UI
    # Generate and display x-t plot, calculate test times
    fig1 = Figure()
    ax1f1 = fig1.add_subplot(111)
    # give the subplot handle to the plotting function
    t_test_perf, term_perf, t_contact2_perf, opt_perf = ecl.exp_tube_xt(us_perf,pdxdt,pGas2,pDriverGas,
                                                                       pGas7.vel,inputDict['L'],ax1f1)
        
    t_test = [t_test_perf]
    term = [term_perf]
    t_contact2 = [t_contact2_perf]
    optimal = [opt_perf]
        
    if not perfect:
        fig2 = Figure()
        ax1f2 = fig2.add_subplot(111)
        rGas2.gamma = rGas2.cp/rGas2.cv # Needed for similarity soln
        rDriverGas.gamma = rDriverGas.cp/rDriverGas.cv
        t_test_real, term_real, t_contact2_real, opt_real = ecl.exp_tube_xt(us_real,rdxdt,rGas2,rDriverGas,
                                                                            rGas7.vel,inputDict['L'],ax1f2)
        t_test.append(t_test_real)
        term.append(term_real)
        t_contact2.append(t_contact2_real)
        optimal.append(opt_real)
   
    # Compile complete output for saving to file
    all_data = {'Test time':t_test,
                'Enthalpy mode':mode,
                'Sound speed ratio':ssr,
                'Perfect shock':(Ms1_perf,Ms2_perf,us1_perf,us2_perf),
                'Termination mode':term,
                'Lengths':inputDict['L'],
                'Contact arrival':t_contact2,
                'Extreme T':(pExtremeT,rExtremeT),
                'Perfect dxdt':pdxdt,
                'Optimal test time':optimal
                }
    
    # Add extra data depending on calculation mode
    if passive:
        all_data['Perfect gases'] = (pTestGas,pGas2,pGas3,pDriverGas,pAccelGas,pGas6,pGas7)
        if not perfect:           
            all_data['Real gases'] = (rTestGas,rGas2,rGas3,rDriverGas,rAccelGas,rGas6,rGas7)
            all_data['Real shock'] = (Ms1_real,Ms2_real,us1_real,us2_real)
            all_data['Real dxdt'] = rdxdt
    elif not passive:
        all_data['Perfect gases'] = (pTestGas,pGas2,pGas3,pPreDetGas,pDriverGas,pAccelGas,pGas6,pGas7)
        all_data['Detonation'] = cj_speed
        if not perfect:           
            all_data['Real gases'] = (rTestGas,rGas2,rGas3,rPreDetGas,rDriverGas,rAccelGas,rGas6,rGas7)
            all_data['Real shock'] = (Ms1_real,Ms2_real,us1_real,us2_real)
            all_data['Real dxdt'] = rdxdt
     
    
    if geom[0] == 'Wedge':
        all_data['Perfect test body'] = pGeomOut
        if not perfect:
            all_data['Real test body'] = rGeomOut
            
    
    # Restore output target
    sys.stdout = old_target
        
    return all_data




def double_contour(ax,d):
    formattingReplacements = {'CO2':r'CO$_2$',
                              'C2H2':r'C$_2$H$_2$'}
    
    
    x = np.asarray(d['x'])*d['xscale']; y = np.asarray(d['y'])*d['yscale']
    data1 = np.asarray(d['data1']); data2 = np.asarray(d['data2'])

    anchor = {'ha': 'left', 'va': 'bottom'}
    
    if d['mainMarker'] is not None:
        ax.plot(x,y,d['mainMarker'],markersize=d['mainSize'],color=d['markerColour'])
    
    # If the selected step size would miss the last entry, include it anyway.
    # However, don't append it if it was already included because this can cause the plotted
    # line and text to be doubled up and appear thicker.
    if len(d['contour1']) % d['step1'] == 1 or d['step1'] == 1:
        cont1Show = d['contour1'][::d['step1']]
    else:
        cont1Show = np.concatenate((d['contour1'][::d['step1']],[d['contour1'][-1]]))
        
    if len(d['contour2']) % d['step2'] == 1 or d['step2'] == 1:
        cont2Show = d['contour2'][::d['step2']]
    else:
        cont2Show = np.concatenate((d['contour2'][::d['step2']],[d['contour2'][-1]]))
    
    if d['topLayer'] == 1:
        z1 = 100; z2 = 1
    elif d['topLayer'] == 2:
        z2 = 100; z1 = 1
    elif d['topLayer'] == 3:
        z1 = d['zLayers'][0]; z2 = d['zLayers'][1]
    
    for cont1 in cont1Show:
        xcont = x[data1==cont1]
        ycont = y[data1==cont1]
        ax.plot(xcont,ycont,color=d['colour1'],linestyle=d['style1'],zorder=z1)
        textx = max(xcont)+d['xOffset1']
        if max(xcont) == min(xcont):
            texty = max(ycont)+d['yOffset1']
        else:
            texty = ycont[xcont.argmax()]+d['yOffset1']
            
        if d['format1'] is not None:
            if cont1 in formattingReplacements:
                cont1 = formattingReplacements[cont1]
            ax.text(textx,texty,d['format1'] % cont1,anchor,
                    color=d['colour1'],rotation=d['angle1'],fontsize=d['fontsize'])
    
    for cont2 in cont2Show:
        xcont = x[data2==cont2]
        ycont = y[data2==cont2]
        ax.plot(xcont,ycont,color=d['colour2'],linestyle=d['style2'],zorder=z2)
        textx = max(xcont)+d['xOffset2']
        if max(xcont) == min(xcont):
            texty = max(ycont)+d['yOffset2']
        else:
            texty = ycont[xcont.argmax()]+d['yOffset2']
            
        if d['format2'] is not None:
            if cont2 in formattingReplacements:
                cont2 = formattingReplacements[cont2]
            ax.text(textx,texty,d['format2'] % cont2,anchor,
                    color=d['colour2'],rotation=d['angle2'],fontsize=d['fontsize'])
    
    if d['cases'] is not None:
        for i,case in enumerate(d['cases']['Names']):
            xcase = d['cases'][d['casex']][i]*d['xscale']
            ycase = d['cases'][d['casey']][i]*d['yscale']
            ax.plot(xcase,ycase,d['caseMarker'],markersize=d['caseSize'])
            ax.text(xcase+d['caseXOffset'],ycase+d['caseYOffset'],case,fontsize=d['caseFontsize'])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(d['xlim'])
    ax.set_ylim(d['ylim'])
    
    ax.set_xlabel(d['xlabel'])
    ax.set_ylabel(d['ylabel'])
    if d['title'] is not None:
        ax.set_title(d['title'])
             
        