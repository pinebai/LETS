# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:08:31 2016

@author: LawsonJoelM
"""
##############################################################################
#** Import external modules
import sys
from PyQt5 import  QtWidgets, uic
from PyQt5.QtGui import QFont
import matplotlib
matplotlib.use('Qt5Agg')
#matplotlib.rcParams['lines.linewidth'] = 0.5 # change default line width in plots
matplotlib.rcdefaults() # restore standard defaults. Comment out if using above line!
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pandas import read_excel
from numpy import linspace,logspace,asarray,sqrt
##############################################################################
#** Import expansion tube libraries
import exp_perf_lib as epl
import exp_real_lib as erl
import exp_common_lib as ecl
import sdtoolbox as sdt

##############################################################################
##** User interface creation
# Load templates for widget layouts
mainTemplate = uic.loadUiType("exp_tube.ui")[0]
mixerTemplate = uic.loadUiType("mixer.ui")[0]


class Main(QtWidgets.QMainWindow, mainTemplate):
    """
    This class defines the main window of the UI.
    The actual layout and the names of the widgets are defined in the file "exp_tube.ui"
    which was generated using QtDesigner.
    
    All functionality of the main windows is encapsulated in the methods of this class.
    These include setup of menus and lists, error-checking inputs, gathering inputs and
    submitting them for calculation, overseeing plot creation and display, and outputting
    results to file.
    
    Upon executing this script from the command line, an instance of this Main class is
    instantiated, and the window containing it will be displayed. This is done at the very
    end of this script.
    """
    # Initialization of UI window
    def __init__(self, parent=None): 
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        
        # Initialize some flags and 'global' (i.e. within the scope of the main UI) variables:
            
        self.existingCombustionMix = None # Used to store mixer settings even after popup closed
        # Flags that indicate 1) if results from a multi-calc are in memory, 2) if it had test body data
        self.activeMulti = [False,False]
        
        # Define the list of dependent variables for multi-plot mode
        # Format = key:value
        # key = 'name to appear in menu'
        # value = [ [pKey,pIndex,pAttr,pConv], [rKey,rIndex,rAttr,rConv], 'unit']
        # pKey = key in all_data dict
        # pIndex = index in corresponding value from all_data dict
        # pAttr = attribute if the selected index is a gas object. 'None' otherwise.
        # pConv = conversion factor (for dimensional consistency) (raw value gets divided by this)
        # Note: use index -1 for freestream gas, since this gas is always the last object in the list
        # -- but there are different length lists for passive vs det
        self.depVarDict = {
                        'Freestream Mach number':[['Perfect gases',-1,'M',1],['Real gases',-1,'M',1],''],
                        'Freestream velocity':[['Perfect gases',-1,'vel',1],['Real gases',-1,'vel',1],'m/s'],
                        'Freestream pressure':[['Perfect gases',-1,'P',1],['Real gases',-1,'P',1e3],'kPa'],
                        'Freestream density':[['Perfect gases',-1,'density',1],['Real gases',-1,'density',1],'kg/m^3'],
                        'Freestream temperature':[['Perfect gases',-1,'T',1],['Real gases',-1,'T',1],'K'],
                        'Freestream total enthalpy':[['Perfect gases',-1,'h0',1e3],['Real gases',-1,'h0',1e6],'MJ/kg'],
                        'Test time':[['Test time',0,None,1e-3],['Test time',1,None,1e-3],'us'],
                        'Primary shock Mach number':[['Perfect shock',0,None,1],['Real shock',0,None,1],''],
                        'Secondary shock Mach number':[['Perfect shock',1,None,1],['Real shock',1,None,1],''],
                        'Primary shock speed':[['Perfect shock',2,None,1],['Real shock',2,None,1],'m/s'],
                        'Secondary shock speed':[['Perfect shock',3,None,1],['Real shock',3,None,1],'m/s'],
                        'Sound speed ratio a3/a2':[['Sound speed ratio',0,None,1],['Sound speed ratio',1,None,1],''],
                        'Oblique shock angle':[['Perfect test body',0,None,1],['Real test body',0,None,1],'deg'],
                        'Detachment angle':[['Perfect test body',1,None,1],['Real test body',1,None,1],'deg']
                          }
        
        # The subset of the above dependent variables that are unique to wedge geometry
        self.wedgeDepVars = ['Oblique shock angle','Detachment angle']
        
        # Formatted names of each gas state, for use in output files. Order follows the convention of
        # numbering the gases 1 thru 7 (with 4a for the post-detonation driver)
        self.gasDesc_passive = ['Initial test gas','Post-shock test gas','Post-exp. driver gas',
                                'Initial driver gas','Initial accel. gas',
                                'Post-shock accel. gas','Final test gas']
        
        self.gasDesc_det = ['Initial test gas','Post-shock test gas','Post-exp. driver gas',
                            'Pre-det. driver','Post-det. driver','Initial accel. gas',
                            'Post-shock accel. gas','Final test gas']
        
        # The available types of geometry to be placed in the freestream (test section)
        self.geomList = ['None','Wedge','Normal shock']
        

        # Bind the widgets in the UI to methods defined in this library        
        self.calc_PushBtn.clicked.connect(self.calculate)  
        self.import_PushBtn.clicked.connect(self.menu_setup)
        self.gasbrowse_PushBtn.clicked.connect(self.input_file)
        self.outbrowse_PushBtn.clicked.connect(self.output_folder)
        self.calcType_Tabs.currentChanged.connect(self.switch_calcMode)
        self.indepVar_ButtonGroup.buttonClicked.connect(self.disable_indep)
        self.plotPerf_CheckBox.clicked.connect(self.multi_plot)
        self.depVar_list.currentItemChanged.connect(self.multi_plot)
        self.abort_PushBtn.clicked.connect(self.abort)
        self.reset_PushBtn.clicked.connect(self.default_setup)
        self.driverMode_ButtonGroup.buttonClicked.connect(self.switch_driverMode)
        self.mixer_PushBtn.clicked.connect(self.open_mixer)
        
        # The calculate button and data inputs are greyed out until gases are imported
        self.inputFrame.setEnabled(False)
        # Some colours need changing manually (they don't grey out even when disabled)
        self.drivunits_GroupBox.setStyleSheet('QGroupBox {color:gray}')
        self.testunits_GroupBox.setStyleSheet('QGroupBox {color:gray}')
        self.accelunits_GroupBox.setStyleSheet('QGroupBox {color:gray}')
        self.calc_PushBtn.setEnabled(False)
        self.abort_PushBtn.setEnabled(False)
        # Set the error message font colour to red
        self.msg_Display.setStyleSheet('QLabel {color:red}')
        importErrorMsg = 'Please import gas data first'
        self.msg_Display.setText(importErrorMsg)
        
        # Bind the editing fields and radio button groups to the top-level checking function 
        # so everything is re-checked everytime something is edited
        self.drivP_Edit.editingFinished.connect(self.all_check)
        self.testP_Edit.editingFinished.connect(self.all_check)
        self.accelP_Edit.editingFinished.connect(self.all_check)
        self.drivT_Edit.editingFinished.connect(self.all_check)
        self.testT_Edit.editingFinished.connect(self.all_check)
        self.accelT_Edit.editingFinished.connect(self.all_check)
        self.drivL_Edit.editingFinished.connect(self.all_check)
        self.testL_Edit.editingFinished.connect(self.all_check)
        self.accelL_Edit.editingFinished.connect(self.all_check)
        
        self.angle_Edit.editingFinished.connect(self.all_check)
        
        self.driver_Tabs.currentChanged.connect(self.all_check)
        self.primShock_Edit.editingFinished.connect(self.all_check)
        
        self.rangeStart_Edit.editingFinished.connect(self.all_check)
        self.rangeStop_Edit.editingFinished.connect(self.all_check)
        self.rangeNum_Edit.editingFinished.connect(self.all_check)
        
        self.indepVar_ButtonGroup.buttonClicked.connect(self.all_check)
        self.rangeType_ButtonGroup.buttonClicked.connect(self.all_check)
        self.drivunits_ButtonGroup.buttonClicked.connect(self.all_check)
        self.shockunits_ButtonGroup.buttonClicked.connect(self.all_check)
        self.testunits_ButtonGroup.buttonClicked.connect(self.all_check)
        self.accelunits_ButtonGroup.buttonClicked.connect(self.all_check)
        
        self.passive_PushBtn.clicked.connect(self.all_check)
        self.detonation_PushBtn.clicked.connect(self.all_check)
        
        self.depVar_list.itemSelectionChanged.connect(self.all_check)
        
        self.drivGas_Dropdown.currentIndexChanged.connect(self.all_check)
        self.testGas_Dropdown.currentIndexChanged.connect(self.all_check)
        self.accelGas_Dropdown.currentIndexChanged.connect(self.all_check)
        
        # Whenever a new geometry is selected for the test body, update the list of 
        # dependent variables in multi-calc mode (i.e. remove/add any dep vars unique
        # to that geometry)
        self.geom_Dropdown.currentIndexChanged.connect(self.setup_depVar)
        # When detonation mixture species changed, erase any stored mixture
        # to force them to open the mixer again (to avoid carrying forward the previous species)
        self.drivGas_Dropdown.currentIndexChanged.connect(self.detGas_changed)
        
        fig = Figure()
        self.add_mpl(fig) # create an empty figure upon initialization
        # otherwise will throw an error when it tries to clear previous figure
        self.clear_output() # get rid of the placeholder output values in .ui file
        
        # Populate various menus and lists
        self.default_setup()
        self.setup_geom()
        self.setup_depVar()
        

 
##############################################################################
#** User interface callbacks
# Populate and enable/disable menus, fields
    def default_setup(self):
        """
        Populate various fields with default values from external file.
        
        Called upon initialization. Can also be called again by pressing
        the "Reset" button.        
        """
        default_filename = 'Data/default.xlsx'
        df = read_excel(open(default_filename,'rb'), sheet_name=0)        
        self.inFile_Edit.setText(df.loc['InputFile'].Value)
        self.outFolder_Edit.setText(df.loc['OutputFolder'].Value)
        self.outName_Edit.setText(df.loc['OutputFile'].Value)
        self.drivP_Edit.setText(str(df.loc['DriverP'].Value))
        self.drivT_Edit.setText(str(df.loc['DriverT'].Value))
        self.drivL_Edit.setText(str(df.loc['DriverL'].Value))
        self.testP_Edit.setText(str(df.loc['DrivenP'].Value))
        self.testT_Edit.setText(str(df.loc['DrivenT'].Value))
        self.testL_Edit.setText(str(df.loc['DrivenL'].Value))
        self.accelP_Edit.setText(str(df.loc['AccelP'].Value))
        self.accelT_Edit.setText(str(df.loc['AccelT'].Value))
        self.accelL_Edit.setText(str(df.loc['AccelL'].Value))
  
        
    def menu_setup(self):
        """
        Callback for pressing the "Import Data" button. Reads in gas names from
        the provided filename and uses these to populate the dropdown menus for
        each section. Also enables the rest of the dialog.
        """
        gas_filename = str(self.inFile_Edit.text())
        try:
            # Read in both lists of gases from the Passive and Detonation sheets
            # Store as pandas dataframes, globally available to UI
            self.passive_df = read_excel(open(gas_filename,'rb'), sheet_name='Passive')
            self.det_df = read_excel(open(gas_filename,'rb'), sheet_name='Detonation')
        except IOError:
            inputErrormsg = 'Incorrect gas data filename.'
            self.msg_Display.setText(inputErrormsg)
            return
        
        # Create a list of the gas names (note: behaviour of map() changed from Py2 to Py3)
        self.passive_gas_names = list(map(str,self.passive_df.index.get_values()))
        self.det_gas_names = list(map(str,self.det_df.index.get_values()))
        
        self.drivGas_Dropdown.clear()
        # Populate the driver gas dropdown depending on active mode
        if self.passive_PushBtn.isChecked():
            self.drivGas_Dropdown.addItems(self.passive_gas_names)
        else:
            self.drivGas_Dropdown.addItems(self.det_gas_names)
        
        # Populate the test and acceleration gas dropdowns with passive gases
        self.testGas_Dropdown.clear()
        self.testGas_Dropdown.addItems(self.passive_gas_names)
        self.accelGas_Dropdown.clear()
        self.accelGas_Dropdown.addItems(self.passive_gas_names)
        
        # Enable the input panel
        self.calc_PushBtn.setEnabled(True)
        self.inputFrame.setEnabled(True)
        # Need to manually change the radio button group title colour
        self.drivunits_GroupBox.setStyleSheet('QGroupBox {color:black}')
        self.testunits_GroupBox.setStyleSheet('QGroupBox {color:black}')
        self.accelunits_GroupBox.setStyleSheet('QGroupBox {color:black}')
        self.msg_Display.setText('')
        
        self.all_check()
        self.enable_lineEdits()
        self.disable_indep()
        return

    def setup_geom(self):
        """
        Populates the test body geometry dropdown menu.
        """
        self.geom_Dropdown.clear()
        self.geom_Dropdown.addItems(self.geomList)
    
    def setup_depVar(self):
        """
        Populates the list of dependent variables for multi-calc mode.
        Makes sure only dependent variables compatible with currently-selected
        test body geometry are displayed.
        """
        # Temporarily block signals since otherwise tries to re-draw plot since
        # removal and addition of list items emits a signal that the currently-selected
        # item was changed
        self.depVar_list.blockSignals(True)
        self.depVar_list.clear()
        for key in self.depVarDict:
            if key in self.wedgeDepVars and self.geom_Dropdown.currentText() == 'Wedge':
                pass
            else:
                self.depVar_list.addItem(key)
        # Items in the list get sorted alphabetically
        self.depVar_list.sortItems()
        self.depVar_list.blockSignals(False)
        self.all_check()
        return
    
    def disable_indep(self):
        """
        Disables the input field for the numerical value of 1 of the 9 variables
        (pressure, temperature, length for driver, test and acceleration section)
        that are available as independent variables in multi-calc mode.
        """
        if self.inputFrame.isEnabled() and self.calcType_Tabs.currentIndex() == 1:
            self.enable_lineEdits()
            # trim the button name to obtain the name of the independent variable
            indepStr = self.indepVar_ButtonGroup.checkedButton().objectName()[5:-6]
            # find the corresponding QLineEdit with matching name
            indepLineEdit = self.findChild(QtWidgets.QLineEdit,str(indepStr+'_Edit'))
            indepLineEdit.setEnabled(False)
        return
    
    def enable_lineEdits(self):
        """
        Enables all 9 of the input fields for the independent variables.
        Used when switching back from multi-calc to single-calc mode.
        """
        # Obtain a list of all 9 line edits in inputFrame
        allLineEdits = self.inputFrame.findChildren(QtWidgets.QLineEdit)
        for lineEdit in allLineEdits:
            lineEdit.setEnabled(True)
        return
    
    def switch_calcMode(self):
        """
        Oversees enabling/disabling of various features when switching between
        single and multi-calc modes.
        """
        if self.calcType_Tabs.currentIndex() == 0 and self.inputFrame.isEnabled(): # single
            self.enable_lineEdits()
            self.all_check()
            self.abort_PushBtn.setEnabled(False)
            self.driver_Tabs.setTabEnabled(1,True)
        if self.calcType_Tabs.currentIndex() == 1: # multi
            self.disable_indep()
            self.all_check()
            self.driver_Tabs.setTabEnabled(1,False)
        return
    
    def switch_driverMode(self):
        """
        Updates driver gas dropdown menu when switching between passive and detonation modes.
        Also toggles the availability of the "Mixer" button (only req'd for detonation).
        """
        if self.passive_PushBtn.isChecked():
            self.drivGas_Dropdown.clear()
            self.drivGas_Dropdown.addItems(self.passive_gas_names)
            self.mixer_PushBtn.setEnabled(False)
        else:
            self.drivGas_Dropdown.clear()
            self.drivGas_Dropdown.addItems(self.det_gas_names)
            self.mixer_PushBtn.setEnabled(True)
        self.all_check()
            
    def detGas_changed(self):
        """
        Erase any existing saved combustion mixture if detonation gas type changed.
        """
        self.existingCombustionMix = None
        self.all_check()
            
    def open_mixer(self):
        """
        Create, display and obtain values from Mixer popup window.
        """
        # Create instance of MixerPopup class (defined in this library)
        self.mixerWindow = MixerPopup(parent=self,
                            selectedGas=str(self.drivGas_Dropdown.currentText()),
                            df=self.det_df,fodInit=self.existingCombustionMix)
        self.mixerWindow.show()
        # Obtain and store values from Mixer window only if it was closed by 
        # clicking the "OK" button (i.e. not the "Cancel" or "X" buttons)
        if self.mixerWindow.exec_():
            self.existingCombustionMix = self.mixerWindow.getSliderValues()
            self.construct_CompStr()
        self.all_check()
            
    def construct_CompStr(self):
        """ 
        Takes the 3 floats representing the mole fractions of the fuel,
        oxidizer, and diluent, and constructs a string that Cantera can interpret
        as a composition for the pre-detonation mixture.
        """
        selectedGas=str(self.drivGas_Dropdown.currentText())
        fuelSym = self.det_df.loc[selectedGas,'Fuel']
        oxSym = self.det_df.loc[selectedGas,'Oxidizer'] 
        dilSym = self.det_df.loc[selectedGas,'Diluent'] 
        
        xfuel,xox,xdil,self.phi,_ = self.existingCombustionMix
        
        self.detCompStr = fuelSym+':'+str(xfuel)+','+oxSym+':'+str(xox)+','+dilSym+':'+str(xdil)
     
    def abort(self):
        """
        Changes the abort flag to True when "Abort" button pressed. This flag
        is checked during each iteration in multi-calc mode. Allows the UI to exit
        cleanly at the end of the current iteration.
        """
        self.abortRequested = True
        return


##############################################################################
#** Check for correct type and value of user-entered variables
    def all_check(self):
        """
        Top-level error checking method. 
        Calls all other error checks according to a defined precedence.
        If any of them return an error code, the appropriate error message (defined here)
        is displayed in the UI, and the "Calculate" button is disabled.
        
        Only once all error checks have been performed with no error codes returned
        will the "Calculate" button be enabled.
        """
        # Define error messages
        typeErrorMsg = 'Please enter numbers only.'
        PErrorMsg = 'Need pressures such that driver > driven > accel > 0.'
        TErrorMsg = 'Need temperatures > 0.'
        LErrorMsg = 'Need lengths > 0.'
        mixerErrorMsg = 'Please set the composition of the pre-detonation mixture.'
        depVarErrormsg = 'Please select a dependent variable first.'
        shockMachErrorMsg = 'The primary shock Mach number exceeds the attainable perfect gas limit for these gases.'
        shockVelErrorMsg = 'The primary shock velocity exceeds the attainable perfect gas limit for these gases.'
        shockNegativeErrorMsg = 'Please enter positive values for the primary shock.'
        angleErrorMsg = 'Need angle between 0 and 90 deg.'
        
        # Define some shorthand booleans for various modes being active
        if self.calcType_Tabs.currentIndex() == 1:
            multimode = True
        else:
            multimode = False
            
        if self.driver_Tabs.currentIndex() == 1:
             shockmode = True
        else:
            shockmode = False            
        # if self.geom_Dropdown.currentText() != 'None'
        if self.geom_Dropdown.currentText() == 'Wedge':
            geommode = True
        else:
            geommode = False
        
        if self.detonation_PushBtn.isChecked():
            detmode = True
        else:
            detmode = False
        
        self.calc_PushBtn.setEnabled(False)
        # Ensure importing data takes precedence for error messages/disabling calc
        if self.inputFrame.isEnabled():
            # Core pressure/temperature/length checks
            if any(e == 1 for e in [self.pressure_error(),self.temperature_error(),self.length_error()]):
                self.msg_Display.setText(typeErrorMsg)
            elif self.pressure_error() == 2:
                self.msg_Display.setText(PErrorMsg)
            elif self.temperature_error() == 2:
                self.msg_Display.setText(TErrorMsg)
            elif self.length_error() == 2:
                self.msg_Display.setText(LErrorMsg)
            
            # Only perform checks for a mode if that mode is active (to avoid errors if an incorrect value
            # remains entered in an inactive tab/section)
            
            # Detonation mode checks
            elif detmode and self.existingCombustionMix is None:
                self.msg_Display.setText(mixerErrorMsg)
                
            # Multi-calc mode checks
            elif multimode and self.depVar_list.currentItem() is None:
                self.msg_Display.setText(depVarErrormsg)
                
            # Shock-specified mode checks
            elif shockmode and self.shock_limit() == 1:
                self.msg_Display.setText(typeErrorMsg)
            elif shockmode and self.shock_limit() == 2:
                self.msg_Display.setText(shockMachErrorMsg)
            elif shockmode and self.shock_limit() == 3:
                self.msg_Display.setText(shockVelErrorMsg)
            elif shockmode and self.shock_limit() == 4:
                self.msg_Display.setText(shockNegativeErrorMsg)
                
            # Test body geometry checks
            elif geommode and self.wedge_error() == 1:
                self.msg_Display.setText(typeErrorMsg)
            elif geommode and self.wedge_error() == 2:
                self.msg_Display.setText(angleErrorMsg)
                
            # If all checks are passed:
            else:
                self.calc_PushBtn.setEnabled(True)
                self.msg_Display.setText('')
            
             
        
    def pressure_error(self):
        """
        Tests driver P > test P > acceleration P.
        Also ensures all pressures are greater than zero.
        
        Returns 0 if there are no errors, or else returns an error code:        
        Error code 1 if incompatible format (e.g. string)
        Error code 2 if above two conditions not met

        CAUTION: code only takes into account the input pressure, which is fine
        for passive mode, but for detonation mode, this is the pre-detonation pressure.
        In theory, this could be less than the test gas pressure while still being valid
        since the post-detonation pressure is higher. However it will still return an error in this case.
        To get around this robustly, would have to calculate the post-detonation state every time this
        is called, which is impractically slow.
        However, most realistic conditions should not encounter this problem.        
        """
        # Convert to kPa first if necessary. Check that a compatible string is entered.
        indepStr = ''
        try: 
            drivP = float(self.drivP_Edit.text())
            testP = float(self.testP_Edit.text())
            accelP = float(self.accelP_Edit.text())
            
            # Lump the shock-input field in here too, just to check for formatting
            if self.driver_Tabs.currentIndex() == 1:
                float(self.primShock_Edit.text())
            
            # Obtain Start and Stop values if multi-calc indep var is any of the pressures
            if self.calcType_Tabs.currentIndex() == 1:
                indepStr = str(self.indepVar_ButtonGroup.checkedButton().objectName()[5:-6])
                if indepStr[-1] == 'P':
                    if self.linear_radio.isChecked():
                        indepStart = float(self.rangeStart_Edit.text())
                        indepStop = float(self.rangeStop_Edit.text())
                    if self.log_radio.isChecked():
                        indepStart = 10**float(self.rangeStart_Edit.text())
                        indepStop = 10**float(self.rangeStop_Edit.text())
                    int(abs(float(self.rangeNum_Edit.text())))
                    
        except ValueError:
            return 1
         
        # Convert units if they weren't entered in kPa
        if self.drivbar_radio.isChecked():
            drivP = ecl.bar2kpa(drivP)
            if indepStr == 'drivP':
                indepStart = ecl.bar2kpa(indepStart)
                indepStop = ecl.bar2kpa(indepStop)
        if self.testtorr_radio.isChecked():
            testP = ecl.torr2kpa(testP)
            if indepStr == 'testP':
                indepStart = ecl.torr2kpa(indepStart)
                indepStop = ecl.torr2kpa(indepStop)
        if self.acceltorr_radio.isChecked():
            accelP = ecl.torr2kpa(accelP)
            if indepStr == 'accelP':
                indepStart = ecl.torr2kpa(indepStart)
                indepStop = ecl.torr2kpa(indepStop)
        
        # Check conditions P4 > P1 > P5 > 0
        if self.calcType_Tabs.currentIndex() == 1 and indepStr[-1] == 'P':
            if (indepStr == 'drivP' and indepStart > testP and indepStop > testP
                and testP > accelP and all(P > 0 for P in [drivP,testP,accelP,indepStart,indepStop])):
                return 0
            elif (indepStr == 'testP' and drivP > indepStart and drivP > indepStop
                  and indepStart > accelP and indepStop > accelP
                  and all(P > 0 for P in [drivP,testP,accelP,indepStart,indepStop])):
                return 0
            elif (indepStr == 'accelP' and drivP > testP and testP > indepStart
                and testP > indepStop and all(P > 0 for P in [drivP,testP,accelP,indepStart,indepStop])):
                return 0
            else:
                return 2
        else:
            if drivP > testP and testP > accelP and all(P > 0 for P in [drivP,testP,accelP]):
                return 0
            else:
                return 2
            
            
    def temperature_error(self):
        """
        Tests whether temperatures are all positive.
        
        Returns 0 if there are no errors, or else returns an error code:        
        Error code 1 if incompatible format (e.g. string)
        Error code 2 if above condition not met
        """
        indepStr = ''
        try:
            drivT = float(self.drivT_Edit.text())
            testT = float(self.testT_Edit.text())
            accelT = float(self.accelT_Edit.text())
            
            # Obtain Start and Stop values if multi-calc indep var is any of the temperatures
            if self.calcType_Tabs.currentIndex() == 1:
                indepStr = str(self.indepVar_ButtonGroup.checkedButton().objectName()[5:-6])
                if indepStr[-1] == 'T':
                    if self.linear_radio.isChecked():
                        indepStart = float(self.rangeStart_Edit.text())
                        indepStop = float(self.rangeStop_Edit.text())
                    if self.log_radio.isChecked():
                        indepStart = 10**float(self.rangeStart_Edit.text())
                        indepStop = 10**float(self.rangeStop_Edit.text())
                    int(abs(float(self.rangeNum_Edit.text())))
            
        except ValueError:
            return 1   
        
        # Check conditions T1,T4,T5 > 0
        if self.calcType_Tabs.currentIndex() == 1 and indepStr[-1] == 'T':
            if all(T > 0 for T in [drivT,testT,accelT,indepStart,indepStop]):
                return 0
            else:
                return 2
        else:
            if all(T > 0 for T in [drivT,testT,accelT]):
                return 0
            else:
                return 2
            
    def length_error(self):
        """
        Tests whether lengths are all positive.
        
        Returns 0 if there are no errors, or else returns an error code:        
        Error code 1 if incompatible format (e.g. string)
        Error code 2 if above condition not met 
        """
        indepStr = ''
        try:
            drivL = float(self.drivL_Edit.text())
            testL = float(self.testL_Edit.text())
            accelL = float(self.accelL_Edit.text())
            
            # Obtain Start and Stop values if multi-calc indep var is any of the lengths
            if self.calcType_Tabs.currentIndex() == 1:
                indepStr = str(self.indepVar_ButtonGroup.checkedButton().objectName()[5:-6])
                if indepStr[-1] == 'L':
                    if self.linear_radio.isChecked():
                        indepStart = float(self.rangeStart_Edit.text())
                        indepStop = float(self.rangeStop_Edit.text())
                    if self.log_radio.isChecked():
                        indepStart = 10**float(self.rangeStart_Edit.text())
                        indepStop = 10**float(self.rangeStop_Edit.text())
                    int(abs(float(self.rangeNum_Edit.text())))
            
        except ValueError:
            return 1
        
        # Check conditions L1,L4,L5 > 0
        if self.calcType_Tabs.currentIndex() == 1 and indepStr[-1] == 'L':            
            if all(L > 0 for L in [drivL,testL,accelL,indepStart,indepStop]):
                return 0
            else:
                return 2
        else:
            if all(L > 0 for L in [drivL,testL,accelL]):
                return 0
            else:
                return 2
    
    def wedge_error(self):
        """
        Tests whether the angle input for test body geometry is between 0 and 90 deg.
        
        Returns 0 if there are no errors, or else returns an error code:        
        Error code 1 if incompatible format (e.g. string)
        Error code 2 if above condition not met 
        """
        try:
            angle = float(self.angle_Edit.text())
        except ValueError:
            return 1
        
        if angle >= 90 or angle <= 0:
            return 2
        
        return 0
        
            
    def shock_limit(self):
        """
        Uses the currently-entered values for the driver and driven gases to determine
        whether the requested primary shock can be achieved.
        
        Returns 0 if there are no errors, or else returns an error code:
        Error code 1 if incompatible format (e.g. string)
        Error code 2 if shock Mach number exceeds limit
        Error code 3 if shock velocity exceeds limit
        Error code 4 if input is negative
        
        Note that checks for correct format of other inputs don't need to be performed
        because their own checks take higher precedence than this check.
        
        Based on the asymptotic limit for shock Mach number in a simple, perfect-gas
        shocktube:
            
        (Ms^2 - 1)/Ms --> a4/a1 * (g1 + 1)/(g4 - 1)
        
        For the sake of reusing code, this check actually gets all inputs and creates
        PerfectGas objects as if it were going to begin a calculation. However,
        anything created here will be overwritten when the Calculate button is actually pressed.
        
        Also for the sake of speed, always checks based on a perfect gas, even if in real gas mode.
        This means that some conditions that have real gas solutions may be blocked, or that other conditions
        that don't have real gas solutions will be passed. This should be a rather narrow band around
        the perfect gas limit, however.
        """ 
        # 
        if self.driver_Tabs.currentIndex() == 0:
            return 0
        
        # Check format of shock input
        try:
            primShock = float(self.primShock_Edit.text())
            allGasStr,allP,allT,allL,allu,rawP,primShock,geom = self.get_single_inputs()
        except ValueError:
            return 1
        
        
        pDriverGas = epl.PerfectGas(allGasStr[0],allP[0],allT[0],allu[0],self.passive_df)            
        pTestGas = epl.PerfectGas(allGasStr[1],allP[1],allT[1],allu[1],self.passive_df)
        
        x = pDriverGas.a/pTestGas.a*(pTestGas.gamma+1)/(pDriverGas.gamma-1)
        
        # Take positive root of the quadratic
        Ms_limit = (x + sqrt(x**2 + 4))/2
        
        us_limit = Ms_limit*pTestGas.a
        
        if self.primmach_radio.isChecked() and Ms_limit < primShock:
            return 2
        elif self.primspeed_radio.isChecked() and us_limit < primShock:
            return 3
        elif primShock < 0:
            return 4
        else:       
            return 0

        
##############################################################################
#** Add and remove matplotlib canvas     
    def add_mpl(self, fig):
        """
        Add matplotlib: prepares a canvas and then draws a provided figure on it.
        """
        self.canvas = FigureCanvas(fig)
        if self.calcType_Tabs.currentIndex() == 0:
            self.single_mpl.addWidget(self.canvas)
        if self.calcType_Tabs.currentIndex() == 1:
            self.multi_mpl.addWidget(self.canvas)
        self.canvas.draw()
        
    def remove_mpl(self):
        """
        Remove matplotlib: deletes the figure and canvas.
        """
        if self.calcType_Tabs.currentIndex() == 0:
            self.single_mpl.removeWidget(self.canvas)
        if self.calcType_Tabs.currentIndex() == 1:
            self.multi_mpl.removeWidget(self.canvas)
        self.canvas.close()
   
##############################################################################
#** Assign inputs to variables, solve
    def get_single_inputs(self):
        """
        Read values from all input fields (line edits and dropdown menus).
        Does not contain any error checks (e.g. for type) since this function is
        usually only called by calculate(), which can only happen once all error
        checks in all_check() have passed, so formats must already be correct.
        
        The other place this gets called is during the shock_limit() check. There
        it is wrapped in a try/except ValueError statement which should catch incorrect
        input formats. Even this is unnecessary, as the precedence of other checks means
        that shock_limit() will never be called unless pressure_error(), temperature_error()
        and length_error() all pass first.
        """
        # Read the gas names, pressures, temperatures and section lengths
        drivGasStr = str(self.drivGas_Dropdown.currentText())
        testGasStr = str(self.testGas_Dropdown.currentText())
        accelGasStr = str(self.accelGas_Dropdown.currentText())        
        
        testP = float(self.testP_Edit.text())
        accelP = float(self.accelP_Edit.text())
        
        # Store original (unconverted) pressure values in separate variables if needed:
        # We keep these so that they can be appended to the output file (so it also displays
        # the nicely rounded original value in bar or torr instead of just the long messy
        # converted value in kPa
        rawDrivP = None
        rawTestP = None
        rawAccelP = None
        
        # Convert units if they weren't entered in kPa        
        if self.testtorr_radio.isChecked():
            rawTestP = testP
            testP = ecl.torr2kpa(testP)
        if self.acceltorr_radio.isChecked():
            rawAccelP = accelP
            accelP = ecl.torr2kpa(accelP)
        
        if self.driver_Tabs.currentIndex() == 0:
            drivP = float(self.drivP_Edit.text())
            if self.drivbar_radio.isChecked():
                rawDrivP = drivP
                drivP = ecl.bar2kpa(drivP)
            primShock = None                
        else:
            drivP = 1
            primShock = float(self.primShock_Edit.text())
            
        
        drivT = float(self.drivT_Edit.text())
        testT = float(self.testT_Edit.text())
        accelT = float(self.accelT_Edit.text())
                
        drivL = float(self.drivL_Edit.text())
        testL = float(self.testL_Edit.text())
        accelL = float(self.accelL_Edit.text())
        
        # Initially zero velocity in all 3 sections
        drivu = 0
        testu = 0
        accelu = 0
        
        # Test body geometry
        geom = [self.geom_Dropdown.currentText()]
        # if self.geom_Dropdown.currentText() != 'None':
        if self.geom_Dropdown.currentText() == 'Wedge':
            geom.append(float(self.angle_Edit.text()))
        
        # Group together for returning output. 
        # Use lists not tuples so can be edited when varying a parameter in multi-calc mode.
        allGasStr = [drivGasStr,testGasStr,accelGasStr]
        allP = [drivP,testP,accelP]
        allT = [drivT,testT,accelT]
        allL = [drivL,testL,accelL]
        allu = [drivu,testu,accelu]
        rawP = [rawDrivP,rawTestP,rawAccelP]
        
        return allGasStr,allP,allT,allL,allu,rawP,primShock,geom
    
    
    def calculate(self):
        """
        Callback for pressing the "Calculate" button. Top-level method for governing entire solution process,
        including input, solving, and output to both UI and file.
        
        - Disables and enables parts of the interface to avoid changes being made mid-calculation.
        - Reads in all formatted parameters from user input fields using get_single_inputs()
        - Passes parameters to the method single_solve(), which orchestrates the actual solving process.
        
        - Decides whether to perform a single or multiple calculation.
        - For multi-calc, reads in and prepares array of independent variables, sets up loop to repeatedly
          call single_solve() on each element in this array. Passes entire solution dataset to multi_plot()
          to create UI plots.
          
        - If the user requests file outputs, also calls single_file_output() or multi_file_output() as appropriate.
        """
        
        ####################################################################
        # Setup
        ####################################################################
        
        # This flag ensures that if an error is raised during calculation, that its displayed
        # error message doesn't get overwritten by other messages at the end of this method
        self.calcError = False
        # Error messages
        outErrormsg = 'Invalid output folder. Output file not generated.'
        
        # Disabling parts of interface
        self.calc_PushBtn.setEnabled(False)
        self.reset_PushBtn.setEnabled(False)
        # Disable the other calculation type tab to prevent switching during calculation
        self.calcType_Tabs.setTabEnabled(1-self.calcType_Tabs.currentIndex(),False)
        self.msg_Display.setText('Calculating, please wait...')
        QtWidgets.QApplication.processEvents()
        
        # Get user inputs
        allGasStr,allP,allT,allL,allu,rawP,primShock,geom = self.get_single_inputs()
        
        # Restructure the 9 possible independent variables into a dict
        # so that they can be referenced using the radio button name (in multicalc mode)
        self.indepDict = {'drivP':[allP[0],'Driver pressure',['kPa','bar'],rawP[0]],
                     'testP':[allP[1],'Driven pressure',['kPa','Torr'],rawP[1]],
                     'accelP':[allP[2],'Acceleration pressure',['kPa','Torr'],rawP[2]],
                     'drivT':[allT[0],'Driver temperature','K'],
                     'testT':[allT[1],'Driven temperature','K'],
                     'accelT':[allT[2],'Acceleration temperature','K'],
                     'drivL':[allL[0],'Driver length','m'],
                     'testL':[allL[1],'Driven length','m'],
                     'accelL':[allL[2],'Acceleration length','m']}
        
        ####################################################################
        # Enter either single or multiple calculation mode:
        
        
        ####################################################################
        # Single calculation
        ####################################################################
        if self.calcType_Tabs.currentIndex() == 0:
            all_data = self.single_solve(allGasStr,allP,allT,allL,allu,primShock,geom)
            if self.out_CheckBox.isChecked():
                try:
                    self.single_file_output(all_data)
                except IOError:
                    self.msg_Display.setText(outErrormsg)
                    self.calcError = True
        
        ####################################################################
        # Multiple calculations
        ####################################################################
        elif self.calcType_Tabs.currentIndex() == 1:
            self.multi_progressBar.setValue(0)
            self.abortRequested = False
            # Abort signal only acted on after an iteration finished, therefore only useful
            # for multiple calculations
            self.abort_PushBtn.setEnabled(True)
            
            # Determine independent variable, use Start/Stop/Num inputs to construct range array
            self.indepStr = str(self.indepVar_ButtonGroup.checkedButton().objectName()[5:-6])
            indepStart = float(self.rangeStart_Edit.text())
            indepStop = float(self.rangeStop_Edit.text())
            indepNum = int(abs(float(self.rangeNum_Edit.text())))
            
            if self.linear_radio.isChecked():                                
                self.indepVarRange = linspace(indepStart,indepStop,num=indepNum)
            elif self.log_radio.isChecked():
                self.indepVarRange = logspace(indepStart,indepStop,num=indepNum)
            
            # This is the master variable that contains all outputs for every data point in indepVarRange            
            self.allRange = []
            
            for i,indepVar in enumerate(self.indepVarRange):
                # Modify the selected independent variable. Convert units for P if needed
                if (self.indepStr == 'drivP') and self.drivbar_radio.isChecked():
                    self.indepDict[self.indepStr][0] = ecl.bar2kpa(indepVar)
                elif (self.indepStr == 'testP') and self.testtorr_radio.isChecked():
                    self.indepDict[self.indepStr][0] = ecl.torr2kpa(indepVar)
                elif (self.indepStr == 'accelP') and self.acceltorr_radio.isChecked():
                    self.indepDict[self.indepStr][0] = ecl.torr2kpa(indepVar)
                else:              
                    self.indepDict[self.indepStr][0] = indepVar
                
                # Reassemble back into lists:
                allP = [self.indepDict['drivP'][0],self.indepDict['testP'][0],self.indepDict['accelP'][0]]
                allT = [self.indepDict['drivT'][0],self.indepDict['testT'][0],self.indepDict['accelT'][0]]
                allL = [self.indepDict['drivL'][0],self.indepDict['testL'][0],self.indepDict['accelL'][0]]
                
                # Call this every loop to allow the UI to update itself (e.g. modify progress bar)
                QtWidgets.QApplication.processEvents()
                all_data = self.single_solve(allGasStr,allP,allT,allL,allu,primShock,geom)
                
                self.allRange.append(all_data)
                progress = float(i+1)/indepNum*100
                self.multi_progressBar.setValue(progress)
                if self.abortRequested:
                    # Truncate indepVarRange so it matches current length of datasets in allRange
                    # otherwise will cause errors when trying to plot arrays of differing length
                    self.indepVarRange = self.indepVarRange[0:i+1]
                    break
            
            # Create plot
            self.fig2 = Figure()
            # A flag indicating if there is an active multi-calculation dataset available
            self.activeMulti[0] = True
            # A flag indicating whether active dataset included geometry data
            if geom[0] == 'None':
                self.activeMulti[1] = False
            else:
                self.activeMulti[1] = True
            
            self.multi_plot()
            
            if self.out_CheckBox.isChecked():
                try:
                    self.multi_file_output()
                except IOError:
                    self.msg_Display.setText(outErrormsg)
                    self.calcError = True
        
        
        ####################################################################
        # Finish up
        ####################################################################
        self.calc_PushBtn.setEnabled(True)
        self.reset_PushBtn.setEnabled(True)
        self.calcType_Tabs.setTabEnabled(1-self.calcType_Tabs.currentIndex(),True)
        if not self.calcError: # Don't overwrite any error messages
            self.msg_Display.setText('Calculations completed successfully!')
            print('Completed.')
            
            
        ####################################################################
        # Optional testbed section
        ####################################################################        
       
        # Use this area to call temporary functions for testing on full datasets. 
        # Typically want to give either all_data or self.allRange as input 
        #(depending on if in single or multi-calc).
        # Generally outputs will just be to console behind UI, but could also point to 
        # temporary widgets on the UI for testing.
        # To avoid clutter, define these testing functions at the end of exp_common_lib.
        
        #ecl.BL_Re_lm(all_data,allL)
        
        return
    
    
    def single_solve(self,allGasStr,allP,allT,allL,allu,primShock,geom):
        """
        This method is the interface between the UI and the core solving methods in the
        other libraries. It solves a single expansion tube problem (i.e. one iteration
        in a multi-calc).
        
        This is the longest method of the Main class. It could probably be broken up,
        but that would require passing around a large number of intermediate variables,
        or else abusing "global" variables via use of "self" to an extreme extent...
        
        Using the formatted inputs passed to it from calculate(), it sets up instances
        of the PerfectGas and/or RealGas classes. It does this differently depending on
        what calculation mode is active (i.e. passive, detonation, shock-specified).
        
        It then passes these objects to the top-level solver functions in the perfect or
        real libraries, i.e. epl.polar_solve() and erl.real_polar_solve().
        
        These get called twice each (for the two shocktube problems that comprise the expansion
        tube problem). Hence the actual solving process only occupies a few lines.
        Much of this method is instead dedicated to either preparing the gas objects, or post-processing
        the solution to obtain additional data such as sound speed ratio, constructing x-t diagrams,
        and formatting the output into data structures to be passed back to calculate().
        """

        print('\n\n***************')
        print('Beginning new calculation...')       
        
        # Unpack:
        drivGasStr,testGasStr,accelGasStr = allGasStr
        drivP,testP,accelP = allP
        drivT,testT,accelT = allT
        drivL,testL,accelL = allL
        drivu,testu,accelu = allu
        
        ####################################################################
        # Create objects for perfect driver,test and acceleration gases
        ####################################################################
        
        # If operating in detonation mode, pDriverGas and rDriverGas represent the
        # post-detonation gas, i.e. state 4a.
        # In this case, we need to create a RealGas object representing the pre-detonation
        # mixture, then use SDT to solve for the post-detonation state. From this, we
        # can create an approximate PerfectGas object with fixed gamma etc.
        try:
            if self.driver_Tabs.currentIndex() == 0 and self.passive_PushBtn.isChecked():
                pDriverGas = epl.PerfectGas(drivGasStr,drivP,drivT,drivu,self.passive_df)
            elif self.driver_Tabs.currentIndex() == 0 and self.detonation_PushBtn.isChecked():
                print('Solving for post-detonation state...')
                preDetState = (self.det_df.loc[drivGasStr,'Cantera'],drivP*1e3,
                               drivT,self.detCompStr,drivGasStr) #mech,P,T,comp,name
                pDriverGas,rDriverGas,pPreDetGas,rPreDetGas,cj_speed = erl.detonation_setup(preDetState)              
            else:
                # If in shock-specified mode, set an arbitrary driver pressure,
                # in case the pressure entry on other tab was invalid/blank
                pDriverGas = epl.PerfectGas(drivGasStr,1.,drivT,drivu,self.passive_df)
            
            pTestGas = epl.PerfectGas(testGasStr,testP,testT,testu,self.passive_df)
            pAccelGas = epl.PerfectGas(accelGasStr,accelP,accelT,accelu,self.passive_df)
        except AttributeError:
            dataErrormsg = 'Gas data file format incorrect.'
            self.msg_Display.setText(dataErrormsg)
            self.calcError = True
            self.clear_output()
            return
        
        ####################################################################
        # Begin the solving process for perfect gas
        ####################################################################
        self.realCalc = False
        if self.driver_Tabs.currentIndex() == 0: # Driver-specified mode
            # First polar intersection: use states 1 & 4 to find 2 & 3
            Ms1_perf,pGas2,pGas3,pdxdt1 = epl.polar_solve(pTestGas,pDriverGas,self.passive_df)
            us1_perf = Ms1_perf*pTestGas.a
            solvedDrivP_perf = pDriverGas.P
        else: # Shock-specified mode
            if self.primmach_radio.isChecked():
                Ms1_perf = primShock
                us1_perf = Ms1_perf*pTestGas.a
            else:
                Ms1_perf = primShock/pTestGas.a
                us1_perf = primShock
            pGas2 = epl.shock_jump(pTestGas,Ms1_perf) # fully defined State 2
            pGas3,pdxdt1 = epl.unsteady_exp(pDriverGas,pGas2.vel,pGas2.P,self.passive_df) # fully defined State 3
            solvedDrivP_perf = epl.polar_P(pGas3.P,pGas3.vel,pDriverGas.vel,pGas3.gamma,pGas3.a) # kPa
            pDriverGas.P = solvedDrivP_perf
        
        
        # Second polar intersection: use states 2 & 5 to find 6 & 7
        Ms2_perf,pGas6,pGas7,pdxdt2 = epl.polar_solve(pAccelGas,pGas2,self.passive_df)
        us2_perf = Ms2_perf*pAccelGas.a
        
        ####################################################################
        # If required, create objects for corresponding real gases and solve real system
        ####################################################################
        if self.modelreal_radio.isChecked():
            self.realCalc = True
            conv = 1e3
            
            # Creation. Note Cantera requires Pa, not kPa
            if self.driver_Tabs.currentIndex() == 0 and self.passive_PushBtn.isChecked():
                # Already created above for detonation
                rDriverGas = erl.RealGas(self.passive_df.loc[drivGasStr,'Cantera'],vel=drivu,T=drivT,
                                         P=drivP*conv,X=self.passive_df.loc[drivGasStr,'X'],gas=drivGasStr,equil=True)
            elif self.driver_Tabs.currentIndex() == 1:
                # For real case, don't use arbitrary pressure. Instead, use solved driver
                # pressure from perfect case above
                rDriverGas = erl.RealGas(self.passive_df.loc[drivGasStr,'Cantera'],vel=drivu,T=drivT,
                                         P=solvedDrivP_perf*conv,X=self.passive_df.loc[drivGasStr,'X'],gas=drivGasStr,equil=True)    
            
            rTestGas = erl.RealGas(self.passive_df.loc[testGasStr,'Cantera'],vel=testu,T=testT,
                                     P=testP*conv,X=self.passive_df.loc[testGasStr,'X'],gas=testGasStr,equil=True)
            rAccelGas = erl.RealGas(self.passive_df.loc[accelGasStr,'Cantera'],vel=accelu,T=accelT,
                                     P=accelP*conv,X=self.passive_df.loc[accelGasStr,'X'],gas=accelGasStr,equil=True)
           
            # Initialize states 2, 3 and 6 as copies of 1, 4, and 5 respectively
            if self.passive_PushBtn.isChecked():
                rGas3 = erl.RealGas(self.passive_df.loc[drivGasStr,'Cantera'],vel=rDriverGas.vel,T=rDriverGas.T,
                                     P=rDriverGas.P,X=rDriverGas.X,gas=rDriverGas.gas,equil=True)
            else:
                rGas3 = erl.RealGas(self.det_df.loc[drivGasStr,'Cantera'],vel=rDriverGas.vel,T=rDriverGas.T,
                                     P=rDriverGas.P,X=rDriverGas.X,gas=rDriverGas.gas,equil=True)   
                        
            rGas2 = erl.RealGas(self.passive_df.loc[testGasStr,'Cantera'],vel=testu,T=testT,
                                     P=testP*conv,X=self.passive_df.loc[testGasStr,'X'],gas=testGasStr,equil=True)
            rGas6 = erl.RealGas(self.passive_df.loc[accelGasStr,'Cantera'],vel=accelu,T=accelT,
                                     P=accelP*conv,X=self.passive_df.loc[accelGasStr,'X'],gas=accelGasStr,equil=True)
            
            print('Gases 1-6 initialized!')
            # Solving            
            if self.driver_Tabs.currentIndex() == 0:
                try:
                    print('Beginning first shocktube problem...')
                    Ms1_real, rdxdt1 = erl.real_polar_solve(rTestGas,rGas2,rGas3,rDriverGas,Ms1_perf)  
                    print('First shocktube problem completed!')
                except RuntimeError as ctError:
                    ctErrormsg = '{0}'.format(ctError)
                    # Format error message by removing * and \n so it fits in UI
                    ctErrormsg = ctErrormsg.translate(None,'*\n') 
                    self.msg_Display.setText(ctErrormsg+' in first shocktube problem')
                    self.calcError = True
                    self.clear_output()
                    return
                except ValueError as valError:
                    self.msg_Display.setText(valError.args[0]+' in first shocktube problem')
                    self.calcError = True
                    self.clear_output()
                    return
                us1_real = Ms1_real*rTestGas.a
                
            else:
                if self.primmach_radio.isChecked():
                    Ms1_real = primShock
                    us1_real = Ms1_real*rTestGas.a
                else:
                    Ms1_real = primShock/rTestGas.a
                    us1_real = primShock
                # Use the perfect driver gas solution as an initial guess for the real solution
                solvedDrivP_real,rdxdt1 = erl.real_polar_solve_inv(rTestGas,rGas2,rGas3,rDriverGas,
                                                                   Ms1_real,solvedDrivP_perf*conv)
                
            # Initialize state 7 as a copy of state 2
            rGas7 = erl.RealGas(self.passive_df.loc[testGasStr,'Cantera'],vel=rGas2.vel,T=rGas2.T,
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
                self.msg_Display.setText(ctErrormsg+' in second shocktube problem')
                self.calcError = True
                self.clear_output()
                return
            except ValueError as valError:
                self.msg_Display.setText(valError.args[0]+' in second shocktube problem')
                self.calcError = True
                self.clear_output()
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
        
        if self.modelreal_radio.isChecked():
            rGas7.reynolds()
            ssr_real = rGas3.a/rGas2.a # sound speed ratio
            ssr.append(ssr_real)
            if ssr_real < 1:
                mode_real = 'High'
            else:
                mode_real = 'Low'
            mode.append(mode_real)
        
        # Identify regions of min/max temperature in facility
        if self.passive_PushBtn.isChecked():
            pTarray = asarray([pTestGas.T,pGas2.T,pGas3.T,pDriverGas.T,pAccelGas.T,pGas6.T,pGas7.T])
            if self.modelreal_radio.isChecked():
                rTarray = asarray([rTestGas.T,rGas2.T,rGas3.T,rDriverGas.T,rAccelGas.T,rGas6.T,rGas7.T])
            gasDesc = self.gasDesc_passive
        elif self.detonation_PushBtn.isChecked():
            pTarray = asarray([pTestGas.T,pGas2.T,pGas3.T,rPreDetGas.T,pDriverGas.T,pAccelGas.T,pGas6.T,pGas7.T])
            if self.modelreal_radio.isChecked():
                rTarray = asarray([rTestGas.T,rGas2.T,rGas3.T,rPreDetGas.T,rDriverGas.T,rAccelGas.T,rGas6.T,rGas7.T])
            gasDesc = self.gasDesc_det
       
        pMinTidx = pTarray.argmin()
        pMaxTidx = pTarray.argmax()        
        pExtremeT = (pTarray[pMinTidx],gasDesc[pMinTidx],pTarray[pMaxTidx],gasDesc[pMaxTidx])
        if self.modelreal_radio.isChecked():          
            rMinTidx = rTarray.argmin()
            rMaxTidx = rTarray.argmax()
            rExtremeT = (rTarray[rMinTidx],gasDesc[rMinTidx],rTarray[rMaxTidx],gasDesc[rMaxTidx])
        
        # Use freestream solution to calculate flow over test body geometry, if required
        # If additional geometries are added, should probably spin this out into own method
        pGeomOut = None
        rGeomOut = None
        geom
        #print('[%s]' % ', '.join(map(str, geom)))
        if geom[0] == 'Wedge':
            pGeomOut = list(epl.wedge(pGas7,geom[1]))
            pGeomOut.extend(geom)            
            if self.modelreal_radio.isChecked():
                print('Solving for flow over wedge...')
                rPostOblGas = erl.RealGas(self.passive_df.loc[testGasStr,'Cantera'],vel=rGas7.vel,T=rGas7.T,
                                     P=rGas7.P,X=rGas7.X,gas=testGasStr,equil=True)
                print('Post-oblique shock gas initialized!')
                rGeomOut = list(erl.wedge(rGas7,geom[1],rPostOblGas))
                rGeomOut.extend([rPostOblGas])
                rGeomOut.extend(geom)
                print('Wedge flow problem complete!')
        elif geom[0] == 'Normal shock':
            pGas8 = epl.shock_jump(pGas7,pGas7.M)
            pGas8.vel = pGas7.vel-pGas8.vel
            pGas8.calc_properties()
            pGeomOut = [pGas8]
            pGeomOut.extend(geom)
            if self.modelreal_radio.isChecked():
                rGas8_eq_dummy = sdt.postshock.PostShock_eq(rGas7.vel,rGas7.P,rGas7.T,rGas7.X,rGas7.ID+'.cti')
                rGas8_fr_dummy = sdt.postshock.PostShock_fr(rGas7.vel,rGas7.P,rGas7.T,rGas7.X,rGas7.ID+'.cti')
                
                rGas8_eq = erl.RealGas(rGas7.ID+'.cti',P=rGas8_eq_dummy.P,
                                   T=rGas8_eq_dummy.T,X=rGas8_eq_dummy.X,
                                   equil=True,vel=0,gas=rGas7.gas)
                
                rGas8_fr = erl.RealGas(rGas7.ID+'.cti',P=rGas8_fr_dummy.P,
                                   T=rGas8_fr_dummy.T,X=rGas8_fr_dummy.X,
                                   equil=False,vel=0,gas=rGas7.gas)
                
                rGas8_eq.vel = rGas7.vel*rGas7.density/rGas8_eq.density
                rGas8_fr.vel =  rGas7.vel*rGas7.density/rGas8_fr.density
                
                rGas8_eq.extra()
                rGas8_eq.reynolds()
                rGas8_fr.extra()
                rGas8_fr.reynolds()
                
                rGeomOut = [rGas8_eq,rGas8_fr]
                rGeomOut.extend(geom)
                rGeomOut
       
        # Compile outputs for use by UI display
        if self.passive_PushBtn.isChecked():
            us_perf = (us1_perf,us2_perf)
            Ms_perf = (Ms1_perf,Ms2_perf)
        else:
            # Include detonation wave speed and Mach number in these tuples
            # Calculate Mach number based on sound speed in pre-detonation gas mix
            # Use perfect or real versions of this gas as appropriate
            us_perf = (us1_perf,us2_perf,cj_speed)
            Ms_perf = (Ms1_perf,Ms2_perf,cj_speed/pPreDetGas.a)  
        pdxdt = (pdxdt1,pdxdt2)
        
        perf_outputs = (pGas2,pDriverGas,pGas7,us_perf,
                        Ms_perf,pdxdt,ssr_perf,mode_perf,pExtremeT,pGeomOut)
        
        real_outputs = ()
        if self.modelreal_radio.isChecked():
            if self.passive_PushBtn.isChecked():
                us_real = (us1_real,us2_real)
                Ms_real = (Ms1_real,Ms2_real)
            else:
                us_real = (us1_real,us2_real,cj_speed)
                Ms_real = (Ms1_real,Ms2_real,cj_speed/rPreDetGas.a)
            rdxdt = (rdxdt1,rdxdt2)
        
            real_outputs = (rGas2,rDriverGas,rGas7,us_real,
                            Ms_real,rdxdt,ssr_real,mode_real,rExtremeT,rGeomOut)
            

        # Generate x-t diagram and display output to UI
        fighandle,t_test,term,t_contact2,optimal = self.display_output(allL,perf_outputs,real_outputs)
       
        # Compile complete output for saving to file
        all_data = {'Test time':t_test,
                    'x-t diagram':fighandle,
                    'Enthalpy mode':mode,
                    'Sound speed ratio':ssr,
                    'Perfect shock':(Ms1_perf,Ms2_perf,us1_perf,us2_perf),
                    'Termination mode':term,
                    'Lengths':allL,
                    'Contact arrival':t_contact2,
                    'Optimal test time':optimal
                    }
        
        # Add extra data depending on calculation mode
        if self.passive_PushBtn.isChecked():
            all_data['Perfect gases'] = (pTestGas,pGas2,pGas3,pDriverGas,pAccelGas,pGas6,pGas7)
            if self.modelreal_radio.isChecked():           
                all_data['Real gases'] = (rTestGas,rGas2,rGas3,rDriverGas,rAccelGas,rGas6,rGas7)
                all_data['Real shock'] = (Ms1_real,Ms2_real,us1_real,us2_real)
        elif self.detonation_PushBtn.isChecked():
            all_data['Perfect gases'] = (pTestGas,pGas2,pGas3,pPreDetGas,pDriverGas,pAccelGas,pGas6,pGas7)
            all_data['Detonation'] = cj_speed
            if self.modelreal_radio.isChecked():           
                all_data['Real gases'] = (rTestGas,rGas2,rGas3,rPreDetGas,rDriverGas,rAccelGas,rGas6,rGas7)
                all_data['Real shock'] = (Ms1_real,Ms2_real,us1_real,us2_real)
         
        
        if geom[0] == 'Wedge':
            all_data['Perfect test body'] = pGeomOut
            if self.modelreal_radio.isChecked():
                all_data['Real test body'] = rGeomOut
            
        return all_data
    
    
    def multi_data_prep(self,depName):
        """
        This method prepares the data contained in self.allRange (i.e. the entire raw data output
        from performing a multi-calc).
        
        Input: - depName: str of selected dependent variable from UI list
        
        Output: - pdepVarRange, rdepVarRange: lists of the dependent variable range.
                    Used both for plotting against indepVarRange, and for making output csv file
                - indepNameUnit, depNameUnit: str of formatted variable names and units
                    Also used both for plot labels and headers in csv file
        """
        
        depUnit = self.depVarDict[depName][2]
        pKey,pIndex,pAttr,pConv = self.depVarDict[depName][0]
        rKey,rIndex,rAttr,rConv = self.depVarDict[depName][1]
        
        pdepVarRange = [] # perfect
        rdepVarRange = [] # real
        for dataset in self.allRange:
            if pAttr is not None:
                if getattr(dataset[pKey][pIndex],pAttr) is not None:
                    pdepVarRange.append(getattr(dataset[pKey][pIndex],pAttr)/pConv)
                else:
                    # Do this step to avoid dividing NoneType by pConv (e.g. in the case of detached shock where beta = None)
                    pdepVarRange.append(None)
                # Only perform next step if real mode select AND if real calculation
                # data is in memory (in case they change the radio button and change dep.
                # variable without clicking Calculate)
                if self.modelreal_radio.isChecked() and self.realCalc:
                    if getattr(dataset[rKey][rIndex],rAttr) is not None:
                        rdepVarRange.append(getattr(dataset[rKey][rIndex],rAttr)/rConv)
                    else:
                        rdepVarRange.append(None)
            else:
                if dataset[pKey][pIndex] is not None:
                    pdepVarRange.append((dataset[pKey][pIndex])/pConv)
                else:
                    pdepVarRange.append(None)
                if self.modelreal_radio.isChecked() and self.realCalc:
                    if dataset[rKey][rIndex] is not None:
                        rdepVarRange.append((dataset[rKey][rIndex])/rConv)
                    else:
                        rdepVarRange.append(None)
        
        # Prepare labels      
        indepNameUnit,_,__ = ecl.stateVar_formatted(self.indepStr,self.indepDict,
                                           (self.drivbar_radio.isChecked(),
                                           self.testtorr_radio.isChecked(),
                                           self.acceltorr_radio.isChecked()) )                 
        
        if depUnit != '': # Don't append units if unitless
            depNameUnit = depName+' ['+depUnit+']'
        else:
            depNameUnit = depName       
        
        return pdepVarRange,rdepVarRange,indepNameUnit,depNameUnit
    
    def multi_plot(self):
        """
        Creates the plot for multi-calc mode. Takes the arrays for the independent variable
        and the currently-selected dependent variable --  generated by multi_data_prep() --
        as well as the formatted labels provided by that method, and plots them against each other.
        
        Contains formatting instructions for perfect and real datasets, as well as whether to use
        linear or logarithmic axes.
        """
        if self.activeMulti[0]:
            self.fig2.clear()             
            self.remove_mpl()
            self.ax1f2 = self.fig2.add_subplot(111)
            depName = self.depVar_list.currentItem().text()
            if depName in self.wedgeDepVars and not self.activeMulti[1]:
                # If the selected dependent variable relies on test body data (e.g. oblique shock angle)
                # but the existing calculation didn't have a test body, don't plot
                # This functionality only exists in case a geometry is selected from the dropdown menu,
                # thereby updating the dependent variable list, while an old calculation is still active.                            
                pass
            else:
                pdepVarRange,rdepVarRange,indepNameUnit,depNameUnit = self.multi_data_prep(depName)
                
                # Plot
                if self.linear_radio.isChecked():
                    if self.modelperf_radio.isChecked() or self.plotPerf_CheckBox.isChecked():
                        self.ax1f2.plot(self.indepVarRange,pdepVarRange,'bo',label='Perfect')
                    if self.modelreal_radio.isChecked() and self.realCalc:
                        self.ax1f2.plot(self.indepVarRange,rdepVarRange,'ro',label='Real')
                        
                elif self.log_radio.isChecked():
                    if self.modelperf_radio.isChecked() or self.plotPerf_CheckBox.isChecked():
                        self.ax1f2.semilogx(self.indepVarRange,pdepVarRange,'bo',label='Perfect')
                    if self.modelreal_radio.isChecked() and self.realCalc:
                        self.ax1f2.semilogx(self.indepVarRange,rdepVarRange,'ro',label='Real')
                
                # Axis labels       
                self.ax1f2.set_xlabel(indepNameUnit)
                self.ax1f2.set_ylabel(depNameUnit)
                
                # Legend: place outside and above axes, expand to fill width
                self.ax1f2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
                
                # Don't scale or offset the y-axis even if the range is narrow
                self.ax1f2.get_yaxis().get_major_formatter().set_useOffset(False)
                
                # Add the plot to the UI
                self.add_mpl(self.fig2)
         
        return
    
        if self.mirels.isChecked():
            mirels = erl.boundary_layer
            print(mirels)
        return
    
    
##############################################################################
#** Output results for displaying in UI 
    def display_output(self,lengths,perf_outputs,real_outputs):
        """
        Displays selected values in the output section of the UI, in single-calc mode.
        Also generates and displays the x-t diagram in single-calc mode.
        
        While this is a long method, most of it is repetitive, since individual calls are made
        in order to edit the value of each label (the widget used to display static output).
        
        For each label: 
            - If in perfect mode, just display the perfect value in the first column (named 'test')
            - If in real mode, display the real value in the first column ('test')
                               display the perfect value in the second column ('perf')
                               calculate then display the % difference in the third column ('diff')
                               
        Most of the labels are always used for the same values and have names (set in the UI file)
        that reflect this. There are 5 rows of variable labels that are used for different purposes
        (e.g. shock-specified or detonation modes).
            
        
        """
        self.clear_output()
        self.msg_Display.setText('')
        # Unpack perfect output (and real, only if we are doing real calculation)
        pGas2,pGas4,pGas7,us_perf,Ms_perf,pdxdt,ssr_perf,mode_perf,pExtremeT,pGeomOut = perf_outputs
        if self.modelreal_radio.isChecked():
            rGas2,rGas4,rGas7,us_real,Ms_real,rdxdt,ssr_real,mode_real,rExtremeT,rGeomOut = real_outputs
            rGas2.gamma = rGas2.cp/rGas2.cv # Needed for similarity soln
            rGas4.gamma = rGas4.cp/rGas4.cv
        
        # Generate and display x-t plot, calculate test times
        fig1 = Figure()
        ax1f1 = fig1.add_subplot(111)
        # give the subplot handle to the plotting function
        t_test_perf, term_perf, t_contact2_perf, opt_perf = ecl.exp_tube_xt(us_perf,pdxdt,pGas2,pGas4,
                                                                            pGas7.vel,lengths,ax1f1)
        
        # Prepare arrays for passing back to 'calculate'
        t_test = [t_test_perf]
        term = [term_perf]
        t_contact2 = [t_contact2_perf]
        optimal = [opt_perf]
        
        if self.modelreal_radio.isChecked():
            fig1 = Figure()
            ax1f1 = fig1.add_subplot(111)
            t_test_real, term_real, t_contact2_real, opt_real = ecl.exp_tube_xt(us_real,rdxdt,rGas2,rGas4,
                                                                                rGas7.vel,lengths,ax1f1)
            t_test.append(t_test_real)
            term.append(term_real)
            t_contact2.append(t_contact2_real)
            optimal.append(opt_real)
            
        
        fig1.subplots_adjust(top=0.80) # scale axes to make room for legend
        
        # Only display values and plot in UI if in single-calc mode
        if self.calcType_Tabs.currentIndex() == 0:
            self.add_mpl(fig1)
            
            # Display values
            if self.modelreal_radio.isChecked(): # Real
                # Column titles
                self.column1_Display.setText('Real')
                self.column2_Display.setText('Perfect')
                self.column3_Display.setText('Error')
                self.column4_Display.setText('Real')
                self.column5_Display.setText('Perfect')
                self.column6_Display.setText('Error')                                         
                
                # Speeds and times                             
                self.testM_Display.setText('%.3g' % rGas7.M) # round to fixed s.f.
                self.testu_Display.setNum(round(rGas7.vel,0)) # round to fixed d.p.
                self.testtime_Display.setNum(round(t_test_real*1e3,0)) # test time in us
                self.perfM_Display.setText('%.3g' % pGas7.M) 
                self.perfu_Display.setNum(round(pGas7.vel,0))
                self.perftime_Display.setNum(round(t_test_perf*1e3,0)) # test time in us
                self.diffM_Display.setText('%.1f' % ((pGas7.M-rGas7.M)/rGas7.M*100.)+' %')
                self.diffu_Display.setText('%.1f' % ((pGas7.vel-rGas7.vel)/rGas7.vel*100.)+' %')
                self.difftime_Display.setText('%.1f' % ((t_test_perf-t_test_real)/t_test_real*100.)+' %')                
                
                # Freestream thermodynamic properties
                self.testT_Display.setNum(round(rGas7.T,0))
                self.testrho_Display.setText('%.3f' % (rGas7.density*1e3))
                self.testP_Display.setText('%.3g' % (rGas7.P/1e3))
                self.testh0_Display.setText('%.4g' % (rGas7.h0/1e6))
                self.perfT_Display.setNum(round(pGas7.T,0))
                self.perfrho_Display.setText('%.3f' % (pGas7.density*1e3))
                self.perfP_Display.setText('%.3g' % pGas7.P)
                self.perfh0_Display.setText('%.4g' % (pGas7.h0/1e3))
                self.diffT_Display.setText('%.1f' % ((pGas7.T-rGas7.T)/rGas7.T*100.)+' %')
                self.diffrho_Display.setText('%.1f' % ((pGas7.density-rGas7.density)/rGas7.density*100.)+' %')
                self.diffP_Display.setText('%.1f' % ((pGas7.P-rGas7.P/1e3)/(rGas7.P/1e3)*100.)+' %')
                self.diffh0_Display.setText('%.1f' % ((pGas7.h0/1e3-rGas7.h0/1e6)/(rGas7.h0/1e6)*100.)+' %')
                
                # Freestream transport properties
                self.testRem_Display.setNum(round(rGas7.Rem*1e-3,0))
                self.testk_Display.setNum(round(rGas7.thermal_conductivity*1e3,0))
                
                # Shocks                
                self.testus1_Display.setNum(round(us_real[0],0))                
                self.testMs1_Display.setText('%.3g' % Ms_real[0])                
                self.perfus1_Display.setNum(round(us_perf[0],0))               
                self.perfMs1_Display.setText('%.3g' % Ms_perf[0])                
                self.diffus1_Display.setText('%.1f' % ((us_perf[0]-us_real[0])/us_real[0]*100.)+' %')                
                self.diffMs1_Display.setText('%.1f' % ((Ms_perf[0]-Ms_real[0])/Ms_real[0]*100.)+' %')
            
                self.testus2_Display.setNum(round(us_real[1],0))
                self.testMs2_Display.setText('%.3g' % Ms_real[1])
                self.perfus2_Display.setNum(round(us_perf[1],0))
                self.perfMs2_Display.setText('%.3g' % Ms_perf[1])
                self.diffus2_Display.setText('%.1f' % ((us_perf[1]-us_real[1])/us_real[1]*100.)+' %')
                self.diffMs2_Display.setText('%.1f' % ((Ms_perf[1]-Ms_real[1])/Ms_real[1]*100.)+' %')
                
                # Operating mode
                self.testEnd_Display.setText(term_real)
                self.perfEnd_Display.setText(term_perf)
                self.ssr_Display.setText('%.2f' % ssr_real) 
                self.mode_Display.setText(mode_real)
                self.perfssr_Display.setText('%.2f' % ssr_perf) 
                self.perfmode_Display.setText(mode_perf)
                
                # Extreme temperatures
                self.testLowT_Display.setNum(round(rExtremeT[0],0))
                self.testHighT_Display.setNum(round(rExtremeT[2],0))
                self.testLowTLoc_Display.setText(rExtremeT[1])
                self.testHighTLoc_Display.setText(rExtremeT[3])
                self.perfLowT_Display.setNum(round(pExtremeT[0],0))
                self.perfHighT_Display.setNum(round(pExtremeT[2],0))
                self.perfLowTLoc_Display.setText(pExtremeT[1])
                self.perfHighTLoc_Display.setText(pExtremeT[3])
                #print('[%s]' % ', '.join(map(str, rGeomOut)))
                #print('[%s]' % ', '.join(map(str, pGeomOut)))
                # Test body geometry
                if rGeomOut is not None:
                    if rGeomOut[-2] == 'Wedge':
                        if rGeomOut[0] is None:
                            self.testOblique_Display.setText('Detached')
                        else:
                            self.testOblique_Display.setNum(round(rGeomOut[0],1))
                        self.testDetach_Display.setNum(round(rGeomOut[1],1))
                    # Terminal display for now
                    elif rGeomOut[-1] == 'Normal shock':
                        print('\nPost-shock conditions (real gas - equilibrium):')
                        print('Pressure = '+str(rGeomOut[0].P*1e-3)+' kPa')
                        print('Temperature = '+str(rGeomOut[0].T)+' K')
                        print('Density = '+str(rGeomOut[0].density)+' kg/m^3')
                        print('Velocity = '+str(rGeomOut[0].vel)+' m/s')
                        print('Mach number = '+str(rGeomOut[0].M))
                        print('Unit Reynolds number = '+str(rGeomOut[0].Rem)+' m^-1')
                        print('Total enthalpy = '+str(rGeomOut[0].h0*1e-6)+' MJ/kg')
                        rGeomOut[0]()
                        
                        print('\n\n\nPost-shock conditions (real gas - frozen):')
                        print('Pressure = '+str(rGeomOut[1].P*1e-3)+' kPa')
                        print('Temperature = '+str(rGeomOut[1].T)+' K')
                        print('Density = '+str(rGeomOut[1].density)+' kg/m^3')
                        print('Velocity = '+str(rGeomOut[1].vel)+' m/s')
                        print('Mach number = '+str(rGeomOut[1].M))
                        print('Unit Reynolds number = '+str(rGeomOut[1].Rem)+' m^-1')
                        print('Total enthalpy = '+str(rGeomOut[1].h0*1e-6)+' MJ/kg')
                        rGeomOut[1]()
                        
                    
                    
                if pGeomOut is not None:
                    if pGeomOut[-2] == 'Wedge':
                        if pGeomOut[0] is None:
                            self.perfOblique_Display.setText('Detached')
                        else:
                            self.perfOblique_Display.setNum(round(pGeomOut[0],1))
                            self.diffOblique_Display.setText('%.1f' % ((pGeomOut[0]-rGeomOut[0])/rGeomOut[0]*100.)+' %')
                        
                        self.perfDetach_Display.setNum(round(pGeomOut[1],1))
                        self.diffDetach_Display.setText('%.1f' % ((pGeomOut[1]-rGeomOut[1])/rGeomOut[1]*100.)+' %')
                        
                    # Terminal display for now
                    elif pGeomOut[-1] == 'Normal shock':
                        print('\n\n\nPost-shock conditions (perfect gas):')
                        print('Pressure = '+str(pGeomOut[0].P)+' kPa')
                        print('Temperature = '+str(pGeomOut[0].T)+' K')
                        print('Density = '+str(pGeomOut[0].density)+' kg/m^3')
                        print('Velocity = '+str(pGeomOut[0].vel)+' m/s')
                        print('Mach number = '+str(pGeomOut[0].M))
                        print('Total enthalpy = '+str(pGeomOut[0].h0*1e-3)+' MJ/kg')
                  
                # Variable slots:
                    # Driver pressure
                if self.driver_Tabs.currentIndex() == 1:
                    self.varTitle_Display.setText('Shock-specified mode')
                    if rGas4.P > 1e6:                    
                        self.varLabel1_Display.setText('Driver pressure (MPa)')
                        rConv = 1e6
                        pConv = 1e3                        
                    else:
                        self.varLabel1_Display.setText('Driver pressure (kPa)')
                        rConv = 1e3
                        pConv = 1   
                        
                    self.testLabel1_Display.setText('%.3g' % (rGas4.P/rConv)) # MPa
                    self.perfLabel1_Display.setText('%.3g' % (pGas4.P/pConv)) 
                    self.diffLabel1_Display.setText('%.1f' % ((pGas4.P/pConv-rGas4.P/rConv)/(rGas4.P/rConv)*100.)+' %')
                        
                    # Post-detonation state
                if self.detonation_PushBtn.isChecked():
                    # No point in showing perfect calc or difference, because it's the same data
                    self.varTitle_Display.setText('Detonation')
                    self.varLabel1_Display.setText('Post-det. pressure (MPa)')
                    self.testLabel1_Display.setText('%.3g' % (rGas4.P/1e6)) # MPa
                    self.varLabel2_Display.setText('Post-det. temperature (K)')
                    self.testLabel2_Display.setNum(round(rGas4.T,0)) # K
                    self.varLabel3_Display.setText('Post-det. velocity (to left) (m/s)')
                    self.testLabel3_Display.setNum(round(abs(rGas4.vel),0)) # m/s
                    self.varLabel4_Display.setText('Detonation wave speed (m/s)')
                    self.testLabel4_Display.setNum(round((us_real[2]),0)) # m/s
                    self.varLabel5_Display.setText('Detonation wave Mach number')
                    self.testLabel5_Display.setText('%.3g' % Ms_real[2])
                     
            else: # Perfect
                # Column titles
                self.column1_Display.setText('Perfect')
                self.column4_Display.setText('Perfect')
                
                # Speeds and times
                self.testM_Display.setText('%.3g' % pGas7.M) 
                self.testu_Display.setNum(round(pGas7.vel,0))
                self.testtime_Display.setNum(round(t_test_perf*1e3,0)) # test time in us
                
                # Freestream thermodynamic properties                            
                self.testT_Display.setNum(round(pGas7.T,0))
                self.testrho_Display.setText('%.3f' % (pGas7.density*1e3))
                self.testP_Display.setText('%.3g' % pGas7.P)
                self.testh0_Display.setText('%.4g' % (pGas7.h0/1e3))
                
                # Shocks
                self.testus1_Display.setNum(round(us_perf[0],0))                
                self.testMs1_Display.setText('%.3g' % Ms_perf[0])
                
                self.testus2_Display.setNum(round(us_perf[1],0))
                self.testMs2_Display.setText('%.3g' % Ms_perf[1])
                
                # Operating mode
                self.testEnd_Display.setText(term_perf) 
                self.ssr_Display.setText('%.2f' % ssr_perf) 
                self.mode_Display.setText(mode_perf)
                
                # Extreme temperatures
                self.testLowT_Display.setNum(round(pExtremeT[0],0))
                self.testHighT_Display.setNum(round(pExtremeT[2],0))
                self.testLowTLoc_Display.setText(pExtremeT[1])
                self.testHighTLoc_Display.setText(pExtremeT[3]) 
                
                # Test body geometry
                if pGeomOut is not None:
                    if pGeomOut[-2] == 'Wedge':                    
                        if pGeomOut[0] is None:
                            self.testOblique_Display.setText('Detached')
                        else:
                            self.testOblique_Display.setNum(round(pGeomOut[0],1))
                        self.testDetach_Display.setNum(round(pGeomOut[1],1))
                        
                    elif pGeomOut[-1] == 'Normal shock':
                        print('Post-shock conditions (perfect gas):')
                        print('Pressure = '+str(pGeomOut[0].P)+' kPa')
                        print('Temperature = '+str(pGeomOut[0].T)+' K')
                        print('Density = '+str(pGeomOut[0].density)+' kg/m^3')
                        print('Velocity = '+str(pGeomOut[0].vel)+' m/s')
                
                # Variable slots:
                    # Driver pressure
                if self.driver_Tabs.currentIndex() == 1:
                    self.varTitle_Display.setText('Shock-specified mode')
                    self.varLabel1_Display.setText('Driver pressure (MPa)')
                    self.testLabel1_Display.setText('%.3g' % (pGas4.P/1e3)) # MPa
                    
                if self.detonation_PushBtn.isChecked():
                    self.varTitle_Display.setText('Detonation')
                    self.varLabel1_Display.setText('Post-det. pressure (MPa)')
                    self.testLabel1_Display.setText('%.3g' % (pGas4.P/1e3)) # MPa                    
                    self.varLabel2_Display.setText('Post-det. temperature (K)')
                    self.testLabel2_Display.setNum(round(pGas4.T,0)) # K
                    self.varLabel3_Display.setText('Post-det velocity (to left) (m/s)')
                    self.testLabel3_Display.setNum(round(abs(pGas4.vel),0)) # m/s
                    self.varLabel4_Display.setText('Detonation wave speed (m/s)')
                    self.testLabel4_Display.setNum(round((us_perf[2]),0)) # m/s
                    self.varLabel5_Display.setText('Detonation wave Mach number')
                    self.testLabel5_Display.setText('%.3g' % Ms_perf[2])
                    
                    
        # Call temporary output to terminal here (for debugging)
          
        return fig1,t_test,term,t_contact2,optimal
   
    def clear_output(self):
        """
        Clears all labels used for displaying output (by setting their text to a blank string).
        There is probably a better way to do this (e.g. adding all labels to some sort of group
        and performing this en masse) but I can't figure out if it's possible.
        
        If any new labels are added to the UI, ensure they are also added to this method.
        """
        self.remove_mpl() # clean up any previous plot
        
        self.column1_Display.setText('')
        self.column2_Display.setText('')
        self.column3_Display.setText('')
        self.column4_Display.setText('')
        self.column5_Display.setText('')
        self.column6_Display.setText('')
        
        self.testM_Display.setText('') 
        self.testu_Display.setText('') 
        self.testtime_Display.setText('') 
        self.perfM_Display.setText('')
        self.perfu_Display.setText('')
        self.perftime_Display.setText('')
        self.diffM_Display.setText('')
        self.diffu_Display.setText('')
        self.difftime_Display.setText('')
        self.testEnd_Display.setText('')
        self.perfEnd_Display.setText('')
        
        self.testT_Display.setText('')
        self.testrho_Display.setText('')
        self.testP_Display.setText('')
        self.testh_Display.setText('')
        self.testh0_Display.setText('')
        self.perfT_Display.setText('')
        self.perfrho_Display.setText('')
        self.perfP_Display.setText('')
        self.perfh_Display.setText('')
        self.perfh0_Display.setText('')
        self.diffT_Display.setText('')
        self.diffrho_Display.setText('')
        self.diffP_Display.setText('')
        self.diffh0_Display.setText('')
        
        self.testRem_Display.setText('')
        self.perfRem_Display.setText('')
        self.diffRem_Display.setText('')
        
        self.testk_Display.setText('')
        self.perfk_Display.setText('')
        self.diffk_Display.setText('')
        
        self.testus1_Display.setText('')
        self.testus2_Display.setText('')
        self.testMs1_Display.setText('')
        self.testMs2_Display.setText('')
        self.perfus1_Display.setText('')
        self.perfus2_Display.setText('')
        self.perfMs1_Display.setText('')
        self.perfMs2_Display.setText('')        
        self.diffus1_Display.setText('')
        self.diffus2_Display.setText('')
        self.diffMs1_Display.setText('')
        self.diffMs2_Display.setText('')
        
        self.ssr_Display.setText('') 
        self.mode_Display.setText('')
        self.perfssr_Display.setText('') 
        self.perfmode_Display.setText('')
        
        self.testLowT_Display.setText('')
        self.testHighT_Display.setText('')
        self.testLowTLoc_Display.setText('')
        self.testHighTLoc_Display.setText('')
        self.perfLowT_Display.setText('')
        self.perfHighT_Display.setText('')
        self.perfLowTLoc_Display.setText('')
        self.perfHighTLoc_Display.setText('')
        
        self.testOblique_Display.setText('')
        self.testDetach_Display.setText('')
        self.perfOblique_Display.setText('')
        self.perfDetach_Display.setText('')      
        self.diffOblique_Display.setText('')
        self.diffDetach_Display.setText('')
        
        self.varTitle_Display.setText('')
        self.varLabel1_Display.setText('')
        self.testLabel1_Display.setText('')
        self.perfLabel1_Display.setText('') 
        self.diffLabel1_Display.setText('')
        self.varLabel2_Display.setText('')
        self.testLabel2_Display.setText('')
        self.perfLabel2_Display.setText('') 
        self.diffLabel2_Display.setText('')
        self.varLabel3_Display.setText('')
        self.testLabel3_Display.setText('')
        self.perfLabel3_Display.setText('') 
        self.diffLabel3_Display.setText('')
        self.varLabel4_Display.setText('')
        self.testLabel4_Display.setText('')
        self.perfLabel4_Display.setText('') 
        self.diffLabel4_Display.setText('')
        self.varLabel5_Display.setText('')
        self.testLabel5_Display.setText('')
        self.perfLabel5_Display.setText('') 
        self.diffLabel5_Display.setText('')
        

##############################################################################
# External file/folder dialogs
    def input_file(self):
        """
        Opens a file dialog native to the operating system.
        The user is restricted to selecting Excel files.
        Dialog returns a string of the full path to the "Gas data filename" line edit.
        """
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setNameFilters(["Excel files (*.xls *.xlsx)"])
        
		 # Only return filename if "OK" pressed
        if dlg.exec_():
            filename = dlg.selectedFiles()[0]
            self.inFile_Edit.setText(filename)
            
        return
    
    def output_folder(self):
        """
        Opens a file dialog native to the operating system.
        The user is restricted to selecting folders (directories).
        Dialog returns a string of the full path to the "Output folder" line edit.
        """
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.Directory)
        
        if dlg.exec_():
            foldername = dlg.selectedFiles()[0]
            self.outFolder_Edit.setText(foldername)            
        return
    
##############################################################################
# Define formats for writing to output files   
    def single_file_output(self,all_data):
        """
        Takes the data provided by calculate() in the all_data variable and uses it
        to generate output files in single-calc mode.
        
        Two files are created: a txt file containing the full facility setup and the
        thermodynamic data of all gas states, and a pdf file containing a copy of the
        x-t diagram displayed in the UI.
        
        The actual writing is largely done by a collection of functions in exp_common_lib,
        this method dictates the order in which sections of data appear in the file.
        """
        # Set directory and filename from user input, open txt file
        outdir = self.outFolder_Edit.text()
        outfn = self.outName_Edit.text()        
        outfile = open(outdir+'/'+outfn+'_data.txt','w')
        
        # Save x-t diagram directly from matplotlib figure to pdf format
        all_data['x-t diagram'].savefig(outdir+'/'+outfn+'_xt.pdf', bbox_inches='tight')
        
        # Set up differently based on calculation mode
        realMode = False
        if self.passive_PushBtn.isChecked():
            gasDesc = self.gasDesc_passive
            pFreestream = all_data['Perfect gases'][6]
            detMode = False
            detMix = None
            if self.modelreal_radio.isChecked():
                rFreestream = all_data['Real gases'][6]              
        else:
            gasDesc = self.gasDesc_det
            pFreestream = all_data['Perfect gases'][7]
            detMode = True
            detMix = (self.detCompStr,self.phi)
            if self.modelreal_radio.isChecked():
                rFreestream = all_data['Real gases'][7]
        
        # Write real data ahead of perfect data
        if self.modelreal_radio.isChecked():
            realMode = True
            outfile.write('Real gas results:\n\n')
            ecl.write_facility(self.indepDict,all_data['Real gases'],all_data['Lengths'],
                               outfile,(realMode,detMode),
                               (self.drivbar_radio.isChecked(),
                               self.testtorr_radio.isChecked(),
                               self.acceltorr_radio.isChecked()),
                               detMix)
            ecl.write_key(rFreestream,all_data['Test time'][1],outfile,True)
            ecl.write_shock(all_data['Real shock'],outfile)
            ecl.write_modes(all_data['Enthalpy mode'][1],all_data['Termination mode'][1],
                        all_data['Sound speed ratio'][1],outfile)
            ecl.write_real_states(all_data['Real gases'],outfile,gasDesc)
            
            # Append real test body flow to end of data if present
            if 'Real test body' in all_data:
                ecl.write_real_wedge(all_data['Real test body'],outfile)
                
            # In real mode, only append perfect data if requested (via checkbox)
            if self.includePerf_CheckBox.isChecked():
                outfile.write('End of real gas data.\n\n')
                outfile.write('\n*************************\n*************************\n')                
            else:
                outfile.close()
                return
        
        # Write perfect data
        outfile.write('Perfect gas results:\n\n')
        ecl.write_facility(self.indepDict,all_data['Perfect gases'],all_data['Lengths'],
                               outfile,(realMode,detMode),
                               (self.drivbar_radio.isChecked(),
                               self.testtorr_radio.isChecked(),
                               self.acceltorr_radio.isChecked()),
                               detMix)
        ecl.write_key(pFreestream,all_data['Test time'][0],outfile,False)
        ecl.write_shock(all_data['Perfect shock'],outfile)
        ecl.write_modes(all_data['Enthalpy mode'][0],all_data['Termination mode'][0],
                        all_data['Sound speed ratio'][0],outfile)
        ecl.write_perf_states(all_data['Perfect gases'],outfile,gasDesc)
        
        # Append perfect test body flow to end of data if present
        if 'Perfect test body' in all_data:
            ecl.write_perf_wedge(all_data['Perfect test body'],outfile)
        
        outfile.close()
        
        return
    
    def multi_file_output(self):
        """
        Creates csv file with independent variable range and corresponding
        dependent variable values (i.e. the same data found in the multi-plot).
        
        Using this csv file, the user can externally recreate all the plots shown
        in multi-calc mode.
        
        Also generates a txt file that lists the facility state (i.e. the other 8 
        candidates for independent variable that were fixed as constants, as well as
        additional parameters such as gas species, and initial fuel mix, if detonation mode)
        """
        import csv 
        outdir = self.outFolder_Edit.text()
        outfn = self.outName_Edit.text()
        
        # Extract all dependent variable ranges for every dependent variable listed
        pOutputDict = {}
        rOutputDict = {}
        for depName in self.depVarDict:
            if depName in self.wedgeDepVars and not self.activeMulti[1]:
                # If the selected dependent variable relies on test body data (e.g. oblique shock angle)
                # but the existing calculation didn't have a test body, don't attempt lookup, and don't include in csv                        
                pass
            else:
                pdepVarRange,rdepVarRange,indepNameUnit,depNameUnit = self.multi_data_prep(depName)
                pOutputDict[depNameUnit] = pdepVarRange
                rOutputDict[depNameUnit] = rdepVarRange
        
        # Sort columns alphabetically by header
        pkeys = sorted(pOutputDict.keys())
        rkeys = sorted(rOutputDict.keys())
        
        # Add independent range after sorting, then place independent key in the first position
        # of the list of alphabetically sorted keys, to ensure independent data in first column of csv
        pOutputDict[indepNameUnit] = self.indepVarRange
        rOutputDict[indepNameUnit] = self.indepVarRange
        
        pkeys.insert(0,indepNameUnit)
        rkeys.insert(0,indepNameUnit)
        
        
        ######################################################################### 
        # Write to csv file (create separate files for real and perfect data, if needed)
        if self.modelperf_radio.isChecked() or self.includePerf_CheckBox.isChecked():
            with open(outdir+'/'+outfn+'_perf_dep_vars.csv', 'w') as csvfile:
               writer = csv.writer(csvfile)
               writer.writerow(pkeys)
               writer.writerows(zip(*[pOutputDict[pkey] for pkey in pkeys]))           
        
        if self.modelreal_radio.isChecked():
            with open(outdir+'/'+outfn+'_real_dep_vars.csv', 'w') as csvfile:
               writer = csv.writer(csvfile)
               writer.writerow(rkeys)
               writer.writerows(zip(*[rOutputDict[rkey] for rkey in rkeys]))
        
        ######################################################################### 
        # Create txt file listing the values of the other 8 (fixed) state variables
             
        txtfile = open(outdir+'/'+outfn+'_constants.txt','w')
        txtfile.write('Facility constants:\n\n')
        
        if self.detonation_PushBtn.isChecked():
            txtfile.write('Driver mode: Detonation\n\n')
            accelGas = self.allRange[0]['Perfect gases'][5]
        else:
            txtfile.write('Driver mode: Passive\n\n')
            accelGas = self.allRange[0]['Perfect gases'][4]
        
        # Can get gas names from any all_data entry in allRange (since always the same)
        # Also can always be from perfect (since same as real)
        drivGasName = getattr(self.allRange[0]['Perfect gases'][3],'gas')
        testGasName = getattr(self.allRange[0]['Perfect gases'][0],'gas')
        accelGasName = accelGas.gas
        
        txtfile.write('Driver gas: '+drivGasName+'\n')
        
        # Add mixture ratio information if detonation mode
        if self.detonation_PushBtn.isChecked():
            txtfile.write('\t Initial mole ratio: '+self.detCompStr+'\n')
            txtfile.write('\t Initial equivalence ratio: '+str(self.phi)+'\n\n')
            
        txtfile.write('Test gas: '+testGasName+'\n')
        txtfile.write('Acceleration gas: '+accelGasName+'\n\n')
        
        # Need to create ordered list, since keys aren't ordered
        orderedStateVars = ['drivP','drivT','drivL','testP','testT','testL','accelP','accelT','accelL']
        orderedStateVars.remove(self.indepStr)
        for i,stateVar in enumerate(orderedStateVars):
            _,__,formattedEntry = ecl.stateVar_formatted(stateVar,self.indepDict,
                                                    (self.drivbar_radio.isChecked(),
                                                    self.testtorr_radio.isChecked(),
                                                    self.acceltorr_radio.isChecked()) )
            
            # Put a blank line between driv, test, and accel states.
            # Do this by checking first letter of state var string and seeing if it matches
            # the previous entry - only works because 'driv', 'test', and 'accel' all
            # start with different characters!
            if (i > 0) and (orderedStateVars[i][0] != orderedStateVars[i-1][0]):
                txtfile.write('\n')
            txtfile.write(formattedEntry+'\n')

        txtfile.close()
        return
    
class MixerPopup(QtWidgets.QDialog,mixerTemplate):    
    """
    This class defines the popup window that is used to mix the desired ratios of
    fuel, oxidizer and diluent in detonation mode.
    
    Its layout is defined in the file "mixer.ui" which was generated using QtDesigner.
    
    The majority of the methods for this class are dedicated to updating the sliders so that
    the rest move correctly if one slider is altered by the user.
    """
    def __init__(self,parent=None,selectedGas=None,df=None,fodInit=None):
        from numpy import ceil,log10
        QtWidgets.QDialog.__init__(self,parent)
        self.setupUi(self)
        self.setModal(True) # Must close dialog before interacting with main window
        
        # Local parameters
        
        # The stoichiometric ratio for the current chemistry, taken from the input file
        self.stoich = float(df.loc[selectedGas,'Stoichiometry'])
        # The maximum value of phi that the slider can move to
        self.phiMax = 5.0
        # Number of increments for F/O/D sliders (e.g. 1000 = increments of 0.1%)
        self.fodInc = 10000
        # Number of decimal places to display for the percentages, based on number of increments
        self.dp = max(int(ceil(log10(self.fodInc/100))),0)        
        initDil = float(df.loc[selectedGas,'Diluent Percent'])
        
        self.fuelSlider.setMaximum(self.fodInc)
        self.oxSlider.setMaximum(self.fodInc)
        self.dilSlider.setMaximum(self.fodInc)
        self.phiSlider.setMaximum(self.fodInc)
        
        self.fuelSlider.setPageStep(int(self.fodInc/10))
        self.oxSlider.setPageStep(int(self.fodInc/10))
        self.dilSlider.setPageStep(int(self.fodInc/10))
        self.phiSlider.setPageStep(int(self.fodInc/10))
        
        self.dilSlider.setSliderPosition(initDil/100*self.fodInc)
        self.updatePosition(0,fodInit=fodInit)
        self.lockSlider()
        
        self.fuelNameLabel.setText(df.loc[selectedGas,'Fuel Name'])
        self.oxNameLabel.setText(df.loc[selectedGas,'Oxidizer Name'])
        self.dilNameLabel.setText(df.loc[selectedGas,'Diluent Name'])
        
        # Differentiate between the slider moving by being clicked and by having its
        # value altered by a function
        self.fuelSlider.valueChanged.connect(lambda: self.updatePosition(1))
        self.oxSlider.valueChanged.connect(lambda: self.updatePosition(2))
        self.dilSlider.valueChanged.connect(lambda: self.updatePosition(3))
        self.phiSlider.valueChanged.connect(lambda: self.updatePosition(4))
        
        # Connect the lock buttons to the lockSlider() method
        self.fuelLock.clicked.connect(self.lockSlider)
        self.oxLock.clicked.connect(self.lockSlider)
        self.dilLock.clicked.connect(self.lockSlider)
        self.phiLock.clicked.connect(self.lockSlider)
        
        # Cause the "Reset to stoichiometric conditions" button to set phi to 1
        self.resetStoich_PushBtn.clicked.connect(lambda: self.calcPosition_indepPhi(target=1))
        
        
    def updatePosition(self,changed,fodInit=None):
        """
        Updates the positions of all sliders when a change is made.
        5 behavioural options depending on the value of the variable "changed":
            0: used for initialization
            1: used when fuel slider moved by user
            2: used when oxidizer slider moved by user
            3: used when diluent slider moved by user
            4: used when phi slider moved by user
        """        
        if changed == 0:            
            if fodInit is None:
                # If no previous F/O/D state saved:
                # Initially set to stoichiometric conditions with no diluent
                self.setStoich()
            else:
                # Reinstate previous settings, including checkbox state
                # Discard phi (otherwise overconstrained since F/O/D all known)
                fuel,ox,dil,_,checked = fodInit
                # Important to block signals from sliders as appropriate, to prevent infinite loops
                # - otherwise changes made by these methods will emit their own signals and trigger
                # the methods again
                self.fuelSlider.blockSignals(True)
                self.oxSlider.blockSignals(True)
                self.dilSlider.blockSignals(True)
                self.phiSlider.blockSignals(True)
                
                self.fuelSlider.setSliderPosition(fuel*self.fodInc)
                self.oxSlider.setSliderPosition(ox*self.fodInc)
                self.dilSlider.setSliderPosition(dil*self.fodInc)
                self.phi_rev(self.calcPhi())
                
                if checked == 1:
                    self.fuelLock.setChecked(True)
                elif checked == 2:
                    self.oxLock.setChecked(True)
                elif checked == 3:
                    self.dilLock.setChecked(True)
                elif checked == 4:
                    self.phiLock.setChecked(True)
                
                # Make sure to unblock signals before returning
                self.fuelSlider.blockSignals(False)
                self.oxSlider.blockSignals(False)
                self.dilSlider.blockSignals(False)
                self.phiSlider.blockSignals(False)
                
                
        elif changed < 4:
            self.phiSlider.blockSignals(True)
            if changed == 1:
                if self.oxLock.isChecked():
                    self.calcPosition_depPhi(self.fuelSlider,self.dilSlider,self.oxSlider)
                elif self.dilLock.isChecked():
                    self.calcPosition_depPhi(self.fuelSlider,self.oxSlider,self.dilSlider)
                elif self.phiLock.isChecked():
                    self.calcPosition_lockedPhi(self.fuelSlider)
                    
            elif changed == 2:
                if self.fuelLock.isChecked():
                    self.calcPosition_depPhi(self.oxSlider,self.dilSlider,self.fuelSlider)
                elif self.dilLock.isChecked():
                    self.calcPosition_depPhi(self.oxSlider,self.fuelSlider,self.dilSlider)
                elif self.phiLock.isChecked():
                    self.calcPosition_lockedPhi(self.oxSlider)
                        
            elif changed == 3:
                if self.fuelLock.isChecked():
                    self.calcPosition_depPhi(self.dilSlider,self.oxSlider,self.fuelSlider)
                elif self.oxLock.isChecked():
                    self.calcPosition_depPhi(self.dilSlider,self.fuelSlider,self.oxSlider)
                elif self.phiLock.isChecked():
                    self.calcPosition_lockedPhi(self.dilSlider)
                    
            # Calculate and set phi
            if not self.phiLock.isChecked():
                self.phi_rev(self.calcPhi())
            self.phiSlider.blockSignals(False)
            
        elif changed == 4:
            self.calcPosition_indepPhi()
            
        self.updateValues()
        return
    
    
    def calcPosition_depPhi(self,indep,dep,fixed):
        """
        Calculates and sets the positions of the sliders when phi is a dependent variable
        (i.e. when phi isn't locked and isn't the slider being moved by the user).
        
        The 3 inputs are handles to the other 3 sliders:
            indep = independent slider (i.e. the one being moved by user)
            dep = dependent slider (unlocked but not being moved by user)
            fixed = fixed slider (the locked slider that cannot move)
        """
        dep.blockSignals(True)
        
        # Ensure the 3 sliders sum to 100% (measured in the fixed increment size, fodInc)
        dep.setSliderPosition(self.fodInc-fixed.value()-indep.value())
        
        # Don't allow the indep slider to keep increasing once dep slider can't decrease further
        if indep.value()+dep.value()+fixed.value()>self.fodInc:
            dep.setSliderPosition(dep.minimum())
            indep.setSliderPosition(self.fodInc-dep.value()-fixed.value())
        
        # Don't allow indep slider to keep changing if the PageStep phi slider value exceeded
        if self.calcPhi() > self.phiMax:
            minOx = self.fuelSlider.value()/(self.stoich*self.phiMax)
            if indep.objectName() == 'oxSlider':                
                indep.setSliderPosition(minOx)                
            elif indep.objectName() == 'dilSlider':
                indep.setSliderPosition(self.fodInc-minOx-self.fuelSlider.value())
            
        dep.blockSignals(False)
        return
    
    def calcPosition_indepPhi(self,target=None):
        """
        Calculates and sets the positions of the F/O/D sliders when phi is the
        independent variable. Respects the 'locked' status of each of the other
        3 sliders.
        
        Usually called when the phi slider is moved, but if the optional keyword
        argument 'target' is given as a value, will instead set phi to this value,
        then determine the other sliders as usual. This optional functionality
        is used to implement the button that resets everything to stoichiometric
        conditions, but it can also be used for phi != 1
        """
        if target is not None:
            self.fuelSlider.blockSignals(True)
            self.oxSlider.blockSignals(True)
            self.dilSlider.blockSignals(True)
            self.phiSlider.blockSignals(True)
            self.phi_rev(target)
            self.fuelSlider.blockSignals(False)
            self.oxSlider.blockSignals(False)
            self.dilSlider.blockSignals(False)
            self.phiSlider.blockSignals(False)
            
        
        if self.fuelLock.isChecked():            
            self.oxSlider.blockSignals(True)
            self.dilSlider.blockSignals(True)
            
            if self.phi_fwd() != 0:
                ox = self.fuelSlider.value()/(self.stoich*self.phi_fwd())
            
            oxMax = float(self.fodInc-self.fuelSlider.value()) # No diluent
            phiMin = (self.fuelSlider.value()/oxMax)/self.stoich
            if self.phi_fwd() < phiMin:
                self.phi_rev(phiMin)
                ox = oxMax
                
            self.oxSlider.setSliderPosition(ox)
            self.dilSlider.setSliderPosition(self.fodInc-self.fuelSlider.value()-ox)
            
            self.oxSlider.blockSignals(False)
            self.dilSlider.blockSignals(False)
        elif self.oxLock.isChecked():
            self.fuelSlider.blockSignals(True)
            self.dilSlider.blockSignals(True)
            
            fuel = self.oxSlider.value()*(self.stoich*self.phi_fwd())
            
            if fuel > self.fodInc-self.oxSlider.value(): # no diluent
                fuel = float(self.fodInc-self.oxSlider.value())
                phiMax = (fuel/self.oxSlider.value())/self.stoich
                self.phi_rev(phiMax)
            
            self.fuelSlider.setSliderPosition(fuel)
            self.dilSlider.setSliderPosition(self.fodInc-self.oxSlider.value()-fuel)
            
            self.fuelSlider.blockSignals(False)
            self.dilSlider.blockSignals(False)
        else:
            self.fuelSlider.blockSignals(True)
            self.oxSlider.blockSignals(True)
            
            ratio = self.stoich*self.phi_fwd()
            fuel = ratio/(ratio+1)*(self.fodInc-self.dilSlider.value())
            self.fuelSlider.setSliderPosition(fuel)
            self.oxSlider.setSliderPosition(self.fodInc-self.dilSlider.value()-fuel)
            
            self.fuelSlider.blockSignals(False)
            self.oxSlider.blockSignals(False)
        
        self.updateValues()
        return
        
    def calcPosition_lockedPhi(self,indep):
        """
        Calculates and sets the positions of the F/O/D sliders when the phi slider is locked.
        In this case, none of the other 3 can be locked (since the problem would be overconstrained).
        
        indep = the slider moved by the user
        """
        self.fuelSlider.blockSignals(True)
        self.oxSlider.blockSignals(True)
        self.dilSlider.blockSignals(True)
        # Define this constant for compactness
        K = self.phi_fwd()*self.stoich
        
        if indep.objectName() == 'fuelSlider':
            fuelMin = self.oxSlider.minimum()*K
            fuelMax = self.fodInc/(1+1/K)
            if indep.value() > fuelMax:
                indep.setSliderPosition(fuelMax)
            if indep.value() < fuelMin:
                indep.setSliderPosition(fuelMin)
            
            ox = indep.value()/K
            self.oxSlider.setSliderPosition(ox)
            self.dilSlider.setSliderPosition(self.fodInc-indep.value()-ox)
            
        elif indep.objectName() == 'oxSlider':          
            oxMax = self.fodInc/(1+K)
            if indep.value() > oxMax:
                indep.setSliderPosition(oxMax)

            fuel= indep.value()*K
            self.fuelSlider.setSliderPosition(fuel)
            self.dilSlider.setSliderPosition(self.fodInc-indep.value()-fuel)
            
        elif indep.objectName() == 'dilSlider':
            fuel = K/(K+1)*(self.fodInc-indep.value())
            self.fuelSlider.setSliderPosition(fuel)
            self.oxSlider.setSliderPosition(self.fodInc-indep.value()-fuel)
            
        self.fuelSlider.blockSignals(False)
        self.oxSlider.blockSignals(False)
        self.dilSlider.blockSignals(False)

    
    def updateValues(self):
        """
        Updates the percentage values printed next to each slider. Takes into account
        the defined size of the slider increment.
        """
        percentConv = 100./self.fodInc
        strformat = '%.'+str(self.dp)+'f'
        self.fuelCurrentVal.setText(strformat % (self.fuelSlider.value()*percentConv)+' %')
        self.oxCurrentVal.setText(strformat % (self.oxSlider.value()*percentConv)+' %')
        self.dilCurrentVal.setText(strformat % (self.dilSlider.value()*percentConv)+' %')   
        self.phiCurrentVal.setText('%.3f' % self.phi_fwd())
        return
        
    def lockSlider(self):
        """
        Disables slider if its corresponding lock button is checked.
        Have to use this "double negative" syntax because there is no 
        setDisabled method, only setEnabled
        """
        self.fuelSlider.setEnabled(not self.fuelLock.isChecked())
        self.oxSlider.setEnabled(not self.oxLock.isChecked())
        self.dilSlider.setEnabled(not self.dilLock.isChecked())
        self.phiSlider.setEnabled(not self.phiLock.isChecked())
        return
    
    def calcPhi(self):
        """
        Calculates phi from the known fuel and oxygen %, and the known stoichiometric ratio
        """
        phi = (float(self.fuelSlider.value())/float(self.oxSlider.value()))/self.stoich
        return phi
    
    def phi_fwd(self):
        """
        Forward transfer function for phi slider.
        Takes raw slider value (integer in range 0-fodInc) and returns actual phi value
        either between 0-1 or 1-max_phi.
        
        Here the phi slider is split into the above two regions at the midpoint.
        Within each region the progression is linear.
        """
        midpoint = self.fodInc/2
        raw = self.phiSlider.value()
        if raw < midpoint:
            phi = raw/float(midpoint) # convert to float to avoid int/int division (would discard decimal portion)
        else:
            phi = 1 + (raw/float(midpoint)-1)*(self.phiMax-1)            
        return phi
    
    def phi_rev(self,phi):
        """ 
        Reverse transfer function for phi slider. Takes an actual phi value and calculates
        the slider position based on the ranges defined in phi_fwd()
        """
        midpoint = self.fodInc/2
        if phi < 1:
            self.phiSlider.setSliderPosition(phi*midpoint)
        else:
            self.phiSlider.setSliderPosition(midpoint*((phi-1)/(self.phiMax-1)+1))
        return
    
    def setStoich(self):
        """ 
        Set sliders to stoichiometric conditions.
        Does not alter current diluent level.
        Used on first initialization, based on diluent level specified in input file.
        """
        self.fuelSlider.blockSignals(True)
        self.oxSlider.blockSignals(True)
        self.phiSlider.blockSignals(True)
        
        fuel = self.stoich/(self.stoich+1)*(self.fodInc-self.dilSlider.value())
        
        self.fuelSlider.setSliderPosition(fuel)
        self.oxSlider.setSliderPosition(self.fodInc-self.fuelSlider.value()-self.dilSlider.value())
        self.phi_rev(1)
        
        self.fuelSlider.blockSignals(False)
        self.oxSlider.blockSignals(False)
        self.phiSlider.blockSignals(False)
        self.updateValues()
        
        return
         
    def getSliderValues(self):
        """
        Reads the raw values of the 4 sliders and converts them to fractions.
        These are returned as variables, along with "checked", which stores the current
        state of the lock checkboxes.
        
        This method is used to pass these data back to the main window (since the local scope
        is lost once the popup window is closed).
        """
        xfuel = self.fuelSlider.value()/float(self.fodInc)
        xox = self.oxSlider.value()/float(self.fodInc)
        xdil = self.dilSlider.value()/float(self.fodInc)
        phi = self.phi_fwd()
        
        if self.fuelLock.isChecked():
            checked = 1
        elif self.oxLock.isChecked():
            checked = 2
        elif self.dilLock.isChecked():
            checked = 3
        elif self.phiLock.isChecked():
            checked = 4
            
        return xfuel,xox,xdil,phi,checked




#########################################################################       
# Execution
#########################################################################   

# The 'if' statement means the GUI will only be created if this is run
# as a top-level script, not if it is imported by something else
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv) # prerequisite to running
    app.setFont(QFont('Cantarell', 9))
    myWindow = Main(None) # create the application instance
    myWindow.showMaximized() # display the UI, maximize it upon initialization
    sys.exit(app.exec_()) # close properly
