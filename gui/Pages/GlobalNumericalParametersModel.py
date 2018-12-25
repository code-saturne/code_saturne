# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

import sys, unittest
from code_saturne.Base.XMLvariables import Model
from code_saturne.Base.XMLengine import *
from code_saturne.Base.XMLmodel import *


class GlobalNumericalParametersModel(Model):

    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        self.case = case
        self._numericalNode = self.case.xmlInitNode('numerical_parameters')
        self._restartTimeStep = self._numericalNode.xmlInitNode('restart_time_step_parameter')


    def defaultValues(self):
        default = {}

        default['restart_time_step']                    = 'off'
        default['potential_state']                      = 'off'
        default['faces_reconstruction']                 = 'on'
        default['multigrid']                            = 'on'
        default['sub_cycles']                           = 1
        default['alpha-pressure_cycles']                = 50
        default['sum_alpha']                            = 1e-8
        default['velocity_update']                      = 'pressure_gradient_increment'
        default['max_number_restart']                   = 1
        default['time_step_splitting']                  = 1.0
        default['pressure_relaxation']                  = 0.6
        default['pressure_symetrisation']               = 'off'
        default['MinPressure']                          = -1.e15
        default['MaxPressure']                          = 1.e15
        default['upwind_scheme']                        = 'on'
        default['stop_no_convergence']                  = 'off'
        default['velocity_predictor_algorithm_std']     = 'coupled_difvitc'
        default['velocity_predictor_algorithm_bubble']  = 'mean_velocity_relative_velocity'
        default['pressure_gradient']                    = 'mass_ponderation'
        default['regulate_bad_cells']                   = 'off'
        return default


    @Variables.undoLocal
    def setRestartTimeStep(self, status):
        """
        Set if restart time step
        only if it 's different of default value
        """
        self.isOnOff(status)
        if status == self.defaultValues()['restart_time_step'] :
            self._numericalNode.xmlRemoveChild('restart_time_step')
        else:
            node = self._numericalNode.xmlInitNode('restart_time_step')
            node['status'] = status


    @Variables.noUndo
    def getRestartTimeStep(self):
        """
        Get if restart time step
        """
        value = self.defaultValues()['restart_time_step']
        node = self._numericalNode.xmlGetNode('restart_time_step')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setPotentielState(self, status):
        """
        Set if potentiel state
        only if it 's different of default value
        """
        self.isOnOff(status)
        if status == self.defaultValues()['potential_state'] :
            self._numericalNode.xmlRemoveChild('potential_state')
        else:
            node = self._numericalNode.xmlInitNode('potential_state')
            node['status'] = status


    @Variables.noUndo
    def getPotentielState(self):
        """
        Get if potentiel state
        """
        value = self.defaultValues()['potential_state']
        node = self._numericalNode.xmlGetNode('potential_state')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setFacesReconstruction(self, status):
        """
        Set maximum internal faces reconstruction status
        only if it 's different of default value
        """
        self.isOnOff(status)
        if status == self.defaultValues()['faces_reconstruction'] :
            self._numericalNode.xmlRemoveChild('max_faces_reconstruction')
        else:
            node = self._numericalNode.xmlInitNode('max_faces_reconstruction')
            node['status'] = status


    @Variables.noUndo
    def getFacesReconstruction(self):
        """
        Get maximum internal faces reconstruction status
        """
        value = self.defaultValues()['faces_reconstruction']
        node = self._numericalNode.xmlGetNode('max_faces_reconstruction')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setRegulateBadCells(self, status):
        """
        Activate bad cells regulations.
        """
        self.isOnOff(status)
        if status == self.defaultValues()['regulate_bad_cells']:
            self._numericalNode.xmlRemoveChild('regulate_bad_cells')
        else:
            node = self._numericalNode.xmlInitNode('regulate_bad_cells')
            node['status'] = status


    @Variables.noUndo
    def getRegulateBadCElls(self):
        """
        Get bad cells regulation activation status.
        """
        value = self.defaultValues()['regulate_bad_cells']
        node  = self._numericalNode.xmlGetNode('regulate_bad_cells')
        if node:
            value = node['status']
        return value


    @Variables.undoLocal
    def setMultigridStatus(self, status):
        """
        Set multigrid status for pressure
        only if it 's different of default value
        """
        self.isOnOff(status)
        if status == self.defaultValues()['multigrid'] :
            self._numericalNode.xmlRemoveChild('pressure_multigrid')
        else:
            node = self._numericalNode.xmlInitNode('pressure_multigrid')
            node['status'] = status


    @Variables.noUndo
    def getMultigridStatus(self):
        """
        Get multigrid status for pressure
        """
        value = self.defaultValues()['multigrid']
        node = self._numericalNode.xmlGetNode('pressure_multigrid')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setVelocityUpdate(self, model):
        """
        Set velocity update methods
        only if it 's different of default value
        """
        self.isInList(model, ('pressure_gradient_increment', 'flow_rate', 'flumas_increment', 'conv_diff_equation'))
        if model == self.defaultValues()['velocity_update']:
            self._numericalNode.xmlRemoveChild('velocity_update')
        else:
            node = self._numericalNode.xmlInitNode('velocity_update')
            node['choice'] = model


    @Variables.noUndo
    def getVelocityUpdate(self):
        """
        Get velocity update methods
        """
        value = self.defaultValues()['velocity_update']
        node = self._numericalNode.xmlGetNode('velocity_update')
        if node:
            value = node['choice']
        return value


    @Variables.undoLocal
    def setAlphaPressureCycles(self, value):
        """
        Set alpha pressure cycles
        only if it 's different of default value
        """
        self.isInt(value)
        if value != self.defaultValues()['alpha-pressure_cycles']:
            self._numericalNode.xmlSetData('alpha-pressure_cycles', value)
        else:
            self._numericalNode.xmlRemoveChild('alpha-pressure_cycles')


    @Variables.noUndo
    def getAlphaPressureCycles(self):
        """
        Get alpha pressure cycles
        """
        value = self._numericalNode.xmlGetInt('alpha-pressure_cycles')
        if value == None:
            value = self.defaultValues()['alpha-pressure_cycles']
        return value


    @Variables.undoLocal
    def setSumAlpha(self, value):
        """
        Set value for 1- sum(alpha)
        only if it 's different of default value
        """
        self.isPositiveFloat(value)

        if value != self.defaultValues()['sum_alpha']:
            self._numericalNode.xmlSetData('sum_alpha', value)
        else:
            self._numericalNode.xmlRemoveChild('sum_alpha')


    @Variables.noUndo
    def getSumAlpha(self):
        """
        Get value for 1- sum(alpha)
        """

        value = self._numericalNode.xmlGetDouble('sum_alpha')
        if value == None:
            value = self.defaultValues()['sum_alpha']
        return value


    @Variables.undoLocal
    def setMaxNumberOfRestart(self, value):
        """
        Set max number of restart
        """
        self.isInt(value)
        if value != self.defaultValues()['max_number_restart']:
            self._restartTimeStep.xmlSetData('max_number_restart', value)
        else:
            self._restartTimeStep.xmlRemoveChild('max_number_restart')


    @Variables.noUndo
    def getMaxNumberOfRestart(self):
        """
        Get max number of restart
        """
        value = self._restartTimeStep.xmlGetInt('max_number_restart')
        if value == None:
            value = self.defaultValues()['max_number_restart']
        return value


    @Variables.undoLocal
    def setTimeSplit(self, value):
        """
        Set Time Step Splitting
        """
        self.isPositiveFloat(value)

        if value != self.defaultValues()['time_step_splitting']:
            self._restartTimeStep.xmlSetData('time_step_splitting', value)
        else:
            self._restartTimeStep.xmlRemoveChild('time_step_splitting')


    @Variables.noUndo
    def getTimeSplit(self):
        """
        Get Time Step Splitting
        """
        value = self._restartTimeStep.xmlGetDouble('time_step_splitting')
        if value == None:
            value = self.defaultValues()['time_step_splitting']
        return value


    @Variables.undoLocal
    def setPressureRelaxation(self, value):
        """
        Set pressure relaxation
        """
        self.isPositiveFloat(value)

        if value != self.defaultValues()['pressure_relaxation']:
            self._restartTimeStep.xmlSetData('pressure_relaxation', value)
        else:
            self._restartTimeStep.xmlRemoveChild('pressure_relaxation')


    @Variables.noUndo
    def getPressureRelaxation(self):
        """
        Get pressure relaxation
        """
        value = self._restartTimeStep.xmlGetDouble('pressure_relaxation')
        if value == None:
            value = self.defaultValues()['pressure_relaxation']
        return value


    @Variables.undoLocal
    def setPressureSymetrisation(self, status):
        """
        Set pressure matrix symetrisation
        """
        self.isOnOff(status)
        if status == self.defaultValues()['pressure_symetrisation']:
            self._numericalNode.xmlRemoveChild('pressure_symetrisation')
        else:
            node = self._numericalNode.xmlInitNode('pressure_symetrisation')
            node['status'] = status


    @Variables.noUndo
    def getPressureSymetrisation(self):
        """
        Get pressure matrix symetrisation status
        """
        value = self.defaultValues()['pressure_symetrisation']
        node = self._numericalNode.xmlGetNode('pressure_symetrisation')
        if node:
            value = node['status']
        return value


    @Variables.undoLocal
    def setMinPressure(self, value):
        """
        Set Min Pressure for clipping
        """
        self.isFloat(value)

        if value != self.defaultValues()['MinPressure']:
            self._numericalNode.xmlSetData('MinPressure', value)
        else:
            self._numericalNode.xmlRemoveChild('MinPressure')


    @Variables.noUndo
    def getMinPressure(self):
        """
        Get Min Pressure for clipping
        """
        value = self._numericalNode.xmlGetDouble('MinPressure')
        if value == None:
            value = self.defaultValues()['MinPressure']
        return value


    @Variables.undoLocal
    def setMaxPressure(self, value):
        """
        Set Max Pressure for clipping
        """
        self.isFloat(value)

        if value != self.defaultValues()['MaxPressure']:
            self._numericalNode.xmlSetData('MaxPressure', value)
        else:
            self._numericalNode.xmlRemoveChild('MaxPressure')


    @Variables.noUndo
    def getMaxPressure(self):
        """
        Get Max Pressure for clipping
        """
        value = self._numericalNode.xmlGetDouble('MaxPressure')
        if value == None:
            value = self.defaultValues()['MaxPressure']
        return value


    @Variables.undoLocal
    def setUpwindScheme(self, status):
        """
        Set status for upwind scheme for alpha and energy
        """
        self.isOnOff(status)
        if status == self.defaultValues()['upwind_scheme'] :
            self._restartTimeStep.xmlRemoveChild('upwind_scheme')
        else:
            node = self._restartTimeStep.xmlInitNode('upwind_scheme')
            node['status'] = status


    @Variables.noUndo
    def getUpwindScheme(self):
        """
        Get status for upwind scheme for alpha and energy
        """
        value = self.defaultValues()['upwind_scheme']
        node = self._restartTimeStep.xmlGetNode('upwind_scheme')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setStopNoConvergence(self, status):
        """
        Set the stop if no convergence status
        """
        self.isOnOff(status)
        if status == self.defaultValues()['stop_no_convergence'] :
            self._numericalNode.xmlRemoveChild('stop_no_convergence')
        else:
            node = self._numericalNode.xmlInitNode('stop_no_convergence')
            node['status'] = status


    @Variables.noUndo
    def getStopNoConvergence(self):
        """
        Get stop if no convergence status
        """
        value = self.defaultValues()['stop_no_convergence']
        node = self._numericalNode.xmlGetNode('stop_no_convergence')
        if node :
            value = node['status']
        return value


    @Variables.undoLocal
    def setVelocityPredictorAlgo(self, value):
        """
        Set velocity predictor algorithme
        """
        self.isInList(value, ('standard_difvit', 'coupled_difvitc', 'mean_velocity_relative_velocity'))
        node = self._numericalNode.xmlInitNode('velocity_predictor_algorithm')
        node['choice'] = value


    @Variables.noUndo
    def getVelocityPredictorAlgo(self):
        """
        Get velocity predictor algorithme
        """
        from code_saturne.Pages.MainFieldsModel import MainFieldsModel

        node = self._numericalNode.xmlGetNode('velocity_predictor_algorithm')
        if node == None:
            model = ""
            if (MainFieldsModel(self.case).getPredefinedFlow() == "boiling_flow") :
                model = self.defaultValues()['velocity_predictor_algorithm_bubble']
            else:
                model = self.defaultValues()['velocity_predictor_algorithm_std']
            self.setVelocityPredictorAlgo(model)
            node = self._numericalNode.xmlGetNode('velocity_predictor_algorithm')
        model = node['choice']

        return model


    @Variables.undoLocal
    def setPressureGradient(self, value):
        """
        Set pressure gradient
        """
        self.isInList(value, ('mass_ponderation', 'standard', 'controlled_mass_pound', 'momentum_pound', 'gravity_momentum'))
        if value == self.defaultValues()['pressure_gradient']:
            self._numericalNode.xmlRemoveChild('pressure_gradient')
        else:
            node = self._numericalNode.xmlInitNode('pressure_gradient')
            node['choice'] = value


    @Variables.noUndo
    def getPressureGradient(self):
        """
        Get pressure gradient
        """
        value = self.defaultValues()['pressure_gradient']
        node = self._numericalNode.xmlGetNode('pressure_gradient')
        if node:
            value = node['choice']
        return value


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class GlobalNumericalParametersTestCase(ModelTest):
    """
    """
    def checkGlobalNumericalParametersInstantiation(self):
        """Check whether the GlobalNumericalParametersModel class could be instantiated"""
        model = None
        model = GlobalNumericalParametersModel(self.case)
        assert model != None, 'Could not instantiate GlobalNumericalParametersModel'


    def checkGetandSetRestartTimeStep(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Restart Time Step"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setRestartTimeStep('on')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <restart_time_step status="on"/>
                 </numerical_parameters>'''
        assert mdl.getRestartTimeStep() == 'on',\
            'Could not get RestartTimeStep'


    def checkGetandSetPotentielState(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Potentiel State"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setPotentielState('on')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <potential_state status="on"/>
                 </numerical_parameters>'''
        assert mdl.getPotentielState() == 'on',\
            'Could not get Potentiel State'


    def checkGetandSetFacesReconstruction(self):
        """Check whether the GlobalNumericalParametersModel class could get and set InternalFacesReconstruction"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setFacesReconstruction('off')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <max_faces_reconstruction status="off"/>
                 </numerical_parameters>'''
        assert mdl.getFacesReconstruction() == 'off',\
            'Could not get InternalFacesReconstruction'


    def checkGetandSetMultigridStatus(self):
        """Check whether the GlobalNumericalParametersModel class could get and set MultigridStatus"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setMultigridStatus('on')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <pressure_multigrid status="on"/>
                 </numerical_parameters>'''
        assert mdl.getMultigridStatus() == 'on',\
            'Could not get MultigridStatus'


    def checkGetandSetVelocityUpdate(self):
        """Check whether the GlobalNumericalParametersModel class could get and set VelocityUpdate"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setVelocityUpdate('flow_rate')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <velocity_update choice="flow_rate"/>
                 </numerical_parameters>'''
        assert mdl.getVelocityUpdate() == 'flow_rate',\
            'Could not get Velocity Update'


    def checkGetandSetAlphaPressureCycles(self):
        """Check whether the GlobalNumericalParametersModel class could get and set AlphaPressureCycles"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setAlphaPressureCycles(3)
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <alpha-pressure_cycles>
                                 3
                         </alpha-pressure_cycles>
                 </numerical_parameters>'''
        assert mdl.getAlphaPressureCycles() == 3,\
            'Could not get Alpha Pressure Cycles'


    def checkGetandSetSumAlpha(self):
        """Check whether the GlobalNumericalParametersModel class could get and set SumAlpha"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setSumAlpha(5.5)
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <sum_alpha>
                                 5.5
                         </sum_alpha>
                 </numerical_parameters>'''
        assert mdl.getSumAlpha() == 5.5,\
            'Could not get SumAlpha'


    def checkGetandSetMaxNumberOfRestart(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Max Number Of Restart"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setMaxNumberOfRestart(60)
        doc = '''<restart_time_step_parameter>
                         <max_number_restart>
                                 60
                         </max_number_restart>
                 </restart_time_step_parameter>'''
        assert mdl.getMaxNumberOfRestart() == 60,\
            'Could not get Max Number Of Restart'


    def checkGetandSetTimeSplit(self):
        """Check whether the GlobalNumericalParametersModel class could get and set TimeSplit"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setTimeSplit(54)
        doc = '''<restart_time_step_parameter>
                         <time_step_splitting>
                                 54
                         </time_step_splitting>
                 </restart_time_step_parameter>'''
        assert mdl.getTimeSplit() == 54,\
            'Could not get TimeSplit'


    def checkGetandSetPressureRelaxation(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Pressure Relaxation"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setPressureRelaxation(118)
        doc = '''<restart_time_step_parameter>
                         <pressure_relaxation>
                                 118
                         </pressure_relaxation>
                 </restart_time_step_parameter>'''
        assert mdl.getPressureRelaxation() == 118,\
            'Could no get Pressure Relaxation'


    def checkGetandSetPressureSymetrisation(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Implicitation Coefficient"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setPressureSymetrisation('on')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <pressure_symetrisation status="on"/>
                 </numerical_parameters>'''
        assert mdl.getPressureSymetrisation() == self.xmlNodeFromString(doc),\
            'Could not set pressure symetrisation'
        assert mdl.getPressureSymetrisation() == 'on',\
            'Could not get pressure symetrisation'


    def checkGetandSetMinMaxPressure(self):
        """Check whether the GlobalNumericalParametersModel class could get and set Min Max Pressure"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setMinPressure(118712)
        mdl.setMaxPressure(191258)
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <MinPressure>
                                 118712
                         </MinPressure>
                         <MaxPressure>
                                 191258
                         </MaxPressure>
                 </numerical_parameters>'''
        assert mdl.getMinPressure() == 118712,\
            'Could not get Min Pressure'
        assert mdl.getMaxPressure() == 191258,\
            'Could not get Max Pressure'


    def checkGetandSetUpwindScheme(self):
        """Check whether the GlobalNumericalParametersModel class could get and set UpwindScheme"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setUpwindScheme('off')
        doc = '''<restart_time_step_parameter>
                        <upwind_scheme status="off"/>
                 </restart_time_step_parameter>'''
        assert mdl.getUpwindScheme() == 'off',\
            'Could not get UpwindScheme'


    def checkGetandSetStopNoConvergence(self):
        """Check whether the GlobalNumericalParametersModel class could get and set StopNoConvergence"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setStopNoConvergence('on')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <stop_no_convergence status="on"/>
                 </numerical_parameters>'''
        assert mdl.getStopNoConvergence() == 'on',\
            'Could not get StopNoConvergence'


    def checkGetandSetVelocityPredictorAlgo(self):
        """Check whether the GlobalNumericalParametersModel class could get and set VelocityPredictorAlgo"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setVelocityPredictorAlgo('coupled_difvitc')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <velocity_predictor_algorithm choice="coupled_difvitc"/>
                 </numerical_parameters>'''
        assert mdl.getVelocityPredictorAlgo() == 'coupled_difvitc',\
            'Could not get VelocityPredictorAlgo'


    def checkGetandSetPressureGradient(self):
        """Check whether the GlobalNumericalParametersModel class could get and set PressureGradient"""
        mdl = GlobalNumericalParametersModel(self.case)
        mdl.setPressureGradient('standard')
        doc = '''<numerical_parameters>
                         <restart_time_step_parameter/>
                         <pressure_gradient choice="standard"/>
                 </numerical_parameters>'''
        assert mdl.getPressureGradient() == 'standard',\
            'Could not get PressureGradient'


def suite():
    testSuite = unittest.makeSuite(GlobalNumericalParametersTestCase, "check")
    return testSuite


def runTest():
    print("GlobalNumericalParametersTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
