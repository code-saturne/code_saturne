# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

"""
This module defines the turbulence model data management.

This module contains the following classes and function:
- TurbulenceModel
- TurbulenceTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.NumericalParamGlobalModel import NumericalParamGlobalModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# Turbulence model class
#-------------------------------------------------------------------------------

class TurbulenceModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc about Turbulence
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_lagr   = self.case.xmlGetNode('lagrangian')
        self.node_coal   = self.node_models.xmlGetChildNode('solid_fuels',
                                                            'model')
        self.node_joule  = self.node_models.xmlGetChildNode('joule_effect',
                                                            'model')
        self.node_gas    = self.node_models.xmlGetChildNode('gas_combustion',
                                                            'model')
        self.node_turb   = self.node_models.xmlInitChildNode('turbulence',
                                                             'model')
        self.node_bc     = self.case.xmlGetNode('boundary_conditions')
        self.node_ana    = self.case.xmlInitNode('analysis_control')
        self.node_prof   = self.node_ana.xmlInitNode('profiles')
        self.node_ava    = self.node_ana.xmlInitNode('time_averages')

        self.__turbModel = ('off',
                            'mixing_length',
                            'k-epsilon',
                            'k-epsilon-PL',
                            'Rij-epsilon',
                            'Rij-SSG',
                            'Rij-EBRSM',
                            'v2f-BL-v2/k',
                            'k-omega-SST',
                            'Spalart-Allmaras',
                            'LES_Smagorinsky',
                            'LES_dynamique',
                            'LES_WALE')

        self.__turbLESmodel = ('LES_Smagorinsky',
                               'LES_dynamique',
                               'LES_WALE')

        self.__allVariables = ['r11',
                               'r22',
                               'r33',
                               'r12',
                               'r13',
                               'r23',
                               'rij',
                               'k',
                               'epsilon',
                               'phi',
                               'alpha',
                               'omega',
                               'nu_tilda']


    def turbulenceModels(self):
        """
        Return all turbulence models known.

        @rtype: C{Tuple}
        @return: list of turbulence models
        """
        return self.__turbModel


    def LESmodels(self):
        """
        Return only LES turbulence models.

        @rtype: C{Tuple}
        @return: list of LES models
        """
        return self.__turbLESmodel


    def RANSmodels(self):
        """
        Return only RANS turbulence models.

        @rtype: C{Tuple}
        @return: list of RANS models
        """
        l = []
        for m in self.__turbModel:
            if m not in  self.__turbLESmodel and \
               m not in ("off", "mixing_length"):
                l.append(m)
        return l


    def defaultTurbulenceValues(self):
        """
        Return in a dictionnary which contains default values.

        @rtype: C{Dictionary}
        @return: default values
        """
        default = {}
        default['turbulence_model']          = "k-epsilon-PL"
        default['length_scale']              = 1.0
        default['turbulent_diffusion_model'] = 'shir'
        default['gravity_terms']             = "on"
        default['reference_velocity']        = 1.0
        default['reference_length_choice']   = 'automatic'
        default['reference_length']          = 1.0

        return default


    def turbulenceModelsList(self):
        """
        Create a tuple with the turbulence models allowed by the calculation
        features (multi-phases model, and reactive flow models).

        @rtype: C{Tuple}
        @return: list of avalaible models
        """
        turbList = self.__turbModel

        if self.node_lagr and self.node_lagr['model'] == 'on':
            turbList = self.RANSmodels()
            turbList.insert(0, "off")

        if self.node_gas and self.node_gas['model'] != 'off':
            turbList = self.RANSmodels()

        if self.node_coal and self.node_coal['model'] != 'off':
            turbList = ('off', 'k-epsilon', 'k-epsilon-PL')

        return turbList

    def wall_function_types(self):
        """
        Return dictionary of wall function types.
        """
        # TODO: switch keys to explicit names (needs to be done in C also).

        turb_model = self.node_turb['model']

        all_types = {'0': "No wall function",
                     '1': "One scale power law",
                     '2': "One scale log law",
                     '3': "Two scales log law",
                     '4': "Scalable wall function",
                     '5': "Two scales Van Driest",
                     '6': "Two scales smooth/rough)",
                     '7': "All y+ (2-scale model)"}

        wall_f_default = ''
        wall_f_types = {'-1': "Default"}

        wall_f_list = []

        if turb_model in ('off', 'v2f-BL-v2/k'):
            wall_f_list = ['0']
            wall_f_default = '0'

        elif turb_model == 'mixing_length':
            wall_f_list = ['2']
            wall_f_default = '2'

        elif turb_model == 'Rij-EBRSM':
            # Combo - piecewise laws (iwallf=2,3) unavailable through the GUI
            wall_f_list = ['0', '7']
            wall_f_default = '7'

        elif turb_model == 'k-omega-SST':
            # Combo - power law (iwallf=1) unavailable through the GUI
            wall_f_list = ['0', '2', '3', '7', '4']
            wall_f_default = '3'

        elif turb_model == 'Spalart-Allmaras':
            wall_f_list = ['2']
            wall_f_default = '2'

        else:
            # Combo - power law (iwallf=1) unavailable through the GUI
            wall_f_list = ['0', '2', '3', '4']
            wall_f_default = '3'

        for o in wall_f_list:
            wall_f_types[o] = all_types[o]

        return wall_f_types, wall_f_default


    def __removeVariablesAndProperties(self, varList, propName):
        """
        Delete variables and property that are useless accordingly to the model.
        """
        for v in self.__allVariables:
            if v not in varList:
                self.node_turb.xmlRemoveChild('variable', name=v)
                for node in self.node_prof.xmlGetNodeList('profile'):
                    node.xmlRemoveChild('var_prop', name=v)
                for node in self.node_ava.xmlGetNodeList('time_average'):
                    node.xmlRemoveChild('var_prop', name=v)
        self.node_turb.xmlRemoveChild('property', name=propName)
        for node in self.node_prof.xmlGetNodeList('profile'):
            node.xmlRemoveChild('var_prop', name=propName)
        for node in self.node_ava.xmlGetNodeList('time_average'):
            node.xmlRemoveChild('var_prop', name=propName)


    @Variables.undoGlobal
    def setTurbulenceModel(self, model_turb):
        """
        Input ITURB
        """
        self.isInList(model_turb, self.turbulenceModelsList())

        model_turb_old = self.node_turb['model']
        self.node_turb['model'] = model_turb

        NumericalParamGlobalModel(self.case).setTimeSchemeOrder(1)

        if model_turb == 'mixing_length':
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__removeVariablesAndProperties([], 'smagorinsky_constant^2')
            self.node_turb.xmlRemoveChild('wall_function')

        elif model_turb in ('k-epsilon', 'k-epsilon-PL'):
            lst = ('k', 'epsilon')
            for v in lst:
                self.setNewVariable(self.node_turb, v, label=v)
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')

        elif model_turb in ('Rij-epsilon', 'Rij-SSG', 'Rij-EBRSM'):
            # Rij is now considered as a tensor (vector of length 6,
            # since it is symmetric)
            lst = ['rij', 'epsilon']
            if model_turb == 'Rij-EBRSM':
                lst.append('alpha')
            for v in lst:
                if v == 'rij':
                    self.setNewVariable(self.node_turb, 'rij', label='Rij', dim='6')
                else:
                    self.setNewVariable(self.node_turb, v, label=v)
                if v == 'alpha':
                    v_n = node.xmlGetNode('variable', name=v)
                    v_n['_convect'] = 'no'
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')

        elif model_turb in self.LESmodels():
            if model_turb == 'LES_dynamique':
                self.setNewProperty(self.node_turb, 'smagorinsky_constant^2')
            else:
                self.__removeVariablesAndProperties([], 'smagorinsky_constant^2')

            if self.node_lagr['model'] != "off":
                lst = ('k', 'epsilon')
                for v in lst:
                    self.setNewVariable(self.node_turb, v, label=v)

            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.node_turb.xmlRemoveChild('wall_function')

            from code_saturne.model.TimeStepModel import TimeStepModel
            TimeStepModel(self.case).setTimePassing(0)
            del TimeStepModel

            NumericalParamGlobalModel(self.case).setTimeSchemeOrder(2)

        elif model_turb == 'v2f-BL-v2/k':
            lst = ('k', 'epsilon', 'phi', 'alpha')
            for v in lst:
                self.setNewVariable(self.node_turb, v, label=v)
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')
            wall_function = 0
            self.setWallFunction(wall_function)

        elif model_turb == 'k-omega-SST':
            lst = ('k', 'omega')
            for v in lst:
                self.setNewVariable(self.node_turb, v, label=v)
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')

        elif model_turb == 'Spalart-Allmaras':
            lst = ('nu_tilda')
            self.setNewVariable(self.node_turb, 'nu_tilda', label='nu_tilda')
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')
            self.node_turb.xmlRemoveChild('wall_function')

        else:
            model_turb = 'off'
            self.node_turb.xmlRemoveChild('variable')
            self.node_turb.xmlRemoveChild('property')
            self.node_turb.xmlRemoveChild('wall_function')
            self.__removeVariablesAndProperties([], 'smagorinsky_constant^2')
            wall_function = 0
            self.setWallFunction(wall_function)

        # If the user changes the turbulence model, reset the default
        # numerical parameters
        if model_turb != model_turb_old:
            from code_saturne.model.NumericalParamEquationModel import NumericalParamEquationModel
            NumericalParamEquationModel(self.case).setSchemeDefaultValues()
            del NumericalParamEquationModel

        DefineUserScalarsModel(self.case).setTurbulentFluxGlobalModel(model_turb)


    def __updateInletsForTurbulence(self):
        """
        Put boundaries conditions if it's necessary
        """
        from code_saturne.model.Boundary import Boundary
        for nodbc in self.node_bc.xmlGetChildNodeList('inlet'):
            model = Boundary('inlet', nodbc['label'], self.case)
            model.getTurbulenceChoice()

        del Boundary

    @Variables.noUndo
    def getTurbulenceModel(self):
        """
        Return the current turbulence model.
        """
        model = self.node_turb['model']
        if model not in self.turbulenceModelsList():
            model = self.defaultTurbulenceValues()['turbulence_model']
            self.setTurbulenceModel(model)
        return model


    @Variables.undoLocal
    def setLengthScale(self, l_scale):
        """
        Input XLOMLG.
        """
        self.isGreater(l_scale, 0.0)
        if self.getTurbulenceModel() == 'mixing_length':
            self.node_turb.xmlSetData('mixing_length_scale', l_scale)


    @Variables.noUndo
    def getLengthScale(self):
        """
        Return XLOMLG.
        """
        l_scale = self.node_turb.xmlGetDouble('mixing_length_scale')
        if l_scale is None:
            l_scale = self.defaultTurbulenceValues()['length_scale']
            self.setLengthScale(l_scale)
        return l_scale


    @Variables.undoLocal
    def setVelocity(self, value):
        """
        Set value of reference velocity into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_turb.xmlSetData('reference_velocity',value)


    @Variables.noUndo
    def getVelocity(self):
        """
        Return the value of reference velocity.
        """
        value = self.node_turb.xmlGetDouble('reference_velocity')
        if value is None:
            value = self.defaultTurbulenceValues()['reference_velocity']
            self.setVelocity(value)

        return value


    @Variables.undoLocal
    def setLengthChoice(self, choice):
        """
        Set the Length choice.
        """
        self.isInList(choice, ['automatic','prescribed'])

        node_init = self.node_turb.xmlInitNode('reference_length')
        node_init['choice'] = choice
        if choice == 'automatic':
            self.node_turb.xmlRemoveChild('reference_length')


    @Variables.noUndo
    def getLengthChoice(self):
        """
        Get the Length choice.
        """
        node_init = self.node_turb.xmlInitNode('reference_length')
        choice = node_init['choice']
        if choice is None:
            choice = self.defaultTurbulenceValues()['reference_length_choice']
            self.setLengthChoice(choice)
        return choice


    @Variables.undoLocal
    def setLength(self, value):
        """
        Set value of reference length into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_turb.xmlSetData('reference_length',value)


    @Variables.noUndo
    def getLength(self):
        """
        Return the value of reference length.
        """
        value = self.node_turb.xmlGetDouble('reference_length')
        if value is None:
            value = self.defaultTurbulenceValues()['reference_length']
            self.setLength(value)

        return value


    @Variables.noUndo
    def getWallFunction(self):
        """
        Return wall function from advanced options.
        """
        wall_function = self.node_turb.xmlGetInt('wall_function')
        if wall_function is None:
            wall_function = '-1' # default

        wall_f_types, wall_f_default = self.wall_function_types()
        if str(wall_function) not in wall_f_types:
            wall_function = '-1'
            self.setWallFunction(wall_function)

        return int(wall_function)


    @Variables.undoLocal
    def setWallFunction(self, wall_function):
        """
        Input wall function for advanced options.
        """
        k = int(wall_function)
        self.isIntInList(k, [-1, 0, 1, 2, 3, 4, 5, 7])
        if wall_function != -1:
            self.node_turb.xmlSetData('wall_function', wall_function)
        else:
            self.node_turb.xmlRemoveChild('wall_function')


    @Variables.noUndo
    def getTurbDiffModel(self):
        """
        Return turbulent diffusion model from advanced options.
        """
        turb_diff_model = self.node_turb.xmlGetString('turbulent_diffusion_model')
        if not turb_diff_model:
            turb_diff_model = self.defaultTurbulenceValues()['turbulent_diffusion_model']
        return turb_diff_model


    @Variables.undoLocal
    def setTurbDiffModel(self, turb_diff_model):
        """
        Input turbulent diffusion model for advanced options.
        """
        if turb_diff_model == self.defaultTurbulenceValues()['turbulent_diffusion_model']:
            self.node_turb.xmlRemoveChild('turbulent_diffusion_model')
        else:
            self.node_turb.xmlSetData('turbulent_diffusion_model', turb_diff_model)


    @Variables.noUndo
    def getRijCoupled(self):
        """
        Return Rij component coupling
        """
        status = 'on'
        node_c = self.node_turb.xmlGetNode('coupled_rij', 'status')
        if node_c:
            status = node_c['status']

        return status


    @Variables.undoLocal
    def setRijCoupled(self, status):
        """
        Set Rij component coupling
        """
        self.isOnOff(status)
        if status == 'off':
            node_c = self.node_turb.xmlInitNode('coupled_rij', 'status')
            node_c ['status'] = status
        else:
            self.node_turb.xmlRemoveChild('coupled_rij')


    @Variables.noUndo
    def getGravity(self):
        """
        Return scale model from advanced options .
        """
        node_gravity = self.node_turb.xmlInitNode('gravity_terms', 'status')
        gravity = node_gravity['status']
        if not gravity:
            gravity = self.defaultTurbulenceValues()['gravity_terms']
            self.setGravity(gravity)

        # force gravity force to off for Spalart-Allmaras model
        if self.getTurbulenceModel() == 'Spalart-Allmaras':
            gravity = 'off'
            self.setGravity(gravity)

        return gravity


    @Variables.undoLocal
    def setGravity(self, gravity):
        """
        Input gravity for advanced options.
        """
        self.isOnOff(gravity)
        node_gravity = self.node_turb.xmlInitNode('gravity_terms', 'status')
        node_gravity ['status'] = gravity


    def getTurbulenceVariable(self):
        """
        Return the turbulence <variable> markup list.
        """
        model = self.getTurbulenceModel()
        nodeList = []

        if model in ('k-epsilon','k-epsilon-PL'):
            nodeList.append(self.node_turb.xmlGetNode('variable', name='k'))
            nodeList.append(self.node_turb.xmlGetNode('variable', name='epsilon'))
        elif model in self.LESmodels():
            if self.node_lagr['model'] != "off":
                nodeList.append(self.node_turb.xmlGetNode('variable', name='k'))
                nodeList.append(self.node_turb.xmlGetNode('variable', name='epsilon'))
        elif model in ('Rij-epsilon', 'Rij-SSG', 'Rij-EBRSM'):
            for var in ('r11', 'r22', 'r33',
                        'r12', 'r13', 'r23', 'epsilon'):
                nodeList.append(self.node_turb.xmlGetNode('variable', name=var))
            if model == 'Rij-EBRSM':
                nodeList.append(self.node_turb.xmlGetNode('variable', name='alpha'))
        elif model == 'v2f-BL-v2/k':
            nodeList.append(self.node_turb.xmlGetNode('variable', name='k'))
            nodeList.append(self.node_turb.xmlGetNode('variable', name='epsilon'))
            nodeList.append(self.node_turb.xmlGetNode('variable', name='phi'))
            nodeList.append(self.node_turb.xmlGetNode('variable', name='alpha'))
        elif model == 'k-omega-SST':
            nodeList.append(self.node_turb.xmlGetNode('variable', name='k'))
            nodeList.append(self.node_turb.xmlGetNode('variable', name='omega'))
        elif model == 'Spalart-Allmaras':
            nodeList.append(self.node_turb.xmlGetNode('variable', name='nu_tilda'))
        return nodeList

#-------------------------------------------------------------------------------
# TurbulenceModel test case
#-------------------------------------------------------------------------------

class TurbulenceModelTestCase(ModelTest):
    """
    """
    def checkTurbulenceInstantiation(self):
        """Check whether the TurbulenceModel class could be instantiated"""
        model = None
        model = TurbulenceModel(self.case)
        assert model != None, 'Could not instantiate TurbulenceModel'

    def checkTurbulenceModelsList(self):
        """Check whether the TurbulenceModelList could be get"""
        from code_saturne.model.LagrangianModel import LagrangianModel
        LagrangianModel(self.case).setLagrangianStatus('on')
        del LagrangianModel
        mdl = TurbulenceModel(self.case)

        l = mdl.RANSmodels()
        l.insert(0, "off")
        assert mdl.turbulenceModelsList() == l, \
               'Could not return turbulence models for particles tracking'

        mdl.node_gas['model'] = 'on'
        assert mdl.turbulenceModelsList() == mdl.RANSmodels(), \
            'Could not return turbulence models for particular physics'

    def checkSetMixingLength(self):
        """Check whether the mixing length turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.setTurbulenceModel('mixing_length')
        mdl.setLengthScale(1)
        doc ='''<turbulence model="mixing_length">
                    <property label="TurbVisc" name="turbulent_viscosity"/>
                    <initialization choice="reference_velocity">
                        <reference_velocity>1</reference_velocity>
                    </initialization>
                    <mixing_length_scale>1</mixing_length_scale>
              </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set the mixing length turbulence model'

    def checkSetkepsilon(self):
        """Check whether the k-epsilon turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('k-epsilon')
        doc ='''<turbulence model="k-epsilon">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1</reference_velocity>
                </initialization>
               </turbulence>'''

        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set the k-epsilon turbulence model'

    def checkSetkepsilonPL(self):
        """Check whether the k-epsilon turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('k-epsilon-PL')
        doc ='''<turbulence model="k-epsilon-PL">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1</reference_velocity>
                </initialization>
              </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set the linear production k-epsilon turbulence model'

    def checkSetRijepsilon(self):
        """Check whether the Rij-epsilon turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.setTurbulenceModel('Rij-epsilon')
        doc ='''<turbulence model="Rij-epsilon">
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <variable label="R11" name="r11"/>
                <variable label="R22" name="r22"/>
                <variable label="R33" name="r33"/>
                <variable label="R12" name="r12"/>
                <variable label="R13" name="r13"/>
                <variable label="R23" name="r23"/>
                <variable label="Dissip" name="epsilon"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1</reference_velocity>
                </initialization>
            </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set the Rij-epsilon turbulence model'

    def checkSetRijepsilonSSG(self):
        """Check whether the Rij-epsilon SSG turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.setTurbulenceModel('Rij-SSG')
        truc = mdl.node_turb
        doc ='''<turbulence model="Rij-SSG">
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1</reference_velocity>
                </initialization>
                <variable label="R11" name="r11"/>
                <variable label="R22" name="r22"/>
                <variable label="R33" name="r33"/>
                <variable label="R12" name="r12"/>
                <variable label="R13" name="r13"/>
                <variable label="R23" name="r23"/>
                <variable label="Dissip" name="epsilon"/>
              </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the Rij-epsilon SSG turbulence model'

    def checkSetRijepsilonEBRSM(self):
        """Check whether the Rij-epsilon EBRSM turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.setTurbulenceModel('Rij-EBRSM')
        truc = mdl.node_turb
        doc ='''<turbulence model="Rij-EBRSM">
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1</reference_velocity>
                </initialization>
                <variable label="R11" name="r11"/>
                <variable label="R22" name="r22"/>
                <variable label="R33" name="r33"/>
                <variable label="R12" name="r12"/>
                <variable label="R13" name="r13"/>
                <variable label="R23" name="r23"/>
                <variable label="Dissip" name="epsilon"/>
                <variable label="alpha" name="alpha"/>
              </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the Rij-epsilon EBRSM turbulence model'

    def checkSetLESSmagorinsky(self):
        """Check whether the classical LES turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.node_turb.xmlRemoveChild('property')
        mdl.node_turb.xmlRemoveChild('initialization')
        mdl.setTurbulenceModel('LES_Smagorinsky')
        truc = mdl.node_turb
        doc ='''<turbulence model="LES_Smagorinsky">
                    <property label="Csdyn2" name="smagorinsky_constant^2"/>
               </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
             'Could not set the LES turbulence model'

    def checkSetLESdynamique(self):
        """Check whether the dynamique LES turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.node_turb.xmlRemoveChild('property')
        mdl.node_turb.xmlRemoveChild('initialization')
        mdl.setTurbulenceModel('LES_dynamique')
        truc = mdl.node_turb
        doc = '''<turbulence model="LES_dynamique">
                 <property label="Csdyn2" name="smagorinsky_constant^2"/>
               </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the dynamique LES turbulence model'

    def checkSetV2F(self):
        """Check whether the v2f phi turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('v2f-BL-v2/k')
        doc = '''<turbulence model="v2f-BL-v2/k">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <variable label="phi" name="phi"/>
                <variable label="alpha" name="alpha"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1.0</reference_velocity>
                </initialization>
              </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the v2f phi turbulence model'

    def checkkOmegaSST(self):
        """Check whether the k-Omega SST turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('k-omega-SST')
        doc = '''<turbulence model="k-omega-SST">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <variable label="omega" name="omega"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1.0</reference_velocity>
                </initialization>
            </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the k_Omega SST turbulence model'

    def checkSpalartAllmaras(self):
        """Check whether the Spalart-Allmaras turbulence model could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('Spalart-Allmaras')
        doc = '''<turbulence model="Spalart-Allmaras">
                <variable label="NuTilda" name="nu_tilda"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                  <reference_velocity>1.0</reference_velocity>
                </initialization>
            </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the Spalart-Allmaras turbulence model'

    def checkGetTurbulenceModel(self):
        """Check whether the turbulence model could be get"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('Rij-epsilon')

        assert mdl.getTurbulenceModel() == 'Rij-epsilon', \
            'Could not get the turbulence model'

    def checkSetLengthScale(self):
        """Check whether the mixing length scale could be set"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('mixing_length')
        mdl.node_turb.xmlRemoveChild('variable')
        mdl.node_turb.xmlRemoveChild('property')
        mdl.node_turb.xmlRemoveChild('initialization')
        mdl.setLengthScale(123.0)
        doc = '''<turbulence model="mixing_length">
                  <mixing_length_scale>123</mixing_length_scale>
                </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
           'Could not set the mixing length scale'

    def checkSetandGetWallFunction(self):
        """Check whether the wall function could be get"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('k-epsilon')
        mdl.setWallFunction(2)

        doc = '''<turbulence model="k-epsilon">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                    <reference_velocity>1</reference_velocity>
                </initialization>
                <wall_function>2</wall_function>
                </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set the wall function '
        assert mdl.getWallFunction() == 2,\
            'Could not get the wall function '

    def checkSetandGetgravity(self):
        """Check whether the gravity could be get"""
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('k-epsilon')
        mdl.setGravity('off')
        doc = '''<turbulence model="k-epsilon">
                <variable label="TurbEner" name="k"/>
                <variable label="Dissip" name="epsilon"/>
                <property label="TurbVisc" name="turbulent_viscosity"/>
                <initialization choice="reference_velocity">
                    <reference_velocity>1</reference_velocity>
                </initialization>
                    <gravity_terms status="off"/>
                </turbulence>'''
        assert mdl.node_turb == self.xmlNodeFromString(doc),\
            'Could not set gravity status '
        assert mdl.getGravity() == "off",\
            'Could not get gravity status '


def suite():
    testSuite = unittest.makeSuite(TurbulenceModelTestCase, "check")
    return testSuite

def runTest():
    print("TurbulenceModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
