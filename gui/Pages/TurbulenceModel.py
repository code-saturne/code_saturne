# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel

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
        self.node_coal   = self.node_models.xmlGetChildNode('solid_fuels',    'model')
        self.node_joule  = self.node_models.xmlGetChildNode('joule_effect',   'model')
        self.node_gas    = self.node_models.xmlGetChildNode('gas_combustion', 'model')
        self.node_turb   = self.node_models.xmlInitChildNode('turbulence',    'model')
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
        default['turbulence_model'] = "k-epsilon-PL"
        default['length_scale']     = 1.0
        default['wall_function']      = 3
        default['gravity_terms']    = "on"

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

        elif model_turb in ('Rij-epsilon', 'Rij-SSG'):
            lst = ('rij', 'epsilon')
            for v in lst:
                # Rij is now considered as a tensor (vector of length 6, since it is symmetric)
                if v == 'rij':
                    self.setNewVariable(self.node_turb, v, label='Rij', dim='6')
                else:
                    self.setNewVariable(self.node_turb, v, label=v)
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')

        elif model_turb == 'Rij-EBRSM':
            lst = ('rij', 'epsilon', 'alpha')
            for v in lst:
                self.setNewVariable(self.node_turb, v, label=v)
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.__updateInletsForTurbulence()
            self.__removeVariablesAndProperties(lst, 'smagorinsky_constant^2')
            wall_function = 0
            self.setWallFunction(wall_function)

        elif model_turb in self.LESmodels():
            if model_turb == 'LES_dynamique':
                self.setNewProperty(self.node_turb, 'smagorinsky_constant^2')
            else:
                self.__removeVariablesAndProperties([], 'smagorinsky_constant^2')
            self.setNewProperty(self.node_turb, 'turbulent_viscosity')
            self.node_turb.xmlRemoveChild('wall_function')

            from code_saturne.Pages.TimeStepModel import TimeStepModel
            TimeStepModel(self.case).setTimePassing(0)
            del TimeStepModel

            NumericalParamGlobalModel(self.case).setTimeSchemeOrder(2)

            from code_saturne.Pages.NumericalParamEquationModel import NumericalParamEquationModel
            NumericalParamEquationModel(self.case).setSchemeDefaultValues()
            del NumericalParamEquationModel

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

        DefineUserScalarsModel(self.case).setTurbulentFluxGlobalModel(model_turb)


    def __updateInletsForTurbulence(self):
        """
        Put boundaries conditions if it's necessary
        """
        from code_saturne.Pages.Boundary import Boundary
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
        if l_scale == None:
            l_scale = self.defaultTurbulenceValues()['length_scale']
            self.setLengthScale(l_scale)
        return l_scale


    @Variables.noUndo
    def getWallFunction(self):
        """
        Return wall function from advanced options.
        """
        wall_function = self.node_turb.xmlGetInt('wall_function')
        if not wall_function:
            wall_function = -1 # for next test
        if wall_function < 0 or wall_function > 5:
            wall_function = self.defaultTurbulenceValues()['wall_function']
            self.setWallFunction(wall_function)
        return wall_function


    @Variables.undoLocal
    def setWallFunction(self, wall_function):
        """
        Input wall function for advanced options.
        """
        self.isIntInList(wall_function, [0, 1, 2, 3, 4, 5])
        self.node_turb.xmlSetData('wall_function', wall_function)


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
        from code_saturne.Pages.LagrangianModel import LagrangianModel
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
