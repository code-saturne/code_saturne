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

"""
This module defines the differents possible outputs : listings, for ensights
chronologic and historic files .... captors .... used variables ...

This module defines the following classes:
- NumericalParamEquationModel
- NumericalParamEquationModelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import *
import code_saturne.Base.Toolbox as Tool
from code_saturne.Base.XMLvariables import Variables, Model
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
from code_saturne.Pages.TurbulenceModel import TurbulenceModel
from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
from code_saturne.Pages.GroundwaterModel import GroundwaterModel

#-------------------------------------------------------------------------------
# NumericalParamEquat model class
#-------------------------------------------------------------------------------

class NumericalParamEquationModel(Model):
    """
    """
    def __init__(self, case):
        """
        initialization of nodes lists
        """
        self.case = case
        self.node_models  = self.case.xmlGetNode('thermophysical_models')
        self.node_vitpre  = self.node_models.xmlGetNode('velocity_pressure')
        self.node_varVP   = self.node_vitpre.xmlGetNodeList('variable')
        self.node_np      = self.case.xmlInitNode('numerical_parameters')
        self.node_anal    = self.case.xmlGetNode('analysis_control')
        self.model        = XMLmodel(self.case)
        self.darcy        = GroundwaterModel(self.case).getGroundwaterModel()

       #variables list
        self.var = []
        if self.darcy == "off":
            for node in self.node_varVP:
                self.var.append(node['name'])
        else:
            for node in self.node_varVP:
                if node['name'] != "velocity":
                    self.var.append(node['name'])

        for node in self._getThermalScalarNode():
            self.var.append(node['name'])
        for node in self._getAdditionalScalarNodes():
            self.var.append(node['name'])

        self.thermo = []
        for node in self._getThermalScalarNode():
            self.thermo.append(node['name'])

        self.UVW = []
        if self.darcy == "off":
            for node in self.node_varVP:
                if node['name'] == 'velocity':
                    self.UVW.append(node['name'])


    def _defaultValues(self, name=""):
        """ Private method: return default values """
        self.default = {}
        self.default['time_step_factor'] = 1.0
        self.default['verbosity'] = 0
        self.default['solver_precision'] = 1e-8
        self.default['solver_precision_pressure'] = 1e-8
        if NumericalParamGlobalModel(self.case).getTimeSchemeOrder() == 2:
            self.default['solver_precision'] = 1e-5
            self.default['solver_precision_pressure'] = 1e-5
        self.default['slope_test'] = 'on'
        self.default['flux_reconstruction'] = 'on'

        self.default['solver_choice'] = 'automatic'
        self.default['preconditioning_choice'] = 'automatic'

        self.default['order_scheme'] = 'automatic'

        if name not in self.var:
            self.default['blending_factor'] = 0.
        else:
            self.default['blending_factor'] = 1.

        from code_saturne.Pages.CompressibleModel import CompressibleModel
        if CompressibleModel(self.case).getCompressibleModel() != 'off':
            self.default['order_scheme'] = 'upwind'
            self.default['blending_factor'] = 0.
        del CompressibleModel

        if TurbulenceModel(self.case).getTurbulenceModel() in \
            ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            if name in self.UVW:
                self.default['slope_test'] = 'off'
            if name == 'pressure':
                self.default['rhs_reconstruction'] = 5
            else:
                self.default['rhs_reconstruction'] = 10
        else:
            if name == 'pressure':
                self.default['rhs_reconstruction'] = 2
            else:
                self.default['rhs_reconstruction'] = 1

        if name in self.thermo:
            for node in self._getThermalScalarNode():
                if node['name'] == name:
                    mdl = ThermalScalarModel(self.case).getThermalScalarModel()
                    if mdl == 'temperature_celsius':
                        self.default['min_value'] = -273.15
                        self.default['max_value'] = 1e+12
                    elif mdl =='temperature_kelvin':
                        self.default['min_value'] = 0
                        self.default['max_value'] = 1e+12
                    elif mdl =='enthalpy':
                        self.default['min_value'] = 0
                        self.default['max_value'] = 1e+12
                    elif mdl =='potential_temperature':
                        self.default['min_value'] = 0
                        self.default['max_value'] = 1e+12
                    elif mdl =='liquid_potential_temperature':
                        self.default['min_value'] = 0
                        self.default['max_value'] = 1e+12
                    elif mdl =='total_energy':
                        self.default['min_value'] = 0
                        self.default['max_value'] = 1e+12
        else:
            self.default['min_value'] = -1e+12
            self.default['max_value'] = 1e+12

        return self.default


    def _getThermalScalarNode(self):
        """ Private method: return node of thermal scalar """
        node_models = self.case.xmlGetNode('thermophysical_models')
        node_scalar = node_models.xmlGetNode('thermal_scalar')
        thermal_scalar_list = node_scalar.xmlGetNodeList('variable', type='thermal')
        return thermal_scalar_list


    def _getPuCoalScalarsNodes(self):
        """ Private method: return list of pulverized coal scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('solid_fuels', 'model')
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
        return nodList


    def _getGasScalarsNodes(self):
        """ Private method: return list of gas combustion scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('gas_combustion', 'model')
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
        return nodList


    def _getMeteoScalarsNodes(self):
        """ Private method: return list of meteo scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('atmospheric_flows', 'model')
        if not node: return []
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
        return nodList


    def _getElectricalScalarsNodes(self):
        """ Private method: return list of electric scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('joule_effect', 'model')
        if not node: return []
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
        return nodList


    def _getCompressibleScalarsNodes(self):
        """ Private method: return list of compressible scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('compressible_model', 'model')
        if not node: return []
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('variable')
        return nodList


    def _getAdditionalScalarNodes(self):
        """ Private method: return list of additional scalar's nodes """
        n = self.case.xmlGetNode('additional_scalars')
        return n.xmlGetNodeList('variable', type='user')


    def _getAleVariablesNodes(self):
        """ Private method: return list of nodes for ALE"""
        nodList = []
        n = self.node_models.xmlGetNode('ale_method')
        if n['status'] == 'on':
            nodList =  n.xmlGetNodeList('variable')
        return nodList


    def _getClippingNodesList(self):
        """ Return list of nodes for class view Scheme"""
        self.var_clip = []

        for part in (self._getThermalScalarNode(),
                     self._getPuCoalScalarsNodes(),
                     self._getGasScalarsNodes(),
                     self._getMeteoScalarsNodes(),
                     self._getElectricalScalarsNodes(),
                     self._getCompressibleScalarsNodes(),
                     self._getAdditionalScalarNodes()):
            self.var_clip.append(part)
        return self.var_clip


    def _getSchemeNodesList(self):
        """ Return list of nodes for class view Scheme"""
        self.var_shem = []

        for part in (self.node_varVP,
                     self.model.getTurbVariable(),
                     self._getThermalScalarNode(),
                     self._getPuCoalScalarsNodes(),
                     self._getGasScalarsNodes(),
                     self._getMeteoScalarsNodes(),
                     self._getElectricalScalarsNodes(),
                     self._getCompressibleScalarsNodes(),
                     self._getAdditionalScalarNodes()):
            self.var_shem.append(part)
        return self.var_shem


    def _getSolverNodesList(self):
        """ Return list of nodes for class view Solver"""
        self.var_solv = []
        for part in (self.node_varVP,
                     self.model.getTurbVariable(),
                     self._getThermalScalarNode(),
                     self._getPuCoalScalarsNodes(),
                     self._getGasScalarsNodes(),
                     self._getMeteoScalarsNodes(),
                     self._getElectricalScalarsNodes(),
                     self._getAdditionalScalarNodes(),
                     self._getCompressibleScalarsNodes(),
                     self._getAleVariablesNodes()):
            self.var_solv.append(part)
        return self.var_solv


    def _getClippingNameNode(self, name):
        """ Private method: return node called with name'name' for solver scheme"""
        for node in self._getClippingNodesList():
            for n in node:
                if n['name'] == name:
                    return n
        raise ValueError("This name does not exist: " + name)


    def _getSchemeNameNode(self, name):
        """ Private method: return node called with name'name' for scheme nodes"""
        for node in self._getSchemeNodesList():
            for n in node:
                if n['name'] == name:
                    return n
        raise ValueError("This name does not exist: " + name)


    def _getSolverNameNode(self, name):
        """ Private method: return node called with name'name' for solver scheme"""
        for node in self._getSolverNodesList():
            for n in node:
                if n['name'] == name:
                    return n
        raise ValueError("This name does not exist: " + name)


    def _isPressure(self, node):
        """ Return : 1 if name of node is 'pressure', 0 if not """
        if node and node['name'] == 'pressure':
            return 1
        else:
            return 0

    @Variables.undoGlobal
    def setSchemeDefaultValues(self):
        """Usefull for TurbulenceModel in case of LES"""
        for name in self.var:
            try:
                self.setBlendingFactor(name, self._defaultValues(name)['blending_factor'])
                self.setScheme(name, self._defaultValues(name)['order_scheme'])
                self.setSlopeTest(name, self._defaultValues(name)['slope_test'])
                self.setFluxReconstruction(name, self._defaultValues(name)['flux_reconstruction'])
                self.setRhsReconstruction(name, self._defaultValues(name)['rhs_reconstruction'])
            except:
                pass


    @Variables.noUndo
    def getClippingList(self):
        """ Return the variables name list for clipping parameters """
        lst = []
        for node in self._getClippingNodesList():
            for n in node:
                if n['type'] != 'model':
                    lst.append(n['name'])
        return lst


    @Variables.noUndo
    def getSchemeList(self):
        """ Return the variables name list for scheme parameters """
        lst = []
        if self.darcy == "off":
            for node in self._getSchemeNodesList():
                for n in node:
                    lst.append(n['name'])
        else:
            for node in self._getSchemeNodesList():
                for n in node:
                    if n['name'] != "velocity":
                        lst.append(n['name'])
        return lst


    @Variables.noUndo
    def getSolverList(self):
        """ Return the variables name list for solver parameters """
        lst = []
        from code_saturne.Pages.CompressibleModel import CompressibleModel
        comp_model = CompressibleModel(self.case).getCompressibleModel()
        del CompressibleModel

        if self.darcy == "off":
            for node in self._getSolverNodesList():
                for n in node:
                    lst.append(n['name'])
        else:
            for node in self._getSolverNodesList():
                for n in node:
                    if n['name'] != "velocity":
                        lst.append(n['name'])

        return lst


    def isScalar(self, name):
        """
        Return : 1 if type of node is 'user' or 'thermal' or 'model',
                 0 if not.  Only used by the view by solver class
        """
        node = self._getSolverNameNode(name)
        if node:
            if node['type'] in ['user', 'thermal', 'model']:
                return 1
            else:
                return 0


# Following methods for dependances of scheme:

    @Variables.noUndo
    def getScheme(self, name):
        """ Return value of order scheme for variable labelled name """
        node = self._getSchemeNameNode(name)
        if self._isPressure(node):
            return None
        value = self._defaultValues(name)['order_scheme']
        n = node.xmlGetNode('order_scheme')
        if n:
            value = n['choice']
        return value


    @Variables.noUndo
    def getBlendingFactor(self, name):
        """ Return value of blending factor for variable labelled name """
        node = self._getSchemeNameNode(name)
        if self._isPressure(node):
            return None
        value = node.xmlGetDouble('blending_factor')
        if value == None:
            value = self._defaultValues(name)['blending_factor']
        return value


    @Variables.noUndo
    def getSlopeTest(self, name):
        """ Return value of slope test for variable labelled name """
        node = self._getSchemeNameNode(name)
        if self._isPressure(node):
            return None
        value = self._defaultValues(name)['slope_test']
        n = node.xmlGetNode('slope_test')
        if n:
            value = n['status']
        return value


    @Variables.noUndo
    def getFluxReconstruction(self, name):
        """ Return value of flux reconstruction for variable labelled name """
        node = self._getSchemeNameNode(name)
        value = self._defaultValues()['flux_reconstruction']
        if node.xmlGetNode('flux_reconstruction'):
            value = node.xmlGetNode('flux_reconstruction')['status']
        return value


    @Variables.noUndo
    def getRhsReconstruction(self, name):
        """ Return value of RHS reconstruction for variable labelled name """
        node = self._getSchemeNameNode(name)
        value = node.xmlGetDouble('rhs_reconstruction')
        if value == None:
            value = self._defaultValues(name)['rhs_reconstruction']
        return value


    @Variables.undoGlobal
    def setBlendingFactor(self, name, value):
        """
        Put value of blending factor for variable labelled name
        only if it 's different of default value
        """
        node = self._getSchemeNameNode(name)
        if self._isPressure(node):
            return
        self.isGreaterOrEqual(value, 0.)
        self.isLowerOrEqual(value, 1.)
        scheme = self.getScheme(name)
        if scheme == self._defaultValues(name)['order_scheme']:
            if scheme == 'upwind':
                node.xmlRemoveChild('blending_factor')
            else:
                node.xmlSetData('blending_factor', value)
        else:
            node.xmlSetData('blending_factor', value)
#            if value != self._defaultValues(name)['blending_factor']:
#                node.xmlSetData('blending_factor', value)
#            else:
#                node.xmlRemoveChild('blending_factor')


    @Variables.undoGlobal
    def setScheme(self, name, value):
        """
        Put value of order scheme for variable or scalar labelled name
        only if it 's different of default value
        """
        self.isInList(value, ('automatic', 'upwind', 'centered', 'solu'))
        node = self._getSchemeNameNode(name)
        if value == self._defaultValues(name)['order_scheme']:
            node.xmlRemoveChild('order_scheme')
            if self.getBlendingFactor(name) == self._defaultValues(name)['blending_factor']\
                or value == 'centered' and self.getBlendingFactor(name) == 0. \
                or value == 'upwind' and self.getBlendingFactor(name) != 0.:
                    node.xmlRemoveChild('blending_factor')
        else:
            n = node.xmlInitNode('order_scheme')
            n['choice'] = value


    @Variables.undoLocal
    def setSlopeTest(self, name, status):
        """ Put status of slope test for variable labelled name """
        self.isOnOff(status)
        node = self._getSchemeNameNode(name)
        if status == self._defaultValues(name)['slope_test']:
            node.xmlRemoveChild('slope_test')
        else:
            n = node.xmlInitNode('slope_test')
            n['status'] = status


    @Variables.undoLocal
    def setFluxReconstruction(self, name, value):
        """ Put status of flux reconstruction for variable labelled name """
        self.isOnOff(value)
        node = self._getSchemeNameNode(name)
        if value == self._defaultValues()['flux_reconstruction']:
            node.xmlRemoveChild('flux_reconstruction')
        else:
            n = node.xmlInitNode('flux_reconstruction')
            n['status']=value


    @Variables.undoLocal
    def setRhsReconstruction(self, name, value):
        """
        Put value of RHS reconstruction
        """
        self.isInt(value)
        node = self._getSchemeNameNode(name)
        node.xmlSetData('rhs_reconstruction', value)


# Following methods for dependances of solver:
    @Variables.undoLocal
    def setSolverPrecision(self, name, value):
        """ Put value of solver precision for variable labelled name """
        # for pressure default value always equal to 1e-8
        self.isPositiveFloat(value)
        node = self._getSolverNameNode(name)
        if self._isPressure(node):
            default = self._defaultValues()['solver_precision_pressure']
        else:
            default = self._defaultValues()['solver_precision']

        if value != default:
            node.xmlSetData('solver_precision', value)
        else:
            node.xmlRemoveChild('solver_precision')


    @Variables.undoLocal
    def setSolverChoice(self, name, value):
        """ Put choice of solver for variable labelled name """
        self.isInList(value, ('multigrid', 'multigrid_k_cycle',
                              'conjugate_gradient',
                              'flexible_conjugate_gradient',
                              'inexact_conjugate_gradient', 'jacobi',
                              'bi_cgstab', 'bi_cgstab2', 'gmres', 'automatic',
                              'gauss_seidel', 'symmetric_gauss_seidel', 'PCR3'))
        node = self._getSolverNameNode(name)

        default = self._defaultValues()['solver_choice']

        if value != default:
            n = node.xmlInitNode('solver_choice')
            n['choice'] = value
        else:
            node.xmlRemoveChild('solver_choice')

    @Variables.undoLocal
    def setPreconditioningChoice(self, name, value):
        """ Put choice of preconditioning for variable labelled name """
        self.isInList(value, ('multigrid', 'multigrid_k_cycle',
                              'none', 'jacobi',
                              'polynomial', 'automatic'))
        node = self._getSolverNameNode(name)

        default = self._defaultValues()['preconditioning_choice']

        if value != default:
            n = node.xmlInitNode('preconditioning_choice')
            n['choice'] = value
        else:
            node.xmlRemoveChild('preconditioning_choice')


    @Variables.noUndo
    def getSolverPrecision(self, name):
        """ Return value of solver precision for variable labelled name """
        node = self._getSolverNameNode(name)

        if self._isPressure(node):
            default = self._defaultValues()['solver_precision_pressure']
        else:
            default = self._defaultValues()['solver_precision']

        value = node.xmlGetDouble('solver_precision')
        if value == None:
            value = default
        return value


    @Variables.noUndo
    def getSolverAllowMultigrid(self, name):
        """ Return value of variable dimension for variable labelled name """
        node = self._getSolverNameNode(name)
        dim = node['dimension']
        if not dim:
            dim = 1
        if int(dim) <= 1:
            return True
        return False


    @Variables.noUndo
    def getSolverChoice(self, name):
        """ Return choice of solver for variable labelled name """
        node = self._getSolverNameNode(name)
        n = node.xmlGetNode('solver_choice')

        if n:
            value = n['choice']
        else:
            value = self._defaultValues()['solver_choice']
        return value


    @Variables.noUndo
    def getPreconditioningChoice(self, name):
        """ Return choice of preconditioning for variable labelled name """
        node = self._getSolverNameNode(name)
        n = node.xmlGetNode('preconditioning_choice')

        if n:
            value = n['choice']
        else:
            value = self._defaultValues()['preconditioning_choice']
        return value


    @Variables.noUndo
    def getScalarTimeStepFactor(self, name):
        """ Return value of time_step_factor for variable labelled name """
        if self.isScalar(name):
            node = self._getSolverNameNode(name)
            value = node.xmlGetDouble('time_step_factor')
            if value == None:
                value = self._defaultValues()['time_step_factor']
            return value
        else:
            raise ValueError("This method runs only with scalar name")


    @Variables.undoLocal
    def setScalarTimeStepFactor(self, name, value):
        """ Put value of time_step_factor for variable labelled name """
        self.isStrictPositiveFloat(value)
        if self.isScalar(name):
            node = self._getSolverNameNode(name)
            if value != self._defaultValues()['time_step_factor']:
                node.xmlSetData('time_step_factor', value)
            else:
                node.xmlRemoveChild('time_step_factor')
        else:
            raise ValueError("This method runs only with scalar name")


    @Variables.noUndo
    def getVerbosity(self, name):
        """ Return value of verbosity for variable labelled name """
        node = self._getSolverNameNode(name)
        value = node.xmlGetInt('verbosity')
        if value == None:
            value = self._defaultValues()['verbosity']
        return value


    @Variables.undoLocal
    def setVerbosity(self, name, value):
        """ Put value of verbosity for variable labelled name """
        self.isInt(value)
        node = self._getSolverNameNode(name)
        if value != self._defaultValues()['verbosity']:
            node.xmlSetData('verbosity', value)
        else:
            node.xmlRemoveChild('verbosity')


    @Variables.noUndo
    def getMinValue(self, name):
        """Get minimal value from an additional_scalar with label scalar_name"""
        self.isInList(name, self.getClippingList())
        node = self._getClippingNameNode(name)
        min_val = node.xmlGetChildDouble('min_value')
        if min_val == None:
            min_val = self._defaultValues(name)['min_value']
            self.setMinValue(name, min_val)

        return min_val


    @Variables.undoLocal
    def setMinValue(self, name, min_value):
        """
        Put minimal value for an additional_scalar with label scalar_name.
        Method also used by ThermalScalarModel
        """
        self.isFloat(min_value)
        self.isInList(name, self.getClippingList())
        node = self._getClippingNameNode(name)
        node.xmlSetData('min_value', min_value)


    @Variables.noUndo
    def getMaxValue(self, name):
        """Get maximal value from an additional_scalar with label scalar_name"""
        self.isInList(name, self.getClippingList())
        node = self._getClippingNameNode(name)
        max_val = node.xmlGetDouble('max_value')
        if max_val == None:
            max_val = self._defaultValues(name)['max_value']
            self.setMaxValue(name, max_val)
        return max_val


    @Variables.undoLocal
    def setMaxValue(self, name, max_value):
        """
        Put maximal value for an additional_scalar with label scalar_name.
        Method also used by ThermalScalarModel
        """
        self.isFloat(max_value)
        self.isInList(name, self.getClippingList())
        node = self._getClippingNameNode(name)
        node.xmlSetData('max_value', max_value)


#-------------------------------------------------------------------------------
# NumericalParamEquat test case
#-------------------------------------------------------------------------------


class NumericalParamEquatTestCase(ModelTest):
    """
    """
    def checkNumericalParamEquatInstantiation(self):
        """
        Check whether the NumericalParamEquationModel class could be instantiated
        """
        model = None
        model = NumericalParamEquationModel(self.case)
        assert model != None, 'Could not instantiate NumericalParamEquationModel'

    def checkSetAndGetScheme(self):
        """
        Check whether the NumericalParamEquationModel class could set and get scheme
        """
        model = NumericalParamEquationModel(self.case)
        model.setScheme('Velocity', 'upwind')
        doc = """<velocity_pressure>
                        <variable label="Pressure" name="pressure"/>
                        <variable label="Velocity" name="velocity"/>
                                <order_scheme choice="upwind"/>
                        </variable>
                        <property label="total_pressure" name="total_pressure"/>
                        <property label="Yplus" name="yplus" support="boundary"/>
                        <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set scheme in NumericalParamEquationModel'
        assert model.getScheme('VelocitW') == 'upwind',\
                'Could not get scheme in NumericalParamEquationModel'

    def checkSetAndGetBlendingFactor(self):
        """
        Check whether the NumericalParamEquationModel class could set and get blending factor
        """
        model = NumericalParamEquationModel(self.case)
        model.setScheme('Velocity', 'centered')
        model.setBlendingFactor('Velocity', 0.5)
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity"/>
                            <blending_factor>0.5</blending_factor>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set blending factor in NumericalParamEquationModel'
        assert model.getBlendingFactor('VelocitW') == 0.5,\
                'Could not get blending factor in NumericalParamEquationModel'

    def checkSetAndGetSlopeTest(self):
        """
        Check whether the NumericalParamEquationModel class could set and get slope test
        """
        model = NumericalParamEquationModel(self.case)
        model.setSlopeTest('Velocity', 'off')
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity"/>
                            <slope_test status="off"/>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set status of slope test in NumericalParamEquationModel'
        assert model.getSlopeTest('VelocitW') == 'off',\
                'Could not get status of slope test in NumericalParamEquationModel'


    def checkSetAndGetFluxReconstruction(self):
        """
        Check whether the NumericalParamEquationModel class could set and get flux reconstruction
        """
        model = NumericalParamEquationModel(self.case)
        model.setFluxReconstruction('Velocity', 'on')
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity"/>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set status of flux reconstruction in NumericalParamEquationModel'
        assert model.getFluxReconstruction('VelocitW') == 'on',\
                'Could not get status of flux reconstruction in NumericalParamEquationModel'

        model.setFluxReconstruction('Velocity', 'off')
        doc2 = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity"/>
                        <flux_reconstruction status="off"/>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                  </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc2),\
                'Could not set status of flux reconstruction in NumericalParamEquationModel'
        assert model.getFluxReconstruction('VelocitW') == 'off',\
                'Could not get status of flux reconstruction in NumericalParamEquationModel'


    def checkSetAndGetSolverPrecision(self):
        """
        Check whether the NumericalParamEquationModel class could set and get solver precision
        """
        model = NumericalParamEquationModel(self.case)

        assert model.getSolverPrecision('pressure') == 1e-8,\
                'Could not get solver precision for pressure in NumericalParamEquationModel'
        from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
        NumericalParamGlobalModel(self.case).setTimeSchemeOrder(2)
        del NumericalParamGlobalModel
        assert model.getSolverPrecision('velocity') == 1e-5

        model.setSolverPrecision('velocity', 2e-6)
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity">
                            <solver_precision>2e-06</solver_precision>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set solver precision in NumericalParamEquationModel'
        assert model.getSolverPrecision('velocity') == 2e-6,\
                'Could not get solver precision in NumericalParamEquationModel'

    def checkSetAndGetScalarTimeStepFactor(self):
        """
        Check whether the NumericalParamEquationModel class could set and get time step factor
        """
        model = NumericalParamEquationModel(self.case)
        from code_saturne.Pages.ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

##        self.failUnlessRaises(ValueError, model.setScalarTimeStepFactor('VelocitU', 25.), \
##           'Could not set time step factor in NumericalParamEquationModel')

        model.setScalarTimeStepFactor('TempC', 52.)
        node_sca = self.case.xmlGetNode('additional_scalars')
        vit = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="Velocity" name="velocity"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                 </velocity_pressure>"""
        sca = """<additional_scalars>
                    <variable label="TempC" name="temperature_celsius" type="thermal">
                            <initial_value zone_id="1">20.0</initial_value>
                            <min_value>-1e+12 </min_value>
                            <max_value>1e+12</max_value>
                            <time_step_factor>52</time_step_factor>
                    </variable>
                 </additional_scalars>"""
        assert model.node_vitpre == self.xmlNodeFromString(vit),\
                'Could not set time step factor in NumericalParamEquationModel'

        assert node_sca == self.xmlNodeFromString(sca),\
                'Could not set time step factor for scalar in NumericalParamEquationModel'

##        self.failUnlessRaises(ValueError, model.getScalarTimeStepFactor('VelocitV'), \
##           'Could not get time step factor in NumericalParamEquationModel')

        assert model.getScalarTimeStepFactor('TempC') == 52.,\
                'Could not get time step factor for scalar in NumericalParamEquationModel'

def suite():
    testSuite = unittest.makeSuite(NumericalParamEquatTestCase, "check")
    return testSuite

def runTest():
    print("NumericalParamEquatTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
