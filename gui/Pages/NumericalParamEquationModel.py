# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module defines the differents possible outputs : listings, for ensights
chronologic and historic files .... captors .... used variables ...

This module defines the following classes:
- NumericalParamEquatModel
- NumericalParamEquatModelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string
import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Variables, Model
from Base.XMLmodel import XMLmodel, ModelTest
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
from Pages.TurbulenceModel import TurbulenceModel

#-------------------------------------------------------------------------------
# NumericalParamEquat model class
#-------------------------------------------------------------------------------

class NumericalParamEquatModel(Model):
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
        self.node_anal      = self.case.xmlGetNode('analysis_control')
        self.model = XMLmodel(self.case)

       #variables list
        self.var = []
        for node in self.node_varVP:
            self.var.append(node['label'])
        for node in self._getThermalScalarNode():
            self.var.append(node['label'])
        for node in self._getAdditionalScalarNodes():
            self.var.append(node['label'])

        self.UVW = []
        for node in self.node_varVP:
            if node['name'] in ('velocity_U', 'velocity_V', 'velocity_W'):
                self.UVW.append(node['label'])


    def _defaultValues(self, label=""):
        """ Private method: return default values """
        self.default = {}
        self.default['time_step_factor'] = 1.0
        self.default['max_iter_number'] = 10000
        self.default['solveur_precision'] = 1e-8
        self.default['solveur_precision_pressure'] = 1e-8
        if NumericalParamGlobalModel(self.case).getTimeSchemeOrder() == 2:
            self.default['solveur_precision'] = 1e-5
        self.default['slope_test'] = 'on'
        self.default['flux_reconstruction'] = 'on'

        if label not in self.var:
            self.default['order_scheme'] = 'upwind'
            self.default['blending_factor'] = 0.
        else:
            self.default['order_scheme'] = 'centered'
            self.default['blending_factor'] = 1.

        if TurbulenceModel(self.case).getTurbulenceModel() in \
            ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            if label in self.UVW:
                self.default['slope_test'] = 'off'

        return self.default


    def _getThermalScalarNode(self):
        """ Private method: return node of thermal scalar """
        node_scalar = self.case.xmlGetNode('additional_scalars')
        thermal_scalar_list = node_scalar.xmlGetNodeList('scalar', type='thermal')
        return thermal_scalar_list


    def _getPuCoalScalarsNodes(self):
        """ Private method: return list of pulverized coal scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('pulverized_coal', 'model')
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('scalar')
        return nodList


    def _getMeteoScalarsNodes(self):
        """ Private method: return list of pulverized coal scalar's nodes """
        nodList = []
        node = self.node_models.xmlGetNode('atmospheric_flows', 'model')
	if not node: return []
        model = node['model']
        if model != 'off':
            nodList = node.xmlGetNodeList('scalar')
        return nodList


    def _getAdditionalScalarNodes(self):
        """ Private method: return list of additional scalar's nodes """
        n = self.case.xmlGetNode('additional_scalars')
        return n.xmlGetNodeList('scalar', type='user')


    def _getAleVariablesNodes(self):
        """ Private method: return list of nodes for ALE"""
        nodList = []
        n = self.node_models.xmlGetNode('ale_method')
        if n['status'] == 'on':
            nodList =  n.xmlGetNodeList('variable')
        return nodList


    def _getSchemeNodesList(self):
        """ Return list of nodes for class view Scheme"""
        self.var_shem = []

        for part in (self.node_varVP,
                     self.model.getTurbVariable(),
                     self._getThermalScalarNode(),
                     self._getPuCoalScalarsNodes(),
                     self._getMeteoScalarsNodes(),
                     self._getAdditionalScalarNodes()):
            self.var_shem.append(part)
        return self.var_shem


    def _getSolveurNodesList(self):
        """ Return list of nodes for class view Solveur"""
        self.var_solv = []
        for part in (self.node_varVP,
                     self.model.getTurbVariable(),
                     self._getThermalScalarNode(),
                     self._getPuCoalScalarsNodes(),
                     self._getMeteoScalarsNodes(),
                     self._getAdditionalScalarNodes(),
                     self._getAleVariablesNodes()):
            self.var_solv.append(part)
        return self.var_solv


    def _getSchemeLabelNode(self, label):
        """ Private method: return node called with label'label' for scheme nodes"""
        for node in self._getSchemeNodesList():
            for n in node:
                if n['label'] == label:
                    if n['name'] == 'pressure':
                        raise ValueError, "This method does not run with pressure"
                    else:
                        return n
        raise ValueError, "This label does not exist: " + label


    def _getSolveurLabelNode(self, label):
        """ Private method: return node called with label'label' for solveur scheme"""
        for node in self._getSolveurNodesList():
            for n in node:
                if n['label'] == label:
                    return n
        raise ValueError, "This label does not exist: " + label


    def _isPressure(self, node):
        """ Return : 1 if name of node is 'pressure', 0 if not """
        if node and node['name'] == 'pressure': 
            return 1
        else: 
            return 0


    def setSchemeDefaultValues(self):
        """Usefull for TurbulenceModel in case of LES"""
        for label in self.var:
            try:
                self.setBlendingFactor(label, self._defaultValues(label)['blending_factor'])
                self.setScheme(label, self._defaultValues(label)['order_scheme'])
                self.setSlopeTest(label, self._defaultValues(label)['slope_test'])
                self.setFluxReconstruction(label, self._defaultValues(label)['flux_reconstruction'])
            except:
                pass


    def getSchemeList(self):
        """ Return the variables label list for scheme parameters """
        list = []
        for node in self._getSchemeNodesList():
            for n in node:
                if not self._isPressure(n):
                    list.append(n['label'])
        return list


    def getSolveurList(self):
        """ Return the variables label list for solveur parameters """
        list = []
        for node in self._getSolveurNodesList():
            for n in node:
                list.append(n['label'])
        return list


    def isScalar(self, label):
        """
        Return : 1 if type of node is 'user' or 'thermal' or 'model', 
                 0 if not.  Only used by the view by solveur class
        """
        node = self._getSolveurLabelNode(label)
        if node:
            if node['type'] in ['user', 'thermal', 'model']:
                return 1
            else: 
                return 0


# Following methods for dependances of scheme:

    def getScheme(self, label):
        """ Return value of order scheme for variable labelled label """
        node = self._getSchemeLabelNode(label)
##        if self._isPressure(node): 
### FIXME: return nothing ?
##            return
##        else:
        value = self._defaultValues(label)['order_scheme']
        n = node.xmlGetNode('order_scheme')
        if n: 
            value = n['choice']
        return value


    def getBlendingFactor(self, label):
        """ Return value of blending factor for variable labelled label """
        node = self._getSchemeLabelNode(label)
##        if self._isPressure(node): 
##            return
##        else:
        value = node.xmlGetDouble('blending_factor')
        if value == None:
            value = self._defaultValues(label)['blending_factor']
        return value


    def getSlopeTest(self, label):
        """ Return value of slope test for variable labelled label """
        node = self._getSchemeLabelNode(label)
##        if self._isPressure(node): 
##            return
##        else:
        value = self._defaultValues(label)['slope_test']
        n = node.xmlGetNode('slope_test')
        if n: 
            value = n['status']
        return value


    def getFluxReconstruction(self, label):
        """ Return value of flux reconstruction for variable labelled label """
        node = self._getSchemeLabelNode(label) 
##        if self._isPressure(node): 
##            return
##        else:
        value = self._defaultValues()['flux_reconstruction']
        if node.xmlGetNode('flux_reconstruction'):
            value = node.xmlGetNode('flux_reconstruction')['status']
        return value


    def setBlendingFactor(self, label, value):
        """
        Put value of blending factor for variable labelled label 
        only if it 's different of default value
        """
        self.isGreaterOrEqual(value, 0.)
        self.isLowerOrEqual(value, 1.)
        node = self._getSchemeLabelNode(label)
        scheme = self.getScheme(label)
##        if scheme == 'upwind': 
##                node.xmlSetData('blending_factor', 0.)
##        else:
##            if value != self._defaultValues(label)['blending_factor']:
##                node.xmlSetData('blending_factor', value)
        if scheme == self._defaultValues(label)['order_scheme']:
            if scheme == 'upwind':
                node.xmlRemoveChild('blending_factor')
            else:
                node.xmlSetData('blending_factor', value)
        else:
            if value != self._defaultValues(label)['blending_factor']:
                node.xmlSetData('blending_factor', value)


    def setScheme(self, label, value):
        """ 
        Put value of order scheme for variable or scalar labelled label 
        only if it 's different of default value
        """
        self.isInList(value, ('upwind', 'centered', 'solu'))
        node = self._getSchemeLabelNode(label)
        if value == self._defaultValues(label)['order_scheme']:
            node.xmlRemoveChild('order_scheme')
            if self.getBlendingFactor(label) == self._defaultValues(label)['blending_factor']\
                or value == 'centered' and self.getBlendingFactor(label) == 0. \
                or value == 'upwind' and self.getBlendingFactor(label) != 0.:
                    node.xmlRemoveChild('blending_factor')
        else:
            n = node.xmlInitNode('order_scheme')
            n['choice'] = value


    def setSlopeTest(self, label, status):
        """ Put status of slope test for variable labelled label """
        self.isOnOff(status)
        node = self._getSchemeLabelNode(label)
        if status == self._defaultValues(label)['slope_test']:
            node.xmlRemoveChild('slope_test')
        else:
            n = node.xmlInitNode('slope_test')
            n['status'] = status


    def setFluxReconstruction(self, label, value):
        """ Put status of flux reconstruction for variable labelled label """
        self.isOnOff(value)
        node = self._getSchemeLabelNode(label)
        if value == self._defaultValues()['flux_reconstruction']:
            node.xmlRemoveChild('flux_reconstruction')
        else:
            n = node.xmlInitNode('flux_reconstruction')
            n['status']=value


# Following methods for dependances of solveur:

    def setMaxIterNumber(self, label, value):
        """ Put number of maximum iterations for variable labelled label """
        self.isInt(value)
        node = self._getSolveurLabelNode(label)
        if value != self._defaultValues()['max_iter_number']:
            node.xmlSetData('max_iter_number', value)
        else:
            node.xmlRemoveChild('max_iter_number')


    def setSolveurPrecision(self, label, value):
        """ Put value of solveur precision for variable labelled label """
        # for pressure default value always equal to 1e-8
        self.isPositiveFloat(value)
        node = self._getSolveurLabelNode(label)
        if self._isPressure(node): 
            default = self._defaultValues()['solveur_precision_pressure']
        else:
            default = self._defaultValues()['solveur_precision']

        if value != default:
            node.xmlSetData('solveur_precision', value)
        else:
            node.xmlRemoveChild('solveur_precision')


    def getMaxIterNumber(self, label):
        """ Return number of maximum iterations for variable labelled label """
        node = self._getSolveurLabelNode(label)
        value = node.xmlGetInt('max_iter_number')
        if value == None:
            value = self._defaultValues()['max_iter_number']
        return value


    def getSolveurPrecision(self, label):
        """ Return value of solveur precision for variable labelled label """
        node = self._getSolveurLabelNode(label)

        if self._isPressure(node): 
            default = self._defaultValues()['solveur_precision_pressure']
        else:
            default = self._defaultValues()['solveur_precision']

        value = node.xmlGetDouble('solveur_precision')
        if value == None:
            value = default
        return value


    def getScalarTimeStepFactor(self, label):
        """ Return value of time_step_factor for variable labelled label """
        if self.isScalar(label):
            node = self._getSolveurLabelNode(label)
            value = node.xmlGetDouble('time_step_factor')
            if value == None:
                value = self._defaultValues()['time_step_factor']
            return value
        else:
            raise ValueError, "This method runs only with scalar label"


    def setScalarTimeStepFactor(self, label, value):
        """ Put value of time_step_factor for variable labelled label """
        self.isStrictPositiveFloat(value)
        if self.isScalar(label):
            node = self._getSolveurLabelNode(label)
            if value != self._defaultValues()['time_step_factor']:
                node.xmlSetData('time_step_factor', value)
            else:
                node.xmlRemoveChild('time_step_factor')
        else:
            raise ValueError, "This method runs only with scalar label"


#-------------------------------------------------------------------------------
# NumericalParamEquat test case
#-------------------------------------------------------------------------------


class NumericalParamEquatTestCase(ModelTest):
    """
    """
    def checkNumericalParamEquatInstantiation(self):
        """
        Check whether the NumericalParamEquatModel class could be instantiated
        """
        model = None
        model = NumericalParamEquatModel(self.case)
        assert model != None, 'Could not instantiate NumericalParamEquatModel'

    def checkSetAndGetScheme(self):
        """
        Check whether the NumericalParamEquatModel class could set and get scheme
        """
        model = NumericalParamEquatModel(self.case)
        model.setScheme('VelocitW', 'upwind')
        doc = """<velocity_pressure>
                        <variable label="Pressure" name="pressure"/>
                        <variable label="VelocitU" name="velocity_U"/>
                        <variable label="VelocitV" name="velocity_V"/>
                        <variable label="VelocitW" name="velocity_W">
                                <order_scheme choice="upwind"/>
                        </variable>
                        <property label="total_pressure" name="total_pressure"/>
                        <property label="Yplus" name="yplus" support="boundary"/>
                        <property label="Efforts" name="effort" support="boundary"/>
                        <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set scheme in NumericalParamEquationModel'
        assert model.getScheme('VelocitW') == 'upwind',\
                'Could not get scheme in NumericalParamEquationModel'

    def checkSetAndGetBlendingFactor(self):
        """
        Check whether the NumericalParamEquatModel class could set and get blending factor
        """
        model = NumericalParamEquatModel(self.case)
        model.setScheme('VelocitW', 'centered')
        model.setBlendingFactor('VelocitW', 0.5)
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W">
                            <blending_factor>0.5</blending_factor>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set blending factor in NumericalParamEquationModel'
        assert model.getBlendingFactor('VelocitW') == 0.5,\
                'Could not get blending factor in NumericalParamEquationModel'

    def checkSetAndGetSlopeTest(self):
        """
        Check whether the NumericalParamEquatModel class could set and get slope test
        """
        model = NumericalParamEquatModel(self.case)
        model.setSlopeTest('VelocitW', 'off')
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W">
                            <slope_test status="off"/>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set status of slope test in NumericalParamEquationModel'
        assert model.getSlopeTest('VelocitW') == 'off',\
                'Could not get status of slope test in NumericalParamEquationModel'


    def checkSetAndGetFluxReconstruction(self):
        """
        Check whether the NumericalParamEquatModel class could set and get flux reconstruction
        """
        model = NumericalParamEquatModel(self.case)
        model.setFluxReconstruction('VelocitW', 'on')
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W">
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set status of flux reconstruction in NumericalParamEquationModel'
        assert model.getFluxReconstruction('VelocitW') == 'on',\
                'Could not get status of flux reconstruction in NumericalParamEquationModel'
                
        model.setFluxReconstruction('VelocitW', 'off')
        doc2 = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W">
                        <flux_reconstruction status="off"/>
                    </variable>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                  </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc2),\
                'Could not set status of flux reconstruction in NumericalParamEquationModel'
        assert model.getFluxReconstruction('VelocitW') == 'off',\
                'Could not get status of flux reconstruction in NumericalParamEquationModel'

    def checkSetAndGetMaxIterNumber(self):
        """
        Check whether the NumericalParamEquatModel class could set and get max of number of iterations
        """
        model = NumericalParamEquatModel(self.case)
        model.setMaxIterNumber('Pressure', 22222)
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure">
                            <max_iter_number>22222</max_iter_number>
                    </variable>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set max of number of iterations in NumericalParamEquationModel'
        assert model.getMaxIterNumber('Pressure') == 22222,\
                'Could not get max of number of iterations in NumericalParamEquationModel'

    def checkSetAndGetSolveurPrecision(self):
        """
        Check whether the NumericalParamEquatModel class could set and get solveur precision
        """
        model = NumericalParamEquatModel(self.case)
        
        assert model.getSolveurPrecision('Pressure') == 1e-8,\
                'Could not get solveur precision for pressure in NumericalParamEquationModel'
        from Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
        NumericalParamGlobalModel(self.case).setTimeSchemeOrder(2)
        del NumericalParamGlobalModel
        assert model.getSolveurPrecision('VelocitU') == 1e-5
        
        model.setSolveurPrecision('VelocitU', 2e-6)
        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U">
                            <solveur_precision>2e-06</solveur_precision>
                    </variable>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        assert model.node_vitpre == self.xmlNodeFromString(doc),\
                'Could not set solveur precision in NumericalParamEquationModel'
        assert model.getSolveurPrecision('VelocitU') == 2e-6,\
                'Could not get solveur precision in NumericalParamEquationModel'

    def checkSetAndGetScalarTimeStepFactor(self):
        """
        Check whether the NumericalParamEquatModel class could set and get time step factor
        """
        model = NumericalParamEquatModel(self.case)
        from Pages.ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel

##        self.failUnlessRaises(ValueError, model.setScalarTimeStepFactor('VelocitU', 25.), \
##           'Could not set time step factor in NumericalParamEquationModel')

        model.setScalarTimeStepFactor('Temp.C', 52.)
        node_sca = self.case.xmlGetNode('additional_scalars')
        vit = """<velocity_pressure>
                    <variable label="Pressure" name="pressure"/>
                    <variable label="VelocitU" name="velocity_U"/>
                    <variable label="VelocitV" name="velocity_V"/>
                    <variable label="VelocitW" name="velocity_W"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Efforts" name="effort" support="boundary"/>
                    <property label="all_variables" name="all_variables" support="boundary"/>
                 </velocity_pressure>"""
        sca = """<additional_scalars>
                    <scalar label="Temp.C" name="temperature_celsius" type="thermal">
                            <initial_value zone="1">20.0</initial_value>
                            <min_value>-1e+12 </min_value>
                            <max_value>1e+12</max_value>
                            <time_step_factor>52</time_step_factor>
                    </scalar>
                 </additional_scalars>"""
        assert model.node_vitpre == self.xmlNodeFromString(vit),\
                'Could not set time step factor in NumericalParamEquationModel'

        assert node_sca == self.xmlNodeFromString(sca),\
                'Could not set time step factor for scalar in NumericalParamEquationModel'

##        self.failUnlessRaises(ValueError, model.getScalarTimeStepFactor('VelocitV'), \
##           'Could not get time step factor in NumericalParamEquationModel')

        assert model.getScalarTimeStepFactor('Temp.C') == 52.,\
                'Could not get time step factor for scalar in NumericalParamEquationModel'

def suite():
    testSuite = unittest.makeSuite(NumericalParamEquatTestCase, "check")
    return testSuite

def runTest():
    print "NumericalParamEquatTestCase"
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
