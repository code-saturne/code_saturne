# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
chronologic and historic files .... probes .... used variables ...

This module defines the following classes:
- NumericalParamGlobalModel
- NumericalParamGlobalModelTestCase

"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Variables, Model
from Base.XMLmodel import ModelTest
from Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel

#-------------------------------------------------------------------------------
# NumericalParamGlobal model class
#-------------------------------------------------------------------------------

class NumericalParamGlobalModel(Model):
    """
    Manage the input/output markups in the xml doc about Turbulence
    """
    def __init__(self, case):
        """
        Xml node declaration and supplementary default value settings.
        """
        self.case = case
        self.node_np = self.case.xmlInitNode('numerical_parameters')


    def _defaultValues(self):
        """
        default values
        """
        self.default = {}
        self.default['multigrid'] = 'on'
        self.default['gradient_transposed'] = 'on'
        self.default['pressure_relaxation'] = 1.
        self.default['density_relaxation'] = 0.95
        self.default['velocity_pressure_coupling'] ='off'
        self.default['wall_pressure_extrapolation'] = 'neumann'
        self.default['gradient_reconstruction'] = 0
        self.default['time_scheme_order'] = 1
        return self.default


    def getMultigrid(self):
        """
        Return status of multigrid algorithm for pressure
        """
        node = self.node_np.xmlInitNode('multigrid', 'status')
        status = node['status']
        if not status:
            status = self._defaultValues()['multigrid']
            self.setMultigrid(status)
        return status


    def getTransposedGradient(self):
        """
        Return status of transposed gradient
        """
        node = self.node_np.xmlInitNode('gradient_transposed', 'status')
        status = node['status']
        if not status:
            status = self._defaultValues()['gradient_transposed']
            self.setTransposedGradient(status)
        return status


    def getVelocityPressureCoupling(self):
        """
        Return status of IPUCOU value is activated or not
        """
        node = self.node_np.xmlInitNode('velocity_pressure_coupling', 'status')
        status = node['status']
        if not status:
            status = self._defaultValues()['velocity_pressure_coupling']
            self.setVelocityPressureCoupling(status)
        return status


    def getWallPressureExtrapolation(self):
        """
        Return EXTRAG value
        """
        value = self.node_np.xmlGetString('wall_pressure_extrapolation')
        if not value:
            value = self._defaultValues()['wall_pressure_extrapolation']
            self.setWallPressureExtrapolation(value)
        else:
            if value == '0':
                value = 'neumann'
            else:
                value = 'extrapolation'

        return value


    def getPressureRelaxation(self):
        """
        Return RELAXP value
        """
        value = self.node_np.xmlGetDouble('pressure_relaxation')
        if value == None:
            value = self._defaultValues()['pressure_relaxation']
            self.setPressureRelaxation(value)
        return value


    def getDensityRelaxation(self):
        """
        Return SRROM value
        """
        value = self.node_np.xmlGetDouble('density_relaxation')
        if value == None:
            value = self._defaultValues()['density_relaxation']
            self.setDensityRelaxation(value)
        return value


    def getGradientReconstruction(self):
        """
        Return IMRGRA value : 0, 1, 2, 3, 4 or 5
        """
        node = self.node_np.xmlInitNode('gradient_reconstruction', 'choice')
        choice = node['choice']
        if not choice:
            choice = self._defaultValues()['gradient_reconstruction']
            self.setGradientReconstruction(choice)
        return choice


    def setMultigrid(self, status):
        """
        Put status of multigrid algorithm for pressure
        """
        self.isOnOff(status)
        node = self.node_np.xmlInitNode('multigrid', 'status')
        node['status'] = status


    def setTransposedGradient(self, status):
        """
        Put status of gradient transposed
        """
        self.isOnOff(status)
        node = self.node_np.xmlInitNode('gradient_transposed', 'status')
        node['status'] = status


    def setVelocityPressureCoupling(self, status):
        """
        Put status of velocity_pressure_coupling
        """
        self.isOnOff(status)
        node_ipucou = self.node_np.xmlInitNode('velocity_pressure_coupling', 'status')
        node_ipucou['status'] = status
        if status == 'on':
            node_Tx = node_ipucou.xmlInitNode('property', name='weight_matrix_X')
            node_Ty = node_ipucou.xmlInitNode('property', name='weight_matrix_Y')
            node_Tz = node_ipucou.xmlInitNode('property', name='weight_matrix_Z')

            for (node, val) in [(node_Tx, 'weight_matrix_X'),
                                (node_Ty, 'weight_matrix_Y'),
                                (node_Tz, 'weight_matrix_Z')]:
                if not node['label']: node['label'] = Tool.dicoLabel(val)
                node.xmlInitChildNode('listing_printing', status="off")
                node.xmlInitChildNode('postprocessing_recording', status="off")
        else:
            for node in node_ipucou.xmlGetNodeList('property'):
                node.xmlRemoveNode()


    def setPressureRelaxation(self, value):
        """
        Put value of pressure_relaxation
        """
        self.isPositiveFloat(value)
        self.node_np.xmlSetData('pressure_relaxation', value)


    def setDensityRelaxation(self, value):
        """
        Put value of density_relaxation
        """
        self.isGreater(value, 0.0)
        self.isLowerOrEqual(value, 1.0)
        self.node_np.xmlSetData('density_relaxation', value)


    def setWallPressureExtrapolation(self, value):
        """
        Put value of wall pressure extrapolation
        """
        self.isInList(value, ('neumann', 'extrapolation'))
        if value == 'neumann':
            value = '0'
        else:
            value = '1'
        self.node_np.xmlSetData('wall_pressure_extrapolation', value)


    def setGradientReconstruction(self, value):
        """
        Put value of gradient_reconstruction
        """
        self.isInt(value)
        self.isInList(value, (0, 1, 2, 3, 4, 5))
        node = self.node_np.xmlInitNode('gradient_reconstruction', 'choice')
        node['choice'] = value


    def getTimeSchemeOrder(self):
        """
        Return time scheme order for NumericalParamEquationModel
        """
        # getTimeSchemeOrder: used only by NumericalParamEquationModel
        node = self.node_np.xmlGetNode('time_scheme_order')
        if node:
            order = self.node_np.xmlGetInt('time_scheme_order')
        else:
            order = self._defaultValues()['time_scheme_order']
        return order


    def setTimeSchemeOrder(self, order):
        """
        Set or remove markup of time scheme order for turbulence (LES)
        """
        # setTimeSchemeOrder : used only by Turbulence Model
        self.isInt(order)
        if order == 2:
            self.node_np.xmlSetData('time_scheme_order', 2)
        else:
            self.node_np.xmlRemoveChild('time_scheme_order')

#-------------------------------------------------------------------------------
# NumericalParamEquat test case
#-------------------------------------------------------------------------------

class NumericalParamGlobalTestCase(ModelTest):
    """
    """
    def checkNumericalParamGlobalInstantiation(self):
        """
        Check whether the NumericalParamEquatModel class could be instantiated
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        assert model != None, 'Could not instantiate NumericalParamGlobalModel'


    def checkSetandGetGradTransp(self):
        """
        Check whether the NumericalParamEquatModel class
        could be set and get gradient transposed
        """
        model = NumericalParamGlobalModel(self.case)
        model.setTransposedGradient('on')
        doc = """<numerical_parameters>
                    <gradient_transposed status="on"/>
                 </numerical_parameters>"""

        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set gradient transposed in NumericalParamGlobalModel'
        assert model.getTransposedGradient() == 'on',\
                'Could not get gradient transposed in NumericalParamGlobalModel'

    def checkSetandGetVelPesCoupl(self):
        """
        Check whether the NumericalParamEquatModel class
        could be set and get velocity pressure coupling
        """
        model = NumericalParamGlobalModel(self.case)
        model.setVelocityPressureCoupling('on')
        doc = '''<numerical_parameters>
                    <velocity_pressure_coupling status="on">
                        <property label="VPsolve1" name="weight_matrix_X">
                            <listing_printing status="off"/>
                            <postprocessing_recording status="off"/>
                        </property>
                        <property label="VPsolve2" name="weight_matrix_Y">
                            <listing_printing status="off"/>
                            <postprocessing_recording status="off"/>
                        </property>
                        <property label="VPsolve3" name="weight_matrix_Z">
                            <listing_printing status="off"/>
                            <postprocessing_recording status="off"/>
                        </property>
                    </velocity_pressure_coupling>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set velocity_pressure_coupling in NumericalParamGlobalModel'
        assert model.getVelocityPressureCoupling() == 'on',\
                'Could not get velocity_pressure_coupling in NumericalParamGlobalModel'

    def checkSetandGetWallPressureExtrapolation(self):
        """
        Check whether the NumericalParamEquatModel class could be set
        and get wall pressure extrapolation
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        model.setWallPressureExtrapolation('extrapolation')

        doc = '''<numerical_parameters>
                    <wall_pressure_extrapolation>1</wall_pressure_extrapolation>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set wall pressure extrapolation in NumericalParamGlobalModel'
        assert model.getWallPressureExtrapolation() == 'extrapolation',\
                'Could not get wall pressure extrapolation in NumericalParamGlobalModel'

    def checkSetandGetPressureRelaxation(self):
        """
        Check whether the NumericalParamEquatModel class could be set
        and get pressure relaxation
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        model.setPressureRelaxation(0.88)
        doc = '''<numerical_parameters>
                    <pressure_relaxation>0.88</pressure_relaxation>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set pressure relaxation in NumericalParamGlobalModel'
        assert model.getPressureRelaxation() == 0.88,\
                'Could not get pressure relaxation in NumericalParamGlobalModel'

    def checkSetandGetDensityRelaxation(self):
        """
        Check whether the NumericalParamEquatModel class could be set
        and get density relaxation
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        model.setDensityRelaxation(0.91)
        doc = '''<numerical_parameters>
                    <density_relaxation>0.91</density_relaxation>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set density relaxation in NumericalParamGlobalModel'
        assert model.getDensityRelaxation() == 0.91,\
                'Could not get density relaxation in NumericalParamGlobalModel'

    def checkSetandGetGradientReconstructionruction(self):
        """
        Check whether the NumericalParamEquatModel class could be set
        and get gradient_reconstruction
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        model.setGradientReconstruction(3)
        doc = '''<numerical_parameters>
                    <gradient_reconstruction choice="3"/>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set gradient_reconstruction in NumericalParamGlobalModel'
        assert model.getGradientReconstruction() == "3",\
                'Could not get gradient_reconstruction in NumericalParamGlobalModel'

    def checkSetandGetsetTimeSchemeOrder(self):
        """
        Check whether the NumericalParamEquatModel class could be set
        and get time scheme order
        """
        model = None
        model = NumericalParamGlobalModel(self.case)
        model.setTimeSchemeOrder(2)
        doc = '''<numerical_parameters>
                    <time_scheme_order>2</time_scheme_order>
                 </numerical_parameters>'''
        assert model.node_np == self.xmlNodeFromString(doc),\
                'Could not set time scheme order in NumericalParamGlobalModel'
        assert model.getTimeSchemeOrder() == 2,\
                'Could not get time scheme order in NumericalParamGlobalModel'


def suite():
    testSuite = unittest.makeSuite(NumericalParamGlobalTestCase, "check")
    return testSuite


def runTest():
    print("NumericalParamGlobalModelTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
