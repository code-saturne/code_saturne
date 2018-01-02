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
This module defines the values of reference.

This module contains the following classes and function:
- SteadyManagementModel
- SteadyManagementTestCase
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
from code_saturne.Base.XMLmodel import XMLmodel, ModelTest
from code_saturne.Pages.TimeStepModel import TimeStepModel

#-------------------------------------------------------------------------------
#  SteadyManagement model class
#-------------------------------------------------------------------------------

class SteadyManagementModel(Model):
    """
    Manage the input/output markups in the xml doc about Pressure
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_anal   = self.case.xmlGetNode('analysis_control')


    def defaultValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        default['status']                 = "off"
        default['iterations']             = 10
        default['relaxation_coefficient'] = 0.7
        default['zero_iteration']         = "off"

        return default


    @Variables.undoGlobal
    def setSteadyFlowManagement(self, steady):
        """
        Set steady flow management balise into xml file.
        """
        self.isOnOff(steady)
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        node['status'] = steady

        mdl_time = TimeStepModel(self.case)

        if steady == 'on':
            mdl_time.node_time.xmlRemoveChild('property', name='courant_number')
            mdl_time.node_time.xmlRemoveChild('property', name='fourier_number')
            mdl_time.node_time.xmlRemoveChild('property', name='local_time_step')
            self.case.xmlRemoveChild('time_average')
            self.getZeroIteration()
            self.getNbIter()
            self.getRelaxCoefficient()
            from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
            model = NumericalParamGlobalModel(self.case).getVelocityPressureAlgorithm()
            if model == 'piso':
                self.case.xmlRemoveChild('velocity_pressure_algo')
            del NumericalParamGlobalModel
        else:
            mdl_time.setTimePassing(0)
            # Treatment of SIMPLE algorithm
            from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
            model = NumericalParamGlobalModel(self.case).getVelocityPressureAlgorithm()
            if model == 'simple':
                self.case.xmlRemoveChild('velocity_pressure_algo')
            del NumericalParamGlobalModel


    @Variables.noUndo
    def getSteadyFlowManagement(self):
        """
        Get status of steady flow management balise fromxml file.
        """
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        status = node['status']
        if not status:
            status = self.defaultValues()['status']
            self.setSteadyFlowManagement(status)
        return status


    @Variables.undoLocal
    def setRelaxCoefficient(self, value):
        """
        Set value of coefficient of relaxation into xml file.
        """
        self.isGreater(value, 0.)
        self.isLowerOrEqual(value, 1.)
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        node.xmlSetData('relaxation_coefficient', value)


    @Variables.undoLocal
    def setNbIter(self, value):
        """
        Set value of iterations number into xml file.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0.)
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        node.xmlSetData('iterations', value)


    @Variables.undoLocal
    def setZeroIteration(self, status):
        """
        Set status of option of zero iteration into xml file.
        """
        self.isOnOff(status)
        node_steady = self.node_anal.xmlInitNode('steady_management', 'status')
        node = node_steady.xmlInitChildNode('zero_iteration', 'status')
        node['status'] = status


    @Variables.noUndo
    def getRelaxCoefficient(self):
        """
        Get value of coefficient of relaxation from xml file.
        """
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        coef = node.xmlGetDouble('relaxation_coefficient')
        if not coef:
            coef = self.defaultValues()['relaxation_coefficient']
            self.setRelaxCoefficient(coef)

        return coef


    @Variables.noUndo
    def getNbIter(self):
        """
        Get value of coefficient of relaxation from xml file.
        """
        node = self.node_anal.xmlInitNode('steady_management', 'status')
        value = node.xmlGetInt('iterations')
        if not value:
            value = self.defaultValues()['iterations']
            self.setNbIter(value)

        return value


    @Variables.noUndo
    def getZeroIteration(self):
        """
        Get status of option of zero iteration from xml file.
        """
        node_steady = self.node_anal.xmlInitNode('steady_management', 'status')
        node = node_steady.xmlGetChildNode('zero_iteration')
        if not node or not node['status']:
            self.setZeroIteration(self.defaultValues()['zero_iteration'])
            node = node_steady.xmlGetChildNode('zero_iteration')

        return node['status']


#-------------------------------------------------------------------------------
# SteadyManagement Model test case
#-------------------------------------------------------------------------------


class SteadyManagementTestCase(ModelTest):
    """
    """
    def checkSteadyManagementInstantiation(self):
        """Check whether the SteadyManagementModel class could be instantiated"""
        model = None
        model = SteadyManagementModel(self.case)
        assert model != None, 'Could not instantiate SteadyManagementModel'

    def checkSetandGetRelaxCoefficient(self):
        """Check whether the SteadyManagementModel class could be set or get relax coefficient """
        mdl = SteadyManagementModel(self.case)
        mdl.setSteadyFlowManagement('on')
        node = mdl.node_anal.xmlInitNode('steady_management', 'status')
        mdl.setRelaxCoefficient(0.45)
        doc = """<steady_management status="on">
                   <zero_iteration status="off"/>
                   <iterations>10</iterations>
                   <relaxation_coefficient>0.45</relaxation_coefficient>
                 </steady_management>"""
        assert node == self.xmlNodeFromString(doc),\
                     'Could not set a relax coefficient'
        coef = mdl.getRelaxCoefficient()
        assert coef == 0.45, 'Could not get a relax coefficient in SteadyManagementModel '

    def checkSetandGeNbIter(self):
        """
        Check whether the SteadyManagementModel class could be
        set or get number of iterations
        """
        mdl = SteadyManagementModel(self.case)
        mdl.setSteadyFlowManagement('on')
        node = mdl.node_anal.xmlInitNode('steady_management', 'status')
        mdl.setNbIter(33)
        doc = """<steady_management status="on">
                    <zero_iteration status="off"/>
                    <iterations>33</iterations>
                    <relaxation_coefficient>0.7</relaxation_coefficient>
                 </steady_management>"""
        assert node == self.xmlNodeFromString(doc),\
                    'Could not set a number of iterations'
        assert mdl.getNbIter() == 33,\
            'Could not get a number of iterations in SteadyManagementModel'

    def checkSetandGetZeroIteration(self):
        """
        Check whether the SteadyManagementModel class could be
        set or get zero iteration status
        """
        mdl = SteadyManagementModel(self.case)
        mdl.setSteadyFlowManagement('on')
        node = mdl.node_anal.xmlInitNode('steady_management', 'status')
        mdl.setZeroIteration('on')
        doc = """<steady_management status="on">
                   <zero_iteration status="on"/>
                   <iterations>10</iterations>
                   <relaxation_coefficient>0.7</relaxation_coefficient>
                 </steady_management>"""
        assert node == self.xmlNodeFromString(doc),\
                    'Could not set a status for zero iteration'
        stat = mdl.getZeroIteration()
        assert stat == 'on', 'Could not get a status for zero iteration in SteadyManagementModel'

def suite():
    testSuite = unittest.makeSuite(SteadyManagementTestCase, "check")
    return testSuite

def runTest():
    print("SteadyManagementTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
