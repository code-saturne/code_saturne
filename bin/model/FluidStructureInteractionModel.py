# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
- FluidStructureInteractionModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------
from code_saturne.model.XMLvariables import  Model, Variables
from code_saturne.model.XMLmodel     import  ModelTest


#-------------------------------------------------------------------------------
# Constants class
#-------------------------------------------------------------------------------
class Constants:
    """
    Define class that manages constants

    This class is especially useful when we need constant string identifier.
    Indeed, if the constant value must be changed, the changes is done at
    one place

    ### example
    const = Constants()

    # define constant 'a'
    const.a = 10

    # raise Exception, "Cannot reassign constant a"
    const.a = 20
    """
    def __getattr__(self, attr):
        """
        Get an attributs
        """
        try:
            return self.__dict__[attr]
        except(KeyError):
            raise AttributeError('A instance has no attribute %s' % attr)


    def __setattr__(self, attr, value):
        """
        Set an attributs
        """
        if attr in list(self.__dict__.keys()):
            raise Exception("Cannot reassign constant %s" % attr)
        else:
            self.__dict__[attr] = value


#-------------------------------------------------------------------------------
# Constant class
#-------------------------------------------------------------------------------
const = Constants()

const.max_iterations_implicitation           = 'max_iterations_implicitation'
const.implicitation_precision                = 'implicitation_precision'
const.displacement_prediction_alpha          = 'displacement_prediction_alpha'
const.displacement_prediction_beta           = 'displacement_prediction_beta'
const.stress_prediction_alpha                = 'stress_prediction_alpha'
const.monitor_point_synchronisation          = 'monitor_point_synchronisation'

#-------------------------------------------------------------------------------
# Mobil Mesh model class
#-------------------------------------------------------------------------------

class FluidStructureInteractionModel(Model):
    """
    Manage the input/output markups in the xml doc about mobil mesh
    """

    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.__node_models = case.xmlGetNode('thermophysical_models')
        self.__node_ale    = self.__node_models.xmlInitChildNode('ale_method')

        self.__defaults = {}
        self.__defaults[const.max_iterations_implicitation]  = 1
        self.__defaults[const.implicitation_precision]  = 1e-05
        self.__defaults[const.displacement_prediction_alpha] = 0.5
        self.__defaults[const.displacement_prediction_beta] = 0
        self.__defaults[const.stress_prediction_alpha] = 2
        self.__defaults[const.monitor_point_synchronisation] = 'off'


    #------------------------------------------------------------------
    # MaxIterations
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setMaxIterations(self, value):
        """
        Set value of maximum of iteration if implicitation into xml file.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 1)
        self.__node_ale.xmlSetData(const.max_iterations_implicitation, value)


    @Variables.noUndo
    def getMaxIterations(self):
        """
        Get value of maximum of iteration if implicitation from xml file.
        """
        return self.__getIntData(const.max_iterations_implicitation,
                                 self.setMaxIterations )


    #------------------------------------------------------------------
    # Precision
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setPrecision(self, value):
        """
        Set value of precision of implicitation into xml file.
        """
        self.isGreater(value, 0.0)
        self.__node_ale.xmlSetData(const.implicitation_precision, value)


    @Variables.noUndo
    def getPrecision(self):
        """
        Get value of precision of implicitation from xml file.
        """
        return self.__getDoubleData(const.implicitation_precision,
                                    self.setPrecision )


    #------------------------------------------------------------------
    # DisplacementPredictionAlpha
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setDisplacementPredictionAlpha(self, value):
        """
        Set value of isplacement prediction alpha into xml file.
        """
        self.__node_ale.xmlSetData(const.displacement_prediction_alpha, value)


    @Variables.noUndo
    def getDisplacementPredictionAlpha(self):
        """
        Get value of displacement prediction alpha from xml file.
        """
        return self.__getDoubleData(const.displacement_prediction_alpha,
                                     self.setDisplacementPredictionAlpha )


    #------------------------------------------------------------------
    # DisplacementPredictionBeta
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setDisplacementPredictionBeta(self, value):
        """
        Set value of isplacement prediction beta into xml file.
        """
        self.__node_ale.xmlSetData(const.displacement_prediction_beta, value)


    @Variables.noUndo
    def getDisplacementPredictionBeta(self):
        """
        Get value of displacement prediction beta from xml file.
        """
        return self.__getDoubleData(const.displacement_prediction_beta,
                                     self.setDisplacementPredictionBeta )


    #------------------------------------------------------------------
    # StressPredictionAlpha
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setStressPredictionAlpha(self, value):
        """
        Set value of stress prediction alpha into xml file.
        """
        self.__node_ale.xmlSetData(const.stress_prediction_alpha, value)


    @Variables.noUndo
    def getStressPredictionAlpha(self):
        """
        Get value of stress prediction alpha from xml file.
        """
        return self.__getDoubleData(const.stress_prediction_alpha,
                                    self.setStressPredictionAlpha )


    #------------------------------------------------------------------
    # Monitor point synchronisation
    #------------------------------------------------------------------
    @Variables.undoLocal
    def setMonitorPointSynchronisation(self, value):
        """
        Set value of monitor point synchronisation into xml file.
        """
        self.__setOnOffXML(const.monitor_point_synchronisation, value)


    @Variables.noUndo
    def getMonitorPointSynchronisation(self):
        """
        Get value of monitor point synchronisation from xml file.
        """
        return self.__getOnOffXML(const.monitor_point_synchronisation,
                                  self.setMonitorPointSynchronisation)


    #------------------------------------------------------------------
    # Helper function
    #------------------------------------------------------------------
    def __getStringData(self, name, setFunction):
        """
        Get string value from xml file.
        """
        value = self.__node_ale.xmlGetString(name)
        return self.__getDefaultDataIfNone(value, name, setFunction)


    def __getDoubleData(self, name, setFunction):
        """
        Get double value from xml file.
        """
        value = self.__node_ale.xmlGetDouble(name)
        return self.__getDefaultDataIfNone(value, name, setFunction)


    def __getIntData(self, name, setFunction):
        """
        Get int value from xml file.
        """
        value = self.__node_ale.xmlGetInt(name)
        return self.__getDefaultDataIfNone(value, name, setFunction)


    def __getDefaultDataIfNone(self, value, name, setFunction):
        """
        Get default value if value is none.
        """
        if value == None or value == "":
            value = self.__defaults[name]
            setFunction(value)
        return value


    def __setOnOffXML(self, name, value):
        """
        Set value of 'on'/'off' xml attribute
        """
        Model().isInList(value, [ 'on', 'off'])
        xmlNode = self.__node_ale.xmlInitNode(name)
        xmlNode['status'] = value


    def __getOnOffXML(self, name, setFunction):
        """
        Get value of 'on'/'off' xml attribut
        """
        node = self.__node_ale.xmlInitNode(name, 'status')
        value = node['status']

        return self.__getDefaultDataIfNone(value, name, setFunction)



    def getNodeALE(self):
        """
        Return the node ALE
        """
        return self.__node_ale


#-------------------------------------------------------------------------------
# FluidStructureInteraction test case
#-------------------------------------------------------------------------------

class FluidStructureInteractionTestCase(ModelTest):
    """
    Test case for FluidStructureInteraction
    """

    def checkFluidStructureInteractionInstantiation(self):
        """
        Check whether the FluidStructureInteraction class could be instantiated
        """
        model = None
        model = FluidStructureInteractionModel(self.case)
        assert model != None, 'Could not instantiate '


    def checkGetandSetPrecision(self):
        """Check whether the FluidStructureInteraction class could be set and get precision"""
        mdl = FluidStructureInteractionModel(self.case)
        mdl.setPrecision(0.001)

        doc = """<ale_method status="off">
                    <implicitation_precision>
                        0.001
                    </implicitation_precision>
                </ale_method>"""

        assert mdl.getNodeALE() == self.xmlNodeFromString(doc), \
            'Could not set fluid structure interaction precision'
        assert mdl.getPrecision() == 0.001, \
            'Could not get fluid structure interaction precision'


    def checkGetandSetMaxIterations(self):
        """Check whether the FluidStructureInteraction class could be set and get max iterations"""
        mdl = FluidStructureInteractionModel(self.case)
        mdl.setMaxIterations(99)

        doc = """<ale_method status="off">
                <max_iterations_implicitation>
                    99
                </max_iterations_implicitation>
                </ale_method>"""

        assert mdl.getNodeALE() == self.xmlNodeFromString(doc), \
            'Could not set fluid structure interaction max iterations'
        assert mdl.getMaxIterations() == 99, \
            'Could not get fluid structure interaction max iteration'


    def checkGetAndSetAdvancedViewFeatures(self):
        """Check whether the FluidStructureInteraction class could be set and get advanced view features"""
        mdl = FluidStructureInteractionModel(self.case)

        mdl.setStressPredictionAlpha(42.0)
        mdl.setDisplacementPredictionBeta(42.42)
        mdl.setDisplacementPredictionAlpha(42.4242)
        mdl.setMonitorPointSynchronisation('on')

        doc = """<ale_method status="off">
                <stress_prediction_alpha>
                42
                </stress_prediction_alpha>
                <displacement_prediction_beta>
                42.42
                </displacement_prediction_beta>
                <displacement_prediction_alpha>
                42.4242
                </displacement_prediction_alpha>
                <monitor_point_synchronisation status="on"/>
                </ale_method>"""

        assert mdl.getNodeALE() == self.xmlNodeFromString(doc), \
            'Could not set fluid structure interaction advanced view features'
        assert mdl.getStressPredictionAlpha() == 42.0, \
            'Could not get fluid structure interaction stress prediction alpha'
        assert mdl.getDisplacementPredictionBeta() == 42.42, \
            'Could not get fluid structure interaction displacement prediction beta'
        assert mdl.getDisplacementPredictionAlpha() == 42.4242, \
            'Could not get fluid structure interaction displacement prediction alpha'
        assert mdl.getMonitorPointSynchronisation() == 'on', \
            'Could not get fluid structure interaction monitor point synchronisation'



def suite():
    """
    Test Suite for FluidStructureInteractionTestCase
    """
    testSuite = unittest.makeSuite(FluidStructureInteractionTestCase, "check")
    return testSuite


def runTest():
    """
    run test
    """
    print("FluidStructureInteractionTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
