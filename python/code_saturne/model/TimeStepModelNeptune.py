# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import *


class TimeStepModel(MainFieldsModel, Variables, Model):

    """
    This class manages the time step in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case = case
        self.__analysisControl = self.case.xmlGetNode('analysis_control')
        self.__timeParameters = self.__analysisControl.xmlInitNode('time_parameters')


    def defaultValues(self):
        default = {}

        default['time_passing']                 = "constant"
        default['calculation_stop']             = "time"
        default['iterations']                   = 10
        default['reference_time_step']          = 0.001
        default['maximum_time']                 = 1.
        default['max_courant_num']              = 1.0
        default['max_fourier_num']              = 10.0
        default['dt_max_increasing_variation']  = 0.1
        default['dt_max_decreasing_variation']  = 0.5
        default['dtdt0_min']                    = 1e-6
        default['dtdt0_max']                    = 1e6
        return default


    @Variables.undoGlobal
    def setTimePassingChoice(self, model) :
        """
        """
        self.isInList(model, ('constant', 'uniform', 'steady'))

        # CFL/Fo
        for f_id in self.getFieldIdList():
            field_name = self.getLabel(f_id)
            for tag in ['courant_number', 'fourier_number']:
                Variables(self.case).setNewVariableProperty("property",
                                                            "",
                                                            self.__timeParameters,
                                                            f_id,
                                                            tag,
                                                            tag+"_"+field_name)

        oldmodel = None
        childNode = self.__timeParameters.xmlGetNode('time_passing')
        if childNode != None :
            oldmodel = childNode['model']

        childNode = self.__timeParameters.xmlInitChildNode('time_passing')
        childNode.xmlSetAttribute(model = model)

        if oldmodel != None and oldmodel != model and oldmodel == "uniform" :
            childNode = self.__timeParameters.xmlGetNode('dtdt0_min')
            if childNode:
                childNode.xmlRemoveNode()
            childNode = self.__timeParameters.xmlGetNode('dtdt0_max')
            if childNode:
                childNode.xmlRemoveNode()
            childNode = self.__timeParameters.xmlGetNode('dt_max_increasing_variation')
            if childNode:
                childNode.xmlRemoveNode()
            childNode = self.__timeParameters.xmlGetNode('dt_max_decreasing_variation')
            if childNode:
                childNode.xmlRemoveNode()
            for fieldId in self.getFieldIdList():
                childNode = self.__timeParameters.xmlGetNode('max_courant_num',  field_id=fieldId)
                if childNode:
                    childNode.xmlRemoveNode()
                childNode = self.__timeParameters.xmlGetNode('max_fourier_num',  field_id=fieldId)
                if childNode:
                    childNode.xmlRemoveNode()


    @Variables.noUndo
    def getTimePassingChoice(self) :
        """
        """
        childNode = self.__timeParameters.xmlGetNode('time_passing')

        if childNode is None :
            model = self.defaultValues()['time_passing']
            self.setTimePassingChoice(model)

        model = self.__timeParameters.xmlGetNode('time_passing')['model']

        return model


    @Variables.undoGlobal
    def setTimeStopChoice(self, model) :
        """
        """
        self.isInList(model, ('time', 'iteration'))

        oldmodel = None
        childNode = self.__timeParameters.xmlGetNode('calculation_stop_type')
        if childNode != None :
            oldmodel = childNode['model']

        childNode = self.__timeParameters.xmlInitChildNode('calculation_stop_type')
        childNode.xmlSetAttribute(model = model)

        if oldmodel != None and oldmodel != model :
            if model == "time" :
                childNode = self.__timeParameters.xmlGetNode('iterations')
                childNode.xmlRemoveNode()
            else :
                childNode = self.__timeParameters.xmlGetNode('maximum_time')
                childNode.xmlRemoveNode()


    @Variables.noUndo
    def getTimeStopChoice(self) :
        """
        """
        childNode = self.__timeParameters.xmlGetNode('calculation_stop_type')

        if childNode is None :
            model = self.defaultValues()['calculation_stop']
            self.setTimeStopChoice(model)

        model = self.__timeParameters.xmlGetNode('calculation_stop_type')['model']

        return model


    @Variables.undoLocal
    def setTimeStep(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('reference_time_step', value)


    @Variables.noUndo
    def getTimeStep(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('reference_time_step')
        if value is None :
            value = self.defaultValues()['reference_time_step']
            self.setTimeStep(value)
        return value


    @Variables.undoLocal
    def setTimeStepsNumber(self, value) :
        """
        """
        self.isPositiveInt(value)

        self.__timeParameters.xmlSetData('iterations', value)


    @Variables.noUndo
    def getTimeStepsNumber(self) :
        """
        """
        value = self.__timeParameters.xmlGetInt('iterations')
        if value is None :
            value = self.defaultValues()['iterations']
            self.setTimeStepsNumber(value)
        return value


    @Variables.undoLocal
    def setMaximumTime(self, value) :
        """
        """
        self.isFloat(value)

        self.__timeParameters.xmlSetData('maximum_time', value)


    @Variables.noUndo
    def getMaximumTime(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('maximum_time')
        if value is None :
            value = self.defaultValues()['maximum_time']
            self.setMaximumTime(value)
        return value


    @Variables.undoLocal
    def setMinDtDt0Variation(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('dtdt0_min', value)


    @Variables.noUndo
    def getMinDtDt0Variation(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('dtdt0_min')
        if value is None :
            value = self.defaultValues()['dtdt0_min']
            self.setMinDtDt0Variation(value)
        return value


    @Variables.undoLocal
    def setMaxDtDt0Variation(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('dtdt0_max', value)


    @Variables.noUndo
    def getMaxDtDt0Variation(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('dtdt0_max')
        if value is None :
            value = self.defaultValues()['dtdt0_max']
            self.setMaxDtDt0Variation(value)
        return value


    @Variables.undoLocal
    def setMaxDtVariationIncreasing(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('dt_max_increasing_variation', value)


    @Variables.noUndo
    def getMaxDtVariationIncreasing(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('dt_max_increasing_variation')
        if value is None :
            value = self.defaultValues()['dt_max_increasing_variation']
            self.setMaxDtVariationIncreasing(value)
        return value


    @Variables.undoLocal
    def setMaxDtVariationDecreasing(self, value) :
        """
        """
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('dt_max_decreasing_variation', value)


    @Variables.noUndo
    def getMaxDtVariationDecreasing(self) :
        """
        """
        value = self.__timeParameters.xmlGetDouble('dt_max_decreasing_variation')
        if value is None :
            value = self.defaultValues()['dt_max_decreasing_variation']
            self.setMaxDtVariationDecreasing(value)
        return value


    @Variables.undoLocal
    def setMaxCourant(self, fieldId, value) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('max_courant_num', value, field_id = fieldId)


    @Variables.noUndo
    def getMaxCourant(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        value = self.__timeParameters.xmlGetDouble('max_courant_num',  field_id=fieldId)
        if value is None :
            value = self.defaultValues()['max_courant_num']
            self.setMaxCourant(fieldId, value)
        return value


    @Variables.undoLocal
    def setMaxFourier(self, fieldId, value) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isPositiveFloat(value)

        self.__timeParameters.xmlSetData('max_fourier_num', value, field_id = fieldId)


    @Variables.noUndo
    def getMaxFourier(self, fieldId) :
        """
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        value = self.__timeParameters.xmlGetDouble('max_fourier_num',  field_id=fieldId)
        if value is None :
            value = self.defaultValues()['max_fourier_num']
            self.setMaxFourier(fieldId, value)
        return value


    @Variables.undoLocal
    def setStopCriterion(self, type, val):
        """
        Compatibility method with classical code_saturne TimeStepModel,
        for compatibility of scripts and parametric setup.
        """

        if type == 'iterations':
            self.setTimeStepsNumber(val)
        elif type == 'maximum_time':
            self.setMaximumTime(val)
        else:
            msg = "neptune_cfd setStopCriterion does not accept " \
                  + str(type) + " option"
            raise ValueError(msg)


    def getTimeParametersNode(self):

        return self.__timeParameters


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class TimeStepTestCase(ModelTest):
    """
    """
    def checkTimeStepInstantiation(self):
        """Check whether the TimeStepModel class could be instantiated"""
        model = None
        model = TimeStepModel(self.case)
        assert model != None, 'Could not instantiate TimeStepModel'


    def checkGetandSetTimePassingChoice(self):
        """Check whether the TimeStepModel class could set and get TimePassingChoice"""
        mdl = TimeStepModel(self.case)
        mdl.setTimePassingChoice('constant')
        doc = '''<time_parameters>
                         <time_passing model="constant"/>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set TimePassingChoice'
        assert mdl.getTimePassingChoice() == 'constant',\
            'Could not get TimePassingChoice'


    def checkGetandSetTimeStep(self):
        """Check whether the TimeStepModel class could set and get TimeStep"""
        mdl = TimeStepModel(self.case)
        mdl.setTimeStep(1)
        doc = '''<time_parameters>
                         <reference_time_step>
                                 1
                         </reference_time_step>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set TimeStep'
        assert mdl.getTimeStep() == 1,\
            'Could not get TimeStep'


    def checkGetandSetTimeStepsNumber(self):
        """Check whether the TimeStepModel class could set and get TimeStepNumber"""
        mdl = TimeStepModel(self.case)
        mdl.setTimeStepsNumber(125455)
        doc = '''<time_parameters>
                         <iterations>
                                 125455
                         </iterations>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set TimeStepNumber'
        assert mdl.getTimeStepsNumber() == 125455,\
            'Could not get TimeStepNumber'


    def checkGetandSetMaximumTime(self):
        """Check whether the TimeStepModel class could set and get MaximumTime"""
        mdl = TimeStepModel(self.case)
        mdl.setMaximumTime(1.5)
        doc = '''<time_parameters>
                         <maximum_time>
                                 1.5
                         </maximum_time>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set MaximumTime'
        assert mdl.getMaximumTime() == 1.5,\
            'Could not get MaximumTime'


    def checkGetandSetMinMaxDtDt0Variation(self):
        """Check whether the TimeStepModel class could set and get MinMaxDtDt0Variation"""
        mdl = TimeStepModel(self.case)
        mdl.setMaxDtDt0Variation(1)
        mdl.setMinDtDt0Variation(0.5)
        doc = '''<time_parameters>
                         <dtdt0_max>
                                 1
                         </dtdt0_max>
                         <dtdt0_min>
                                 0.5
                         </dtdt0_min>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set MaxDtDt0Variation'
        assert mdl.getMaxDtDt0Variation() == 1,\
            'Could not get MaxDtDt0Variation'
        assert mdl.getMinDtDt0Variation() == 0.5,\
            'Could not get MaxDt0Variation'


    def checkGetandSetMaxDtVariationIncreasingDecreasing(self):
        """Check whether the TimeStepModel class could set and get MaxDtVariationIncreasingDecreasing"""
        mdl = TimeStepModel(self.case)
        mdl.setMaxDtVariationIncreasing(1)
        mdl.setMaxDtVariationDecreasing(0.5)
        doc = '''<time_parameters>
                         <dt_max_increasing_variation>
                                 1
                         </dt_max_increasing_variation>
                         <dt_max_decreasing_variation>
                                 0.5
                         </dt_max_decreasing_variation>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set MaxDtVariationIncreasingDecreasing'
        assert mdl.getMaxDtVariationIncreasing() == 1,\
            'Could not get MaxDtVariationIncreasing'
        assert mdl.getMaxDtVariationDecreasing() == 0.5,\
            'Could not get MaxDtVariationDecreasing'


    def checkGetandSetMaxCourant(self):
        """Check whether the TimeStepModel class could set and get MaxCourant"""
        MainFieldsModel(self.case).addField()
        mdl = TimeStepModel(self.case)
        mdl.setMaxCourant('1',10)
        doc = '''<time_parameters>
                         <max_courant_num field_id="1">
                                 10
                         </max_courant_num>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set MaxCourant'
        assert mdl.getMaxCourant('1') == 10,\
            'Could not get MaxCourant'


    def checkGetandSetMaxFourier(self):
        """Check whether the TimeStepModel class could set and get MaxFourier"""
        MainFieldsModel(self.case).addField()
        mdl = TimeStepModel(self.case)
        mdl.setMaxFourier('1',1)
        doc = '''<time_parameters>
                         <max_fourier_num field_id="1">
                                 1
                         </max_fourier_num>
                 </time_parameters>'''
        assert mdl.getTimeParametersNode() == self.xmlNodeFromString(doc),\
            'Could not set MaxFourier'
        assert mdl.getMaxFourier('1') == 1,\
            'Could not get MaxFourier'


def suite():
    testSuite = unittest.makeSuite(TimeStepTestCase, "check")
    return testSuite


def runTest():
    print("TimeStepTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

