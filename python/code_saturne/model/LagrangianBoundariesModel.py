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

"""
This module defines the lagrangian two phase flow modelling management.

This module contains the following classes and function:
- LagrangianBoundariesModel
- LagrangianBoundariesTestCase
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import sys, unittest, logging


#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.LagrangianModel import LagrangianModel
from code_saturne.model.CoalCombustionModel import CoalCombustionModel


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------


logging.basicConfig()
log = logging.getLogger("LagrangianBoundariesModel")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# lagrangian model class
#-------------------------------------------------------------------------------


class LagrangianBoundariesModel(Model):
    """
    Manage the input/output markups in the xml doc about Lagrangian module.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_boundaries = self.case.root().xmlInitNode('boundary_conditions')
        self.default = self.defaultParticlesBoundaryValues()


    def defaultParticlesBoundaryValues(self):
        """
        Return a dictionnary which contains default values.
        """
        default = {}
        default['particles'] = "inlet"
        default['n_is'] = 0
        default['number'] = 10
        default['frequency'] = 1
        default['statistical_groups'] = 0.
        default['statistical_weight_choice'] = "prescribed"
        default['statistical_weight'] = 1.0
        default['mass_flow_rate'] = 0.
        default['density'] = 1000.
        default['velocity_choice'] = "fluid"
        default['velocity_norm'] = 0.
        default['velocity_value'] = 0.
        default['temperature_choice'] = "fluid"
        default['temperature'] = 20.
        default['specific_heat'] = 1400.
        default['emissivity'] = 0.9
        default['diameter'] = 1.0e-5
        default['diameter_standard_deviation'] = 0.
        default['fouling_index'] = 1.
        default['coal_number'] = 1
        default['coal_temperature'] = 800.
        return default

    def getFoulingStatus(self):
        """
        Return fouling status
        """
        return LagrangianModel(self.case).getCoalFouling()


    @Variables.undoGlobal
    def setBoundaryChoice(self, nature, labelbc, value):
        """
        Update value for the boundary condition. Here we defined the xml nodes
        'self.node_boundary' and 'self.node_particles' used in many functions.
        """
        if nature == "inlet":
            self.isInList(value, ["inlet", "bounce", "outlet"])
        elif nature == "outlet":
            self.isInList(value, ["outlet"])
        elif nature == "free_inlet_outlet":
            self.isInList(value, ["inlet", "outlet"])
        elif nature == "imposed__outlet":
            self.isInList(value, ["outlet"])
        elif nature == "symmetry":
            self.isInList(value, ["part_symmetry", "bounce"])
        elif nature == "wall":
            l = [ "inlet", "bounce", "deposit1", "deposit2"]
            if LagrangianModel(self.case).getCoalFouling() == "on":
                l.append("fouling")
            self.isInList(value, l)
        self.node_boundary = self.node_boundaries.xmlInitChildNode(nature, label=labelbc, field_id='none')
        self.node_particles = self.node_boundary.xmlInitChildNode('particles', 'choice')
        self.node_particles['choice'] = value
        self.setCurrentBoundaryNode(nature, labelbc)


    @Variables.noUndo
    def getBoundaryChoice(self, nature, labelbc):
        """
        Return value for the boundary condition.
        """
        default = { "wall" : "deposit1", "inlet" : "inlet",
                    "outlet" : "outlet", "free_inlet_outlet" : "outlet",
                    "imposed_p_outlet" : "outlet",
                    "symmetry" : "part_symmetry"}
        self.setCurrentBoundaryNode(nature, labelbc)
        if self.node_particles:
            val = self.node_particles['choice']
            if val is None or val == "":
                val = default[nature]
                self.setBoundaryChoice(nature, labelbc, val)
        return val


    @Variables.undoLocal
    def setCurrentBoundaryNode(self, nature, labelbc):
        """
        Update the current boundary node.
        """
        self.node_boundary = self.node_boundaries.xmlInitChildNode(nature, label=labelbc, field_id='none')
        self.node_particles = self.node_boundary.xmlInitChildNode('particles', 'choice')


    def newSetNode(self):
        """
        Add a new 'set' node with child nodes.
        """
        node_set = self.node_particles.xmlAddChild('class')
        node_set.xmlSetData('number', self.default['number'])
        node_set.xmlSetData('frequency', self.default['frequency'])
        node_set.xmlSetData('statistical_groups', self.default['statistical_groups'])
        node_set.xmlSetData('mass_flow_rate', self.default['mass_flow_rate'])
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") == 'off':
            node_set.xmlSetData('density', self.default['density'])
            node_set.xmlInitChildNode('temperature',
                                      choice=self.default['temperature_choice'])
            node_set.xmlSetData('temperature',
                                self.default['temperature'])

        node_set.xmlInitChildNode('statistical_weight',
                                  choice=self.default['statistical_weight_choice'])
        node_set.xmlSetData('statistical_weight',
                            self.default['statistical_weight'])

        node_set.xmlInitChildNode('velocity',
                                  choice=self.default['velocity_choice'])

        node_set.xmlInitChildNode('diameter')
        node_set.xmlSetData('diameter', self.default['diameter'])
        node_set.xmlSetData('diameter_standard_deviation',
                            self.default['diameter_standard_deviation'])
        node_set.xmlSetData('fouling_index', self.default['fouling_index'])


    @Variables.undoGlobal
    def setNumberOfSetsValue(self, value):
        """
        Update the number of sets. Create or delete nodes if necessary.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node_list = self.node_particles.xmlGetChildNodeList('class')
        nnodes = len(node_list)
        if value > nnodes:
            for i in range(value-nnodes):
                self.newSetNode()
        else:
            for i in range(nnodes - value):
                node_to_delete = node_list.pop()
                node_to_delete.xmlRemoveNode()
            # redefine self.node_set


    @Variables.noUndo
    def getNumberOfSetsValue(self):
        """
        Return the number of injection sets.
        """
        node_list = self.node_particles.xmlGetChildNodeList('class')
        value = len(node_list)
        if value is None:
            value = self.defaultParticlesBoundaryValues()['n_is']
            self.setNumberOfSetsValue(value)
        return value


    @Variables.undoLocal
    def setCurrentSetNode(self, iset):
        """
        Update the current set node.
        """
        choice = self.node_particles['choice']
        self.isInList(choice, ["inlet"])
        self.isInt(iset)
        self.node_set = None
        nodes_list = self.node_particles.xmlGetChildNodeList('class')
        if nodes_list:
            nnodes = len(nodes_list)
            self.isLowerOrEqual(iset, nnodes)
            self.node_set = nodes_list[iset-1]


    @Variables.undoLocal
    def setNumberOfParticulesInSetValue(self, label, iset, value):
        """
        Update the number of particles in a set.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('number', value)


    @Variables.noUndo
    def getNumberOfParticulesInSetValue(self, label, iset):
        """
        Return the number of particles in a set.
        """
        value = self.node_set.xmlGetInt('number')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['number']
            self.setNumberOfParticulesInZoneValue(label, iset,value)
        return value


    @Variables.undoLocal
    def setInjectionFrequencyValue(self, label, iset, value):
        """
        Update the injection frequency.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('frequency', value)


    @Variables.noUndo
    def getInjectionFrequencyValue(self, label, iset):
        """
        Return the injection frequency.
        """
        value = self.node_set.xmlGetInt('frequency')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['frequency']
            self.setInjectionFrequencyValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setParticleGroupNumberValue(self, label, iset, value):
        """
        Update the group number of the particle.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('statistical_groups', value)


    @Variables.noUndo
    def getParticleGroupNumberValue(self, label, iset):
        """
        Return the group number of the particle.
        """
        value = self.node_set.xmlGetInt('statistical_groups')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['statistical_groups']
            self.setParticleGroupNumberValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setMassFlowRateValue(self, label, iset, value):
        """
        Update the mass flow rate value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('mass_flow_rate', value)


    @Variables.noUndo
    def getMassFlowRateValue(self, label, iset):
        """
        Return the mass flow rate value.
        """
        value = self.node_set.xmlGetDouble('mass_flow_rate')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['mass_flow_rate']
            self.setMassFlowRateValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setStatisticalWeightChoice(self, label, iset, value):
        """
        Update the condition on statistical weight.
        """
        self.isInList(value, ["rate", "prescribed"])
        node = self.node_set.xmlInitChildNode('statistical_weight', 'choice')
        node['choice'] = value


    @Variables.noUndo
    def getStatisticalWeightChoice(self, label, iset):
        """
        Return the condition on statistical weight.
        """
        val = None
        node = self.node_set.xmlGetChildNode('statistical_weight', 'choice')
        if node:
            val = node['choice']
        if val is None or val == "":
            val = self.defaultParticlesBoundaryValues()['statistical_weight_choice']
        return val


    @Variables.undoLocal
    def setStatisticalWeightValue(self, label, iset, value):
        """
        Update the statistical weight value.
        """
        self.isFloat(value)
        self.isGreater(value, 0)
        self.node_set.xmlSetData('statistical_weight', value)


    @Variables.noUndo
    def getStatisticalWeightValue(self, label, iset):
        """
        Return the statistical weight value.
        """
        value = self.node_set.xmlGetDouble('statistical_weight')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['statistical_weight']
            self.setStatisticalWeightValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setDensityValue(self, label, iset, value):
        """
        Update the density value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('density', value)


    @Variables.noUndo
    def getDensityValue(self, label, iset):
        """
        Return the density value.
        """
        value = self.node_set.xmlGetDouble('density')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['density']
            self.setDensityValue(label, iset, value)
        return value

    @Variables.undoLocal
    def setFoulingIndexValue(self, label, iset, value):
        """
        Update the fouling index value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('fouling_index', value)


    @Variables.noUndo
    def getFoulingIndexValue(self, label, iset):
        """
        Return the fouling index value.
        """
        value = self.node_set.xmlGetDouble('fouling_index')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['fouling_index']
            self.setFoulingIndexValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setVelocityChoice(self, label, iset, choice):
        """
        Update the condition on velocity.
        """
        self.isInList(choice, ["fluid", "components", "norm"])
        node_velocity = self.node_set.xmlInitChildNode('velocity', 'choice')
        node_velocity['choice'] = choice
        if choice in ["fluid", "norm"]:
            node_velocity.xmlRemoveChild('velocity_x')
            node_velocity.xmlRemoveChild('velocity_y')
            node_velocity.xmlRemoveChild('velocity_z')
        elif choice in ["fluid", "components"]:
            node_velocity.xmlRemoveChild('norm')


    @Variables.noUndo
    def getVelocityChoice(self, label, iset):
        """
        Return the condition on velocity.
        """
        node = self.node_set.xmlInitChildNode('velocity', 'choice')
        if node:
            val = node['choice']
            if val is None:
                val = self.defaultParticlesBoundaryValues()['velocity_choice']
                self.setVelocityChoice(val)
        return val


    @Variables.undoLocal
    def setVelocityNormValue(self, label, iset, value):
        """
        Update the velocity norm.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0.)
        node_velocity = self.node_set.xmlInitChildNode('velocity', choice="norm")
        choice = node_velocity['choice']
        self.isInList(choice, ["norm"])
        node_velocity.xmlSetData('norm', value)


    @Variables.noUndo
    def getVelocityNormValue(self, label, iset):
        """
        Return the velocity norm.
        """
        node_velocity = self.node_set.xmlInitChildNode('velocity', choice="norm")
        value = node_velocity.xmlGetDouble('norm')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['velocity_norm']
            self.setVelocityNormValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setVelocityDirectionValue(self, label, iset, idir, value):
        """
        Update the velocity value in the given direction.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0.)
        node_velocity = self.node_set.xmlInitChildNode('velocity', choice="components")
        choice = node_velocity['choice']
        self.isInList(choice, ["components"])
        node_velocity.xmlSetData('velocity_' + idir, value)


    @Variables.noUndo
    def getVelocityDirectionValue(self, label, iset, idir):
        """
        Return the velocity value in the given direction.
        """
        node_velocity = self.node_set.xmlInitChildNode('velocity', choice="components")
        value = self.node_set.xmlGetDouble('velocity_' + idir)
        if value is None:
            value = self.defaultParticlesBoundaryValues()['velocity_value']
            self.setVelocityDirectionValue(label, iset, idir, value)
        return value


    @Variables.undoLocal
    def setTemperatureChoice(self, label, iset, value):
        """
        Update the condition on temperature.
        """
        self.isInList(value, ["prescribed", "fluid"])
        node = self.node_set.xmlInitChildNode('temperature', 'choice')
        node['choice'] = value


    @Variables.noUndo
    def getTemperatureChoice(self, label, iset):
        """
        Return the condition on temperature.
        """
        val = None
        node = self.node_set.xmlGetChildNode('temperature', 'choice')
        if node:
            val = node['choice']
        if val is None:
            val = self.defaultParticlesBoundaryValues()['temperature_choice']
        return val


    @Variables.undoLocal
    def setTemperatureValue(self, label, iset, value):
        """
        Update the temperature value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('temperature', value)


    @Variables.noUndo
    def getTemperatureValue(self, label, iset):
        """
        Return the temperature value.
        """
        value = self.node_set.xmlGetDouble('temperature')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['temperature']
            self.setTemperatureValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setSpecificHeatValue(self, label, iset, value):
        """
        Update the specific heat value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('specific_heat', value)


    @Variables.noUndo
    def getSpecificHeatValue(self, label, iset):
        """
        Return the specific heat value.
        """
        value = self.node_set.xmlGetDouble('specific_heat')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['specific_heat']
            self.setSpecificHeatValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setEmissivityValue(self, label, iset, value):
        """
        Update the emissivity value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('emissivity', value)


    @Variables.noUndo
    def getEmissivityValue(self, label, iset):
        """
        Return the emissivity value.
        """
        value = self.node_set.xmlGetDouble('emissivity')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['emissivity']
            self.setEmissivityValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setDiameterValue(self, label, iset, value):
        """
        Update the particle diameter value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('diameter', value)


    @Variables.noUndo
    def getDiameterValue(self, label, iset):
        """
        Return the particle diameter value.
        """
        value = self.node_set.xmlGetDouble('diameter')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['diameter']
            self.setDiameterValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setDiameterVarianceValue(self, label, iset, value):
        """
        Update the particle diameter variance value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('diameter_standard_deviation', value)


    @Variables.noUndo
    def getDiameterVarianceValue(self, label, iset):
        """
        Return the particle diameter variance value.
        """
        value = self.node_set.xmlGetDouble('diameter_standard_deviation')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['diameter_standard_deviation']
            self.setDiameterVarianceValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setCoalNumberValue(self, label, iset, value):
        """
        Update the coal number of the particle.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('coal_number', value)


    @Variables.noUndo
    def getCoalNumberValue(self, label, iset):
        """
        Return the coal number of the particle.
        """
        value = self.node_set.xmlGetInt('coal_number')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['coal_number']
            self.setCoalNumberValue(label, iset, value)
        return value


    @Variables.undoLocal
    def setCoalTemperatureValue(self, label, iset, value):
        """
        Update the coal temperature.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_set.xmlSetData('coal_temperature', value)


    @Variables.noUndo
    def getCoalTemperatureValue(self, label, iset):
        """
        Return the coal temperature.
        """
        value = self.node_set.xmlGetDouble('coal_temperature')
        if value is None:
            value = self.defaultParticlesBoundaryValues()['coal_temperature']
            self.setCoalTemperatureValue(label, iset, value)
        return value


#-------------------------------------------------------------------------------
# LagrangianBoundaries test case
#-------------------------------------------------------------------------------


class LagrangianBoundariesTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        self.case = Case()
        XMLinit(self.case).initialize()


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case


    def checkLagrangianBoundariesInstantiation(self):
        """
        Check whether the LagrangianBoundariesModel class could be instantiated
        """
        model = None
        model = LagrangianBoundariesModel(self.case)

        assert model != None, 'Could not instantiate LagrangianBoundariesModel'


    def checkLagrangianBoundariesDefaultValues(self):
        """
        Check the default values
        """
        model = LagrangianBoundariesModel(self.case)
        doc = """"""

        assert model.node_output == self.xmlNodeFromString(doc),\
               'Could not get default values for model'


def suite():
    testSuite = unittest.makeSuite(LagrangianBoundariesTestCase, "check")
    return testSuite


def runTest():
    print("LagrangianBoundariesTestCase TODO*********.")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
