# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
This module defines the XML data model in which the user defines the physical
options of the treated case.

This module contains the following classe:
- XMLinit
- XMLinitTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.XMLvariables import Variables
from Base import Toolbox

from Pages.LocalizationModel import Zone, LocalizationModel
from Pages.OutputControlModel import OutputControlModel
from Pages.MobileMeshModel import MobileMeshModel
from Pages.TurbulenceModel import TurbulenceModel
from Pages.InitializationModel import InitializationModel
from Pages.TimeStepModel import TimeStepModel
from Pages.SteadyManagementModel import SteadyManagementModel
from Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from Pages.CoalCombustionModel import CoalCombustionModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.ElectricalModelsModel import ElectricalModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.ThermalRadiationModel import ThermalRadiationModel

#-------------------------------------------------------------------------------
# class XMLinit
#-------------------------------------------------------------------------------

class XMLinit(Variables):
    """
    This class initializes the XML contents of the case.
    """
    def __init__(self, case):
        """
        """
        self.case = case


    def initialize(self):
        """
        Verify that all Heading exist only once in the XMLDocument and
        create the missing heading.
        """
        msg = self.__initHeading()
        if msg:
            return msg

        self.__backwardCompatibility()

        # Initialization (order is important, see turbulenceModelsList method)

        self.node_models = self.case.xmlInitNode('thermophysical_models')
        node = self.node_models.xmlInitNode('velocity_pressure')
        for tag in ('pressure',
                    'velocity_U',
                    'velocity_V',
                    'velocity_W'):
            self.setNewVariable(node, tag)
        self.setNewProperty(node, 'total_pressure')
        n = self.setNewProperty(node, 'yplus')
        n['support'] = 'boundary'
        n['label'] = 'Yplus'
        n = self.setNewProperty(node, 'effort')
        n['support'] = 'boundary'
        n['label'] = 'Efforts'
        n = self.setNewProperty(node, 'all_variables')
        n['support'] = 'boundary'
        OutputControlModel(self.case).addDefaultWriter()
        OutputControlModel(self.case).addDefaultMesh()
        MobileMeshModel(self.case).getMethod()
        TurbulenceModel(self.case).getTurbulenceModel()

        # First Volume Zone definition for all cells -> initialization

        zones = LocalizationModel("VolumicZone", self.case).getZones()
        iok = 0
        for zone in zones:
            if zone.getLabel() == 'all_cells':
                iok = 1
        if iok == 0:
            zone = Zone("VolumicZone", label = 'all_cells', localization = 'all[]')
            LocalizationModel("VolumicZone", self.case).addZone(zone)
            zone = LocalizationModel("VolumicZone", self.case).getCodeNumberOfZoneLabel('all_cells')
            InitializationModel(self.case).getInitialTurbulenceChoice(zone)

        # Time step

        TimeStepModel(self.case).getTimeStep()
        TimeStepModel(self.case).getIterationsNumber()
        TimeStepModel(self.case).getTimePassing()

        # Thermodynamics definitinon

        m = FluidCharacteristicsModel(self.case)
        for tag in ('density',
                    'molecular_viscosity',
                    'specific_heat',
                    'thermal_conductivity'):
            m.getInitialValue(tag)

        # Calculation features

        SteadyManagementModel(self.case).getSteadyFlowManagement()
        ThermalScalarModel(self.case).getThermalScalarModel()
        CoalCombustionModel(self.case).getCoalCombustionModel()
        GasCombustionModel(self.case).getGasCombustionModel()
        ElectricalModel(self.case).getElectricalModel()
        ThermalRadiationModel(self.case).getRadiativeModel()

        return msg


    def __initHeading(self):
        """
        Create if necessary headings from the root element of the case.
        """
        msg = ""
        tagList = ('solution_domain',
                   'thermophysical_models',
                   'numerical_parameters',
                   'physical_properties',
                   'additional_scalars',
                   'boundary_conditions',
                   'analysis_control',
                   'calculation_management')

        for tag in tagList:
            nodeList = self.case.root().xmlInitChildNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. " \
                      "There is more than one occurence of the tag: \n\n" + tag +  \
                      "\n\nThe application will finish. Sorry."

        for tag in tagList:
            nodeList = self.case.xmlInitNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. " \
                      "There is more than one occurence of the tag: \n\n" + tag +  \
                      "\n\nThe application will finish. Sorry."

        return msg


    def __backwardCompatibility(self):
        """
        Change XML in order to ensure backward compatibility.
        """
        for node in self.case.xmlGetNodeList('initial_value', 'zone'):
            node['zone_id'] = node['zone']

        for varNode in self.case.xmlGetNodeList('variable'):
            value = varNode.xmlGetDouble('solveur_precision')
            if value:
                varNode.xmlSetData('solver_precision', value)
                varNode.xmlRemoveChild('solveur_precision')

        XMLSolutionDomainNode = self.case.xmlInitNode('solution_domain')
        self.__XMLVolumicConditionsNode = XMLSolutionDomainNode.xmlInitNode('volumic_conditions')
        for node in self.__XMLVolumicConditionsNode.xmlGetNodeList('zone'):
            if node['id'] == None:
                node['id'] = node['name']

        oldnode = self.case.xmlGetNode('calcul_management')
        if oldnode:
            newnode = self.case.xmlInitNode('calculation_management')
            newnode.xmlChildsCopy(oldnode)
            oldnode.xmlRemoveNode()


#-------------------------------------------------------------------------------
# XMLinit test case
#-------------------------------------------------------------------------------


class XMLinitTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from Base import XMLengine
        Toolbox.GuiParam.lang = 'en'
        self.doc = XMLengine.XMLDocument("")
        self.case = XMLengine.Case(None)


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case
        del self.doc


    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()


    def checkXMLinitInstantiation(self):
        """
        Check whether the Case class could be instantiated
        """
        xmldoc = None
        xmldoc = XMLinit(self.case)
        assert xmldoc != None, 'Could not instantiate XMLinit'


    def checkInitHeading(self):
        """
        Check whether the headings markups could be initialized
        """
        doc = \
        '<Code_Saturne_GUI case="" study="" version="1.0">'\
        '<solution_domain/>'\
        '<thermophysical_models>'\
                '<velocity_pressure>'\
                        '<variable label="Pressure" name="pressure"/>'\
                        '<variable label="VelocityX" name="velocity_U"/>'\
                        '<variable label="VelocityY" name="velocity_V"/>'\
                        '<variable label="VelocityZ" name="velocity_W"/>'\
                        '<property label="total_pressure" name="total_pressure"/>'\
                '</velocity_pressure>'\
                '<turbulence model="k-epsilon">'\
                        '<variable label="TurbEner" name="turb_k"/>'\
                        '<variable label="Dissip" name="turb_eps"/>'\
                        '<property label="TurbVisc" name="turb_viscosity"/>'\
                        '<initialization choice="reference_velocity">'\
                                '<reference_velocity>1.0</reference_velocity>'\
                        '</initialization>'\
                '</turbulence>'\
                '<initialization>'\
                        '<zone name="1">0</zone>'\
                '</initialization>'\
                '<thermal_scalar model="off"/>'\
                '<gas_combustion model="off"/>'\
                '<pulverized_coal model="off"/>'\
                '<joule_effect model="off"/>'\
                '<radiative_transfer model="off"/>'\
        '</thermophysical_models>'\
        '<numerical_parameters/>'\
        '<physical_properties>'\
                '<fluid_properties>'\
                        '<property choice="constant" label="Density" name="density">'\
                                '<listing_printing status="off"/>'\
                                '<postprocessing_recording status="off"/>'\
                                '<initial_value>1.17862</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="LamVisc" name="molecular_viscosity">'\
                                '<listing_printing status="off"/>'\
                                '<postprocessing_recording status="off"/>'\
                                '<initial_value>1.83e-05</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="SpecHeat" name="specific_heat">'\
                                '<listing_printing status="off"/>'\
                                '<postprocessing_recording status="off"/>'\
                                '<initial_value>1017.24</initial_value>'\
                        '</property>'\
                        '<property choice="constant" label="ThermalCond" name="thermal_conductivity">'\
                                '<listing_printing status="off"/>'\
                                '<postprocessing_recording status="off"/>'\
                                '<initial_value>0.02495</initial_value>'\
                        '</property>'\
                '</fluid_properties>'\
        '</physical_properties>'\
        '<additional_scalars/>'\
        '<boundary_conditions/>'\
        '<analysis_control>'\
                '<time_parameters>'\
                        '<time_step_ref>0.1</time_step_ref>'\
                        '<property label="CourantNb" name="courant_number">'\
                        '<property label="FourierNb" name="fourier_number">'\
                '</time_parameters>'\
        '</analysis_control>'\
        '<calculation_management/>'\
        '</Code_Saturne_GUI>'

        XMLinit(self.case).initialize()

        assert self.case.root() == self.xmlNodeFromString(doc), \
               'Could not use the constructor of the XMLinit class'


def suite():
    testSuite = unittest.makeSuite(XMLinitTestCase, "check")
    return testSuite


def runTest():
    print("XMLinitTestCase to be completed...")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of XMLinit
#-------------------------------------------------------------------------------
