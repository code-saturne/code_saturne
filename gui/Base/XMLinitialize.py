# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
from Pages.ElectricalModel import ElectricalModel
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
        self.setNewVariable(node, 'pressure')
        self.setNewVariable(node, 'velocity', dim = '3')
        self.setNewProperty(node, 'total_pressure')
        n = self.setNewProperty(node, 'yplus')
        n['support'] = 'boundary'
        n['label'] = 'Yplus'
        n = self.setNewProperty(node, 'effort')
        n['support'] = 'boundary'
        n['label'] = 'Efforts'
        if not node.xmlGetChildNode('property', name='effort_tangential'):
            n = self.setNewProperty(node, 'effort_tangential')
            n['label'] = 'Efforts, tangential'
            n['support'] = 'boundary'
            n.xmlInitNode('postprocessing_recording')['status']= "off"
        if not node.xmlGetChildNode('property', name='effort_normal'):
            n = self.setNewProperty(node, 'effort_normal')
            n['label'] = 'Efforts, normal'
            n['support'] = 'boundary'
            n.xmlInitNode('postprocessing_recording')['status']= "off"

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
            zone = Zone("VolumicZone", self.case, label = 'all_cells', localization = 'all[]')
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

        # Reference values
        XMLThermoPhysicalNode = self.case.xmlInitNode('thermophysical_models')
        self.__XMLVelocityPressureNode = XMLThermoPhysicalNode.xmlInitNode('velocity_pressure')
        self.__RefValuesNode = XMLThermoPhysicalNode.xmlInitNode('reference_values')

        nodeP = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="pressure")
        if nodeP:
            value = nodeP.xmlGetDouble('reference_pressure')
            if value:
                self.__RefValuesNode.xmlSetData('pressure', value)
                nodeP.xmlRemoveChild('reference_pressure')

        nodeTurb = XMLThermoPhysicalNode.xmlInitNode('turbulence', 'model')
        nodeInit = nodeTurb.xmlGetNode('initialization')
        if nodeInit:
            value = nodeInit.xmlGetDouble('reference_velocity')
            if value:
                self.__RefValuesNode.xmlSetData('velocity', value)
                nodeInit.xmlRemoveChild('reference_velocity')

            value = nodeInit.xmlGetDouble('reference_length')
            if value:
                self.__RefValuesNode.xmlSetData('length', value)
                nodeInit.xmlRemoveChild('reference_length')

        for node in self.case.xmlGetNodeList('scalar'):
            value = node.xmlGetDouble('initial_value', zone_id="1")
            if value != None:
                formula = node['label'] + " = " + str(value) + ";"
                n = node.xmlInitChildNode('formula', zone_id="1")
                n.xmlSetTextNode(formula)
                node.xmlRemoveChild('initial_value', zone_id="1")

        # solver
        XMLNumParameterNode = self.case.xmlInitNode('numerical_parameters')
        node = XMLNumParameterNode.xmlGetNode('multigrid')
        if node:
            if node['status'] == "off":
                if nodeP:
                    nodeP.xmlInitNode('solver_choice', choice='conjugate_gradient')
            node.xmlRemoveNode()

        # hydrostatique pressure
        XMLPhysicalPropNode = self.case.xmlInitNode('physical_properties')
        node = XMLPhysicalPropNode.xmlGetNode('hydrostatic_pressure')
        if node:
            stat = node['status']
            XMLNumParameterNode.xmlInitNode('hydrostatic_pressure', status=stat)
            node.xmlRemoveNode()

        # properties
        nodeF = XMLPhysicalPropNode.xmlInitNode('fluid_properties')
        for prop in ['density', 'molecular_viscosity', 'specific_heat',
                     'thermal_conductivity', 'volumic_viscosity']:
            node = nodeF.xmlGetNode('property', name=prop)
            if node:
                if node['choice'] == 'user_law':
                    node['choice'] = 'variable'

        # Profiles
        compt = 0
        for node in self.case.xmlGetNodeList('profile'):
            if node:
                n = node.xmlGetNode('formula')
                f=n.xmlGetTextNode().replace("t", "s")
                n.xmlSetTextNode(f)
            if node:
                n = node.xmlGetNode("output_type")
                if n == None:
                    freq = node.xmlGetInt("output_frequency")
                    if freq == -1:
                        node.xmlSetData('output_type', "end")
                    else:
                        node.xmlSetData('output_type', "frequency")
            nodeInit = node.xmlGetNode('x1')
            if nodeInit:
                node.xmlRemoveNode()
                compt = compt + 1
        if compt != 0:
            print("Profiles have been removed from your files due to  incompatibility")
            print("You must re-create them")

        # thermal scalar
        for phys in ['solid_fuels', 'gas_combustion', 'joule_effect', 'atmospheric_flows', 'compressible_model']:
            node = XMLThermoPhysicalNode.xmlInitNode(phys, 'model')
            mdl = node['model']
            if mdl and mdl != 'off':
                if phys != 'atmospheric_flows' and phys != 'compressible_model':
                    n = node.xmlGetNode('scalar', name="Enthalpy")
                    if n:
                        n.xmlRemoveNode()
                    ThermalScalarModel(self.case).setThermalModel('enthalpy')
                elif phys == 'atmospheric_flows':
                    if (mdl == "dry"):
                        n = node.xmlGetNode('scalar', name="potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('potential_temperature')
                    else:
                        n = node.xmlGetNode('scalar', name="liquid_potential_temperature")
                        if n:
                            n.xmlRemoveNode()
                        ThermalScalarModel(self.case).setThermalModel('liquid_potential_temperature')
                else:
                    n = node.xmlGetNode('scalar', name="EnergieT")
                    if n:
                        n.xmlRemoveNode()
                    ThermalScalarModel(self.case).setThermalModel('total_energy')

        node = self.case.xmlGetNode('additional_scalars')
        n = node.xmlGetNode('scalar', type='thermal')
        if n:
            nth = XMLThermoPhysicalNode.xmlGetNode('thermal_scalar')
            nthvar = nth.xmlInitNode('scalar', 'type')
            nthvar['type']  = "thermal"
            nthvar['name']  = n['name']
            nthvar['label'] = n['label']
            nthvar.xmlChildsCopy(n)
            n.xmlRemoveNode()

        # update velocity node
        nodeV = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_U")
        if nodeV:
            nodeV['name'] = 'velocity'
            nodeV['dimension'] = '3'

            nodeTmp = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_V")
            if nodeTmp:
                nodeTmp.xmlRemoveNode()
            nodeTmp = self.__XMLVelocityPressureNode.xmlGetNode('variable', name="velocity_W")
            if nodeTmp:
                nodeTmp.xmlRemoveNode()



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
                '<solid_fuels model="off"/>'\
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
