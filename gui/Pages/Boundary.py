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

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox      import GuiParam
from code_saturne.Base.XMLvariables import Model, Variables
from code_saturne.Base.XMLmodel     import ModelTest
from code_saturne.Base.XMLengine    import *

from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.Pages.ThermalScalarModel     import ThermalScalarModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("Boundary")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class Boundary(object) :
    """
    Abstract class
    """
    def __new__(cls, nature , label, case) :
        """
        Factory
        """
        if nature == 'inlet':
            return InletBoundary.__new__(InletBoundary, label, case)
        elif nature == 'coal_inlet':
            from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
            return CoalInletBoundary.__new__(CoalInletBoundary, label, case)
        elif nature == 'compressible_outlet':
            from code_saturne.Pages.CompressibleModel import CompressibleModel
            return CompressibleOutletBoundary.__new__(CompressibleOutletBoundary, label, case)
        elif nature == 'outlet':
            return OutletBoundary.__new__(OutletBoundary, label, case)
        elif nature == 'symmetry':
            return SymmetryBoundary.__new__(SymmetryBoundary, label, case)
        elif nature == 'wall':
            return WallBoundary.__new__(WallBoundary, label, case)
        elif nature == 'radiative_wall':
            from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
            Model().isNotInList(ThermalRadiationModel(case).getRadiativeModel(), ("off",))
            return RadiativeWallBoundary.__new__(RadiativeWallBoundary, label, case)
        elif nature == 'coupling_mobile_boundary':
            return CouplingMobilBoundary.__new__(CouplingMobilBoundary, label, case)
        elif nature == 'meteo_inlet' or nature == 'meteo_outlet':
            return MeteoBoundary.__new__(MeteoBoundary, label, case)
        elif nature == 'joule_inlet' or nature == 'joule_outlet' or nature == 'joule_wall':
            return JouleBoundary.__new__(JouleBoundary, label, case)
        elif nature == 'free_inlet_outlet':
            return InletOutletBoundary.__new__(InletOutletBoundary, label, case)
        elif nature == 'imposed_p_outlet':
            return ImposedPressureOutletBoundary.__new__(ImposedPressureOutletBoundary, label, case)
        elif nature == 'free_surface':
            return FreeSurfaceBoundary.__new__(FreeSurfaceBoundary, label, case)
        elif nature == 'groundwater':
            return GroundwaterBoundary.__new__(GroundwaterBoundary, label, case)
        else :
            raise ValueError("Unknown boundary nature: " + nature)


    def __init__(self, nature, label, case) :
        """
        """
        self._label = label
        self._nature = nature
        self.case = case
        self._XMLBoundaryConditionsNode = self.case.xmlGetNode('boundary_conditions')
        self._thermalLabelsList = ('temperature',
                                   'enthalpy',
                                   'total_energy')

        self.sca_model = DefineUserScalarsModel(self.case)

        # Create nodes
        if nature not in ["coal_inlet",
                          "compressible_outlet",
                          "radiative_wall",
                          "coupling_mobile_boundary",
                          "meteo_inlet",
                          "meteo_outlet",
                          "joule_inlet",
                          "joule_outlet",
                          "joule_wall"]:
            self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode(nature, label = label, field_id='none')

        else:
            if nature == "coal_inlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('inlet', label = label, field_id='none')

            elif nature == "compressible_outlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('outlet', label = label, field_id='none')

            elif nature in ["radiative_wall",
                            "coupling_mobile_boundary"]:
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('wall', label = label, field_id='none')

            elif nature == "meteo_inlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('inlet', label = label, field_id='none')

            elif nature == "meteo_outlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('outlet', label = label, field_id='none')

            elif nature == "joule_inlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('inlet', label = label, field_id='none')

            elif nature == "joule_outlet":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('outlet', label = label, field_id='none')

            elif nature == "joule_wall":
                self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('wall', label = label)

        self._initBoundary()
        self._initALEBoundary()


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node (vitual method)
        """
        pass

    def _initALEBoundary(self):
        """
        Initialize the possible choices.
        """
        self.__ALEChoices = ["fixed_boundary",
                             "sliding_boundary",
                             "internal_coupling",
                             "external_coupling",
                             "fixed_velocity",
                             "fixed_displacement"]

        self._defaultValues = {}
        self._defaultValues['ale_choice'] = self.__ALEChoices[0]

        formula_velocity = 'ale_formula_' + "fixed_velocity"
        formula_displacement = 'ale_formula_' + "fixed_displacement"
        self._defaultValues[formula_velocity] = 'mesh_velocity_U = 0;\nmesh_velocity_V = 0;\nmesh_velocity_W = 0;'
        self._defaultValues[formula_displacement] = 'mesh_x = 0;\nmesh_y = 0;\nmesh_z = 0;'


    @Variables.noUndo
    def getALEChoice(self):
        """
        Get the choice ALE
        """
        node = self.boundNode.xmlInitNode('ale', 'choice')
        choice = node['choice']

        # Create a defaut choice if it does not exist
        if not node['choice']:
            choice = self._defaultValues['ale_choice']
            self.setALEChoice(choice)

        return choice


    @Variables.undoGlobal
    def setALEChoice(self, value):
        """
        Set the ALE according to choice
        """
        Model().isInList(value, self.__ALEChoices)
        node = self.boundNode.xmlInitNode('ale')

        # if something has changed
        if node['choice'] != value:
            node['choice'] = value
            if  value in ["fixed_velocity", "fixed_displacement"]:
                self.setFormula(self._getDefaultFormula())
            else:
                node.xmlRemoveChild('formula')


    @Variables.noUndo
    def getFormula(self):
        """
        Get the formula from the xml
        """
        node = self.boundNode.xmlInitChildNode('ale')
        value = node.xmlGetChildString('formula')

        if not value:
            value = self._getDefaultFormula()
            self.setFormula(value)

        return value


    @Variables.undoLocal
    def setFormula(self, value):
        """
        Set the formula into the xml
        """
        node = self.boundNode.xmlInitChildNode('ale')
        node.xmlSetData('formula', value)


    def _getDefaultFormula(self):
        """
        Return the default value for the formula
        """
        if  self.getALEChoice() in ["fixed_velocity", "fixed_displacement"]:
            aleChoice = self.getALEChoice()
            return self._defaultValues[ 'ale_formula_' + aleChoice ]
        else:
            return ''


    def updateScalarTypeAndName(self, scalarNode, scalarName):
        """
        Check and update type of scalar scalarName for boundary conditions for wall
        """
        if self.sca_model.getMeteoScalarsNameList() != None:
            if scalarName in self.sca_model.getMeteoScalarsNameList():
                scalarNode['type'] = self.sca_model.getMeteoScalarType(scalarName)
            else :
                scalarNode['type'] = self.sca_model.getScalarTypeByName(scalarName)
        elif self.sca_model.getElectricalScalarsNameList() != None:
            if scalarName in self.sca_model.getElectricalScalarsNameList():
                scalarNode['type'] = self.sca_model.getElectricalScalarType(scalarName)
            else :
                scalarNode['type'] = self.sca_model.getScalarTypeByName(scalarName)
        else:
            scalarNode['type'] = self.sca_model.getScalarTypeByName(scalarName)


    @Variables.noUndo
    def getLabel(self):
        """
        Return the label
        """
        return self._label


    @Variables.noUndo
    def getNature(self):
        """
        Return the nature
        """
        return self._nature


    def delete(self):
        """
        Delete Boundary
        """
        self.boundNode.xmlRemoveNode()

#-------------------------------------------------------------------------------
# Inlet boundary
#-------------------------------------------------------------------------------

class InletBoundary(Boundary):
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self.__velocityChoices   = ['norm',
                                    'flow1',
                                    'flow2',
                                    'norm_formula',
                                    'flow1_formula',
                                    'flow2_formula']
        self.__directionChoices  = ['normal',
                                    'coordinates',
                                    'formula']
        self.__directionTags     = ['direction_x',
                                    'direction_y',
                                    'direction_z',
                                    'direction_formula']
        self.__turbulenceChoices = ['hydraulic_diameter',
                                    'turbulent_intensity',
                                    'formula']
        self._scalarChoicesList  = ['dirichlet',
                                    'neumann',
                                    'exchange_coefficient',
                                    'dirichlet_formula',
                                    'neumann_formula',
                                    'exchange_coefficient_formula']

        self.typeList            = ['imposed_inlet',
                                    'subsonic_inlet_PH']
        self.typeListGasComb     = ['oxydant',
                                    'fuel',
                                    'unburned',
                                    'burned']

        self.__mappingDirectionTags = ['translation_x',
                                       'translation_y',
                                       'translation_z']

        self.th_model = ThermalScalarModel(self.case)

        # Initialize nodes if necessary

        self.getVelocityChoice()
        self.getDirectionChoice()
        self.getTurbulenceChoice()

        for name in self.sca_model.getScalarNameList():
            if self.getScalarChoice(name) == 'dirichlet':
                self.getScalarValue(name,'dirichlet')

        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        from code_saturne.Pages.CompressibleModel import CompressibleModel
        if CoalCombustionModel(self.case).getCoalCombustionModel() == "off":
            self.boundNode.xmlRemoveChild('coal')
            if CompressibleModel(self.case).getCompressibleModel() == "off":
                self.boundNode.xmlRemoveChild('temperature')


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['velocityChoice']      = 'norm'
        dico['directionChoice']     = 'normal'
        dico['turbulenceChoice']    = 'hydraulic_diameter'
        dico['hydraulic_diameter']  = 1
        dico['turbulent_intensity'] = 2
        dico['velocity']            = 0.0
        dico['flow1']               = 1
        dico['flow2']               = 1
        dico['norm']                = 1
        dico['flow1_formula']       = "q_m = 1;"
        dico['flow2_formula']       = "q_v = 1;"
        dico['norm_formula']        = "u_norm = 1;"
        dico['direction_x']         = 0.0
        dico['direction_y']         = 0.0
        dico['direction_z']         = 0.0
        dico['direction_formula']   = "dir_x = 1;\ndir_y = 0;\ndir_z = 0;\n"
        dico['scalar']              = 0.0
        dico['scalarChoice']        = 'dirichlet'
        dico['pressure']            = 101300.
        dico['density']             = 0.0
        dico['temperature']         = 273.
        dico['energy']              = 0.0
        dico['total_pressure']      = 101300.
        dico['enthalpy']            = 0.0
        dico['status']              = 'off'
        dico['compressible_type']   = 'imposed_inlet'
        dico['fraction']            = 0.0
        dico['pressure']            = 0.0
        dico['pressure_flux']       = 0.0
        dico['hydraulicHeadChoice'] = 'dirichlet'

        from code_saturne.Pages.GasCombustionModel import GasCombustionModel
        model = GasCombustionModel(self.case).getGasCombustionModel()
        del GasCombustionModel
        if model == 'lwp' or model == 'ebu':
            dico['gas_type'] = 'unburned'
        elif model == 'd3p':
            dico['gas_type'] = 'oxydant'

        from code_saturne.Pages.ReferenceValuesModel import ReferenceValuesModel
        dico['temperatureGas'] = ReferenceValuesModel(self.case).getTemperature()

        return dico


    def __initChoiceForVelocityAndDirection(self):
        """
        Get the choice of velocity and direction.
        """
        node = self.boundNode.xmlInitNode('velocity_pressure', 'choice', 'direction')
        choice = node['choice']
        if not choice:
            choice = self.__defaultValues()['velocityChoice']
            self.setVelocityChoice(choice)
        dir = node['direction']
        if not dir:
            dir = self.__defaultValues()['directionChoice']
            self.setDirectionChoice(dir)
        return choice, dir


    def __getscalarList(self):
        """
        return list of scalars name
        """
        scalar_list =  self.sca_model.getScalarNameList()
        if self.sca_model.getMeteoScalarsNameList() != None:
            for sca in self.sca_model.getMeteoScalarsNameList():
                scalar_list.append(sca)
        if len(self.sca_model.getThermalScalarName()) > 0:
            scalar_list.append(self.sca_model.getThermalScalarName()[0])

        return scalar_list


    def __deleteScalarNodes(self, name, tag):
        """
        Delete nodes of scalars
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(tag, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        for tt in self._scalarChoicesList:
            if tt != tag:
                scalarNode.xmlRemoveChild(tt)


    @Variables.noUndo
    def getVelocityChoice(self):
        """
        Get the choice of velocity.
        """
        choice, dir = self.__initChoiceForVelocityAndDirection()
        return choice


    @Variables.noUndo
    def getDirectionChoice(self):
        """
        Get the choice of direction.
        """
        choice, dir = self.__initChoiceForVelocityAndDirection()
        return dir


    @Variables.noUndo
    def getVelocity(self):
        """
        Get value of velocity beyond choice.
        """
        choice = self.getVelocityChoice()
        Model().isInList(choice, self.__velocityChoices)

        XMLVelocityNode = self.boundNode.xmlGetNode('velocity_pressure')

        if choice in ('norm', 'flow1', 'flow2'):
            value = XMLVelocityNode.xmlGetChildDouble(choice)
        elif choice in ('norm_formula', 'flow1_formula', 'flow2_formula'):
            value = XMLVelocityNode.xmlGetChildString(choice)
        if value == None:
            value = self.__defaultValues()[choice]
            self.setVelocity(value)

        return value


    @Variables.undoGlobal
    def setVelocity(self, value):
        """
        Set value of velocity.
        """
        choice = self.getVelocityChoice()
        Model().isInList(choice, self.__velocityChoices)

        if choice in ('norm', 'flow1', 'flow2'):
            Model().isFloat(value)

        XMLVelocityNode = self.boundNode.xmlInitNode('velocity_pressure')
        XMLVelocityNode.xmlSetData(choice, value)


    @Variables.noUndo
    def getDirection(self, component):
        """
        Get the component velocity
        """
        Model().isInList(component, self.__directionTags)

        XMLVelocityNode = self.boundNode.xmlGetNode('velocity_pressure')

        if XMLVelocityNode['direction'] == 'coordinates':
            Model().isInList(component, ('direction_x', 'direction_y', 'direction_z'))
            value = XMLVelocityNode.xmlGetChildDouble(component)
        elif XMLVelocityNode['direction'] == 'formula':
            Model().isInList(component, ('direction_formula',))
            value = XMLVelocityNode.xmlGetChildString(component)

        if value == None :
            value = self.__defaultValues()[component]
            self.setDirection(component, value)
        return value


    @Variables.undoLocal
    def setDirection(self, component, value):
        """
        Set the component velocity for fieldLabel
        """
        Model().isInList(component, self.__directionTags)
        if component != 'direction_formula':
            Model().isFloat(value)

        XMLVelocityNode = self.boundNode.xmlInitNode('velocity_pressure')
        XMLVelocityNode.xmlSetData(component, value)


    @Variables.undoGlobal
    def setVelocityChoice(self, value):
        """
        Set the velocity definition according to choice
        """
        Model().isInList(value, self.__velocityChoices)

        # Check if value is a new velocity choice value
        XMLVelocityNode = self.boundNode.xmlInitNode('velocity_pressure')
        if XMLVelocityNode['choice'] != None :
            if XMLVelocityNode['choice'] == value:
                return

        # Update velocity choice
        XMLVelocityNode['choice'] = value
        self.getVelocity()

        for tag in self.__velocityChoices:
            if tag != value:
                XMLVelocityNode.xmlRemoveChild(tag)


    @Variables.undoGlobal
    def setDirectionChoice(self, value):
        """
        Set the direction of the flow definition according to choice.
        """
        Model().isInList(value, self.__directionChoices)

        # Check if value is a new direction choice
        XMLVelocityNode = self.boundNode.xmlInitNode('velocity_pressure')
        if XMLVelocityNode['direction'] != None :
            if XMLVelocityNode['direction'] == value:
                return

        # Update direction choice
        XMLVelocityNode['direction'] = value

        if value == 'coordinates':
            self.getDirection('direction_x')
            self.getDirection('direction_y')
            self.getDirection('direction_z')
            XMLVelocityNode.xmlRemoveChild('direction_formula')

        elif value == 'formula':
            self.getDirection('direction_formula')
            for tag in ('direction_x', 'direction_y', 'direction_z'):
                XMLVelocityNode.xmlRemoveChild(tag)


    @Variables.noUndo
    def getTurbulenceChoice(self):
        """
        Get the turbulence choice
        """
        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')

        choice = XMLTurbulenceNode['choice']
        if choice not in self.__turbulenceChoices :
            choice = self.__defaultValues()['turbulenceChoice']
            self.setTurbulenceChoice(choice)

        return choice


    @Variables.undoGlobal
    def setTurbulenceChoice(self, value):
        """
        Set the choice turbulence
        """
        Model().isInList(value, self.__turbulenceChoices)

        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')

        if XMLTurbulenceNode['choice'] != None:
            if XMLTurbulenceNode['choice'] == value:
                return

        XMLTurbulenceNode['choice'] = value

        # Update values
        if value == 'hydraulic_diameter' :
            self.getHydraulicDiameter()
            XMLTurbulenceNode.xmlRemoveChild('turbulent_intensity')

        elif value == 'turbulent_intensity' :
            self.getHydraulicDiameter()
            self.getTurbulentIntensity()

        elif value == 'formula' :
            self.getTurbFormula()


    @Variables.noUndo
    def getHydraulicDiameter(self):
        """
        Get hydraulic diameter
        """
        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'],  self.__turbulenceChoices)
        value = XMLTurbulenceNode.xmlGetDouble('hydraulic_diameter')
        if value == None :
            value = self.__defaultValues()['hydraulic_diameter']
            self.setHydraulicDiameter(value)
        return value


    @Variables.undoLocal
    def setHydraulicDiameter(self, value):
        """
        Set hydraulic diameter
        """
        Model().isStrictPositiveFloat(value)

        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], self.__turbulenceChoices)
        XMLTurbulenceNode.xmlSetData('hydraulic_diameter', value)


    @Variables.undoLocal
    def setTurbFormula(self, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')
        if not XMLTurbulenceNode:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = XMLTurbulenceNode.xmlInitChildNode('formula')
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getTurbFormula(self):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')

        formula = XMLTurbulenceNode.xmlGetString('formula')
        return formula


    @Variables.noUndo
    def getDefaultTurbFormula(self, turb_model):
        """
        Get defaut turbulence formula
        """
        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            formula = """k = 0.;
epsilon = 0.;"""

        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            formula = """r11 = 0.;
r22 = 0.;
r33 = 0.;
r12 = 0.;
r13 = 0.;
r23 = 0.;
epsilon = 0.;"""

        elif turb_model == 'Rij-EBRSM':
            formula = """r11 = 0.;
r22 = 0.;
r33 = 0.;
r12 = 0.;
r13 = 0.;
r23 = 0.;
epsilon = 0.;
alpha = 0.;"""

        elif turb_model == 'v2f-BL-v2/k':
            formula = """k = 0.;
epsilon = 0.;
phi = 0.;
alpha = 0.;"""

        elif turb_model == 'k-omega-SST':
            formula = """k = 0.;
omega = 0.;"""

        elif turb_model == 'Spalart-Allmaras':
            formula = """nu_tilda = 0.;"""

        return formula


    @Variables.noUndo
    def getTurbulentIntensity(self):
        """
        Get turbulent intensity
        """
        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], ('turbulent_intensity',))
        value = XMLTurbulenceNode.xmlGetDouble('turbulent_intensity')
        if value == None :
            value = self.__defaultValues()['turbulent_intensity']
            self.setTurbulentIntensity(value)

        return value


    @Variables.undoLocal
    def setTurbulentIntensity(self, value):
        """
        Set turbulent intensity
        """
        Model().isStrictPositiveFloat(value)

        XMLTurbulenceNode = self.boundNode.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], ('turbulent_intensity',))
        XMLTurbulenceNode.xmlSetData('turbulent_intensity', value)


    @Variables.noUndo
    def getScalarChoice(self, scalarName):
        """
        Get scalar choice
        """
        Model().isInList(scalarName, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarChoice']
            self.setScalarChoice(scalarName, choice)

        return choice


    @Variables.undoGlobal
    def setScalarChoice(self, scalarName, choice) :
        """
        Set scalar choice
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        self.__deleteScalarNodes(scalarName, choice)


    @Variables.noUndo
    def getScalarValue(self, scalarName, choice) :
        """
        Get scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        value = scalarNode.xmlGetChildDouble(choice)
        if value == None :
            value = self.__defaultValues()['scalar']
            self.setScalarValue(scalarName, choice, value)
        return value


    @Variables.undoGlobal
    def setScalarValue(self, scalarName, choice, value):
        """
        Set scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isFloat(value)
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        scalarNode.xmlSetData(choice, value)


    @Variables.noUndo
    def getDefaultScalarFormula(self, scalarName, scalar_model):
        """
        Get defaut scalar formula
        """
        if scalar_model == 'dirichlet_formula':
            formula = scalarName+""" = 0;\n"""
        elif scalar_model == 'neumann_formula':
            formula = """flux = 0;\n"""
        elif scalar_model == 'exchange_coefficient_formula':
            formula = scalarName+""" = 0;\nhc = 0;\n"""

        return formula


    @Variables.noUndo
    def getScalarFormula(self, scalarName, choice):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)
        formula = scalarNode.xmlGetChildString(choice)

        if not formula:
            formula = self.getDefaultScalarFormula(scalarName, choice)

        return formula


    @Variables.undoLocal
    def setScalarFormula(self, scalarName, choice, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        n = scalarNode.xmlSetData(choice, formula)


    @Variables.noUndo
    def getInletType(self):
        """
        Return type for velocities's boundary conditions for inlet gas combustion.
        """
        node = self.boundNode.xmlGetNode('velocity_pressure')
        n = node.xmlInitNode('compressible_type')
        type = n['choice']
        if type == None:
            type = self.__defaultValues()['compressible_type']
            self.setInletType(type)
        return type


    @Variables.noUndo
    def getThermoStatus(self, var):
        """
        Return status of var for the initialisation
        """
        node = self.boundNode.xmlGetNode('velocity_pressure')
        n = node.xmlInitNode(var, 'status')
        status = n['status']
        if status == "":
            status = self.__defaultValues()['status']
            self.setThermoStatus(var, status)
        return status


    @Variables.undoLocal
    def setThermoStatus(self, var, status):
        """
        Put status of var for the initialisation
        """
        node = self.boundNode.xmlGetNode('velocity_pressure')
        if status == 'off':
            node.xmlRemoveChild(var)
        n = node.xmlInitNode(var, 'status')
        n['status'] = status


    @Variables.undoLocal
    def setInletType(self, type):
        """
        Set type of inlet.
        """
        Model().isInList(type, self.typeList)
        node = self.boundNode.xmlGetNode('velocity_pressure')
        n= node.xmlInitNode('compressible_type')
        old_model = n['choice']
        if old_model != type and type == "subsonic_inlet_PH":
            for nn in ('pressure', 'temperature', 'energy', 'density'):
                node.xmlRemoveChild(nn)
            self.setThermoStatus('total_pressure', 'on')
            self.setThermoStatus('enthalpy', 'on')
        elif old_model != type and type != "subsonic_inlet_PH":
            node.xmlRemoveChild('total_pressure')
            node.xmlRemoveChild('enthalpy')

        n['choice'] = type


    @Variables.noUndo
    def getThermoValue(self, var):
        """
        Return value of the variable
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        value = n.xmlGetChildDouble(var, 'status')
        if value == None:
            value = self.__defaultValues()[var]
            self.setThermoValue(var, value)

        return value


    @Variables.undoLocal
    def setThermoValue(self, var, value):
        """
        Set value of the variable
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('velocity_pressure')
        node.xmlSetData(var, value)


    @Variables.noUndo
    def getCheckedBoxList(self):
        """
        Public method.
        return a list of selected variable for initialisation
        """
        box_list = []
        status = self.getThermoStatus('pressure')
        if status == 'on':
            box_list.append('Pressure')
        status = self.getThermoStatus('density')
        if status == 'on':
            box_list.append('Density')
        status = self.getThermoStatus('temperature')
        if status == 'on':
            box_list.append('Temperature')
        status = self.getThermoStatus('energy')
        if status == 'on':
            box_list.append('Energy')
        return box_list


    @Variables.noUndo
    def getListValue(self):
        """
        Public method.
        return a list of value for initialisation
        """
        value_list = []
        status = self.getThermoStatus('pressure')
        if status == 'on':
            value_list.append(self.getThermoValue('pressure'))
        status = self.getThermoStatus('density')
        if status == 'on':
            value_list.append(self.getThermoValue('density'))
        status = self.getThermoStatus('temperature')
        if status == 'on':
            value_list.append(self.getThermoValue('temperature'))
        status = self.getThermoStatus('energy')
        if status == 'on':
            value_list.append(self.getThermoValue('energy'))

        return value_list


    def deleteCompressibleInlet(self):
        """
        Delete all information of compressible model in boundary conditions.
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        n.xmlRemoveChild('compressible_type')
        n.xmlRemoveChild('density')
        n.xmlRemoveChild('pressure')
        n.xmlRemoveChild('temperature')
        n.xmlRemoveChild('energy')
        n.xmlRemoveChild('total_pressure')
        n.xmlRemoveChild('enthalpy')


    @Variables.noUndo
    def getInletGasCombustionType(self):
        """
        Return type for velocities's boundary conditions for inlet gas combustion.
        """
        node = self.boundNode.xmlGetNode('velocity_pressure')
        n = node.xmlInitNode('gas_type')
        type = n['choice']

        if type == None:
            type = self.__defaultValues()['gas_type']
            self.setInletGasCombustionType(type)

        return type


    @Variables.undoLocal
    def setInletGasCombustionType(self, type):
        """
        Set type for velocities's boundary conditions for inlet gas combustion.
        """
        Model().isInList(type, self.typeListGasComb)

        node = self.boundNode.xmlGetNode('velocity_pressure')
        n= node.xmlInitNode('gas_type')
        n['choice'] = type


    @Variables.noUndo
    def getGasCombustionTemperature(self):
        """
        Return value of the temperature for inlet gas combustion.
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        temperature = n.xmlGetChildDouble('temperature')
        if temperature == None:
            temperature = self.__defaultValues()['temperatureGas']
            self.setGasCombustionTemperature(temperature)

        return temperature


    @Variables.undoLocal
    def setGasCombustionTemperature(self, value):
        """
        Set value of the temperature for inlet gas combustion.
        """
        Model().isFloat(value)
        self.boundNode.xmlInitNode('velocity_pressure').xmlSetData('temperature',value)


    @Variables.noUndo
    def getMeanMixtureFraction(self):
        """
        Return value of the mean mixture fraction
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        fraction = n.xmlGetChildDouble('fraction')
        if fraction == None:
            fraction = self.__defaultValues()['fraction']
            self.setMeanMixtureFraction(fraction)
        return fraction


    @Variables.undoLocal
    def setMeanMixtureFraction(self, value):
        """
        Set value of the mean mixture fraction
        """
        Model().isGreaterOrEqual(value, 0.0)
        Model().isLowerOrEqual(value, 1.0)
        self.boundNode.xmlInitNode('velocity_pressure').xmlSetData('fraction',value)


    def deleteGas(self):
        """
        Delete all information of coal combustion in boundary conditions.
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        n.xmlRemoveChild('gas_type')
        n.xmlRemoveChild('fraction')
        n.xmlRemoveChild('temperature')


    @Variables.noUndo
    def getHydraulicHeadValue(self):
        """
        Return value of the hydraulic head
        """
        hydraulic_head = self.boundNode.xmlGetDouble('dirichlet', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head']
            self.setHydraulicHeadValue(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadValue(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('dirichlet', name='hydraulic_head')
        self.boundNode.xmlSetData('dirichlet', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadFlux(self):
        """
        Return value of the Hydraulic head flux
        """
        hydraulic_head = self.boundNode.xmlGetDouble('neumann', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head_flux']
            self.setHydraulicHeadFlux(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadFlux(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('neumann', name='hydraulic_head')
        self.boundNode.xmlSetData('neumann', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadChoice(self):
        """
        Get hydraulic head choice
        """
        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['hydraulicHeadChoice']
            self.setHydraulicHeadChoice(choice)
        return choice


    @Variables.undoGlobal
    def setHydraulicHeadChoice(self, choice) :
        """
        Set hydraulic head choice
        """
        Model().isInList(choice, ['dirichlet', 'neumann', 'dirichlet_formula'])

        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        if scalarNode['choice'] == choice:
            return

        old_choice = scalarNode['choice']
        scalarNode['choice'] = choice
        # delete nodes
        if old_choice == 'dirichlet':
            self.boundNode.xmlRemoveChild(old_choice)


    @Variables.noUndo
    def getHydraulicHeadFormula(self):
        """
        Public method.
        Return the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')
        formula = scalarNode.xmlGetChildString('formula')

        if not formula:
            formula = 'H = 0.;'

        return formula


    @Variables.undoLocal
    def setHydraulicHeadFormula(self, formula):
        """
        Public method.
        Set the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')

        n = scalarNode.xmlSetData('formula', formula)


    @Variables.noUndo
    def getConvectiveInletStatus(self):
        """
        Get the convective inlet status
        """
        status = None
        n = self.boundNode.xmlGetNode('convective_inlet')
        if n:
            status = n['status']
        if not status:
            status = 'off'
        return status


    @Variables.undoLocal
    def setConvectiveInletStatus(self, status):
        """
        Set the convective inlet status
        """
        Model().isInList(status, ('on', 'off'))

        n = self.boundNode.xmlInitNode('convective_inlet', status=status)
        if status == 'off':
            self.boundNode.xmlRemoveChild('convective_inlet')
        else:
            n['status'] = status


    @Variables.noUndo
    def getMappedInletStatus(self):
        """
        Get the mapped inlet status
        """
        status = None
        n = self.boundNode.xmlGetNode('mapped_inlet')
        if n:
            status = n['status']
        if not status:
            status = 'off'
        return status


    @Variables.undoLocal
    def setMappedInletStatus(self, status):
        """
        Set the mapped inlet status
        """
        Model().isInList(status, ('on', 'off'))

        n = self.boundNode.xmlInitNode('mapped_inlet', status=status)
        if status == 'off':
            self.boundNode.xmlRemoveChild('mapped_inlet')
        else:
            n['status'] = status


    @Variables.noUndo
    def getMappedInletTranslation(self, component):
        """
        Get the translation component
        """
        Model().isInList(component, self.__mappingDirectionTags)

        value = None
        n = self.boundNode.xmlInitNode('mapped_inlet', status='on')
        if n:
            value = n.xmlGetChildDouble(component)
        if not value:
            value = 0.0
        return value


    @Variables.undoLocal
    def setMappedInletTranslation(self, component, value):
        """
        Set the transation component
        """
        Model().isInList(component, self.__mappingDirectionTags)
        Model().isFloat(value)

        self.setMappedInletStatus('on')
        n = self.boundNode.xmlInitNode('mapped_inlet', status='on')
        n.xmlSetData(component, value)


#-------------------------------------------------------------------------------
# Atmospheric flow inlet/outlet boundary.
#-------------------------------------------------------------------------------

class MeteoBoundary(Boundary) :
    """
    Atmospheric flow inlet/outlet boundary.
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['meteo_data'] = 'off'
        dico['meteo_automatic'] = 'off'

        return dico


    @Variables.noUndo
    def getMeteoDataStatus(self):
        """
        Return if one reads the meteorological data.
        """
        node = self.boundNode.xmlInitNode('velocity_pressure').xmlInitNode('meteo_data')
        if node['status'] == None:
            self.setMeteoDataStatus (self.__defaultValues()['meteo_data'])
        return node['status']


    @Variables.undoLocal
    def setMeteoDataStatus(self, status):
        """
        """
        Model().isOnOff(status)
        self.boundNode.xmlInitNode('velocity_pressure').xmlInitNode('meteo_data')['status'] = status
        if status == 'on':
            for n in self.boundNode.xmlGetNodeList('scalar'):
                if n:
                    n.xmlRemoveNode()


    @Variables.noUndo
    def getAutomaticNatureStatus(self):
        """
        The boundary could be set to an inlet or an outlet automaticaly.
        """
        node = self.boundNode.xmlInitNode('velocity_pressure').xmlInitNode('meteo_automatic')
        if node['status'] == None:
            self.setAutomaticNatureStatus(self.__defaultValues()['meteo_automatic'])
        return node['status']


    @Variables.undoLocal
    def setAutomaticNatureStatus(self, status):
        """
        The boundary could be set to an inlet or an outlet automaticaly.
        """
        Model().isOnOff(status)
        self.boundNode.xmlInitNode('velocity_pressure').xmlInitNode('meteo_automatic')['status'] = status

#-------------------------------------------------------------------------------
# Joule effects inlet/outlet/wall boundary.
#-------------------------------------------------------------------------------

class JouleBoundary(Boundary) :
    """
    Joule effects inlet/outlet/wall boundary.
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['scalarChoice']  = 'dirichlet'
        dico['PotVectChoice'] = 'dirichlet_formula'
        dico['scalar']        = 0.

        return dico


    def __getscalarList(self):
        """
        return list of scalars name
        """
        scalar_list = []
        self.sca_model = DefineUserScalarsModel(self.case)
        for sca in self.sca_model.getElectricalScalarsNameList():
            scalar_list.append(sca)

        return scalar_list


    def scalarChoicesList(self):
        """
        Return choice list for boundary conditions
        """
        scalarChoicesList = ['dirichlet', 'neumann', 'exchange_coefficient',
                             'dirichlet_formula', 'neumann_formula', 'exchange_coefficient_formula',
                             'dirichlet_implicit', 'neumann_implicit']
        return scalarChoicesList


    def __deleteScalarNodes(self, name, tag):
        """
        Delete nodes of scalars
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(tag, self.scalarChoicesList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        for tt in self.scalarChoicesList():
            if tt != tag:
                scalarNode.xmlRemoveChild(tt)


    def getPotentialVectorChoice(self, name):
        """
        Get potential vector choice
        """
        Model().isInList(name, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['PotVectChoice']
            self.setPotentialVectorChoice(name, choice)

        return choice


    def setPotentialVectorChoice(self, name, choice):
        """
        Set potential vector choice
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, self.scalarChoicesList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        if choice == 'dirichlet_formula' or choice == 'neumann_formula':
            self.__deleteScalarNodes(name, 'dirichlet')
        else:
            self.__deleteScalarNodes(name, choice)


    def getElecScalarChoice(self, name):
        """
        Get scalar choice
        """
        Model().isInList(name, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarChoice']
            self.setElecScalarChoice(name, choice)

        return choice


    def setElecScalarChoice(self, name, choice) :
        """
        Set scalar choice
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, self.scalarChoicesList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        if choice == 'dirichlet_formula' or choice == 'neumann_formula':
            self.__deleteScalarNodes(name, 'dirichlet')
        else:
            self.__deleteScalarNodes(name, choice)


    def getElecScalarValue(self, scalarName, choice) :
        """
        Get scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self.scalarChoicesList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        value = scalarNode.xmlGetChildDouble(choice)
        if value == None :
            value = self.__defaultValues()['scalar']
            self.setElecScalarValue(scalarName, choice, value)
        return value


    def setElecScalarValue(self, scalarName, choice, value):
        """
        Set scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isFloat(value)
        Model().isInList(choice, self.scalarChoicesList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        scalarNode.xmlSetData(choice, value)


    def getElecScalarFormula(self, scalarName, choice):
        """
        Public method.
        Return the formula for a scalar variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)
        formula = scalarNode.xmlGetChildString(choice)
        return formula


    def setElecScalarFormula(self, scalarName, choice, formula):
        """
        Public method.
        Set the formula for a scalar variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        n = scalarNode.xmlSetData(choice, formula)

#-------------------------------------------------------------------------------
# free inlet outlet boundary
#-------------------------------------------------------------------------------

class InletOutletBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        pass


    @Variables.noUndo
    def getHeadLossesFormula(self):
        """
        Public method.
        Return the formula for an external head loss.
        """
        XMLHeadLossNode = self.boundNode.xmlInitNode('headLoss')

        formula = XMLHeadLossNode.xmlGetString('formula')
        return formula


    @Variables.undoLocal
    def setHeadLossesFormula(self, formula):
        """
        Public method.
        Set the formula for an external head loss.
        """
        XMLHeadLossNode = self.boundNode.xmlInitNode('headLoss')
        if not XMLHeadLossNode:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = XMLHeadLossNode.xmlInitChildNode('formula')
        n.xmlSetTextNode(formula)

#-------------------------------------------------------------------------------
# imposed pressure outlet boundary
#-------------------------------------------------------------------------------

class ImposedPressureOutletBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        pass


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['pressure']    = 101300.

        return dico


    @Variables.noUndo
    def getPressureValue(self):
        """
        Return value of the pressure
        """
        pressure = self.boundNode.xmlGetDouble('dirichlet', name='pressure')
        if pressure == None:
            pressure = self.__defaultValues()['pressure']
            self.setPressureValue(pressure)

        return pressure


    @Variables.undoLocal
    def setPressureValue(self, value):
        """
        Set value of the pressure
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('dirichlet', name='pressure')
        self.boundNode.xmlSetData('dirichlet', value, name='pressure')


#-------------------------------------------------------------------------------
# Coal flow inlet boundary
#-------------------------------------------------------------------------------

class CoalInletBoundary(InletBoundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        InletBoundary._initBoundary(self)
        self.typeList = ['oxydantFlow', 'coalFlow']

        self.__updateCoalInfo()

        # Initialize nodes if necessary
        tpe = self.getInletType()
        self.setInletType(tpe)


    def __updateCoalInfo(self):
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel = CoalCombustionModel(self.case)
        self.coalNumber = CoalCombustionModel.getCoalNumber()
        log.debug("__updateCoalInfo coalNumber: %i " % self.coalNumber)
        self.coalClassesNumber = []
        for c in range(0, self.coalNumber):
            self.coalClassesNumber.append(CoalCombustionModel.getClassNumber(str(c+1)))
            log.debug("__updateCoalInfo number of classes: %i " % self.coalClassesNumber[c])


    def __deleteCoalNodes(self):
        """
        Delete all nodes udes for coal. Private method
        """
        node = self.boundNode.xmlGetNode('velocity_pressure')
        for n in node.xmlGetChildNodeList('coal'):
            n.xmlRemoveNode()


    def __getClassCoalRatio(self, coal, classe):
        """
        Return ratio for classe for coal. Private method
        """
        Model().isInt(coal)
        Model().isInt(classe)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        num = '%2.2i' % (coal+1)
        nc = n.xmlInitNode('coal', name="coal"+num)
        num = '%2.2i' % (classe+1)
        nratio = nc.xmlGetChildNode('ratio', name="class"+num)
        if nratio:
            ratio = nc.xmlGetChildDouble('ratio', name="class"+num)
        else:
            #self.__updateCoalInfo()
            if self.coalClassesNumber[coal] > 1:
                if classe == 0:
                    ratio = 100.
                    self.__setClassCoalRatio(ratio, coal, classe)
                else:
                    ratio = self.__defaultValues()['ratio']
                    self.__setClassCoalRatio(ratio, coal, classe)
            else:
                ratio = 100.
                self.__setClassCoalRatio(ratio, coal, classe)

        return ratio


    def __setClassCoalRatio(self, value, coal, classe):
        """
        Put value of ratio when several classes for coal. Private method
        """
        Model().isFloat(value)
        Model().isLowerOrEqual(value, 100.)
        Model().isInt(coal)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        num = '%2.2i' % (coal+1)
        nc = n.xmlInitNode('coal', name="coal"+ num)

        num = '%2.2i' % (classe+1)
        nc.xmlSetData('ratio', value, name="class"+ num)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['flow'] = 1.0
        dico['ratio'] = 0.0
        dico['oxydant'] = 1
        from code_saturne.Pages.ReferenceValuesModel import ReferenceValuesModel
        dico['temperature'] = ReferenceValuesModel(self.case).getTemperature()

        return dico


    @Variables.noUndo
    def getInletType(self):
        """
        Return type (oxydant or oxydant+coal) for velocities's boundary conditions for inlet coal flow.
        """
        if self.boundNode.xmlGetNode('velocity_pressure').xmlGetChildNodeList('coal'):
            type = "coalFlow"
        else:
            type = "oxydantFlow"
        return type


    @Variables.undoGlobal
    def setInletType(self, type):
        """
        Set type (oxydant or oxydant+coal) for velocities's boundary conditions for inlet coal flow.
        """
        Model().isInList(type, self.typeList)

        self.getOxydantTemperature()

        if type == "oxydantFlow":
            self.__deleteCoalNodes()
        elif type == "coalFlow":
            for coal_idx in range(0, self.coalNumber):
                self.getCoalFlow(coal_idx)
                self.getCoalTemperature(coal_idx)
                self.getCoalRatios(coal_idx)


    @Variables.noUndo
    def getCoalFlow(self, coal_idx):
        """
        Return value of flow for coal
        """
        Model().isInt(coal_idx)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        num = '%2.2i' % (coal_idx+1)
        n2 = n.xmlInitNode('coal', name = "coal"+num)
        flow = n2.xmlGetDouble('flow1')
        if flow == None:
            flow = self.__defaultValues()['flow']
            self.setCoalFlow(flow, coal_idx)

        return flow


    @Variables.undoLocal
    def setCoalFlow(self, value, coal):
        """
        Put value of flow for coal
        """
        Model().isFloat(value)
        Model().isInt(coal)

        num = '%2.2i' % (coal+1)
        n = self.boundNode.xmlGetNode('velocity_pressure')
        n.xmlInitNode('coal', name = "coal"+num).xmlSetData('flow1', value)


    @Variables.noUndo
    def getOxydantTemperature(self):
        """
        Return value of the temperature for oxydant for coal choice
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        temperature = n.xmlGetChildDouble('temperature')
        if temperature == None:
            temperature = self.__defaultValues()['temperature']
            self.setOxydantTemperature(temperature)

        return temperature


    @Variables.undoLocal
    def setOxydantNumber(self, value):
        """
        Set value of the oxydant number.
        """
        Model().isInt(value)
        self.boundNode.xmlInitNode('velocity_pressure').xmlSetData('oxydant',value)


    @Variables.noUndo
    def getOxydantNumber(self):
        """
        Return value of oxydant number.
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        oxydant = n.xmlGetInt('oxydant')
        if oxydant == None:
            oxydant = self.__defaultValues()['oxydant']
            self.setOxydantNumber(oxydant)

        return oxydant


    @Variables.undoLocal
    def setOxydantTemperature(self, value):
        """
        Set value of the temperature for oxydant for coal choice
        """
        Model().isFloat(value)
        self.boundNode.xmlInitNode('velocity_pressure').xmlSetData('temperature',value)


    @Variables.noUndo
    def getCoalTemperature(self, coal):
        """
        Return value of temperature for coal for coal choice
        """
        Model().isInt(coal)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        num = '%2.2i' % (coal+1)
        nt = n.xmlInitNode('coal', name="coal"+ num)
        temperature = nt.xmlGetChildDouble('temperature')
        if temperature == None:
            temperature = self.__defaultValues()['temperature']
            self.setCoalTemperature(temperature, coal)

        return temperature


    @Variables.undoLocal
    def setCoalTemperature(self, value, coal_idx):
        """
        Put value of temperature for coal for coal choice
        """
        Model().isFloat(value)
        Model().isInt(coal_idx)

        num = '%2.2i' % (coal_idx+1)
        n = self.boundNode.xmlInitNode('velocity_pressure')
        n.xmlInitNode('coal', name="coal"+ num).xmlSetData('temperature',value)


    @Variables.noUndo
    def getCoalRatios(self, coal_idx):
        """
        Put list of values of classe's ratio for one coal
        """
        Model().isInt(coal_idx)
        Model().isLowerOrEqual(coal_idx, self.coalNumber-1)

        lst = []
        for i in range(0, self.coalClassesNumber[coal_idx]):
            lst.append(self.__getClassCoalRatio(coal_idx, i))

        # check if sum of ratios of coal mass is equal to 100%
        som = 0.
        for i in range(0, self.coalClassesNumber[coal_idx]):
            som += lst[i]
        Model().isFloatEqual(som, 100.)

        return lst


    @Variables.undoLocal
    def setCoalRatios(self, coal, lst):
        """
        Put list of values of classe's ratio for one coal
        """
        Model().isInt(coal)
        Model().isIntEqual(len(lst), self.coalClassesNumber[coal])
        som = 0.
        for i in range(0, self.coalClassesNumber[coal]):
            som += lst[i]
        Model().isFloatEqual(som, 100.)

        n = self.boundNode.xmlInitNode('velocity_pressure')
        num = '%2.2i' % (coal+1)
        nc = n.xmlInitNode('coal', name="coal"+ num)

        for i in range(0, len(lst)):
            num = '%2.2i' % (i+1)
            nc.xmlSetData('ratio', lst[i], name="class"+ num)


    def deleteCoalFlow(self, coal, nbCoals):
        """
        Delete coal with name = coal.
        Usefull only for CoalCombustionView.
        """
        Model().isInt(coal)
        Model().isInt(nbCoals)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        if n:
            num = '%2.2i' % (coal+1)
            n2 = n.xmlGetNode('coal', name="coal"+ num)
            # delete coal
            if n2:
                n2.xmlRemoveNode()
                # rename other coals
                for icoal in range(coal+1, nbCoals):
                    self.__renameCoalFlow(icoal)


    def __renameCoalFlow(self, coal):
        """
        coaln become coaln-1.
        Usefull only for CoalCombustionView.
       """
        Model().isInt(coal)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        if n:
            num = '%2.2i' % (coal+1)
            newNum = '%2.2i' % coal
            n2 = n.xmlGetNode('coal', name="coal"+ num)
            if n2:
                n2['name'] = "coal"+newNum


    def updateCoalRatios(self, coal):
        """
        Delete ratio with name = classe. Usefull only for CoalCombustionView.
        """
        Model().isInt(coal)

        n = self.boundNode.xmlGetNode('velocity_pressure')
        if n:
            num = '%2.2i' % (coal)
            n2 = n.xmlGetNode('coal', name="coal"+num)
            if n2:
                n2.xmlRemoveChild('ratio')

            coal_idx = coal-1
            self.getCoalRatios(coal_idx)


    def deleteCoals(self):
        """
        Delete all information of coal combustion in boundary conditions.
        """
        n = self.boundNode.xmlGetNode('velocity_pressure')
        n.xmlRemoveChild('oxydant')
        n.xmlRemoveChild('coal')
        n.xmlRemoveChild('temperature')

#-------------------------------------------------------------------------------
# Compressible flow outlet boundary
#-------------------------------------------------------------------------------

class CompressibleOutletBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self.typeList = ['supersonic_outlet', 'subsonic_outlet']

        # Initialize nodes if necessary
        type = self.getOutletType()
        self.setOutletType(type)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['pressure']    = 101300.
        dico['compressible_type']   = 'supersonic_outlet'

        return dico


    @Variables.noUndo
    def getOutletType(self):
        """
        Return type of boundary conditions for outlet
        """
        node = self.boundNode.xmlInitNode('compressible_type')
        type = node['choice']
        if type == None:
            type = self.__defaultValues()['compressible_type']
            self.setOutletType(type)
        return type


    @Variables.undoLocal
    def setOutletType(self, type):
        """
        Set type of boundary conditions for outlet
        """
        Model().isInList(type, self.typeList)

        node = self.boundNode.xmlInitNode('compressible_type')
        node['choice'] = type
        if type == 'supersonic_outlet':
            self.boundNode.xmlRemoveChild('dirichlet', name='pressure')


    @Variables.noUndo
    def getPressureValue(self):
        """
        Return value of the pressure
        """
        pressure = self.boundNode.xmlGetDouble('dirichlet', name='pressure')
        if pressure == None:
            pressure = self.__defaultValues()['pressure']
            self.setPressureValue(pressure)

        return pressure


    @Variables.undoLocal
    def setPressureValue(self, value):
        """
        Set value of the pressure
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('dirichlet', name='pressure')
        self.boundNode.xmlSetData('dirichlet', value, name='pressure')


    def deleteCompressibleOutlet(self):
        """
        Delete all information of compressible model in boundary conditions.
        """
        n = self.boundNode
        n.xmlRemoveChild('compressible_type')
        n.xmlRemoveChild('pressure')


#-------------------------------------------------------------------------------
# Outlet boundary
#-------------------------------------------------------------------------------

class OutletBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self._scalarChoicesList = ['dirichlet', 'neumann',
                                   'dirichlet_formula', 'neumann_formula']
        self.getReferencePressure()

        for name in self.sca_model.getScalarNameList():
            if self.getScalarChoice(name) == 'dirichlet' or \
               self.getScalarChoice(name) == 'neumann':
                self.getScalarValue(name, self.getScalarChoice(name))


    def __getscalarList(self):
        """
        return list of scalars name
        """
        scalar_list =  self.sca_model.getScalarNameList()
        if self.sca_model.getMeteoScalarsNameList() != None:
            for sca in self.sca_model.getMeteoScalarsNameList():
                scalar_list.append(sca)
        if len(self.sca_model.getThermalScalarName()) > 0:
            scalar_list.append(self.sca_model.getThermalScalarName()[0])

        return scalar_list


    def __deleteScalarNodes(self, name, tag):
        """
        Delete nodes of scalars
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(tag, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        for tt in self._scalarChoicesList:
            if tt != tag:
                scalarNode.xmlRemoveChild(tt)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        from code_saturne.Pages.ReferenceValuesModel import ReferenceValuesModel
        dico['reference_pressure'] = ReferenceValuesModel(self.case).getPressure()
        dico['scalarChoice']  = 'neumann'
        dico['scalar']        = 0.
        dico['hydraulic_head']      = 0.
        dico['hydraulic_head_flux'] = 0.
        dico['hydraulicHeadChoice'] = 'dirichlet'

        return dico


    @Variables.noUndo
    def getScalarChoice(self, name):
        """
        Get scalar choice
        """
        Model().isInList(name, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarChoice']
            self.setScalarChoice(name, choice)

        return choice


    @Variables.undoLocal
    def setScalarChoice(self, name, choice) :
        """
        Set scalar choice
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        if choice == 'dirichlet_formula' or choice == 'neumann_formula':
            self.__deleteScalarNodes(name, 'dirichlet')
        else:
            self.getScalarValue(name, choice)
            self.__deleteScalarNodes(name, choice)


    @Variables.noUndo
    def getScalarValue(self, name, choice) :
        """
        Get variableName variable
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = self.getScalarChoice(name)

        value = scalarNode.xmlGetChildDouble(choice)
        if value == None :
            value = self.__defaultValues()['scalar']
            self.setScalarValue(name, choice, value)
        return value


    @Variables.undoLocal
    def setScalarValue(self, name, choice, value) :
        """
        Set variableName variable
        """
        Model().isInList(name, self.__getscalarList())
        Model().isFloat(value)
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = self.getScalarChoice(name)
        scalarNode.xmlSetData(choice, value)


    @Variables.undoLocal
    def setScalarFormula(self, name, choice, formula):
        """
        Public method.
        Set the formula for a scalar variable.
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        n = scalarNode.xmlSetData(choice, formula)


    @Variables.noUndo
    def getDefaultScalarFormula(self, scalarName, scalar_model):
        """
        Get defaut scalar formula
        """
        if scalar_model == 'dirichlet_formula':
            formula = scalarName+""" = 0;\n"""
        elif scalar_model == 'neumann_formula':
            formula = """flux = 0;\n"""
        elif scalar_model == 'exchange_coefficient_formula':
            formula = scalarName+""" = 0;\nhc = 0;\n"""

        return formula


    @Variables.noUndo
    def getScalarFormula(self, name, choice):
        """
        Public method.
        Return the formula for a scalar variable.
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        formula = scalarNode.xmlGetChildString(choice)

        if not formula:
            formula = self.getDefaultScalarFormula(name, choice)

        return formula


    @Variables.noUndo
    def getPressureChoice(self) :
        """
        Return if the value of pressure exist or not of boundary conditions for outlet.
        """
        choice = "off"
        if self.boundNode.xmlGetChildNode('dirichlet', name='reference_pressure'):
            choice = "on"

        return choice


    @Variables.undoGlobal
    def setPressureChoice(self, choice) :
        """
        Set balise of pressure beyond the choice for boundary conditions for outlet
        """
        Model().isOnOff(choice)
        if choice == 'off':
            self.setReferencePressure(self, 0.0)
        else:
            if node.xmlGetDouble('dirichlet', name='pressure') == None:
                self.setReferencePressure(self.__defaultValues()['pressure'])


    @Variables.noUndo
    def getReferencePressure(self) :
        """
        Get reference pressure
        """
        pressure = self.boundNode.xmlGetDouble('dirichlet', name='pressure')
        if pressure == None:
            return 0

        return pressure


    @Variables.undoLocal
    def setReferencePressure(self, value) :
        """
        Set reference pressure
        """
        Model().isPositiveFloat(value)

        node = self.boundNode.xmlInitNode('dirichlet', name='pressure')
        if value == 0:
            node.xmlRemoveNode()
        else:
            self.boundNode.xmlSetData('dirichlet', value, name='pressure')


    @Variables.noUndo
    def getHydraulicHeadValue(self):
        """
        Return value of the hydraulic head
        """
        hydraulic_head = self.boundNode.xmlGetDouble('dirichlet', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head']
            self.setHydraulicHeadValue(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadValue(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('dirichlet', name='hydraulic_head')
        self.boundNode.xmlSetData('dirichlet', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadFlux(self):
        """
        Return value of the hydraulic head
        """
        hydraulic_head = self.boundNode.xmlGetDouble('neumann', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head_flux']
            self.setHydraulicHeadFlux(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadFlux(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('neumann', name='hydraulic_head')
        self.boundNode.xmlSetData('neumann', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadChoice(self):
        """
        Get hydraulic head choice
        """
        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['hydraulicHeadChoice']
            self.setHydraulicHeadChoice(choice)
        return choice


    @Variables.undoGlobal
    def setHydraulicHeadChoice(self, choice) :
        """
        Set hydraulic head choice
        """
        Model().isInList(choice, ['dirichlet', 'neumann', 'dirichlet_formula'])

        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        if scalarNode['choice'] == choice:
            return

        old_choice = scalarNode['choice']
        scalarNode['choice'] = choice
        # delete nodes
        if old_choice == 'dirichlet':
            self.boundNode.xmlRemoveChild(old_choice)


    @Variables.noUndo
    def getHydraulicHeadFormula(self):
        """
        Public method.
        Return the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')
        formula = scalarNode.xmlGetChildString('formula')

        if not formula:
            formula = 'H = 0.;'

        return formula


    @Variables.undoLocal
    def setHydraulicHeadFormula(self, formula):
        """
        Public method.
        Set the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')

        n = scalarNode.xmlSetData('formula', formula)


#-------------------------------------------------------------------------------
# Symmetry boundary
#-------------------------------------------------------------------------------

class SymmetryBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        pass


#-------------------------------------------------------------------------------
# free surface boundary
#-------------------------------------------------------------------------------

class FreeSurfaceBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        pass


    def deleteFreeSurface(self):
        """
        Delete all information of free surface in boundary conditions.
        """
        self.boundNode.xmlRemoveNode()


#-------------------------------------------------------------------------------
# GroundwaterBoundary boundary
#-------------------------------------------------------------------------------

class GroundwaterBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self._scalarChoicesList  = ['dirichlet',
                                    'neumann',
                                    'exchange_coefficient',
                                    'dirichlet_formula',
                                    'neumann_formula',
                                    'exchange_coefficient_formula']
        pass


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['hydraulic_head'] = 0.
        dico['hydraulic_head_flux'] = 0.
        dico['hydraulicHeadChoice'] = 'dirichlet'
        dico['scalar']              = 0.0
        dico['scalarChoice']        = 'dirichlet'

        return dico


    @Variables.noUndo
    def getHydraulicHeadValue(self):
        """
        Return value of the hydraulic head
        """
        hydraulic_head = self.boundNode.xmlGetDouble('dirichlet', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head']
            self.setHydraulicHeadValue(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadValue(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('dirichlet', name='hydraulic_head')
        self.boundNode.xmlSetData('dirichlet', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadFlux(self):
        """
        Return value of the hydraulic head
        """
        hydraulic_head = self.boundNode.xmlGetDouble('neumann', name='hydraulic_head')
        if hydraulic_head == None:
            hydraulic_head = self.__defaultValues()['hydraulic_head_flux']
            self.setHydraulicHeadFlux(hydraulic_head)

        return hydraulic_head


    @Variables.undoLocal
    def setHydraulicHeadFlux(self, value):
        """
        Set value of the hydraulic head
        """
        Model().isFloat(value)
        node = self.boundNode.xmlInitNode('neumann', name='hydraulic_head')
        self.boundNode.xmlSetData('neumann', value, name='hydraulic_head')


    @Variables.noUndo
    def getHydraulicHeadChoice(self):
        """
        Get hydraulic head choice
        """
        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['hydraulicHeadChoice']
            self.setHydraulicHeadChoice(choice)
        return choice


    @Variables.undoGlobal
    def setHydraulicHeadChoice(self, choice) :
        """
        Set hydraulic head choice
        """
        Model().isInList(choice, ['dirichlet', 'neumann', 'dirichlet_formula'])

        scalarNode = self.boundNode.xmlInitNode('hydraulicHead')

        if scalarNode['choice'] == choice:
            return

        old_choice = scalarNode['choice']
        scalarNode['choice'] = choice
        # delete nodes
        if old_choice == 'dirichlet':
            self.boundNode.xmlRemoveChild(old_choice)


    @Variables.noUndo
    def getHydraulicHeadFormula(self):
        """
        Public method.
        Return the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')
        formula = scalarNode.xmlGetChildString('formula')

        if not formula:
            formula = 'H = 0.;'

        return formula


    @Variables.undoLocal
    def setHydraulicHeadFormula(self, formula):
        """
        Public method.
        Set the formula for hydraulic head variable.
        """
        scalarNode = self.boundNode.xmlInitNode('dirichlet_formula', name='hydraulicHead')

        n = scalarNode.xmlSetData('formula', formula)


    def __getscalarList(self):
        """
        return list of scalars name
        """
        scalar_list =  self.sca_model.getScalarNameList()
        if self.sca_model.getMeteoScalarsNameList() != None:
            for sca in self.sca_model.getMeteoScalarsNameList():
                scalar_list.append(sca)
        if len(self.sca_model.getThermalScalarName()) > 0:
            scalar_list.append(self.sca_model.getThermalScalarName()[0])

        return scalar_list


    def __deleteScalarNodes(self, name, tag):
        """
        Delete nodes of scalars
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(tag, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        for tt in self._scalarChoicesList:
            if tt != tag:
                scalarNode.xmlRemoveChild(tt)


    @Variables.noUndo
    def getScalarChoice(self, scalarName):
        """
        Get scalar choice
        """
        Model().isInList(scalarName, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarChoice']
            self.setScalarChoice(scalarName, choice)

        return choice


    @Variables.undoGlobal
    def setScalarChoice(self, scalarName, choice) :
        """
        Set scalar choice
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        self.__deleteScalarNodes(scalarName, choice)


    @Variables.noUndo
    def getScalarValue(self, scalarName, choice) :
        """
        Get scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        value = scalarNode.xmlGetChildDouble(choice)
        if value == None :
            value = self.__defaultValues()['scalar']
            self.setScalarValue(scalarName, choice, value)
        return value


    @Variables.undoGlobal
    def setScalarValue(self, scalarName, choice, value):
        """
        Set scalar value
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isFloat(value)
        Model().isInList(choice, self._scalarChoicesList)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, scalarName)

        scalarNode.xmlSetData(choice, value)


    @Variables.noUndo
    def getDefaultScalarFormula(self, scalarName, scalar_model):
        """
        Get defaut scalar formula
        """
        if scalar_model == 'dirichlet_formula':
            formula = scalarName+""" = 0;\n"""
        elif scalar_model == 'neumann_formula':
            formula = """flux = 0;\n"""
        elif scalar_model == 'exchange_coefficient_formula':
            formula = scalarName+""" = 0;\nhc = 0;\n"""

        return formula


    @Variables.noUndo
    def getScalarFormula(self, scalarName, choice):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)
        formula = scalarNode.xmlGetChildString(choice)

        if not formula:
            formula = self.getDefaultScalarFormula(scalarName, choice)

        return formula


    @Variables.undoLocal
    def setScalarFormula(self, scalarName, choice, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        Model().isInList(scalarName, self.__getscalarList())
        Model().isInList(choice, self._scalarChoicesList)
        scalarNode = self.boundNode.xmlInitNode('scalar', name=scalarName)

        n = scalarNode.xmlSetData(choice, formula)


#-------------------------------------------------------------------------------
# Wall boundary
#-------------------------------------------------------------------------------

class WallBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self._fluxChoices=['temperature', 'flux']
        self._scalarChoices = ['dirichlet', 'neumann', 'exchange_coefficient',
                               'dirichlet_formula', 'neumann_formula',
                               'exchange_coefficient_formula']
        self._scalarFormula = ['dirichlet_formula', 'neumann_formula',
                               'exchange_coefficient_formula']

        # Initialize nodes if necessary
        self.getVelocityChoice()
        self.getRoughnessChoice()

        # Scalars
        for name in self.sca_model.getScalarNameList():
            self.getScalarChoice(name)


    def __deleteVelocities(self, node):
        """
        Delete nodes of velocity
        """
        node.xmlRemoveChild('dirichlet', name='velocity', component="0")
        node.xmlRemoveChild('dirichlet', name='velocity', component="1")
        node.xmlRemoveChild('dirichlet', name='velocity', component="2")


    def __getscalarList(self):
        """
        return list of scalars name
        """
        scalar_list =  self.sca_model.getScalarNameList()
        if self.sca_model.getMeteoScalarsNameList() != None:
            for sca in self.sca_model.getMeteoScalarsNameList():
                scalar_list.append(sca)
        if len(self.sca_model.getThermalScalarName()) > 0:
            scalar_list.append(self.sca_model.getThermalScalarName()[0])

        return scalar_list


    def __deleteScalarNodes(self, name, tag):
        """
        Delete nodes of scalars
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(tag, self._scalarChoices)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        if tag == 'exchange_coefficient':
            for tt in self._scalarChoices:
                scalarNode.xmlRemoveChild(tt)
        else:
            for tt in self._scalarChoices:
                if tt != tag:
                    scalarNode.xmlRemoveChild(tt)


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['velocityChoice'] = "off"
        dico['velocityValue']  = 0.
        dico['scalarChoice']   = "neumann"
        dico['scalarValue']   = 0.
        dico['roughness'] = 0.01
        dico['flux'] = 0
        return dico


    @Variables.noUndo
    def getVelocityChoice(self):
        """
        Get the velocity choice
        """
        node = self.boundNode.xmlInitNode('velocity_pressure', 'choice')
        choice = node['choice']
        if not choice:
            choice = self.__defaultValues()['velocityChoice']
            self.setVelocityChoice(choice)
        return node['choice']


    @Variables.undoGlobal
    def setVelocityChoice(self, choice):
        """
        Set the velocity choice
        """
        Model().isOnOff(choice)
        #
        # Check if value is a new velocity choice value
        XMLVelocityNode = self.boundNode.xmlInitNode('velocity_pressure')

        if XMLVelocityNode['choice'] == choice:
            return

        # Update velocity choice
        XMLVelocityNode['choice'] = choice
        #
        # Velocities updating
        if choice == 'on':
            # Create norm node if necessary
            self.getVelocities()
        else:
            # Delete 'flow1', 'flow2' and 'direction' nodes
            self.__deleteVelocities(XMLVelocityNode)


    @Variables.noUndo
    def getVelocities(self):
        """
        Set the velocity definition according to choice
        """
        node = self.boundNode.xmlInitNode('velocity_pressure')
        Model().isInList(node['choice'],('on',))

        n = node.xmlGetChildNode('dirichlet', name='velocity', component="0")
        if n:
            u = node.xmlGetChildDouble('dirichlet', name='velocity', component="0")
        else:
            u = self.__defaultValues()['velocityValue']
        n = node.xmlGetChildNode('dirichlet', name='velocity', component="1")
        if n:
            v = node.xmlGetChildDouble('dirichlet', name='velocity', component="1")
        else:
            v = self.__defaultValues()['velocityValue']
        n = node.xmlGetChildNode('dirichlet', name='velocity', component="2")
        if n:
            w = node.xmlGetChildDouble('dirichlet', name='velocity', component="2")
        else:
            w = self.__defaultValues()['velocityValue']
        self.setVelocities(u, v, w)

        return u, v, w


    @Variables.undoLocal
    def setVelocities(self, u, v, w):
        """
        Set the velocity definition according to choice
        """
        Model().isFloat(u)
        Model().isFloat(v)
        Model().isFloat(w)

        node = self.boundNode.xmlInitNode('velocity_pressure')
        Model().isInList(node['choice'],('on',))

        node.xmlSetData('dirichlet', u, name='velocity', component="0")
        node.xmlSetData('dirichlet', v, name='velocity', component="1")
        node.xmlSetData('dirichlet', w, name='velocity', component="2")


    @Variables.undoLocal
    def setVelocityComponent(self, val, component):
        """
        Set the value of component of the velocity - Method for the view
        """
        Model().isFloat(val)
        Model().isInList(component, ('0', '1', '2'))
        node = self.boundNode.xmlInitNode('velocity_pressure')
        Model().isInList(node['choice'], ('on', 'off'))

        node.xmlSetData('dirichlet', val, name="velocity", component=component)


    @Variables.noUndo
    def getRoughnessChoice(self):
        """
        Return if the value of roughness height exist or not of boundary conditions for wall.
        """
        choice = "off"
        node = self.boundNode.xmlInitNode('velocity_pressure')
        if node.xmlGetChildNode('roughness'):
            choice = "on"

        return choice


    @Variables.undoGlobal
    def setRoughnessChoice(self, choice):
        """
        Update balise of roughness beyond the choice for boundary conditions for wall.
        """
        Model().isOnOff(choice)
        node = self.boundNode.xmlInitNode('velocity_pressure')
        if choice == 'off':
            self.setRoughness(0.0)
        else:
            if node.xmlGetDouble('roughness') == None:
                self.setRoughness(self.__defaultValues()['roughness'])


    @Variables.noUndo
    def getRoughness(self):
        """
        Get the value of roughness height if it's exist of boundary conditions for wall.
        """
        node = self.boundNode.xmlInitNode('velocity_pressure')

        val = node.xmlGetDouble('roughness')
        if val == None:
            return 0

        return val


    @Variables.undoLocal
    def setRoughness(self, value):
        """
        Put value of roughness height in xmlfile
        """
        Model().isGreaterOrEqual(value, 0.)

        node = self.boundNode.xmlInitNode('velocity_pressure', 'choice')
        if value == 0.:
            node.xmlRemoveChild('roughness')
        else:
            node.xmlSetData('roughness', value)


    @Variables.noUndo
    def getScalarChoice(self, name):
        """
        Get scalar choice
        """
        Model().isInList(name, self.__getscalarList())

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update type and name of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        choice = scalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarChoice']
            self.setScalarChoice(name, choice)

        return choice


    @Variables.undoGlobal
    def setScalarChoice(self, name, choice) :
        """
        Set scalar choice
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, self._scalarChoices)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        if scalarNode['choice'] == choice:
            return

        scalarNode['choice'] = choice
        if choice in self._scalarFormula:
            self.__deleteScalarNodes(name, choice)
        elif choice == 'exchange_coefficient':
            self.getScalarValue(name, 'dirichlet')
            self.getScalarValue(name, 'exchange_coefficient')
            self.__deleteScalarNodes(name, choice)
        else:
            value = self.getScalarValue(name, choice)
            self.__deleteScalarNodes(name, choice)


    @Variables.noUndo
    def getScalarValue(self, name, choice) :
        """
        Get variableName variable
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(self.getScalarChoice(name), ('dirichlet', 'neumann', 'exchange_coefficient'))

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        value = scalarNode.xmlGetChildDouble(choice)
        if value == None :
            value = self.__defaultValues()['scalarValue']
            self.setScalarValue(name, choice, value)
        return value


    @Variables.undoLocal
    def setScalarValue(self, name, choice, value) :
        """
        Set variableName variable
        """
        Model().isInList(name, self.__getscalarList())
        Model().isFloat(value)
        Model().isInList(choice, self._scalarChoices)

        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        #update name and type of scalar
        self.updateScalarTypeAndName(scalarNode, name)

        scalarNode.xmlSetData(choice, value)


    @Variables.undoLocal
    def setScalarFormula(self, name, choice, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula', 'exchange_coefficient_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)

        scalarNode.xmlSetData(choice, formula)


    @Variables.noUndo
    def getDefaultScalarFormula(self, scalarName, scalar_model):
        """
        Get defaut scalar formula
        """
        if scalar_model == 'dirichlet_formula':
            formula = scalarName+""" = 0;\n"""
        elif scalar_model == 'neumann_formula':
            formula = """flux = 0;\n"""
        elif scalar_model == 'exchange_coefficient_formula':
            formula = scalarName+""" = 0;\nhc = 0;\n"""

        return formula


    @Variables.noUndo
    def getScalarFormula(self, name, choice):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        Model().isInList(name, self.__getscalarList())
        Model().isInList(choice, ('dirichlet_formula', 'neumann_formula', 'exchange_coefficient_formula'))
        scalarNode = self.boundNode.xmlInitNode('scalar', name=name)
        formula = scalarNode.xmlGetChildString(choice)

        if not formula:
            formula = self.getDefaultScalarFormula(name, choice)

        return formula


#-------------------------------------------------------------------------------
# Radiative wall boundary
#-------------------------------------------------------------------------------

class RadiativeWallBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the radiative boundary node
        """
        self._radiativeChoices = ['itpimp', 'ipgrno', 'ifgrno']

        self.head_list = ['emissivity',
                          'wall_thermal_conductivity',
                          'thickness', 'flux',
                          'external_temperature_profile',
                          'internal_temperature_profile']


    def _getListValRay(self, choice):
        """
        Return list of radiative variables according to choice's
        """
        Model().isInList(choice, self._radiativeChoices)
        lst = []
        if choice == 'itpimp':
            lst = ('emissivity', 'internal_temperature_profile')
        elif choice == 'ipgrno':
            lst = ('emissivity', 'wall_thermal_conductivity', 'thickness',
                    'external_temperature_profile',
                    'internal_temperature_profile')
        elif choice == 'ifgrno':
            lst = ('emissivity', 'flux', 'internal_temperature_profile')

        return lst


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['emissivity'] = 0.8
        dico['wall_thermal_conductivity'] = 3.0
        dico['thickness'] = 0.10
        dico['flux'] = 0.
        dico['external_temperature_profile'] = 300.
        dico['internal_temperature_profile'] = 300.
        dico['choice_condition'] = 'itpimp'
        return dico


    @Variables.noUndo
    def getRadiativeChoice(self):
        """
        Return variables according to choice of type of condition for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        choice = nod_ray_cond['choice']
        if not choice:
            choice = self.__defaultValues()['choice_condition']
            self.setRadiativeChoice(choice)
        return choice


    @Variables.undoLocal
    def setRadiativeChoice(self, choice):
        """
        Put variables according to choice of type of condition for the radiative wall
        """
        Model().isInList(choice, self._radiativeChoices)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond['choice'] = choice
        lst = self._getListValRay(choice)
        for i in lst:
            if not nod_ray_cond.xmlGetChildNode(i):
                nod_ray_cond.xmlSetData(i, self.__defaultValues()[i])


    @Variables.noUndo
    def getEmissivity(self):
        """
        Return value of emissivity for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        val = nod_ray_cond.xmlGetChildDouble('emissivity')
        if not val:
            val = self.__defaultValues()['emissivity']
            self.setEmissivity(val)

        return val


    @Variables.undoLocal
    def setEmissivity(self, val):
        """
        Put value of emissivity for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)
        Model().isLowerOrEqual(val, 1.)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData('emissivity', val)


    @Variables.noUndo
    def getThermalConductivity(self):
        """
        Return value of thermal conductivity for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        val = nod_ray_cond.xmlGetChildDouble('wall_thermal_conductivity')
        if not val:
            val = self.__defaultValues()['wall_thermal_conductivity']
            self.setThermalConductivity(val)

        return val


    @Variables.undoLocal
    def setThermalConductivity(self, val):
        """
        Put value of thermal conductivity for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData('wall_thermal_conductivity', val)


    @Variables.noUndo
    def getThickness(self):
        """
        Return value of thickness for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        val = nod_ray_cond.xmlGetChildDouble('thickness')
        if not val:
            val = self.__defaultValues()['thickness']
            self.setThickness(val)

        return val


    @Variables.undoLocal
    def setThickness(self, val):
        """
        Put value of thickness for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData('thickness', val)


    @Variables.noUndo
    def getExternalTemperatureProfile(self):
        """
        Return value of external temperature profile for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        val = nod_ray_cond.xmlGetChildDouble('external_temperature_profile')
        if not val:
            val = self.__defaultValues()['external_temperature_profile']
            self.setExternalTemperatureProfile(val)

        return val


    @Variables.undoLocal
    def setExternalTemperatureProfile(self, val):
        """
        Put value of external temperature profile for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData('external_temperature_profile',val)


    @Variables.noUndo
    def getInternalTemperatureProfile(self):
        """
        Return value of internal temperature profile for the radiative wall
        """
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        val = nod_ray_cond.xmlGetChildDouble('internal_temperature_profile')
        if not val:
            val = self.__defaultValues()['internal_temperature_profile']
            self.setInternalTemperatureProfile(val)

        return val


    @Variables.undoLocal
    def setInternalTemperatureProfile(self, val):
        """
        Put value of internal temperature profile for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData('internal_temperature_profile',val)


    @Variables.noUndo
    def getFlux(self):
        """
        Return value of flux for the radiative wall
        """
        val = self.getValRay('flux')

        return val


    @Variables.undoGlobal
    def setFlux(self, val):
        """
        Put value of flux for the radiative wall
        """
        Model().isGreaterOrEqual(val, 0.)

        self.setValRay(val, 'flux')


    @Variables.noUndo
    def getValRay(self, rayvar):
        """
        Return value of radiative variable named 'var' for the radiative wall
        """
        Model().isInList(rayvar, self.head_list)

        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        if rayvar == "output_zone":
            val = nod_ray_cond.xmlGetInt(rayvar)
        else:
            val = nod_ray_cond.xmlGetDouble(rayvar)
        if not val:
            val = self.__defaultValues()[rayvar]
            self.setValRay(val, rayvar)

        return val


    @Variables.undoLocal
    def setValRay(self, val, rayvar):
        """
        Put value of radiative variable named 'rayvar' for the radiative wall
        """
        Model().isInList(rayvar, self.head_list)
        if rayvar == "output_zone":
            Model().isInt(val)
        else:
            Model().isFloat(val)
        nod_ray_cond = self.boundNode.xmlInitChildNode('radiative_data')
        nod_ray_cond.xmlSetData(rayvar, val)

#-------------------------------------------------------------------------------
# CouplingMobilBoundary
#-------------------------------------------------------------------------------

class CouplingMobilBoundary(Boundary) :
    """
    Boundary class for coupling mobil wall.
    """
    def __new__(cls, label, case) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize default values
        """
        # Default values
        self._defaultValues = {}
        self._defaultValues['initial_displacement_X'    ] = 0
        self._defaultValues['initial_displacement_Y'    ] = 0
        self._defaultValues['initial_displacement_Z'    ] = 0
        self._defaultValues['initial_velocity_X'        ] = 0
        self._defaultValues['initial_velocity_Y'        ] = 0
        self._defaultValues['initial_velocity_Z'        ] = 0
        self._defaultValues['equilibrium_displacement_X'] = 0
        self._defaultValues['equilibrium_displacement_Y'] = 0
        self._defaultValues['equilibrium_displacement_Z'] = 0

        defaultMatrix = '%(t)s11=0;\n%(t)s22=0;\n%(t)s33=0;\n'
        defaultMatrix += '%(t)s12=0;\n%(t)s13=0;\n%(t)s23=0;\n'
        defaultMatrix += '%(t)s21=0;\n%(t)s31=0;\n%(t)s32=0;'
        defaultFluidMatrix = "fx = fluid_fx;\nfy = fluid_fy;\nfz = fluid_fz;"

        self._defaultValues['mass_matrix_formula'       ] = defaultMatrix % {'t':'m'}
        self._defaultValues['damping_matrix_formula'    ] = defaultMatrix % {'t':'c'}
        self._defaultValues['stiffness_matrix_formula'  ] = defaultMatrix % {'t':'k'}
        self._defaultValues['fluid_force_matrix_formula'] = defaultFluidMatrix
        self._defaultValues['DDLX_choice'               ] = 'off'
        self._defaultValues['DDLY_choice'               ] = 'off'
        self._defaultValues['DDLZ_choice'               ] = 'off'


    # Accessors
    #----------

    def _setData(self, nodeName, subNodeName, value):
        """
        Set value into xml file.
        """
        aleNode = self.boundNode.xmlGetNode('ale')
        node = aleNode.xmlInitChildNode(nodeName)
        node.xmlSetData(subNodeName, value)


    def _getDoubleData(self, nodeName, subNodeName, setter ):
        """
        Get value from xml file.
        """
        aleNode = self.boundNode.xmlGetNode('ale')
        node = aleNode.xmlInitChildNode(nodeName)
        value = node.xmlGetChildDouble(subNodeName)
        if not value:
            value = self._defaultValues[nodeName + '_' + subNodeName]
            setter(value)
        return value


    def _getStringData(self, nodeName, subNodeName, setter ):
        """
        Get value from xml file.
        """
        aleNode = self.boundNode.xmlGetNode('ale')
        node = aleNode.xmlInitChildNode(nodeName)
        value = node.xmlGetChildString(subNodeName)
        if not value:
            value = self._defaultValues[nodeName + '_' + subNodeName]
            setter(value)
        return value


    # InitialDisplacement
    #--------------------

    @Variables.undoLocal
    def setInitialDisplacementX(self, value ):
        """
        Set value of initial displacement X into xml file.
        """
        self._setData('initial_displacement', 'X', value)


    @Variables.noUndo
    def getInitialDisplacementX(self ):
        """
        Get value of initial displacement X from xml file.
        """
        return self._getDoubleData('initial_displacement', 'X', self.setInitialDisplacementX)


    @Variables.undoLocal
    def setInitialDisplacementY(self, value ):
        """
        Set value of initial displacement Y into xml file.
        """
        self._setData('initial_displacement', 'Y', value)


    @Variables.noUndo
    def getInitialDisplacementY(self ):
        """
        Get value of initial displacement Y from xml file.
        """
        return self._getDoubleData('initial_displacement', 'Y', self.setInitialDisplacementY)


    @Variables.undoLocal
    def setInitialDisplacementZ(self, value ):
        """
        Set value of initial displacement Z into xml file.
        """
        self._setData('initial_displacement', 'Z', value)


    @Variables.noUndo
    def getInitialDisplacementZ(self ):
        """
        Get value of initial displacement Z from xml file.
        """
        return self._getDoubleData('initial_displacement', 'Z', self.setInitialDisplacementZ)


    # EquilibriumDisplacement
    #------------------------

    @Variables.undoLocal
    def setEquilibriumDisplacementX(self, value):
        """
        Set value of equilibrium displacement X into xml file.
        """
        self._setData('equilibrium_displacement', 'X', value)


    @Variables.noUndo
    def getEquilibriumDisplacementX(self):
        """
        Get value of equilibrium displacement X from xml file.
        """
        return self._getDoubleData('equilibrium_displacement', 'X', self.setEquilibriumDisplacementX)


    @Variables.undoLocal
    def setEquilibriumDisplacementY(self, value):
        """
        Set value of equilibrium displacement Y into xml file.
        """
        self._setData('equilibrium_displacement', 'Y', value)


    @Variables.noUndo
    def getEquilibriumDisplacementY(self):
        """
        Get value of equilibrium displacement Y from xml file.
        """
        return self._getDoubleData('equilibrium_displacement', 'Y', self.setEquilibriumDisplacementY)


    @Variables.undoLocal
    def setEquilibriumDisplacementZ(self, value):
        """
        Set value of equilibrium displacement Z into xml file.
        """
        self._setData('equilibrium_displacement', 'Z', value)


    @Variables.noUndo
    def getEquilibriumDisplacementZ(self):
        """
        Get value of equilibrium displacement X from xml file.
        """
        return self._getDoubleData('equilibrium_displacement', 'Z', self.setEquilibriumDisplacementZ)


    # InitialDisplacement
    #--------------------

    @Variables.undoLocal
    def setInitialVelocityX(self, value):
        """
        Set value of initial velocity X into xml file.
        """
        self._setData('initial_velocity', 'X', value)


    @Variables.noUndo
    def getInitialVelocityX(self):
        """
        Get value of initial velocity X from xml file.
        """
        return self._getDoubleData('initial_velocity', 'X', self.setInitialVelocityX)


    @Variables.undoLocal
    def setInitialVelocityY(self, value):
        """
        Set value of initial velocity Y into xml file.
        """
        self._setData('initial_velocity', 'Y', value)


    @Variables.noUndo
    def getInitialVelocityY(self):
        """
        Get value of initial velocity Y from xml file.
        """
        return self._getDoubleData('initial_velocity', 'Y', self.setInitialVelocityY)


    @Variables.undoLocal
    def setInitialVelocityZ(self, value):
        """
        Set value of initial velocity Z into xml file.
        """
        self._setData('initial_velocity', 'Z', value)


    @Variables.noUndo
    def getInitialVelocityZ(self):
        """
        Get value of initial velocity Z from xml file.
        """
        return self._getDoubleData('initial_velocity', 'Z', self.setInitialVelocityZ)


    # Matrix
    #-------

    @Variables.undoLocal
    def setMassMatrix(self, value):
        """
        Set values of massMatrix into xml file.
        """
        self._setData('mass_matrix', 'formula', value)


    @Variables.noUndo
    def getMassMatrix(self):
        """
        Get values of massMatrix from xml file.
        """
        return self._getStringData('mass_matrix', 'formula', self.setMassMatrix)


    @Variables.undoLocal
    def setStiffnessMatrix(self, value):
        """
        Set values of stiffnessMatrix into xml file.
        """
        self._setData('stiffness_matrix', 'formula', value)


    @Variables.noUndo
    def getStiffnessMatrix(self):
        """
        Get values of stiffnessMatrix from xml file.
        """
        return self._getStringData('stiffness_matrix', 'formula', self.setStiffnessMatrix)


    @Variables.undoLocal
    def setDampingMatrix(self, value):
        """
        Set values of dampingMatrix into xml file.
        """
        self._setData('damping_matrix', 'formula', value)


    @Variables.noUndo
    def getDampingMatrix(self):
        """
        Get values of dampingMatrix from xml file.
        """
        return self._getStringData('damping_matrix', 'formula', self.setDampingMatrix)


    @Variables.undoLocal
    def setFluidForceMatrix(self, value):
        """
        Set values of fluid force matrix into xml file.
        """
        self._setData('fluid_force_matrix', 'formula', value)


    @Variables.noUndo
    def getFluidForceMatrix(self):
        """
        Get values of fluid force matrix from xml file.
        """
        return self._getStringData('fluid_force_matrix', 'formula', self.setFluidForceMatrix)

    # DOF
    #----

    def _setChoice(self, nodeName, value):
        """
        Set the choice
        """
        Model().isInList(value, ['on', 'off'])
        aleNode = self.boundNode.xmlGetNode('ale')
        node = aleNode.xmlInitNode(nodeName)
        node['choice'] = value


    def _getChoice(self, nodeName, setter):
        """
        Get the choice
        """
        aleNode = self.boundNode.xmlGetNode('ale')
        node = aleNode.xmlInitNode(nodeName, 'choice')
        choice = node['choice']

        # Create a defaut choice if it does not exist.
        if not node['choice']:
            choice = self._defaultValues[nodeName + '_choice']
            setter(choice)

        return choice


    @Variables.undoLocal
    def setDDLX(self, value):
        """
        Set the DDLX to xml
        """
        self._setChoice('DDLX', value )


    @Variables.noUndo
    def getDDLX(self):
        """
        Get DDLX from xml
        """
        return self._getChoice('DDLX', self.setDDLX)


    @Variables.undoLocal
    def setDDLY(self, value):
        """
        Set the DDLY to xml
        """
        self._setChoice('DDLY', value )


    @Variables.noUndo
    def getDDLY(self):
        """
        Get DDLY from xml
        """
        return self._getChoice('DDLY', self.setDDLY)


    @Variables.undoLocal
    def setDDLZ(self, value):
        """
        Set the DDLZ to xml
        """
        self._setChoice('DDLZ', value )


    @Variables.noUndo
    def getDDLZ(self):
        """
        Get DDLZ from xml
        """
        return self._getChoice('DDLZ', self.setDDLZ)


#-------------------------------------------------------------------------------
# InletBoundaryModel test case for inlet boundaries conditions
#-------------------------------------------------------------------------------


class InletBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkInletBoundaryInstantiation(self):
        """
        Check whether the InletBoundary class could be instantiated
        """
        model = None
        model = Boundary("inlet", "entree1", self.case)
        assert model != None, 'Could not instantiate InletBoundary'


    def checkSetAndGetVelocityChoice(self):
        """Check whether the velocity choice could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        model.setVelocityChoice('flow1+direction')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="flow1+direction">
                            <flow1>1</flow1>
                            <direction_x>0</direction_x>
                            <direction_y>0</direction_y>
                            <direction_z>0</direction_z>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set velocity choice for inlet boundary'

        assert model.getVelocityChoice() == "flow1+direction",\
           'Could not get velocity choice for inlet boundary'


    def checkSetAndGetFlowAndDirection(self):
        """Check whether the mass or volumic flow could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        model.setVelocityChoice('flow1+direction')
        node =  model._XMLBoundaryConditionsNode

        model.setFlow('flow1', 3.5)
##        model.setFlow('flow2', 3.5)
        model.setDirection('direction_z', 2.0)

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="flow1+direction">
                            <flow1>3.5</flow1>
                            <direction_x>0</direction_x>
                            <direction_y>0</direction_y>
                            <direction_z>2.0</direction_z>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set mass flow or volumic flow and directions for inlet boundary'

        assert model.getVelocityChoice() == "flow1+direction",\
           'Could not get mass flow or volumic flow and directions for inlet boundary'


    def checkSetAndGetNorm(self):
        """Check whether the velocity norm could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('flow1+direction')
        model.setVelocityChoice('norm')
        #il faut explicietement supprimer direction si on ne dirige pas la vitesse
        model.deleteDirectionNodes()

        model.setNorm(999.99)
        doc = '''<boundary_conditions>
                        <inlet label="entree1">
                            <velocity_pressure choice="norm">
                                    <norm>999.99</norm>
                            </velocity_pressure>
                            <turbulence choice="hydraulic_diameter">
                                <hydraulic_diameter>1</hydraulic_diameter>
                            </turbulence>
                        </inlet>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set velocity norm for inlet boundary'

        assert model.getVelocityChoice() == "norm",\
           'Could not get velocity norm for inlet boundary'


    def checkSetAndGetTurbulenceChoice(self):
        """Check whether the turbulence choice could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        model.setTurbulenceChoice('turbulent_intensity')

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                        </velocity_pressure>
                        <turbulence choice="turbulent_intensity">
                            <hydraulic_diameter>1</hydraulic_diameter>
                            <turbulent_intensity>2</turbulent_intensity>
                        </turbulence>
                    </inlet>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set turbulence choice for inlet boundary'

        assert model.getTurbulenceChoice() == "turbulent_intensity",\
           'Could not get turbulence choice for inlet boundary'


    def checkSetAndGetHydraulicDiameterAndTurbulentIntensity(self):
        """
        Check whether the hydraulic_diameter and turbulent_intensity could be
        set and get for inlet boundary.
        """
        model = Boundary("inlet", "entree1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        model.setTurbulenceChoice('turbulent_intensity')

        model.setHydraulicDiameter(120.)
        model.setTurbulentIntensity(0.005)

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                        </velocity_pressure>
                        <turbulence choice="turbulent_intensity">
                            <hydraulic_diameter>120.</hydraulic_diameter>
                            <turbulent_intensity>0.005</turbulent_intensity>
                        </turbulence>
                    </inlet>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set hydraulic_diameter and turbulent_intensity for inlet boundary'

        assert model.getHydraulicDiameter() == 120,\
           'Could not get hydraulic_diameter for turbulence for inlet boundary'
        assert model.getTurbulentIntensity() == 0.005,\
           'Could not get turbulent_intensity for inlet boundary'


    def checkSetAndGetThermalScalar(self):
        """Check whether the thermal scalar could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        model.setTurbulenceChoice('hydraulic_diameter')
        model.th_model.setThermalModel('temperature_celsius')
        model.setScalar("temperature", 15.)

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                        <scalar choice="dirichlet" label="TempC" name="temperature" type="thermal">
                            <dirichlet>15</dirichlet>
                        </scalar>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set thermal scalar value for inlet boundary'

        assert model.getScalar("temperature") == 15,\
           'Could not get thermal scalar value for inlet boundary'


    def checkSetAndGetScalar(self):
        """Check whether the user scalar could be set and get for inlet boundary."""
        model = Boundary("inlet", "entree1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        model.setTurbulenceChoice('hydraulic_diameter')
        model.sca_model.addUserScalar('1', 'sca1')
        model.sca_model.addUserScalar('1', 'sca2')
        model.setScalar('sca1', 11.)
        model.setScalar('sca2', 22.)

        doc = '''<boundary_conditions>
                    <inlet label="entree1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                        <scalar choice="dirichlet" label="sca1" name="sca1" type="user">
                            <dirichlet>11</dirichlet>
                        </scalar>
                        <scalar choice="dirichlet" label="sca2" name="sca1" type="user">
                            <dirichlet>22</dirichlet>
                        </scalar>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set user scalar value for inlet boundary'

        assert model.getScalar('sca2') == 22,\
           'Could not get user scalar value for inlet boundary'


def suite():
    testSuite = unittest.makeSuite(InletBoundaryTestCase, "check")
    return testSuite


def runTest():
    print("InletBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# CoalInletBoundaryModel test case for coal inlet boundaries conditions
#-------------------------------------------------------------------------------


class CoalInletBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkCoalInletBoundaryInstantiation(self):
        """
        Check whether the CoalInletBoundary class could be instantiated
        """
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel(self.case).setCoalCombustionModel('homogeneous_fuel')
        model = None
        model = Boundary("coal_inlet", "entree1", self.case)
        assert model != None, 'Could not instantiate CoalInletBoundary'


    def checkSetAndgetInletType(self):
        """Check whether the type of inlet coal could be set and get for coal inlet boundary."""
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel(self.case).setCoalCombustionModel('homogeneous_fuel')

        model = Boundary("inlet", "charb1", self.case)
        coal_model = Boundary("coal_inlet", "charb1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        coal_model.setInletType('coalFlow')

        doc = '''<boundary_conditions>
                    <inlet label="charb1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                            <temperature>1273.15</temperature>
                            <coal name="coal01">
                                <flow1>1</flow1>
                                <temperature>1273.15</temperature>
                                <ratio name="class01">100</ratio>
                            </coal>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set coalFlow type for coal inlet boundary'

        assert coal_model.getInletType() == "coalFlow",\
           'Could not get coalFlow type for coal inlet boundary'

        coal_model.setInletType('oxydantFlow')
        doc1 = '''<boundary_conditions>
                    <inlet label="charb1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                            <temperature>1273.15</temperature>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                </boundary_conditions>'''
        assert node == self.xmlNodeFromString(doc1),\
           'Could not set oxydantFlow type for coal inlet boundary'

        assert coal_model.getInletType() == "oxydantFlow",\
           'Could not get oxydantFlow type for coal inlet boundary'


    def checkSetAndGetOxydantAndCoalTemperature(self):
        """Check whether the temperature of oxydant and coal could be set and get for coal inlet boundary."""
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel(self.case).setCoalCombustionModel('homogeneous_fuel')

        model = Boundary("inlet", "charb1", self.case)
        coal_model = Boundary("coal_inlet", "charb1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('norm')
        coal_model.setInletType('coalFlow')
        coal_model.setOxydantTemperature(500.)
        coal_model.setCoalTemperature(999.99, 0)

        doc = '''<boundary_conditions>
                    <inlet label="charb1">
                        <velocity_pressure choice="norm">
                            <norm>1</norm>
                            <temperature>500.</temperature>
                            <coal name="coal01">
                                <flow1>1</flow1>
                                <temperature>999.99</temperature>
                                <ratio name="class01">100</ratio>
                            </coal>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set oxydant temperature for coal inlet boundary'

        assert coal_model.getOxydantTemperature() == 500.,\
           'Could not get oxydant temperature for coal inlet boundary'

        assert coal_model.getCoalTemperature(0) == 999.99,\
           'Could not get coal temperature for coal inlet boundary'


    def checkSetAndGetCoalFlow(self):
        """Check whether the flow of inlet coal could be set and get for coal inlet boundary."""
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel(self.case).setCoalCombustionModel('homogeneous_fuel')

        model = Boundary("inlet", "charb1", self.case)
        coal_model = Boundary("coal_inlet", "charb1", self.case)
        node =  model._XMLBoundaryConditionsNode
        model.setVelocityChoice('flow1')
        coal_model.setInletType('coalFlow')
        coal_model.setCoalFlow(123.5, 0)

        doc = '''<boundary_conditions>
                    <inlet label="charb1">
                        <velocity_pressure choice="flow1">
                            <temperature>1273.15</temperature>
                            <flow1>1</flow1>
                            <coal name="coal01">
                                <flow1>123.5</flow1>
                                <temperature>1273.15</temperature>
                                <ratio name="class01">100</ratio>
                            </coal>
                        </velocity_pressure>
                        <turbulence choice="hydraulic_diameter">
                            <hydraulic_diameter>1</hydraulic_diameter>
                        </turbulence>
                    </inlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set flow for coal inlet boundary'

        assert coal_model.getCoalFlow(0) == 123.5,\
           'Could not get flow for coal inlet boundary'

def suite2():
    testSuite = unittest.makeSuite(CoalInletBoundaryTestCase, "check")
    return testSuite


def runTest2():
    print("CoalInletBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite2())

#-------------------------------------------------------------------------------
# WallBoundaryModel test case for wall boundaries conditions
#-------------------------------------------------------------------------------


class WallBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkWallBoundaryInstantiation(self):
        """
        Check whether the WallBoundary class could be instantiated WallBoundary
        """
        model = None
        model = Boundary("wall", "paroi", self.case)
        assert model != None, 'Could not instantiate '


    def checkSetAndGetVelocityChoice(self):
        """Check whether the velocity choice could be set and get for wall boundary."""
        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('on')
        model = Boundary("wall", "fenetre", self.case)
        model.setVelocityChoice('off')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="on">
                            <dirichlet name="velocity_U">0.0</dirichlet>
                            <dirichlet name="velocity_V">0.0</dirichlet>
                            <dirichlet name="velocity_W">0.0</dirichlet>
                        </velocity_pressure>
                    </wall>
                    <wall label="fenetre">
                        <velocity_pressure choice="off"/>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set choice of velocity for wall boundary'

        assert model.getVelocityChoice() == 'off',\
           'Could not get set choice of velocity for wall boundary'

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')

        doc2 = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="fenetre">
                        <velocity_pressure choice="off"/>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc2),\
           'Could not set choice of velocity for wall boundary'


    def checkSetAndGetVelocityValue(self):
        """Check whether the velocity value could be set and get for wall boundary."""
        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('on')
        model.setVelocities(1., 2., 3.)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="on">
                            <dirichlet name="velocity_U">1</dirichlet>
                            <dirichlet name="velocity_V">2</dirichlet>
                            <dirichlet name="velocity_W">3</dirichlet>
                        </velocity_pressure>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set values of velocity for wall boundary'

        assert model.getVelocities() == (1, 2 ,3),\
           'Could not get set values of velocity for wall boundary'


    def checkSetAndGetRoughnessChoiceAndValue(self):
        """Check whether the roughness could be set and get for wall boundary."""
        model = Boundary("wall", "mur", self.case)
        model.setRoughnessChoice('on')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off">
                            <roughness>0.01</roughness>
                        </velocity_pressure>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not update roughness choice for wall boundary'

        model.setRoughness(15.33)
        doc1 = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off">
                            <roughness>15.33</roughness>
                        </velocity_pressure>
                    </wall>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc1),\
           'Could not set roughness value for wall boundary'

        assert model.getRoughness() == 15.33,\
           'Could not get roughness value for wall boundary'


    def checkSetAndGetScalarChoice(self):
        """Check whether the scalar choice could be set and get for wall boundary."""
        model = Boundary("wall", "mur", self.case)
        model.sca_model.addUserScalar('1', 'sca1')
        model.sca_model.addUserScalar('1', 'sca2')
        model.setScalarChoice('sca2', 'dirichlet')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                        <scalar choice="dirichlet" label="sca2" name="sca2" type="user">
                            <dirichlet>0.0</dirichlet>
                        </scalar>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set scalar choice for wall boundary'

        assert model.getScalarChoice('sca2') == 'dirichlet',\
           'Could not get scalar choice for wall boundary'

        assert model.getScalarChoice('sca1') == 'neumann',\
           'Could not get scalar choice for wall boundary'

        model.setScalarChoice('sca1', 'exchange_coefficient')
        doc1 = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                        <scalar choice="exchange_coefficient" label="sca1" name="sca1" type="user">
                            <dirichlet>0.0</dirichlet>
                            <exchange_coefficient>0.0</exchange_coefficient>
                        </scalar>
                        <scalar choice="dirichlet" label="sca2" name="sca2" type="user">
                            <dirichlet>0.0</dirichlet>
                        </scalar>
                    </wall>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc1),\
           'Could not set scalar choice for wall boundary'

        assert model.getScalarChoice('sca1') == 'exchange_coefficient',\
           'Could not get scalar choice for wall boundary'


    def checkSetAndGetScalarImposedValueAndExchangeCoefficient(self):
        """Check whether the scalar values could be set and get for wall boundary."""
        model = Boundary("wall", "mur", self.case)
        model.sca_model.addUserScalar('1', 'sca1')
        model.sca_model.addUserScalar('1', 'sca2')
        model.setScalarChoice('sca1', 'exchange_coefficient')
        model.setScalarChoice('sca2', 'dirichlet')
        model.setScalarImposedValue('sca1', 130.)
        model.setScalarExchangeCoefficient('sca1', 0.130)
        model.setScalarImposedValue('sca2', 55.)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                        <scalar choice="exchange_coefficient" label="sca1" name="sca1" type="user">
                            <dirichlet>130.0</dirichlet>
                            <exchange_coefficient>0.130</exchange_coefficient>
                        </scalar>
                        <scalar choice="dirichlet" label="sca2" name="sca2" type="user">
                            <dirichlet>55.</dirichlet>
                        </scalar>
                    </wall>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set scalar imposed value, flux and exchange_coefficient for wall boundary'

        assert model.getScalarImposedValue('sca2') == 55.,\
           'Could not get scalar imposed value for wall boundary'

        assert model.getScalarImposedValue('sca1') == 130.,\
            'Could not get scalar imposed value for wall boundary'

        assert model.getScalarExchangeCoefficient('sca1') == 0.130,\
            'Could not get scalar imposed value for wall boundary'


def suite3():
    testSuite = unittest.makeSuite(WallBoundaryTestCase, "check")
    return testSuite


def runTest3():
    print("WallBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite3())


#-------------------------------------------------------------------------------
# RadiativeWallBoundaryModel test case for radiative boundaries conditions
#-------------------------------------------------------------------------------


class RadiativeWallBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkRadiativeWallBoundaryInstantiation(self):
        """
        Check whether the RadiativeWallBoundary class could be instantiated RadiativeWallBoundary
        """
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = None
        model = Boundary("radiative_wall", "paroi", self.case)
        assert model != None, 'Could not instantiate '


    def checkSetAndgetRadiativeChoice(self):
        """Check whether the type of condition could be set and get for radiative wall boundary."""
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ipgrno')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ipgrno">
                            <emissivity>0.8</emissivity>
                            <wall_thermal_conductivity>3</wall_thermal_conductivity>
                            <thickness>0.1</thickness>
                            <external_temperature_profile>300</external_temperature_profile>
                            <internal_temperature_profile>300</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set type of condition for radiative wall boundary'

        assert model.getRadiativeChoice() == "ipgrno",\
           'Could not get type of condition for radiative wall boundary'


    def checkSetAndGetEmissivity(self):
        """Check whether the emissivity could be set and get for radiative wall boundary."""
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ipgrno')
        model.setEmissivity(0.22)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ipgrno">
                            <emissivity>0.22</emissivity>
                            <wall_thermal_conductivity>3</wall_thermal_conductivity>
                            <thickness>0.1</thickness>
                            <external_temperature_profile>300</external_temperature_profile>
                            <internal_temperature_profile>300</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set emissivity for radiative wall boundary'

        assert model.getEmissivity() == 0.22,\
           'Could not get emissivity for radiative wall boundary'


    def checkSetAndGetThermalConductivity(self):
        """Check whether the thermal conductivity could be set and get for radiative wall boundary."""
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ipgrno')
        model.setThermalConductivity(5.6)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ipgrno">
                            <emissivity>0.8</emissivity>
                            <wall_thermal_conductivity>5.6</wall_thermal_conductivity>
                            <thickness>0.1</thickness>
                            <external_temperature_profile>300</external_temperature_profile>
                            <internal_temperature_profile>300</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set thermal conductivity for radiative wall boundary'

        assert model.getThermalConductivity() == 5.6,\
           'Could not get thermal conductivity for radiative wall boundary'


    def checkSetAndGetThickness(self):
        """Check whether the thickness could be set and get for radiative wall boundary."""
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ipgrno')
        model.setThickness(2.)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ipgrno">
                            <emissivity>0.8</emissivity>
                            <wall_thermal_conductivity>3</wall_thermal_conductivity>
                            <thickness>2.</thickness>
                            <external_temperature_profile>300</external_temperature_profile>
                            <internal_temperature_profile>300</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set thickness for radiative wall boundary'

        assert model.getThickness() == 2.0,\
           'Could not get thickness for radiative wall boundary'


    def checkSetAndGetExternalAndInternalTemperatureProfile(self):
        """
        Check whether the external and internal temperature profile
        could be set and get for radiative wall boundary.
        """
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ipgrno')
        model.setExternalTemperatureProfile(55.55)
        model.setInternalTemperatureProfile(987.)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ipgrno">
                            <emissivity>0.8</emissivity>
                            <wall_thermal_conductivity>3</wall_thermal_conductivity>
                            <thickness>0.1</thickness>
                            <external_temperature_profile>55.55</external_temperature_profile>
                            <internal_temperature_profile>987.</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set external or internal temperature profile for radiative wall boundary'

        assert model.getExternalTemperatureProfile() == 55.55,\
           'Could not get external temperature profile for radiative wall boundary'

        assert model.getInternalTemperatureProfile() == 987.,\
           'Could not get internal temperature profile for radiative wall boundary'



    def checkSetAndGetFlux(self):
        """Check whether the flux could be set and get for radiative wall boundary."""
        from code_saturne.Pages.ThermalRadiationModel import ThermalRadiationModel
        ThermalRadiationModel(self.case).setRadiativeModel('dom')

        model = Boundary("wall", "mur", self.case)
        model.setVelocityChoice('off')
        model = Boundary("radiative_wall", "radiateur", self.case)
        model.setRadiativeChoice('ifgrno')
        model.setFlux(5.65)
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <wall label="mur">
                        <velocity_pressure choice="off"/>
                    </wall>
                    <wall label="radiateur">
                        <radiative_data choice="ifgrno">
                            <emissivity>0.8</emissivity>
                            <flux>5.65</flux>
                            <internal_temperature_profile>300.</internal_temperature_profile>
                            <output_zone>1</output_zone>
                        </radiative_data>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set flux for radiative wall boundary'

        assert model.getFlux() == 5.65,\
           'Could not get flux for radiative wall boundary'

def suite4():
    testSuite = unittest.makeSuite(RadiativeWallBoundaryTestCase, "check")
    return testSuite


def runTest4():
    print("RadiativeWallBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite4())

#-------------------------------------------------------------------------------
# OutletBoundaryModel test case for outlet boundaries conditions
#-------------------------------------------------------------------------------

class OutletBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkOutletBoundaryInstantiation(self):
        """
        Check whether the OutletBoundary class could be instantiated OutletBoundary
        """
        model = None
        model = Boundary("outlet", "sortie", self.case)
        assert model != None, 'Could not instantiate '


    def checkSetAndGetPressure(self):
        """Check whether the reference pressure could be set and get for outlet boundary."""
        model = Boundary("outlet", "sortie", self.case)
        model.setReferencePressure(111333.)

        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <outlet label="sortie">
                        <dirichlet name="pressure">111333</dirichlet>
                    </outlet>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set reference pressure for outlet boundary'

        assert model.getReferencePressure() == 111333,\
           'Could not get reference pressure for outlet boundary'


    def checkSetAndGetScalarChoiceAndValue(self):
        """Check whether the scalar choice and value could be set and get for outlet boundary."""
        model = Boundary("outlet", "sortie", self.case)
        model.sca_model.addUserScalar('1', 'sca1')
        model.sca_model.addUserScalar('1', 'sca2')
        model.setScalarChoice('sca1', 'dirichlet')
        model.setScalar('sca1', 10.10)
        model.setScalarChoice('sca2', 'neumann')
        node =  model._XMLBoundaryConditionsNode

        doc = '''<boundary_conditions>
                    <outlet label="sortie">
                        <scalar choice="dirichlet" label="sca1" name="sca1" type="user">
                            <dirichlet>10.1</dirichlet>
                        </scalar>
                        <scalar choice="neumann" label="sca2" name="sca2" type="user">
                            <neumann>0</neumann>
                        </scalar>
                    </outlet>
            </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set choice and value at scalars for outlet boundary'

        assert model.getScalarChoice('sca1') == 'dirichlet',\
           'Could not get choice of scalar for outlet boundary'

        assert model.getScalar('sca2') == 0,\
           'Could not get choice of scalar for outlet boundary'


def suite5():
    testSuite = unittest.makeSuite(OutletBoundaryTestCase, "check")
    return testSuite


def runTest5():
    print("OutletBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite5())


#-------------------------------------------------------------------------------
# CouplingMobilBoundary test case for coupling mobil wall boundary
#-------------------------------------------------------------------------------

class CouplingMobilBoundaryTestCase(ModelTest):
    """
    Unittest.
    """
    def checkCouplingMobilBoundaryInstantiation(self):
        """
        Check whether the MobilBoundary class could be instantiated
        """
        model = None
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        assert model != None, 'Could not instantiate '


    def checkSetAndGetInitialDisplacement(self):
        """Check whether coupling mobil wall boundary could be set and get initial displacement."""
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        model.boundNode.xmlInitNode('ale')
        model.setInitialDisplacementX(1)
        model.setInitialDisplacementY(2)
        model.setInitialDisplacementZ(3)
        node =  model._XMLBoundaryConditionsNode
        doc = '''<boundary_conditions>
                    <wall label="Wall_1">
                        <ale>
                            <initial_displacement>
                                <X>1</X>
                                <Y>2</Y>
                                <Z>3</Z>
                            </initial_displacement>
                        </ale>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set initial displacement for coupling mobil wall boundary'
        assert model.getInitialDisplacementX() == 1,\
           'Could not get initial displacement X for coupling mobil wall boundary'
        assert model.getInitialDisplacementY() == 2,\
           'Could not get initial displacement Y for coupling mobil wall boundary'
        assert model.getInitialDisplacementZ() == 3,\
           'Could not get initial displacement Z for coupling mobil wall boundary'


    def checkSetAndGetEquilibriumDisplacement(self):
        """Check whether coupling mobil wall boundary could be set and get equilibrium displacement."""
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        model.boundNode.xmlInitNode('ale')
        model.setEquilibriumDisplacementX(1)
        model.setEquilibriumDisplacementY(2)
        model.setEquilibriumDisplacementZ(3)

        node =  model._XMLBoundaryConditionsNode
        doc = '''<boundary_conditions>
                    <wall label="Wall_1">
                        <ale>
                            <equilibrium_displacement>
                                <X>1</X>
                                <Y>2</Y>
                                <Z>3</Z>
                            </equilibrium_displacement>
                        </ale>
                    </wall>
                 </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set equilibrium displacement for coupling mobil wall boundary'
        assert model.getEquilibriumDisplacementX() == 1,\
           'Could not get equilibrium displacement X for coupling mobil wall boundary'
        assert model.getEquilibriumDisplacementY() == 2,\
           'Could not get equilibrium displacement Y for coupling mobil wall boundary'
        assert model.getEquilibriumDisplacementZ() == 3,\
           'Could not get equilibrium displacement Z for coupling mobil wall boundary'


    def checkSetAndGetInitialVelocity(self):
        """Check whether coupling mobil wall boundary could be set and get initial velocity."""
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        model.boundNode.xmlInitNode('ale')
        model.setInitialVelocityX(1)
        model.setInitialVelocityY(2)
        model.setInitialVelocityZ(3)

        node =  model._XMLBoundaryConditionsNode
        doc = '''<boundary_conditions>
                    <wall label="Wall_1">
                        <ale>
                            <initial_velocity>
                                <X>1</X>
                                <Y>2</Y>
                                <Z>3</Z>
                            </initial_velocity>
                        </ale>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set initial velocity for coupling mobil wall boundary'
        assert model.getInitialVelocityX() == 1,\
          'Could not get initial velocity X for coupling mobil wall boundary'
        assert model.getInitialVelocityY() == 2,\
          'Could not get initial velocity Y for coupling mobil wall boundary'
        assert model.getInitialVelocityZ() == 3,\
          'Could not get initial velocity Z for coupling mobil wall boundary'


    def checkSetAndGetMatrix(self):
        """Check whether coupling mobil wall boundary could be set and get matrix."""
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        model.boundNode.xmlInitNode('ale')
        model.setMassMatrix('MassMatrix')
        model.setStiffnessMatrix('StiffnessMatrix')
        model.setDampingMatrix('DampingMatrix')
        model.setFluidForceMatrix('FluidForceMatrix')

        node =  model._XMLBoundaryConditionsNode
        doc = '''<boundary_conditions>
                    <wall label="Wall_1">
                        <ale>
                            <mass_matrix>
                                <formula>MassMatrix</formula>
                            </mass_matrix>
                            <stiffness_matrix>
                                <formula>StiffnessMatrix</formula>
                            </stiffness_matrix>
                            <damping_matrix>
                                <formula>DampingMatrix</formula>
                            </damping_matrix>
                            <fluid_force_matrix>
                                <formula>FluidForceMatrix</formula>
                            </fluid_force_matrix>
                        </ale>
                    </wall>
                </boundary_conditions>'''


        assert node == self.xmlNodeFromString(doc),\
           'Could not set matrix for coupling mobil wall boundary'
        assert model.getMassMatrix() == 'MassMatrix',\
           'Could not get mass matrix for coupling mobil wall boundary'
        assert model.getStiffnessMatrix() == 'StiffnessMatrix',\
           'Could not get stiffness matrix for coupling mobil wall boundary'
        assert model.getDampingMatrix() == 'DampingMatrix',\
           'Could not get damping matrix for coupling mobil wall boundary'
        assert model.getFluidForceMatrix() == 'FluidForceMatrix',\
           'Could not get fluid force matrix for coupling mobil wall boundary'


    def checkSetAndGetDDL(self):
        """Check whether coupling mobil wall boundary could be set and get DDL."""
        model = Boundary("coupling_mobile_boundary", "Wall_1", self.case)
        model.boundNode.xmlInitNode('ale')
        model.setDDLX('on')
        model.setDDLY('on')
        model.setDDLZ('on')

        node =  model._XMLBoundaryConditionsNode
        doc = '''<boundary_conditions>
                    <wall label="Wall_1">
                        <ale>
                            <DDLX choice="on"/>
                            <DDLY choice="on"/>
                            <DDLZ choice="on"/>
                        </ale>
                    </wall>
                </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not set DDL for coupling mobil wall boundary'
        assert model.getDDLX() == 'on',\
           'Could not get DDL X for coupling mobil wall boundary'
        assert model.getDDLY() == 'on',\
           'Could not get DDL Y for coupling mobil wall boundary'
        assert model.getDDLZ() == 'on',\
           'Could not get DDL Z for coupling mobil wall boundary'


def suite7():
    testSuite = unittest.makeSuite(CouplingMobilBoundaryTestCase, "check")
    return testSuite


def runTest7():
    print("CouplingMobilBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite7())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
