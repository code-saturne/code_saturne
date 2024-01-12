# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.Common import GuiParam
from code_saturne.model.MainFieldsModel import *
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModel
from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModelsDescription
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.TimeTablesModel import TimeTablesModel
from code_saturne.model.SpeciesModel import SpeciesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("Boundary")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class Boundary(Model) :
    """
    Abstract class
    """
    def __new__(cls, nature , label, case, fieldId = None) :
        """
        Factory
        """
        if nature == 'inlet':
            return InletBoundary.__new__(InletBoundary, label, case, fieldId)
        elif nature == 'outlet':
            return OutletBoundary.__new__(OutletBoundary, label, case, fieldId)
        elif nature == 'symmetry':
            return SymmetryBoundary.__new__(SymmetryBoundary, label, case, fieldId)
        elif nature == 'wall':
            return WallBoundary.__new__(WallBoundary, label, case, fieldId)
        else :
            raise ValueError("Unknown boundary nature: " + nature)


    def __init__(self, nature, label, case, fieldId = None) :
        """
        """
        self._label = label
        self._nature = nature
        self.case = case
        self._XMLBoundaryConditionsNode = self.case.xmlGetNode('boundary_conditions')
        self._XMLBoundaryNodes = []
        self.mainFieldsModel = MainFieldsModel(case)
        self._fieldId = fieldId
        self.boundNode = None
        self.spm = SpeciesModel(self.case)

        # Create nodes
        if nature == "inlet" :
            for field in self.mainFieldsModel.getFieldIdList():
                self._XMLBoundaryNodes.append(self._XMLBoundaryConditionsNode.xmlInitNode(nature, field_id = field, label = label))

        elif nature == "outlet" :
            self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('outlet', field_id = "none", label = label)
            for field in self.mainFieldsModel.getFieldIdList():
                self._XMLBoundaryNodes.append(self._XMLBoundaryConditionsNode.xmlInitNode(nature, field_id = field, label = label))
        elif nature == "wall" :
            self.boundNode = self._XMLBoundaryConditionsNode.xmlInitNode('wall', field_id = "none", label = label)

        self._initBoundary()


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node (vitual method)
        """
        pass


    def getLabel(self):
        """
        Return the label
        """
        return self._label


    def getNature(self):
        """
        Return the nature
        """
        return self._nature


    def delete(self):
        """
        Delete Boundary
        """
        for node in self._XMLBoundaryNodes :
            node.xmlRemoveNode()
        if self.boundNode != None :
            self.boundNode.xmlRemoveNode()


    @Variables.noUndo
    def getDefaultEnthalpyFormula(self, energy_choice='dirichlet'):
        """
        Get default enthalpy formula
        """

        if energy_choice == 'dirichlet_formula':
            formula = "enthalpy = 50000.0;"
        elif energy_choice == 'flux_formula':
            formula = "flux = 1000.;"
        elif energy_choice in ['temperature_formula', 'timp_K_formula']:
            formula = "temperature = 293.15;"
        else:
            formula = None

        return formula


    @Variables.noUndo
    def getDefaultScalarFormula(self, choice, scalar_name):
        """
        Get default scalar formula.
        """

        if choice == 'dirichlet_formula':
            formula = "%s = 1.0;" % str(scalar_name)
        elif choice == 'flux_formula':
            formula = "flux = 0.;"
        else:
            formula = None

        return formula


    @Variables.noUndo
    def getScalarFormulaComponents(self, fieldId, Scalar):
        """
        Get scalar formula components
        """
        choice = self.getScalarChoice(fieldId, Scalar)
        if choice[-7:] != 'formula':
            raise Exception("Error: Not a formula")

        sname = self.spm.getScalarLabelByName(Scalar)
        if choice == 'dirichlet_formula':
            req = [(sname, sname)]
        elif choice == 'flux_formula':
            req = [('flux', 'Surface flux for scalar \"%s\"' % str(sname))]

        exp = self.getScalarValue(fieldId, Scalar)

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration'),
               ('surface', 'Boundary zone surface')]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        # Time Tables variables
        sym += TimeTablesModel(self.case).getTableVariablesListAll()

        return exp, req, sym


#-------------------------------------------------------------------------------
# InletBoundary class
#-------------------------------------------------------------------------------

class InletBoundary(Boundary):
    """
    """
    def __new__(cls, label, case, fieldId) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self.__velocityChoices = ['norm', 'flow1', 'norm_formula', 'flow1_formula']
        self.__directionChoices = ['normal', 'formula']
        self.__directionTags = ['direction_formula']
        self.__turbulenceChoices = ['hydraulic_diameter', 'turbulent_intensity', 'formula']
        self.__enthalpyChoices = ['flux', 'enthalpy', 'timp_K', 'hsat_P']

        # Initialize nodes if necessary
        for field in self.mainFieldsModel.getFieldIdList():
            self.getVelocityChoice(field)
            self.getDirectionChoice(field)


    def __defaultValues(self, fieldId):
        """
        Default values
        """
        dico = {}
        dico['velocityChoice'] = 'norm'
        dico['directionChoice'] = 'normal'
        dico['turbulenceChoice'] = 'hydraulic_diameter'
        dico['hydraulic_diameter'] = 1
        dico['turbulent_intensity'] = 2
        dico['velocity'] = 0.0
        dico['flow1'] = 1
        dico['norm'] = 1
        dico['flow1_formula'] = "q_m = 1;"
        dico['norm_formula'] = "u_norm = 1;"
        dico['direction_formula'] = "dir_x = 1;\ndir_y = 0;\ndir_z = 0;\n"
        dico['enthalpy'] = 0.
        dico['EnthalpyModel'] = 'flux'
        dico['fraction'] = 0.
        dico['diameter'] = 1.e-3
        dico['noncondensable'] = 0.
        dico['scalar'] = 0.
        dico['scalarModel'] = 'dirichlet'

        # Check turbulence only if a real phase if provided
        if fieldId == "none":
            turb_model = None
        else:
            turb_model =  TurbulenceModel(self.case).getTurbulenceModel(fieldId)

        if turb_model in TurbulenceModelsDescription.dispersedTurbulenceModels:
            dico['turbulenceChoice'] = 'formula'
        return dico


    def __initChoiceForVelocityAndDirection(self, fieldId):
        """
        Get the choice of velocity and direction.
        """
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        childnode = node.xmlInitNode('velocity', 'choice', 'direction')

        choice = childnode['choice']
        if not choice:
            choice = self.__defaultValues(fieldId)['velocityChoice']
            self.setVelocityChoice(fieldId, choice)
        dir = childnode['direction']
        if not dir:
            dir = self.__defaultValues(fieldId)['directionChoice']
            self.setDirectionChoice(fieldId, dir)
        return choice, dir


    @Variables.noUndo
    def getVelocityChoice(self, fieldId):
        """
        Get the choice of velocity.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        choice, dir = self.__initChoiceForVelocityAndDirection(fieldId)
        return choice


    @Variables.noUndo
    def getDirectionChoice(self, fieldId):
        """
        Get the choice of direction.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        choice, dir = self.__initChoiceForVelocityAndDirection(fieldId)
        return dir


    @Variables.noUndo
    def getVelocity(self, fieldId):
        """
        Get value of velocity beyond choice.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        choice = self.getVelocityChoice(fieldId)
        Model().isInList(choice, self.__velocityChoices)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlGetNode('velocity')

        if choice in ('norm', 'flow1'):
            value = XMLVelocityNode.xmlGetChildDouble(choice)
        elif choice in ('norm_formula', 'flow1_formula'):
            value = XMLVelocityNode.xmlGetChildString(choice)
        if value in [None, ""]:
            value = self.__defaultValues(fieldId)[choice]
            self.setVelocity(fieldId, value)

        return value


    @Variables.undoLocal
    def setVelocity(self, fieldId, value):
        """
        Set value of velocity.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        choice = self.getVelocityChoice(fieldId)
        Model().isInList(choice, self.__velocityChoices)

        if choice in ('norm', 'flow1'):
            Model().isFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlInitNode('velocity')
        XMLVelocityNode.xmlSetData(choice, value)


    @Variables.noUndo
    def getDirection(self, fieldId, component):
        """
        Get the component velocity
        """
        Model().isInList(component, self.__directionTags)
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlGetNode('velocity')

        Model().isInList(component, ('direction_formula',))
        value = XMLVelocityNode.xmlGetChildString(component)

        if value is None :
            value = self.__defaultValues(fieldId)[component]
            self.setDirection(fieldId, component, value)
        return value


    @Variables.undoLocal
    def setDirection(self, fieldId, component, value):
        """
        Set the component velocity for fieldLabel
        """
        Model().isInList(component, self.__directionTags)
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlInitNode('velocity')
        XMLVelocityNode.xmlSetData(component, value)


    @Variables.undoLocal
    def setVelocityChoice(self, fieldId, value):
        """
        Set the velocity definition according to choice
        """
        Model().isInList(value, self.__velocityChoices)
        self.mainFieldsModel.isFieldIdValid(fieldId)

        # Check if value is a new velocity choice value
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlInitNode('velocity')
        if XMLVelocityNode['choice'] != None :
            if XMLVelocityNode['choice'] == value:
                return

        # Update velocity choice
        XMLVelocityNode['choice'] = value
        self.getVelocity(fieldId)

        for tag in self.__velocityChoices:
            if tag != value:
                XMLVelocityNode.xmlRemoveChild(tag)


    @Variables.undoLocal
    def setDirectionChoice(self, fieldId, value):
        """
        Set the direction of the flow definition according to choice.
        """
        Model().isInList(value, self.__directionChoices)
        self.mainFieldsModel.isFieldIdValid(fieldId)

        # Check if value is a new direction choice
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLVelocityNode = node.xmlInitNode('velocity')
        if XMLVelocityNode['direction'] != None :
            if XMLVelocityNode['direction'] == value:
                return

        # Update direction choice
        XMLVelocityNode['direction'] = value

        if value == 'formula':
            self.getDirection(fieldId, 'direction_formula')
        else:
            XMLVelocityNode.xmlRemoveChild('direction_formula')


    @Variables.noUndo
    def getTurbulenceChoice(self, fieldId):
        """
        Get the turbulence choice
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')

        choice = XMLTurbulenceNode['choice']
        if choice not in self.__turbulenceChoices :
            choice = self.__defaultValues(fieldId)['turbulenceChoice']
            self.setTurbulenceChoice(fieldId, choice)

        return choice


    @Variables.undoLocal
    def setTurbulenceChoice(self, fieldId, value):
        """
        Set the choice turbulence
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isInList(value, self.__turbulenceChoices)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')

        if XMLTurbulenceNode['choice'] != None:
            if XMLTurbulenceNode['choice'] == value:
                return

        XMLTurbulenceNode['choice'] = value

        # Update values
        if value == 'hydraulic_diameter' :
            self.getHydraulicDiameter(fieldId)
            XMLTurbulenceNode.xmlRemoveChild('turbulent_intensity')
            XMLTurbulenceNode.xmlRemoveChild('formula')

        elif value == 'turbulent_intensity' :
            self.getHydraulicDiameter(fieldId)
            self.getTurbulentIntensity(fieldId)
            XMLTurbulenceNode.xmlRemoveChild('hydraulic_diameter')
            XMLTurbulenceNode.xmlRemoveChild('formula')

        elif value == 'formula' :
            self.getTurbFormula(fieldId)
            XMLTurbulenceNode.xmlRemoveChild('turbulent_intensity')
            XMLTurbulenceNode.xmlRemoveChild('hydraulic_diameter')


    @Variables.noUndo
    def getHydraulicDiameter(self, fieldId):
        """
        Get hydraulic diameter
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'],  self.__turbulenceChoices)
        value = XMLTurbulenceNode.xmlGetDouble('hydraulic_diameter')
        if value is None :
            value = self.__defaultValues(fieldId)['hydraulic_diameter']
            self.setHydraulicDiameter(fieldId, value)
        return value


    @Variables.undoLocal
    def setHydraulicDiameter(self, fieldId, value):
        """
        Set hydraulic diameter
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isStrictPositiveFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], self.__turbulenceChoices)
        XMLTurbulenceNode.xmlSetData('hydraulic_diameter', value)


    @Variables.undoLocal
    def setTurbFormula(self, fieldId, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')

        if not XMLTurbulenceNode:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        Model().isInList(XMLTurbulenceNode['choice'], self.__turbulenceChoices)
        n = XMLTurbulenceNode.xmlInitChildNode('formula')
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getTurbFormula(self, fieldId):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')

        formula = XMLTurbulenceNode.xmlGetString('formula')

        return formula


    @Variables.noUndo
    def getTurbFormulaComponents(self, fieldId, turbModel):

        exp = self.getTurbFormula(fieldId)
        sym = [('x','face center coordinate'),
               ('y','face center coordinate'),
               ('z','face center coordinate'),
               ('t','time'),
               ('dt','time step'),
               ('iter','number of time step'),
               ('surface', 'Boundary zone surface')]

        for (name, val) in NotebookModel(self.case).getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        # Time Tables variables
        sym += TimeTablesModel(self.case).getTableVariablesListAll()

        if turbModel in ('k-epsilon', 'k-epsilon_linear_production'):
            req = [('k', "turbulent energy"),
                   ('eps', "turbulent dissipation")]

        elif turbModel in ('rij-epsilon_ssg', 'rij-epsilon_ebrsm'):
            req = [('R11', "Reynolds stress R11"),
                   ('R22', "Reynolds stress R22"),
                   ('R33', "Reynolds stress R33"),
                   ('R12', "Reynolds stress R12"),
                   ('R23', "Reynolds stress R13"),
                   ('R13', "Reynolds stress R23"),
                   ('eps', "turbulent dissipation")]

        elif turbModel in ('q2-q12-tchen', 'q2-q12'):
            req = [('q2', "turbulent kinetic energy"),
                   ('q12', "covariance")]

        elif turbModel in ('r2-q12'):
            req = [('R11', "Reynolds stress R11"),
                   ('R22', "Reynolds stress R22"),
                   ('R33', "Reynolds stress R33"),
                   ('R12', "Reynolds stress R12"),
                   ('R23', "Reynolds stress R13"),
                   ('R13', "Reynolds stress R23"),
                   ('q12', "covariance")]

        elif turbModel in ('r2-r12-tchen'):
            req = [('R11', "Reynolds stress R11"),
                   ('R22', "Reynolds stress R22"),
                   ('R33', "Reynolds stress R33"),
                   ('R12', "Reynolds stress R12"),
                   ('R23', "Reynolds stress R13"),
                   ('R13', "Reynolds stress R23"),
                   ('R12-11', "Reynolds stress R11"),
                   ('R12-22', "Reynolds stress R22"),
                   ('R12-33', "Reynolds stress R33"),
                   ('R12-12', "Reynolds stress R12"),
                   ('R12-13', "Reynolds stress R13"),
                   ('R12-23', "Reynolds stress R23")]
        else:
            req = []

        return exp, req, sym

    @Variables.noUndo
    def getDefaultTurbFormula(self, turb_model):
        """
        Get defaut turbulence formula
        """
        if turb_model in ('k-epsilon', 'k-epsilon_linear_production'):
            formula = """k = 0.0001;
eps = 0.001;"""

        elif turb_model in ('rij-epsilon_ssg', 'rij-epsilon_ebrsm'):
            formula = """R11 = 5e-05;
R22 = 5e-05;
R33 = 5e-05;
R12 = 5e-05;
R13 = 5e-05;
R23 = 5e-05;
eps = 0.001;"""

        elif turb_model in ('q2-q12-tchen', 'q2-q12'):
            formula = """q2 = 5.e-05;
q12 = 0.0001;"""

        elif turb_model == 'r2-q12':
            formula = """R11 = 5e-05;
R22 = 5e-05;
R33 = 5e-05;
R12 = 5e-05;
R13 = 5e-05;
R23 = 5e-05;
q12 = 0.0001;"""

        elif turb_model == 'r2-r12-tchen':
            formula = """R11 = 5e-05;
R22 = 5e-05;
R33 = 5e-05;
R12 = 5e-05;
R13 = 5e-05;
R23 = 5e-05;
R12-11 = 5e-05;
R12-22 = 5e-05;
R12-33 = 5e-05;
R12-12 = 5e-05;
R12-13 = 5e-05;
R12-23 = 5e-05;"""

        return formula


    @Variables.noUndo
    def getTurbulentIntensity(self, fieldId):
        """
        Get turbulent intensity
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], ('turbulent_intensity',))
        value = XMLTurbulenceNode.xmlGetDouble('turbulent_intensity')
        if value is None :
            value = self.__defaultValues(fieldId)['turbulent_intensity']
            self.setTurbulentIntensity(fieldId, value)

        return value


    @Variables.undoLocal
    def setTurbulentIntensity(self, fieldId, value):
        """
        Set turbulent intensity
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isStrictPositiveFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLTurbulenceNode = node.xmlInitNode('turbulence')
        Model().isInList(XMLTurbulenceNode['choice'], ('turbulent_intensity',))
        XMLTurbulenceNode.xmlSetData('turbulent_intensity', value)


    @Variables.noUndo
    def getEnthalpyChoice(self, fieldId):
        """
        Get the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLEnergyNode = node.xmlInitNode('variable', 'choice', name='enthalpy')

        choice = XMLEnergyNode['choice']
        if not choice:
            choice = self.__defaultValues(fieldId)['EnthalpyModel']
            self.setEnthalpyChoice(fieldId, choice)
        return choice


    @Variables.undoLocal
    def setEnthalpyChoice(self, fieldId, value):
        """
        Set the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLEnergyNode = node.xmlInitNode('variable', 'choice', name='enthalpy')

        XMLEnergyNode['choice'] = value

        if value == 'hsat_P' :
            Childnode = XMLEnergyNode.xmlGetChildNode('value')
            if Childnode != None :
                Childnode.xmlRemoveNode()


    @Variables.noUndo
    def getEnthalpy(self, fieldId) :
        """
        Get energy value for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet",
                                                          field_id = fieldId,
                                                          label = self._label)
        XMLEnergyNode = node.xmlGetChildNode('variable', 'choice', name='enthalpy')

        choice = self.getEnthalpyChoice(fieldId)

        if "formula" == choice[-7:]:
            value = XMLEnergyNode.xmlGetChildString(choice)
        else:
            value = XMLEnergyNode.xmlGetChildDouble('value')

        if value in [None, ""]:
            if choice[-7:] == 'formula':
                value = self.getDefaultEnthalpyFormula(choice)
            else:
                value = self.__defaultValues(fieldId)['enthalpy']
            self.setEnthalpy(fieldId, value)

        return value


    @Variables.undoLocal
    def setEnthalpy(self, fieldId, value) :
        """
        Set energy value for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        _c = self.getEnthalpyChoice(fieldId) \
                if "_formula" in self.getEnthalpyChoice(fieldId) else 'value'

        if _c == 'value':
            Model().isFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet",
                                                          field_id = fieldId,
                                                          label = self._label)

        XMLEnergyNode = node.xmlGetChildNode('variable', 'choice', name='enthalpy')

        XMLEnergyNode.xmlSetData(_c, str(value))


    @Variables.noUndo
    def getFraction(self, fieldId):
        """
        Get fraction for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', choice='dirichlet', name='volume_fraction')

        Childnode = XMLFractionNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues(fieldId)['fraction']
            self.setFraction(fieldId, value)

        value = XMLFractionNode.xmlGetChildDouble('value')
        return value


    @Variables.undoLocal
    def setFraction(self, fieldId, value):
        """
        Set fraction for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isPositiveFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', choice='dirichlet', name='volume_fraction')
        XMLFractionNode.xmlSetData('value', str(value))
        # TODO
        # Check (sum(fractions) <= 1)


    @Variables.noUndo
    def getDiameter(self, fieldId) :
        """
        Get diameter variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLDiameterNode = node.xmlInitNode('variable', choice='dirichlet', name='diameter')

        Childnode = XMLDiameterNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues(fieldId)['diameter']
            self.setDiameter(fieldId, value)

        value = XMLDiameterNode.xmlGetChildDouble('value')
        return value


    @Variables.undoLocal
    def setDiameter(self, fieldId, value) :
        """
        Set diameter variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isStrictPositiveFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLDiameterNode = node.xmlInitNode('variable', choice='dirichlet', name='diameter')
        XMLDiameterNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getNonCondensableValue(self, fieldId, NonCondensable):
        """
        Get non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLNonCondensableNode = node.xmlInitNode('variable', choice='dirichlet', name=NonCondensable)

        Childnode = XMLNonCondensableNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues(fieldId)['noncondensable']
            self.setNonCondensableValue(fieldId, NonCondensable, value)

        value = XMLNonCondensableNode.xmlGetChildDouble('value')
        return value


    @Variables.undoLocal
    def setNonCondensableValue(self, fieldId, NonCondensable, value):
        """
        Set non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isPositiveFloat(value)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("inlet", field_id = fieldId, label = self._label)
        XMLNonCondensableNode = node.xmlInitNode('variable', choice='dirichlet', name=NonCondensable)
        XMLNonCondensableNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getScalarChoice(self, fieldId, scalar):
        """
        Get the boundary condition used for the scalar.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("inlet", field_id = fieldId, label = self._label)

        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=scalar)

        choice = XMLScalarNode['choice']
        if not choice:
            choice = self.__defaultValues(fieldId)['scalarModel']
            self.setScalarChoice(fieldId, scalar, choice)

        return choice


    @Variables.undoLocal
    def setScalarChoice(self, fieldId, scalar, value):
        """
        Set boundary condition type for a scalar.
        """

        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("inlet", field_id = fieldId, label = self._label)

        XMLScalarNode = node.xmlInitNode('variable',
                                         'choice',
                                         name=scalar)

        XMLScalarNode['choice'] = value


    @Variables.noUndo
    def getScalarValue(self, fieldId, Scalar):
        """
        Get non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("inlet",
                                                           field_id = fieldId,
                                                           label = self._label)
        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=Scalar)

        choice = self.getScalarChoice(fieldId, Scalar)

        if choice[-7:] == 'formula':
            value = XMLScalarNode.xmlGetChildString('value')
        else:
            value = XMLScalarNode.xmlGetChildDouble('value')

        if value in ['', None]:
            if choice[-7:] == 'formula':
                value = self.getDefaultScalarFormula(choice,
                                                     self.spm.getScalarLabelByName(Scalar))
            else:
                value = self.__defaultValues(fieldId)['scalar']

            self.setScalarValue(fieldId, Scalar, value)

        return value


    @Variables.undoLocal
    def setScalarValue(self, fieldId, Scalar, value):
        """
        Set non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("inlet", field_id = fieldId, label = self._label)
        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=Scalar)
        XMLScalarNode.xmlSetData('value', str(value))


    def getXMLBoundaryConditionNode(self):

        return self._XMLBoundaryConditionsNode


#-------------------------------------------------------------------------------
# Outlet boundary
#-------------------------------------------------------------------------------

class OutletBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case, fieldId) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self.getReferencePressure()
        self.__enthalpyChoices = ['flux', 'dirichlet', 'timp_K', 'hsat_P']


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['reference_pressure'] = 101325.0
        dico['enthalpy']           = 0.
        dico['EnthalpyModel']      = 'flux'
        dico['FractionModel']      = 'dirichlet'
        dico['PressureModel']      = 'dpdndtau'
        dico['fraction']           = 0.
        dico['noncondensable']     = 0.
        dico['scalar']             = 0.
        dico['scalarModel']        = 'dirichlet'
        return dico


    @Variables.undoLocal
    def setPressureChoice(self, value):
        """
        Set reference boundary condition choice
        """

        Model().isInList(value, ['dirichlet', 'dpdndtau'])

        node = self.boundNode.xmlInitNode('variable', 'choice', name='pressure')
        node['choice'] = value


    @Variables.noUndo
    def getPressureChoice(self):
        """
        Get pressure boundary condition choice
        """

        node = self.boundNode.xmlInitNode('variable', 'choice', name="pressure")
        choice = node['choice']

        if not choice:
            choice = self.__defaultValues()['PressureModel']
            self.setPressureChoice(choice)

        return choice


    @Variables.noUndo
    def getReferencePressure(self) :
        """
        Get reference pressure
        """
        node = self.boundNode.xmlGetChildNode('variable', 'choice', name='pressure')
        if not node:
            node = self.boundNode.xmlInitChildNode('variable', 'choice', name='pressure')

        ValNode = node.xmlGetChildNode('value')
        if not ValNode:
            value = self.__defaultValues()['reference_pressure']
            self.setReferencePressure(value)

        value = node.xmlGetDouble('value')

        return value


    @Variables.undoLocal
    def setReferencePressure(self, value) :
        """
        Set reference pressure
        """
        Model().isFloat(value)

        node = self.boundNode.xmlGetChildNode('variable', 'choice', name='pressure')
        if not node:
            node = self.boundNode.xmlInitChildNode('variable', 'choice', name='pressure')
        node.xmlSetData('value', str(value))


    @Variables.noUndo
    def getEnthalpyChoice(self, fieldId):
        """
        Get the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet",
                                                          field_id = fieldId,
                                                          label = self._label)
        XMLEnergyNode = node.xmlInitNode('variable', 'choice', name='enthalpy')

        choice = XMLEnergyNode['choice']
        if not choice:
            choice = self.__defaultValues()['EnthalpyModel']
            self.setEnthalpyChoice(fieldId, choice)
        return choice


    @Variables.undoLocal
    def setEnthalpyChoice(self, fieldId, value):
        """
        Set the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLEnergyNode = node.xmlInitNode('variable', 'choice', name='enthalpy')

        XMLEnergyNode['choice'] = value

        if value == 'hsat_P' :
            Childnode = XMLEnergyNode.xmlGetChildNode('value')
            if Childnode != None :
                Childnode.xmlRemoveNode()


    @Variables.noUndo
    def getEnthalpy(self, fieldId) :
        """
        Get energy value for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet",
                                                          field_id = fieldId,
                                                          label = self._label)
        XMLEnergyNode = node.xmlGetChildNode('variable', 'choice', name='enthalpy')

        choice = self.getEnthalpyChoice(fieldId)

        if "formula" == choice[-7:]:
            value = XMLEnergyNode.xmlGetChildString(choice)
        else:
            value = XMLEnergyNode.xmlGetChildDouble('value')

        if value in [None, ""]:
            if choice[-7:] == 'formula':
                value = self.getDefaultEnthalpyFormula(choice)
            else:
                value = self.__defaultValues()['enthalpy']
            self.setEnthalpy(fieldId, value)

        return value


    @Variables.undoLocal
    def setEnthalpy(self, fieldId, value) :
        """
        Set energy value for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        if '_formula' not in self.getEnthalpyChoice(fieldId):
            Model().isFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLEnergyNode = node.xmlGetChildNode('variable', 'choice', name='enthalpy')
        XMLEnergyNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getFractionChoice(self, fieldId):
        """
        Get the fraction choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', 'choice', name='volume_fraction')

        choice = XMLFractionNode['choice']
        if not choice:
            choice = self.__defaultValues()['FractionModel']
            self.setFractionChoice(fieldId, choice)
        return choice


    @Variables.undoLocal
    def setFractionChoice(self, fieldId, value):
        """
        Set the fraction choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isInList(value, ['dirichlet','automatic'])

        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', 'choice', name='volume_fraction')

        XMLFractionNode['choice'] = value

        if value == 'automatic' :
            Childnode = XMLFractionNode.xmlGetChildNode('value')
            if Childnode != None :
                Childnode.xmlRemoveNode()


    @Variables.noUndo
    def getFraction(self, fieldId):
        """
        Get fraction for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', choice='dirichlet', name='volume_fraction')

        Childnode = XMLFractionNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues()['fraction']
            self.setFraction(fieldId, value)

        value = XMLFractionNode.xmlGetChildDouble('value')
        return value


    @Variables.undoLocal
    def setFraction(self, fieldId, value):
        """
        Set fraction for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isPositiveFloat(value)

        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLFractionNode = node.xmlInitNode('variable', choice='dirichlet', name='volume_fraction')
        XMLFractionNode.xmlSetData('value', str(value))
        # TODO
        # Check (sum(fractions) <= 1)


    @Variables.noUndo
    def getNonCondensableValue(self, fieldId, NonCondensable):
        """
        Get non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLNonCondensableNode = node.xmlInitNode('variable', choice='dirichlet', name=NonCondensable)

        Childnode = XMLNonCondensableNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues()['noncondensable']
            self.setNonCondensableValue(fieldId, NonCondensable, value)

        value = XMLNonCondensableNode.xmlGetChildDouble('value')
        return value


    @Variables.undoLocal
    def setNonCondensableValue(self, fieldId, NonCondensable, value):
        """
        Set non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isPositiveFloat(value)
        node = self._XMLBoundaryConditionsNode.xmlGetNode("outlet", field_id = fieldId, label = self._label)
        XMLNonCondensableNode = node.xmlInitNode('variable', choice='dirichlet', name=NonCondensable)
        XMLNonCondensableNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getScalarChoice(self, fieldId, scalar):
        """
        Get the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("outlet", field_id = fieldId, label = self._label)
        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=scalar)

        choice = XMLScalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarModel']
            self.setScalarChoice(fieldId, scalar, choice)

        return choice


    @Variables.undoLocal
    def setScalarChoice(self, fieldId, scalar, value):
        """
        Set the enthalpy choice for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)

        node = self._XMLBoundaryConditionsNode.xmlInitNode("outlet", field_id = fieldId, label = self._label)

        XMLScalarNode = node.xmlInitNode('variable',
                                         'choice',
                                         name=scalar)

        XMLScalarNode['choice'] = value


    @Variables.noUndo
    def getScalarValue(self, fieldId, Scalar):
        """
        Get non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("outlet",
                                                           field_id = fieldId,
                                                           label = self._label)
        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=Scalar)

        choice = self.getScalarChoice(fieldId, Scalar)

        if choice[-7:] == 'formula':
            value = XMLScalarNode.xmlGetChildString('value')
        else:
            value = XMLScalarNode.xmlGetChildDouble('value')

        if value in ['', None]:
            if choice[-7:] == 'formula':
                value = self.getDefaultScalarFormula(choice,
                                                     self.spm.getScalarLabelByName(Scalar))
            else:
                value = self.__defaultValues()['scalar']

            self.setScalarValue(fieldId, Scalar, value)

        return value


    @Variables.undoLocal
    def setScalarValue(self, fieldId, Scalar, value):
        """
        Set non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)
        node = self._XMLBoundaryConditionsNode.xmlInitNode("outlet", field_id = fieldId, label = self._label)
        XMLScalarNode = node.xmlInitNode('variable', 'choice', name=Scalar)
        XMLScalarNode.xmlSetData('value', str(value))


    def getXMLBoundaryConditionNode(self):

        return self._XMLBoundaryConditionsNode


#-------------------------------------------------------------------------------
# Symmetry boundary
#-------------------------------------------------------------------------------

class SymmetryBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case, fieldId = None) :
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
# Wall boundary
#-------------------------------------------------------------------------------

class WallBoundary(Boundary) :
    """
    """
    def __new__(cls, label, case, fieldId = None) :
        """
        Constructor
        """
        return object.__new__(cls)


    def _initBoundary(self):
        """
        Initialize the boundary, add nodes in the boundary node
        """
        self._fluxChoices=['temperature', 'flux']
        self._wallModel=['adherence', 'friction', 'du2_dn', 'dvr_dn', 'droplet_friction']


    def __defaultValues(self):
        """
        Default values
        """
        dico = {}
        dico['EnthalpyModel'] = 'flux'
        dico['flux']          = 0.0
        dico['scalarModel']   = 'flux'
        dico['scalar']        = 0.0
        return dico


    @Variables.noUndo
    def getEnthalpyChoice(self, fieldId):
        """
        Get the enthalpy choice for field
        """
        XMLEnergyNode = self.boundNode.xmlInitNode('variable', 'choice', name='enthalpy')

        choice = XMLEnergyNode['choice']
        if not choice:
            choice = self.__defaultValues()['EnthalpyModel']
            self.setEnthalpyChoice(fieldId, choice)
        return choice


    @Variables.undoLocal
    def setEnthalpyChoice(self, fieldId, value):
        """
        Set the enthalpy choice for field
        """

        XMLEnergyNode = self.boundNode.xmlInitNode('variable', 'choice', name='enthalpy')

        XMLEnergyNode['choice'] = value


    @Variables.noUndo
    def getEnthalpy(self, fieldId) :
        """
        Get energy value for field
        """
        XMLEnergyNode = self.boundNode.xmlGetChildNode('variable', 'choice', name='enthalpy')


        Childnode = XMLEnergyNode.xmlGetChildNode('value')
        if Childnode is None :
            value = self.__defaultValues()['flux']
            self.setEnthalpy(fieldId, value)

        choice = self.getEnthalpyChoice("none")
        value = XMLEnergyNode.xmlGetChildString("value")
        if value in [None, ""]:
            if choice[-7:] == 'formula':
                value = self.getDefaultEnthalpyFormula(choice)
            else:
                value = self.__defaultValues()['flux']
            self.setEnthalpy('none', value)

        if choice != "syrthes_coupling" and choice[-7:] != 'formula':
            value = float(value)

        return value


    @Variables.undoLocal
    def setEnthalpy(self, fieldId, value) :
        """
        Set energy value for field
        """
        if self.getEnthalpyChoice(fieldId) in ['temperature','flux']:
            Model().isFloat(value)

        XMLEnergyNode = self.boundNode.xmlGetChildNode('variable', 'choice', name='enthalpy')
        XMLEnergyNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getScalarChoice(self, fieldId, scalar):
        """
        Get the boundary condition used for the scalar.
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)

        XMLScalarNode = self.boundNode.xmlInitNode('variable',
                                                   'choice',
                                                   name=scalar)

        choice = XMLScalarNode['choice']
        if not choice:
            choice = self.__defaultValues()['scalarModel']
            self.setScalarChoice(fieldId, scalar, choice)

        return choice


    @Variables.undoLocal
    def setScalarChoice(self, fieldId, scalar, value):
        """
        Set boundary condition type for a scalar.
        """

        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)

        XMLScalarNode = self.boundNode.xmlInitNode('variable',
                                                   'choice',
                                                   name=scalar)

        XMLScalarNode['choice'] = value


    @Variables.noUndo
    def getScalarValue(self, fieldId, Scalar):
        """
        Get non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)

        XMLScalarNode = self.boundNode.xmlInitNode('variable',
                                                   'choice',
                                                   name=Scalar)

        choice = self.getScalarChoice(fieldId, Scalar)

        if choice[-7:] == 'formula':
            value = XMLScalarNode.xmlGetChildString('value')
        else:
            value = XMLScalarNode.xmlGetChildDouble('value')

        if value in ['', None]:
            if choice[-7:] == 'formula':
                value = self.getDefaultScalarFormula(choice,
                                                     self.spm.getScalarLabelByName(Scalar))
            else:
                value = self.__defaultValues()['scalar']

            self.setScalarValue(fieldId, Scalar, value)

        return value


    @Variables.undoLocal
    def setScalarValue(self, fieldId, Scalar, value):
        """
        Set non condensable variable for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId, strict_check=False)

        XMLScalarNode = self.boundNode.xmlInitNode('variable', 'choice', name=Scalar)

        XMLScalarNode.xmlSetData('value', str(value))


    @Variables.noUndo
    def getWallModel(self, fieldId):
        """
        Get the wall model for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        XMLNode = self.boundNode.xmlInitChildNode('wall_model', field_id = fieldId)

        choice = XMLNode['model']
        field = MainFieldsModel(self.case).getFieldFromId(fieldId)

        if not choice:
            if field.flow_type == "continuous":
                choice = 'friction'
            elif field.phase == "gas":
                choice = 'dvr_dn'
            elif field.phase == "liquid":
                choice = 'droplet_friction'
            elif field.phase == "solid":
                choice = 'du2_dn'
            self.setWallModel(fieldId, choice)
        return choice


    @Variables.undoLocal
    def setWallModel(self, fieldId, mdl):
        """
        Set the wall model for field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        Model().isInList(mdl, self._wallModel)
        XMLNode = self.boundNode.xmlInitChildNode('wall_model', field_id = fieldId)
        XMLNode['model'] = mdl


    @Variables.undoLocal
    def setConjugateHeatTransferCoupling(self, syrthes_name):
        print("toto")


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class InletBoundaryTestCase(ModelTest):
    """
    """
    def checkInletBoundaryInstantiation(self):
        """Check whether the InletBoundaryModel class could be instantiated"""
        model = None
        model = Boundary('inlet','BC_1',self.case)
        assert model != None, 'Could not instantiate InletBoundary'


    def checkGetandSetVelocityChoice(self):
        """Check whether the InletBoundaryModel class could set and get VelocityChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setVelocityChoice('1','flow1')
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="flow1" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <flow1>
                                                 1
                                         </flow1>
                                 </velocity>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set VelocityChoice'
        assert mdl.getVelocityChoice('1') == 'flow1',\
            'Could not get VelocityChoice'


    def checkGetandSetVelocity(self):
        """Check whether the InletBoundaryModel class could set and get Velocity"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setVelocity('1',-2.15)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 -2.15
                                         </norm>
                                 </velocity>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Velocity'
        assert mdl.getVelocity('1') == -2.15,\
            'Could not get Velocity'


    def checkGetandSetDirectionChoice(self):
        """Check whether the InletBoundaryModel class could set and get DirectionChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setDirectionChoice('1','normal')
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="normal">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set DirectionChoice'
        assert mdl.getDirectionChoice('1') == 'normal',\
            'Could not get DirectionChoice'


    def checkGetandSetDirection(self):
        """Check whether the InletBoundaryModel class could set and get Direction"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setDirection('1','direction_x',0.65)
        mdl.setDirection('1','direction_y',-0.0001)
        mdl.setDirection('1','direction_z',15654000.)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0.65
                                         </direction_x>
                                         <direction_y>
                                                 -0.0001
                                         </direction_y>
                                         <direction_z>
                                                 1.5654e+07
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Direction'
        assert mdl.getDirection('1','direction_x') == 0.65,\
            'Could not get Direction x'
        assert mdl.getDirection('1','direction_y') == -0.0001,\
            'Could not get Direction y'
        assert mdl.getDirection('1','direction_z') == 15654000,\
            'Could not get Direction z'


    def checkGetandSetTurbulenceChoice(self):
        """Check whether the InletBoundaryModel class could set and get TurbulenceChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setTurbulenceChoice('1','turbulent_intensity')
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <turbulence choice="turbulent_intensity">
                                         <hydraulic_diameter>
                                                 1
                                         </hydraulic_diameter>
                                         <turbulent_intensity>
                                                 2
                                         </turbulent_intensity>
                                 </turbulence>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set TurbulenceChoice'
        assert mdl.getTurbulenceChoice('1') == 'turbulent_intensity',\
            'Could not get TurbulenceChoice'


    def checkGetandSetHydraulicDiameter(self):
        """Check whether the InletBoundaryModel class could set and get HydraulicDiameter"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setTurbulenceChoice('1','turbulent_intensity')
        mdl.setHydraulicDiameter('1',0.2345)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <turbulence choice="turbulent_intensity">
                                         <hydraulic_diameter>
                                                 0.2345
                                         </hydraulic_diameter>
                                         <turbulent_intensity>
                                                 2
                                         </turbulent_intensity>
                                 </turbulence>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set HydraulicDiameter'
        assert mdl.getHydraulicDiameter('1') == 0.2345,\
            'Could not get HydraulicDiameter'


    def checkGetandSetTurbulentIntensity(self):
        """Check whether the InletBoundaryModel class could set and get TurbulentIntensity"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setTurbulenceChoice('1','turbulent_intensity')
        mdl.setTurbulentIntensity('1',45.26)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <turbulence choice="turbulent_intensity">
                                         <hydraulic_diameter>
                                                 1
                                         </hydraulic_diameter>
                                         <turbulent_intensity>
                                                 45.26
                                         </turbulent_intensity>
                                 </turbulence>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set TurbulentIntensity'
        assert mdl.getTurbulentIntensity('1') == 45.26,\
            'Could not get TurbulentIntensity'


    def checkGetandSetEnthalpyChoice(self):
        """Check whether the InletBoundaryModel class could set and get EnthalpyChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setEnthalpyChoice('1','timp_K')
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <variable choice="timp_K" name="enthalpy"/>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set EnthalpyChoice'
        assert mdl.getEnthalpyChoice('1') == 'timp_K',\
            'Could not get EnthalpyChoice'


    def checkGetandSetEnthalpy(self):
        """Check whether the InletBoundaryModel class could set and get Enthalpy"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setEnthalpyChoice('1','flux')
        mdl.setEnthalpy('1',875.2)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <variable choice="flux" name="enthalpy">
                                         <value>
                                                 875.2
                                         </value>
                                 </variable>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Enthalpy'
        assert mdl.getEnthalpy('1') == 875.2,\
            'Could not get Enthalpy'


    def checkGetandSetFraction(self):
        """Check whether the InletBoundaryModel class could set and get Fraction"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setFraction('1',12.5)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <variable choice="dirichlet" name="volume_fraction">
                                         <value>
                                                 12.5
                                         </value>
                                 </variable>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Fraction'
        assert mdl.getFraction('1') == 12.5,\
            'Could not get Fraction'


    def checkGetandSetDiameter(self):
        """Check whether the InletBoundaryModel class could set and get Diameter"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setDiameter('1',0.0045)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <variable choice="dirichlet" name="diameter">
                                         <value>
                                                 0.0045
                                         </value>
                                 </variable>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Diameter'
        assert mdl.getDiameter('1') == 0.0045,\
            'Could not get Diameter'


    def checkGetandSetNonCondensableValue(self):
        """Check whether the InletBoundaryModel class could set and get NonCondensableValue"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField("2", "field2", 'dispersed', 'gas', 'on', 'on', 'off', 2)
        from code_saturne.model.NonCondensableModel import NonCondensableModel
        NonCondensableModel(self.case).addNonCondensable()
        mdl = Boundary('inlet','BC_1',self.case)
        mdl.setNonCondensableValue('1','H2',16.5)
        doc = '''<boundary_conditions>
                         <inlet field_id="1" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                                 <variable choice="dirichlet" name="H2">
                                         <value>
                                                 16.5
                                         </value>
                                 </variable>
                         </inlet>
                         <inlet field_id="2" label="BC_1">
                                 <velocity choice="norm" direction="coordinates">
                                         <direction_x>
                                                 0
                                         </direction_x>
                                         <direction_y>
                                                 0
                                         </direction_y>
                                         <direction_z>
                                                 0
                                         </direction_z>
                                         <norm>
                                                 1
                                         </norm>
                                 </velocity>
                         </inlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set NonCondensableValue'
        assert mdl.getNonCondensableValue('1','H2') == 16.5,\
            'Could not get NonCondensableValue'


def suite1():
    testSuite = unittest.makeSuite(InletBoundaryTestCase, "check")
    return testSuite


def runTest1():
    print("InletBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite1())

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class OutletBoundaryTestCase(ModelTest):
    """
    """
    def checkOutletBoundaryInstantiation(self):
        """Check whether the OutletBoundaryModel class could be instantiated"""
        model = None
        model = Boundary('outlet','',self.case)
        assert model != None, 'Could not instantiate OutletBoundary'


    def checkGetandSetReferencePressure(self):
        """Check whether the OutletBoundaryModel class could set and get ReferencePressure"""
        mdl = Boundary('outlet','',self.case)
        mdl.setReferencePressure(131313.)
        doc = '''<outlet field_id="none" label="">
                         <dirichlet name="pressure">
                                 131313
                         </dirichlet>
                 </outlet>'''
        assert mdl.boundNode == self.xmlNodeFromString(doc),\
            'Could not set ReferencePressure'
        assert mdl.getReferencePressure() == 131313.,\
            'Could not get ReferencePressure'


    def checkGetandSetEnthalpyChoice(self):
        """Check whether the OutletBoundaryModel class could set and get EnthalpyChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('outlet','',self.case)
        mdl.setEnthalpyChoice('1','timp_K')
        doc = '''<boundary_conditions>
                         <outlet field_id="none" label="">
                                 <dirichlet name="pressure">
                                         101325
                                 </dirichlet>
                         </outlet>
                         <outlet field_id="1" label="">
                                 <variable choice="timp_K" name="enthalpy"/>
                         </outlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set EnthalpyChoice'
        assert mdl.getEnthalpyChoice('1') == 'timp_K',\
            'Could not set EnthalpyChoice'


    def checkGetandSetFraction(self):
        """Check whether the OutletBoundaryModel class could set and get Fraction"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('outlet','',self.case)
        mdl.setFraction('1',12.5)
        doc = '''<boundary_conditions>
                         <outlet field_id="none" label="">
                                 <dirichlet name="pressure">
                                         101325
                                 </dirichlet>
                         </outlet>
                         <outlet field_id="1" label="">
                                 <variable choice="dirichlet" name="volume_fraction">
                                         <value>
                                                 12.5
                                         </value>
                                 </variable>
                         </outlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Fraction'
        assert mdl.getFraction('1') == 12.5,\
            'Could not set Fraction'


    def checkGetandSetEnthalpy(self):
        """Check whether the OutletBoundaryModel class could set and get Enthalpy"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('outlet','',self.case)
        mdl.setEnthalpyChoice('1','flux')
        mdl.setEnthalpy('1',45.26)
        doc = '''<boundary_conditions>
                         <outlet field_id="none" label="">
                                 <dirichlet name="pressure">
                                         101325
                                 </dirichlet>
                         </outlet>
                         <outlet field_id="1" label="">
                                 <variable choice="flux" name="enthalpy">
                                         <value>
                                                 45.26
                                         </value>
                                 </variable>
                         </outlet>
                 </boundary_conditions>'''
        assert mdl.getXMLBoundaryConditionNode() == self.xmlNodeFromString(doc),\
            'Could not set Enthalpy'
        assert mdl.getEnthalpy('1') == 45.26,\
            'Could not set Enthalpy'


def suite2():
    testSuite = unittest.makeSuite(OutletBoundaryTestCase, "check")
    return testSuite


def runTest2():
    print("OutletBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite2())


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class SymmetryBoundaryTestCase(ModelTest):
    """
    """
    def checkSymmetryBoundaryInstantiation(self):
        """Check whether the SymmetryBoundaryModel class could be instantiated"""
        MainFieldsModel(self.case).addField()
        model = None
        model = Boundary('symmetry','',self.case)
        assert model != None, 'Could not instantiate SymmetryBoundary'


def suite3():
    testSuite = unittest.makeSuite(SymmetryBoundaryTestCase, "check")
    return testSuite


def runTest3():
    print("SymmetryBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite3())


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class WallBoundaryTestCase(ModelTest):
    """
    """
    def checkWallBoundaryInstantiation(self):
        """Check whether the WallBoundaryModel class could be instantiated"""
        model = None
        model = Boundary('wall','',self.case)
        assert model != None, 'Could not instantiate WallBoundary'


    def checkGetandSetEnthalpyChoice(self):
        """Check whether the WallBoudaryModel class could set and get EnthalpyChoice"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('wall','',self.case)
        mdl.setEnthalpyChoice('none','flux')
        doc = '''<wall field_id="none" label="">
                         <variable choice="flux" name="enthalpy"/>
                 </wall>'''
        assert mdl.boundNode == self.xmlNodeFromString(doc),\
            'Could not set EnthalpyChoice'
        assert mdl.getEnthalpyChoice('none') == 'flux',\
            'Could not get EnthalpyChoice'


    def checkGetandSetEnthalpy(self):
        """Check whether the WallBoudaryModel class could set and get Enthalpy"""
        MainFieldsModel(self.case).addField()
        mdl = Boundary('wall','',self.case)
        mdl.setEnthalpyChoice('none','flux')
        mdl.setEnthalpy('none',45.23)
        doc = '''<wall field_id="none" label="">
                         <variable choice="flux" name="enthalpy">
                                 <value>
                                         45.23
                                 </value>
                         </variable>
                 </wall>'''
        assert mdl.boundNode == self.xmlNodeFromString(doc),\
            'Could not set Enthalpy'
        assert mdl.getEnthalpy('none') == 45.23,\
            'Could not get Enthalpy'


def suite4():
    testSuite = unittest.makeSuite(WallBoundaryTestCase, "check")
    return testSuite


def runTest4():
    print("WallBoundaryTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite4())
