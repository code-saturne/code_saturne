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
This module initialize model dynamics variables and model scalars

This module contents the following classes:
- InitializationModel
- InitializationTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLmodel import XMLmodel, ModelTest
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.TurbulenceModel import TurbulenceModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.CompressibleModel import CompressibleModel
from code_saturne.model.ElectricalModel import ElectricalModel
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Variables and Scalar model initialization modelling class
#-------------------------------------------------------------------------------

class InitializationModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.models = self.case.xmlGetNode('thermophysical_models')
        self.node_userscalar = self.case.xmlGetNode('additional_scalars')
        self.node_veloce     = self.models.xmlGetNode('velocity_pressure')
        self.node_turb       = self.models.xmlGetNode('turbulence', 'model')
        self.node_therm      = self.models.xmlGetNode('thermal_scalar', 'model')
        self.node_prop       = self.case.xmlGetNode('physical_properties')
        self.node_fluid      = self.node_prop.xmlInitNode('fluid_properties')
        if CompressibleModel(self.case).getCompressibleModel() != 'off':
            self.node_comp = self.models.xmlGetNode('compressible_model', 'model')

        self.turb = TurbulenceModel(self.case)
        self.therm = ThermalScalarModel(self.case)
        self.turbulenceModes = ('formula',
                                'reference_value')
        self.node_scalartherm = self.models.xmlGetNode('thermal_scalar')
        self.notebook = NotebookModel(self.case)


    def __defaultValues(self):
        """
        Private method.
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['init_mode_turb'] = "reference_value"
        default['velocity']       = 0.0
        default['status']         = 'off'

        return default


    def __verifyZone(self, zone):
        """Private method.
        Verify if zone exists and raise ValueError if not.
        """
        self.isInt(int(zone))
        self.isInList(zone, LocalizationModel('VolumicZone', self.case).getCodeNumbersList())


    @Variables.noUndo
    def getDefaultVelocityFormula(self):
        formula = """velocity[0] = 0.;
velocity[1] = 0.;
velocity[2] = 0.;"""
        return formula


    @Variables.noUndo
    def getDefaultThermalFormula(self):
        name = ''
        formula = ''
        if self.therm.getThermalScalarModel() == "enthalpy":
            name = 'enthalpy'
            formula = """cp = 1017.24;\n""" + name + """ = 300. * cp;"""
        elif self.therm.getThermalScalarModel() == "total_energy":
            name = 'total_energy'
            formula = name + """ = 0.;"""
        elif self.therm.getThermalScalarModel() != "off":
            name = 'temperature'
            formula = name + """ = 300.;"""
        return formula


    @Variables.noUndo
    def getDefaultTurbFormula(self, turb_model):
        self.isInList(turb_model,self.turb.turbulenceModels())
        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            formula = """cmu = 0.09;
_k = 1.5*(0.02*uref)^2;
k = _k;
epsilon = _k^1.5*cmu/almax;"""
        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            formula = """trii   = (0.02*uref)^2;
cmu = 0.09;
r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.;
r23 = 0.;
k = 0.5*(r11+r22+r33);
epsilon = k^1.5*cmu/almax;"""
        elif turb_model == 'Rij-EBRSM':
            formula = """trii   = (0.02*uref)^2;
cmu = 0.09;
r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.;
r23 = 0.;
k = 0.5*(r11+r22+r33);
epsilon = k^1.5*cmu/almax;
alpha = 1.;"""
        elif turb_model == 'v2f-BL-v2/k':
            formula = """cmu = 0.22;
_k = 1.5*(0.02*uref)^2;
k = _k;
epsilon = _k^1.5*cmu/almax;
phi = 2./3.;
alpha = 1.;"""
        elif turb_model == 'k-omega-SST':
            formula = """k = 1.5*(0.02*uref)^2;
omega = k^0.5/almax;"""
        elif turb_model == 'Spalart-Allmaras':
            formula = """nu_tilda = (cmu * k)/eps;"""
        return formula


    @Variables.noUndo
    def getInitialTurbulenceChoice(self, zone):
        """
        Public method.
        Return the turbulence initialization choice.
        """
        self.__verifyZone(zone)
        node_init = self.node_turb.xmlGetNode('initialization', zone_id = zone)
        choice = ''
        if not node_init:
            choice = self.__defaultValues()['init_mode_turb']
            self.setInitialTurbulenceChoice(zone, choice)
            node_init = self.node_turb.xmlGetNode('initialization', zone_id = zone)
        choice = node_init['choice']

        return choice


    @Variables.undoLocal
    def setInitialTurbulenceChoice(self, zone, init_mode):
        """
        Public method.
        Set the initialization mode in the attribute choice.
        """
        self.__verifyZone(zone)
        self.isInList(init_mode, self.turbulenceModes)

        node_init = self.node_turb.xmlInitNode('initialization', zone_id = zone)
        node_init['choice'] = init_mode


    def getTurbFormulaComponents(self, zone, turb_model):
        """
        public method.
        Return formula components (exp, req and symbols)
        """

        exp = self.getTurbFormula(zone, turb_model)

        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            req = [('k', "turbulent energy"),
                   ('epsilon', "turbulent dissipation")]
        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            req = [('r11', "Reynolds stress R11"),
                   ('r22', "Reynolds stress R22"),
                   ('r33', "Reynolds stress R33"),
                   ('r12', "Reynolds stress R12"),
                   ('r23', "Reynolds stress R23"),
                   ('r13', "Reynolds stress R13"),
                   ('epsilon', "turbulent dissipation")]
        elif turb_model == 'Rij-EBRSM':
            req = [('r11', "Reynolds stress R11"),
                   ('r22', "Reynolds stress R22"),
                   ('r33', "Reynolds stress R33"),
                   ('r12', "Reynolds stress R12"),
                   ('r23', "Reynolds stress R23"),
                   ('r13', "Reynolds stress R13"),
                   ('epsilon', "turbulent dissipation"),
                   ('alpha', "alpha")]
        elif turb_model == 'v2f-BL-v2/k':
            req = [('k', "turbulent energy"),
                   ('epsilon', "turbulent dissipation"),
                   ('phi', "variable phi in v2f model"),
                   ('alpha', "variable alpha in v2f model")]
        elif turb_model == 'k-omega-SST':
            req = [('k', "turbulent energy"),
                   ('omega', "specific dissipation rate")]
        elif turb_model == 'Spalart-Allmaras':
            req = [('nu_tilda', "nusa")]

        sym = [('rho0', 'density (reference value)'),
               ('mu0', 'viscosity (reference value)'),
               ('cp0', 'specific heat (reference value)'),
               ('x','cell center coordinate'),
               ('y','cell center coordinate'),
               ('z','cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume'),
               ('uref','reference velocity'),
               ('almax','reference length')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.noUndo
    def getTurbFormula(self, zone, turb_model):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        node = self.node_turb.xmlInitNode('initialization', zone_id=zone)

        formula = node.xmlGetString('formula')
        if not formula:
            formula = self.getDefaultTurbFormula(turb_model)
            self.setTurbFormula(zone, formula)
        return formula


    @Variables.undoLocal
    def setTurbFormula(self, zone, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        node = self.node_turb.xmlGetNode('initialization', zone_id=zone)
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula')
        if formula != None:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()



    def getVelocityFormulaComponents(self, zone):
        exp = self.getVelocityFormula(zone)
        if not exp:
            exp = self.getDefaultVelocityFormula()

        req = [('velocity[0]', "velocity"),
               ('velocity[1]', "velocity"),
               ('velocity[2]', "velocity")]
        sym = [('uref', 'reference velocity'),
               ('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setVelocityFormula(self, zone, formula):
        """
        Public method.
        Set the formula for the velocity.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlInitNode('initialization')
        n = node.xmlInitChildNode('formula', zone_id=zone)
        if formula != None:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getVelocityFormula(self, zone):
        """
        Public method.
        Return the formula for the velocity.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlInitNode('initialization')

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    def getThermalFormulaComponents(self, zone):
        exp = self.getThermalFormula(zone)
        if not exp:
            exp = self.getDefaultThermalFormula()

        if self.therm.getThermalScalarModel() == "enthalpy":
            req = [('enthalpy', 'enthalpy')]
        elif self.therm.getThermalScalarModel() == "total_energy":
            req = [('total_energy', 'total energy')]
        else:
            req = [('temperature', 'temperature')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        # Known fields
        knf = []

        return exp, req, sym, knf


    @Variables.undoLocal
    def setThermalFormula(self, zone, formula):
        """
        Public method.
        Set the formula for thermal scalars.
        """
        self.__verifyZone(zone)
        node = self.node_scalartherm.xmlGetNode('variable')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        if formula:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getThermalFormula(self, zone):
        """
        Public method.
        Return the formula for thermal scalars.
        """
        self.__verifyZone(zone)
        node = self.node_scalartherm.xmlGetNode('variable')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.noUndo
    def getDensityStatus(self, zone):
        """
        Return status of Density for the initialisation
        """
        node = self.node_fluid.xmlGetNode('property', name = 'density')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        status = n['status']
        if not status:
            status = self.__defaultValues()['status']
            self.setDensityStatus(zone, status)
        return status


    @Variables.undoLocal
    def setDensityStatus(self, zone, status):
        """
        Put status of Density for the initialisation
        """
        self.isOnOff(status)
        node = self.node_fluid.xmlGetNode('property', name = 'density')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    @Variables.noUndo
    def getTemperatureStatus(self, zone):
        """
        Return status of Temperature for the initialisation
        """
        node = self.node_comp.xmlGetNode('variable', name = 'temperature')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        status = n['status']
        if not status:
            status = self.__defaultValues()['status']
            self.setTemperatureStatus(zone, status)
        return status


    @Variables.undoLocal
    def setTemperatureStatus(self, zone, status):
        """
        Put status of Temperature for the initialisation
        """
        self.isOnOff(status)
        node = self.node_comp.xmlGetNode('variable', name = 'temperature')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    @Variables.noUndo
    def getEnergyStatus(self, zone):
        """
        Return status of total energy for the initialisation
        """
        node = self.node_scalartherm.xmlGetNode('variable', name = 'total_energy')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        status = n['status']
        if not status:
            status = self.__defaultValues()['status']
            self.setEnergyStatus(zone, status)
        return status


    @Variables.undoLocal
    def setEnergyStatus(self, zone, status):
        """
        Put status of Energy for the initialisation
        """
        self.isOnOff(status)
        node = self.node_scalartherm.xmlGetNode('variable', name = 'total_energy')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    @Variables.noUndo
    def getPressureStatus(self, zone):
        """
        Return status of pressure for the initialisation
        """
        node = self.node_veloce.xmlGetNode('variable', name = 'pressure')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        status = n['status']
        if not status:
            status = self.__defaultValues()['status']
            self.setPressureStatus(zone, status)
        return status


    @Variables.undoLocal
    def setPressureStatus(self, zone, status):
        """
        Put status of pressure for the initialisation
        """
        self.isOnOff(status)
        node = self.node_veloce.xmlGetNode('variable', name = 'pressure')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    def getPressureFormulaComponents(self, zone):
        exp = self.getPressureFormula(zone)
        if not exp:
            exp = """p0 = 0.;
g = 9.81;
ro = 1.17862;
pressure = p0 + g * ro * z;\n"""
        req = [('pressure', 'pressure')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setPressureFormula(self, zone, formula):
        """
        Public method.
        Set the formula for Pressure.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlGetNode('variable', name = 'pressure')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id = zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getPressureFormula(self, zone):
        """
        Public method.
        Return the formula for pressure.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlGetNode('variable', name = 'pressure')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    def getHydraulicHeadFormulaComponents(self, zone):

        exp = self.getHydraulicHeadFormula(zone)
        if not exp:
            exp = """H = z;\n"""
        req = [('H', 'hydraulic head')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setHydraulicHeadFormula(self, zone, formula):
        """
        Public method.
        Set the formula for hydraulic head.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlGetNode('variable', name = 'hydraulic_head')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id = zone)
        if formula != None:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getHydraulicHeadFormula(self, zone):
        """
        Public method.
        Return the formula for hydraulic head.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlGetNode('variable', name = 'hydraulic_head')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    def getDensityFormulaComponents(self, zone):
        exp = self.getDensityFormula(zone)
        if not exp:
            exp = """density = 0;\n"""

        req = [('density', 'density')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym



    @Variables.undoLocal
    def setDensityFormula(self, zone, formula):
        """
        Public method.
        Set the formula for density.
        """
        self.__verifyZone(zone)
        node = self.node_fluid.xmlGetNode('property', name = 'density')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlGetNode('formula', zone_id = zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getDensityFormula(self, zone):
        """
        Public method.
        Return the formula for density.
        """
        self.__verifyZone(zone)
        node = self.node_fluid.xmlGetNode('property', name = 'density')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    def getTemperatureFormulaComponents(self, zone):
        exp = self.getTemperatureFormula(zone)
        if not exp:
            exp = """temperature = 0;\n"""

        req = [('temperature', 'temperature')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setTemperatureFormula(self, zone, formula):
        """
        Public method.
        Set the formula for temperature.
        """
        self.__verifyZone(zone)
        node = self.node_comp.xmlGetNode('variable', name = 'temperature')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlGetNode('formula', zone_id = zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getTemperatureFormula(self, zone):
        """
        Public method.
        Return the formula for temperature.
        """
        self.__verifyZone(zone)
        node = self.node_comp.xmlGetNode('variable', name = 'temperature')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    def getEnergyFormulaComponents(self, zone):
        exp = self.getEnergyFormula(zone)
        if not exp:
            exp = """total_energy = 0;\n"""
        req = [('total_energy', 'Energy')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setEnergyFormula(self, zone, formula):
        """
        Public method.
        Set the formula for totale energy.
        """
        self.__verifyZone(zone)
        node = self.node_scalartherm.xmlGetNode('variable', name = 'total_energy')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlGetNode('formula', zone_id = zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getEnergyFormula(self, zone):
        """
        Public method.
        Return the formula for energy.
        """
        self.__verifyZone(zone)
        node = self.node_scalartherm.xmlGetNode('variable', name = 'total_energy')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.noUndo
    def getCheckedBoxList(self,zone):
        """
        Public method.
        return a list of selected variable for initialisation
        """
        box_list = []
        status = self.getPressureStatus(zone)
        if status == 'on':
            box_list.append('Pressure')
        status = self.getDensityStatus(zone)
        if status == 'on':
            box_list.append('Density')
        status = self.getTemperatureStatus(zone)
        if status == 'on':
            box_list.append('Temperature')
        status = self.getEnergyStatus(zone)
        if status == 'on':
            box_list.append('Energy')
        return box_list


    def getSpeciesFormulaComponents(self, zone, species):
        exp = self.getSpeciesFormula(zone, species)
        name = DefineUserScalarsModel(self.case).getScalarName(species)
        if not exp:
            exp = str(name)+""" = 0;\n"""

        req = [(str(name), str(name))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        # Known fields
        knf = [(str(name), str(name)), ('rho', 'density')]

        return exp, req, sym, knf

    @Variables.undoLocal
    def setSpeciesFormula(self, zone, species, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_userscalar.xmlGetNode('variable', name = str(species))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        if formula:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getSpeciesFormula(self, zone, species):
        """
        Public method.
        Return the formula for a scalar species variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarNameList())
        node = self.node_userscalar.xmlGetNode('variable', name = str(species))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)

        return formula


    def getMeteoFormulaComponents(self, zone, scalar):
        exp = self.getMeteoFormula(zone, scalar)
        if not exp:
            exp = str(scalar)+""" = 0;\n"""
        req = [(str(scalar), str(scalar))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setMeteoFormula(self, zone, scalar, formula):
        """
        Public method.
        Set the formula for a meteo variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getMeteoScalarsNameList())
        node_atmo = self.models.xmlGetNode('atmospheric_flows')
        node = node_atmo.xmlGetNode('variable', name = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        if formula != None:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getMeteoFormula(self, zone, scalar):
        """
        Public method.
        Return the formula for a meteo variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getMeteoScalarsNameList())
        node_atmo = self.models.xmlGetNode('atmospheric_flows')
        node = node_atmo.xmlGetNode('variable', name = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + " should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)

        return formula


    def getCombustionFormulaComponents(self, zone, scalar):
        exp = self.getCombustionFormula(zone, scalar)
        if not exp:
            exp = str(scalar)+""" = 0;\n"""
        req = [(str(scalar), str(scalar))]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    @Variables.undoLocal
    def setCombustionFormula(self, zone, scalar, formula):
        """
        Set the formula for a gas combustion variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getGasCombScalarsNameList())
        node_gas = self.models.xmlGetNode('gas_combustion')
        node = node_gas.xmlGetNode('variable', name = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be present"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        if formula != None:
            n.xmlSetTextNode(formula)
        else:
            n.xmlRemoveNode()


    @Variables.noUndo
    def getCombustionFormula(self, zone, scalar):
        """
        Public method.
        Return the formula for a gas combustion variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getGasCombScalarsNameList())
        node_gas = self.models.xmlGetNode('gas_combustion')
        node = node_gas.xmlGetNode('variable', name = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + " should be present"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)

        return formula

#-------------------------------------------------------------------------------
# InitializationModel test case
#-------------------------------------------------------------------------------


class InitializationTestCase(ModelTest):
    """
    Unittest.
    """
    def checkInitializationInstantiation(self):
        """Check whether the InitializationModel class could be instantiated."""
        model = None
        model = InitializationModel(self.case)
        assert model != None, 'Could not instantiate InitializationModel'


def suite():
    testSuite = unittest.makeSuite(InitializationTestCase, "check")
    return testSuite


def runTest():
    print("InitializationTestCase - OK !!!!")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
