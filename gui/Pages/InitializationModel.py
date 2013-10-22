# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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

import sys, string, unittest
from math import pow

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLmodel import XMLmodel, ModelTest
from Base.XMLvariables import Model, Variables
from Pages.TurbulenceModel import TurbulenceModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.LocalizationModel import LocalizationModel
from Pages.CompressibleModel import CompressibleModel
from Pages.ElectricalModel import ElectricalModel

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
        if CompressibleModel(self.case).getCompressibleModel() != 'off':
            self.node_comp = self.models.xmlGetNode('compressible_model', 'model')

        self.VelocityList  = ('velocity_U', 'velocity_V', 'velocity_W')
        self.Turb_var_List = ('turb_k', 'turb_eps',
                              'component_R11', 'component_R22', 'component_R33',
                              'component_R12', 'component_R13', 'component_R23',
                              'turb_phi', 'turb_al', 'turb_omega', 'turb_nusa',
                              'turb_alpha')

        self.turb = TurbulenceModel(self.case)
        self.turbulenceModes = ('formula',
                                'reference_value')
        self.node_scalartherm = self.models.xmlGetNode('thermal_scalar')

        self.thermalModel = ('TempC',
                             'TempK',
                             'Enthalpy',
                             'PotTemp',
                             'LiqPotTemp')


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
    def getDefaultTurbFormula(self, turb_model):
        self.isInList(turb_model,self.turb.turbulenceModels())
        if turb_model in ('k-epsilon', 'k-epsilon-PL'):
            formula = """cmu = 0.42;
k = 1.5*(0.02*uref)^2;
eps = k^1.5*cmu/almax;"""
        elif turb_model in ('Rij-epsilon', 'Rij-SSG'):
            formula = """trii   = (0.02*uref)^2;
cmu = 0.42;
r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.;
r23 = 0.;
k = 0.5*(r11+r22+r33);
eps = k^1.5*cmu/almax;"""
        elif turb_model == 'Rij-EBRSM':
            formula = """trii   = (0.02*uref)^2;
cmu = 0.42;
r11 = trii;
r22 = trii;
r33 = trii;
r12 = 0.;
r13 = 0.;
r23 = 0.;
k = 0.5*(r11+r22+r33);
eps = k^1.5*cmu/almax;
alpha = 1.;"""
        elif turb_model == 'v2f-BL-v2/k':
            formula = """cmu = 0.22;
k = 1.5*(0.02*uref)^2;
eps = k^1.5*cmu/almax;
phi = 2./3.;
al = 0.;"""
        elif turb_model == 'k-omega-SST':
            formula = """k = 1.5*(0.02*uref)^2;
omega = k^0.5/almax;"""
        elif turb_model == 'Spalart-Allmaras':
            formula = """nusa = (cmu * k)/eps;"""
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
        n.xmlSetTextNode(formula)


    @Variables.undoLocal
    def setVelocityFormula(self, zone, formula):
        """
        Public method.
        Set the formula for the velocity.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlInitNode('initialization')
        n = node.xmlInitChildNode('formula', zone_id=zone)
        n.xmlSetTextNode(formula)


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


    @Variables.undoLocal
    def setThermalFormula(self, zone, scalar, formula):
        """
        Public method.
        Set the formula for tharmal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, self.thermalModel)
        node = self.node_scalartherm.xmlGetNode('scalar', label = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getThermalFormula(self, zone, scalar):
        """
        Public method.
        Return the formula for thermal scalars.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, self.thermalModel)
        node = self.node_scalartherm.xmlGetNode('scalar', label = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.noUndo
    def getDensityStatus(self, zone):
        """
        Return status of Density for the initialisation
        """
        node = self.node_comp.xmlGetNode('property', name = 'Rho')
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
        node = self.node_comp.xmlGetNode('property', name = 'Rho')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    @Variables.noUndo
    def getTemperatureStatus(self, zone):
        """
        Return status of Temperature for the initialisation
        """
        node = self.node_comp.xmlGetNode('scalar', name = 'TempK')
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
        node = self.node_comp.xmlGetNode('scalar', name = 'TempK')
        n = node.xmlInitNode('formula', 'status', zone_id = zone)
        n['status'] = status


    @Variables.noUndo
    def getEnergyStatus(self, zone):
        """
        Return status of total energy for the initialisation
        """
        node = self.node_comp.xmlGetNode('scalar', name = 'EnergieT')
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
        node = self.node_comp.xmlGetNode('scalar', name = 'EnergieT')
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


    @Variables.undoLocal
    def setPressureFormula(self, zone, formula):
        """
        Public method.
        Set the formula for Pressure.
        """
        self.__verifyZone(zone)
        node = self.node_veloce.xmlGetNode('variable', name = 'pressure')

        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlGetNode('formula', zone_id = zone)
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
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.undoLocal
    def setDensityFormula(self, zone, formula):
        """
        Public method.
        Set the formula for density.
        """
        self.__verifyZone(zone)
        node = self.node_comp.xmlGetNode('property', name = 'Rho')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
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
        node = self.node_comp.xmlGetNode('property', name = 'Rho')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.undoLocal
    def setTemperatureFormula(self, zone, formula):
        """
        Public method.
        Set the formula for temperature.
        """
        self.__verifyZone(zone)
        node = self.node_comp.xmlGetNode('scalar', name = 'TempK')


        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
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
        node = self.node_comp.xmlGetNode('scalar', name = 'TempK')


        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        return formula


    @Variables.undoLocal
    def setEnergyFormula(self, zone, formula):
        """
        Public method.
        Set the formula for totale energy.
        """
        self.__verifyZone(zone)
        node = self.node_comp.xmlGetNode('scalar', name = 'EnergieT')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
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
        node = self.node_comp.xmlGetNode('scalar', name = 'EnergieT')
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
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


    @Variables.undoLocal
    def setSpeciesFormula(self, zone, species, formula):
        """
        Public method.
        Set the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarLabelsList())
        node = self.node_userscalar.xmlGetNode('scalar', label = str(species))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getSpeciesFormula(self, zone, species):
        """
        Public method.
        Return the formula for a turbulent variable.
        """
        self.__verifyZone(zone)
        self.isInList(species, DefineUserScalarsModel(self.case).getUserScalarLabelsList())
        node = self.node_userscalar.xmlGetNode('scalar', label = str(species))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        if not formula:
            formula = str(species)+""" = 0;\n"""
            self.setSpeciesFormula(zone, species, formula)

        return formula


    @Variables.undoLocal
    def setMeteoFormula(self, zone, scalar, formula):
        """
        Public method.
        Set the formula for a meteo variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getMeteoScalarsList())
        node_atmo = self.models.xmlGetNode('atmospheric_flows')
        node = node_atmo.xmlGetNode('scalar', label = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)
        n = node.xmlInitChildNode('formula', zone_id=zone)
        n.xmlSetTextNode(formula)


    @Variables.noUndo
    def getMeteoFormula(self, zone, scalar):
        """
        Public method.
        Return the formula for a meteo variable.
        """
        self.__verifyZone(zone)
        self.isInList(scalar, DefineUserScalarsModel( self.case).getMeteoScalarsList())
        node_atmo = self.models.xmlGetNode('atmospheric_flows')
        node = node_atmo.xmlGetNode('scalar', label = str(scalar))
        if not node:
            msg = "There is an error: this node " + str(node) + "should be existed"
            raise ValueError(msg)

        formula = node.xmlGetString('formula', zone_id=zone)
        if not formula:
            formula = str(scalar)+""" = 0;\n"""
            self.setMeteoFormula(zone, scalar, formula)

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
