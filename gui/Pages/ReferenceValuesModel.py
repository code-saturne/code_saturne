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
- ReferenceValuesModel
- ReferenceValuesTestCase
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
from code_saturne.Base.XMLmodel import ModelTest
from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
from code_saturne.Pages.GasCombustionModel import GasCombustionModel
from code_saturne.Pages.ElectricalModel import ElectricalModel
from code_saturne.Pages.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.Pages.CompressibleModel import CompressibleModel
from code_saturne.Pages.FluidCharacteristicsModel import FluidCharacteristicsModel

#-------------------------------------------------------------------------------
# Reference values model class
#-------------------------------------------------------------------------------

class ReferenceValuesModel(Model):
    """
    Manage the input/output markups in the xml doc about Pressure
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.node_models    = self.case.xmlGetNode('thermophysical_models')
        self.node_reference = self.node_models.xmlInitNode('reference_values')
        self.node_veloce    = self.node_models.xmlGetNode('velocity_pressure')
        self.node_coal      = self.node_models.xmlGetNode('solid_fuels', 'model')
        self.node_gas       = self.node_models.xmlGetNode('gas_combustion',  'model')
        self.node_joule     = self.node_models.xmlGetNode('joule_effect',  'model')
        self.node_atmo      = self.node_models.xmlGetNode('atmospheric_flows',  'model')
        self.node_comp      = self.node_models.xmlGetNode('compressible',  'model')


    def defaultValues(self):
        """
        Return reference values by default
        """
        default = {}
        default['reference_pressure']    = 1.01325e+5
        default['reference_velocity']    = 1.0
        default['length_choice']         = 'automatic'
        default['reference_length']      = 1.0
        default['reference_temperature'] = 293.15
        default['fuel_temperature']      = 436.
        default['oxydant_temperature']   = 353.
        if (self.getParticularPhysical() == "atmo" or
            self.getParticularPhysical() == "gas"):
            default['reference_temperature'] = 293.15
        if (self.getParticularPhysical() == "off" and
            FluidCharacteristicsModel(self.case).getMaterials() == "user_material"):
            default['reference_temperature'] = 293.15
        # molar mass for dry air
        default['reference_mass_molar'] = 28.966e-3

        return default


    @Variables.undoLocal
    def setPressure(self, value):
        """
        Set value of reference pressure into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_reference.xmlSetData('pressure',value)


    @Variables.noUndo
    def getPressure(self):
        """
        Return the value of reference pressure.
        """
        value = self.node_reference.xmlGetDouble('pressure')
        if value == None:
            value = self.defaultValues()['reference_pressure']
            self.setPressure(value)

        return value


    @Variables.undoLocal
    def setVelocity(self, value):
        """
        Set value of reference velocity into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_reference.xmlSetData('velocity',value)


    @Variables.noUndo
    def getVelocity(self):
        """
        Return the value of reference velocity.
        """
        value = self.node_reference.xmlGetDouble('velocity')
        if value == None:
            value = self.defaultValues()['reference_velocity']
            self.setVelocity(value)

        return value


    @Variables.undoLocal
    def setLengthChoice(self, choice):
        """
        Set the Length choice.
        """
        self.isInList(choice, ['automatic','prescribed'])

        node_init = self.node_reference.xmlInitNode('length')
        node_init['choice'] = choice
        if choice == 'automatic':
            self.node_reference.xmlRemoveChild('length')


    @Variables.noUndo
    def getLengthChoice(self):
        """
        Get the Length choice.
        """
        node_init = self.node_reference.xmlInitNode('length')
        choice = node_init['choice']
        if choice == None:
            choice = self.defaultValues()['length_choice']
            self.setLengthChoice(choice)
        return choice


    @Variables.undoLocal
    def setLength(self, value):
        """
        Set value of reference length into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_reference.xmlSetData('length',value)


    @Variables.noUndo
    def getLength(self):
        """
        Return the value of reference length.
        """
        value = self.node_reference.xmlGetDouble('length')
        if value == None:
            value = self.defaultValues()['reference_length']
            self.setLength(value)

        return value


    @Variables.undoLocal
    def setTemperature(self, value):
        """
        Set reference temperature.
        """
        self.isGreater(value, 0.0)
        self.node_reference.xmlSetData('temperature', value)


    @Variables.noUndo
    def getTemperature(self):
        """
        Get reference temperature.
        """
        value = self.node_reference.xmlGetDouble('temperature')
        if not value :
            value = self.defaultValues()['reference_temperature']
            self.setTemperature(value)
        return value


    @Variables.undoLocal
    def setTempOxydant(self, value):
        """
        Set reference temperature for Oxydant.
        """
        self.isGreater(value, 0.0)
        self.node_reference.xmlSetData('oxydant_temperature', value)


    @Variables.noUndo
    def getTempOxydant(self):
        """
        Get reference temperaturefor Oxydant.
        """
        value = self.node_reference.xmlGetDouble('oxydant_temperature')
        if not value :
            value = self.defaultValues()['oxydant_temperature']
            self.setTempOxydant(value)
        return value


    @Variables.undoLocal
    def setTempFuel(self, value):
        """
        Set reference temperature.
        """
        self.isGreater(value, 0.0)
        self.node_reference.xmlSetData('fuel_temperature', value)


    @Variables.noUndo
    def getTempFuel(self):
        """
        Get reference temperature.
        """
        value = self.node_reference.xmlGetDouble('fuel_temperature')
        if not value :
            value = self.defaultValues()['fuel_temperature']
            self.setTempFuel(value)
        return value


    @Variables.undoLocal
    def setMassemol(self, value):
        """
        Set reference molar mass.
        """
        self.isGreater(value, 0.0)
        self.node_reference.xmlSetData('mass_molar', value)


    @Variables.noUndo
    def getMassemol(self):
        """
        Get reference molar mass.
        """
        value = self.node_reference.xmlGetDouble('mass_molar')
        if not value :
            value = self.defaultValues()['reference_mass_molar']
            self.setMassemol(value)
        return value


    @Variables.noUndo
    def getParticularPhysical(self):
        """
        Get model for set temperature for relative model
        """
        model = 'off'

        coalModel = CoalCombustionModel(self.case).getCoalCombustionModel()
        gasModel = GasCombustionModel(self.case).getGasCombustionModel()
        jouleModel = ElectricalModel(self.case).getElectricalModel()
        atmoModel = AtmosphericFlowsModel(self.case).getAtmosphericFlowsModel()
        compModel = CompressibleModel(self.case).getCompressibleModel()

        if coalModel != 'off':
            model = "coal"
        elif gasModel != 'off':
            model = "gas"
        elif jouleModel != 'off':
            model = "joule"
        elif atmoModel != 'off':
            model = "atmo"
        elif compModel != 'off':
            model = "comp"

        return model


#-------------------------------------------------------------------------------
# ReferenceValuesModel test case
#-------------------------------------------------------------------------------

class ReferenceValuesTestCase(ModelTest):
    """
    """
    def checkReferenceValuesInstantiation(self):
        """Check whether the ReferenceValuesModel class could be instantiated"""
        model = None
        model = ReferenceValuesModel(self.case)
        assert model != None, 'Could not instantiate ReferenceValuesModel'

    def checkGetandSetPressure(self):
        """Check whether the ReferenceValuesModel class could be set and get Pressure"""
        mdl = ReferenceValuesModel(self.case)
        mdl.setPressure(13e+5)

        doc = """<velocity_pressure>
                    <variable label="Pressure" name="pressure">
                            <reference_pressure>1.3e+06</reference_pressure>
                    </variable>
                    <variable label="Velocity" name="velocity"/>
                    <property label="total_pressure" name="total_pressure"/>
                    <property label="Yplus" name="yplus" support="boundary"/>
                    <property label="Stress" name="stress" support="boundary"/>
                </velocity_pressure>"""
        assert mdl.node_veloce == self.xmlNodeFromString(doc),\
            'Could not set pressure ReferenceValuesModel'
        assert mdl.getPressure() == 13e+5,\
            'Could not get pressure ReferenceValuesModel'

    def checkGetandSetTemperature(self):
        """Check whether the ReferenceValuesModel class could be set and get Temperature"""
        mdl = ReferenceValuesModel(self.case)
        from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel
        CoalCombustionModel(self.case).setCoalCombustionModel('homogeneous_fuel')
        del CoalCombustionModel
        mdl.setTemperature(55.5)

        doc = """<solid_fuels model="homogeneous_fuel">
                    <variable label="Enthalpy" name="enthalpy" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Np_01" name="n_p_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ch01" name="x_p_coal_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ck_01" name="x_p_char_" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Xp_Ent01" name="x_p_h_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_mv1_01" name="fr_mv1_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_mv2_01" name="fr_mv2_01" type="model"><flux_reconstruction status="off"/></variable>
                    <variable label="Fr_HET" name="het_fraction" type="model"><flux_reconstruction status="off"/></variable>
                    <property label="Temp_GAZ" name="t_gas"/>
                    <property label="ROM_GAZ" name="rho_gas"/>
                    <property label="YM_CHx1m" name="ym_chx1m"/>
                    <property label="YM_CHx2m" name="ym_chx2m"/>
                    <property label="YM_CO" name="ym_co"/>
                    <property label="YM_O2" name="ym_o2"/>
                    <property label="YM_CO2" name="ym_co2"/>
                    <property label="YM_H2O" name="ym_h2o"/>
                    <property label="YM_N2" name="ym_n2"/>
                    <property label="XM" name="xm"/>
                    <property label="Temp_CP01" name="t_coal01"/>
                    <property label="Frm_CP01" name="w_solid_coal01"/>
                    <property label="Rho_CP01" name="rho_coal01"/>
                    <property label="Dia_CK01" name="diameter_coal01"/>
                    <property label="Ga_DCH01" name="dissapear_rate_coal01"/>
                    <property label="Ga_DV101" name="m_transfer_v1_coal01"/>
                    <property label="Ga_DV201" name="m_transfer_v2_coal01"/>
                    <property label="Ga_HET01" name="het_ts_coal01"/>
                    <property label="ntLuminance_4PI" name="ntLuminance_4PI"/>
                    <reference_temperature>55.5</reference_temperature>
                 </solid_fuels>"""
        assert mdl.node_coal == self.xmlNodeFromString(doc),\
            'Could not set temperature ReferenceValuesModel'
        assert mdl.getTemperature() == 55.5,\
            'Could not get temperature ReferenceValuesModel'

    def checkGetandSsetMassemol(self):
        """Check whether the ReferenceValuesModel class could be set and get Molar mass"""
        mdl = ReferenceValuesModel(self.case)
        from code_saturne.Pages.GasCombustionModel import GasCombustionModel
        GasCombustionModel(self.case).setGasCombustionModel('ebu')
        del GasCombustionModel
        mdl.setMassemol(50.8e-3)
        doc = """<gas_combustion model="ebu">
                    <reference_mass_molar>0.0508</reference_mass_molar>
                 </gas_combustion>"""
        assert mdl.node_gas == self.xmlNodeFromString(doc),\
            'Could not set molar mass ReferenceValuesModel'
        assert mdl.getMassemol() == 50.8e-3,\
            'Could not get molar mass ReferenceValuesModel'


def suite():
    testSuite = unittest.makeSuite(ReferenceValuesTestCase, "check")
    return testSuite

def runTest():
    print("ReferenceValuesTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
