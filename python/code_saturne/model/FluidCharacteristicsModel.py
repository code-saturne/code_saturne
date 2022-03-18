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
This module defines the Page for the physical properties of the fluid.
These properties can be reference value or initial value

This module contens the following classes and function:
- FluidCaracteristicsModel
- FluidCaracteristicsModelTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import XMLmodel, ModelTest

from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

from code_saturne.model.EosWrapper import eosWrapper

# Perfect gases to exclude from EOS choices

eos_excl = ["Argon", "Nitrogen", "Hydrogen", "Oxygen", "Helium", "Air"]

#-------------------------------------------------------------------------------
# Coolprop
#-------------------------------------------------------------------------------

coolprop_fluids = []
coolprop_warn = False

def load_coolprop_fluids(config):
    """
    Try to load list of coolprop fluids.
    """
    global coolprop_fluids
    global coolprop_warn

    try:
        import sys
        sys.path.insert(0, config.libs['coolprop'].flags['pythonpath'])
        import CoolProp
        sys.path.pop(0)
        self.coolprop_fluids = []
        for f in CoolProp.__fluids__:
            coolprop_fluids.append(f)
        coolprop_fluids.sort()

    except Exception:  # CoolProp might be available but not its Python bindings

        if config.libs['coolprop'].have != "gui_only":
            """
            import traceback
            exc_info = sys.exc_info()
            bt = traceback.format_exception(*exc_info)
            for l in bt:
            print(l)
            del exc_info
            print("Warning: CoolProp Python bindings not available or usable")
            print("         list of fluids based on CoolProp 5.1.1")
            """
            pass
        else:
            coolprop_warn = True

        coolprop_fluids = ['1-Butene', 'Acetone', 'Air', 'Ammonia', 'Argon',
                           'Benzene', 'CarbonDioxide', 'CarbonMonoxide',
                           'CarbonylSulfide', 'CycloHexane', 'CycloPropane',
                           'Cyclopentane', 'D4', 'D5', 'D6', 'Deuterium',
                           'DimethylCarbonate', 'DimethylEther', 'Ethane',
                           'Ethanol', 'EthylBenzene', 'Ethylene', 'Fluorine',
                           'HFE143m', 'HeavyWater', 'Helium', 'Hydrogen',
                           'HydrogenSulfide', 'IsoButane', 'IsoButene',
                           'Isohexane', 'Isopentane', 'Krypton', 'MD2M', 'MD3M',
                           'MD4M', 'MDM', 'MM', 'Methane', 'Methanol',
                           'MethylLinoleate', 'MethylLinolenate', 'MethylOleate',
                           'MethylPalmitate', 'MethylStearate', 'Neon',
                           'Neopentane', 'Nitrogen', 'NitrousOxide', 'Novec649',
                           'OrthoDeuterium', 'OrthoHydrogen', 'Oxygen',
                           'ParaDeuterium', 'ParaHydrogen', 'Propylene',
                           'Propyne', 'R11', 'R113', 'R114', 'R115', 'R116',
                           'R12', 'R123', 'R1233zd(E)', 'R1234yf', 'R1234ze(E)',
                           'R1234ze(Z)', 'R124', 'R125', 'R13', 'R134a', 'R13I1',
                           'R14', 'R141b', 'R142b', 'R143a', 'R152A', 'R161',
                           'R21', 'R218', 'R22', 'R227EA', 'R23', 'R236EA',
                           'R236FA', 'R245fa', 'R32', 'R365MFC', 'R404A',
                           'R407C', 'R41', 'R410A', 'R507A', 'RC318', 'SES36',
                           'SulfurDioxide', 'SulfurHexafluoride', 'Toluene',
                           'Water', 'Xenon', 'cis-2-Butene', 'm-Xylene',
                           'n-Butane', 'n-Decane', 'n-Dodecane', 'n-Heptane',
                           'n-Hexane', 'n-Nonane', 'n-Octane', 'n-Pentane',
                           'n-Propane', 'n-Undecane', 'o-Xylene', 'p-Xylene',
                           'trans-2-Butene']

#-------------------------------------------------------------------------------
# Properties ref values
#-------------------------------------------------------------------------------

class propertiesIdInfo(object):
    """
    Class which handles physical properties id information.
    """

    def __init__(self):
        """
        Constructor
        """

        self._info = {}

    def add_property(self, name, output, ref_tag, unit=''):
        """
        Add a new property.
        """

        self._info[name] = {"name":name,
                            "output":output,
                            "ref_tag":ref_tag,
                            "unit":unit}

    def get_property(self, name):
        """
        Returns the dictionnary related to a given property based on its name.
        """

        prop = None
        if name in self._info.keys():
            prop = self._info[name]

        return prop

    def get_property_list(self):
        """
        Returns a list containing the names of all properties.
        """

        return list(self._info.keys())

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class FluidCharacteristicsModel(Variables, Model):
    """
    Class to manipulate Molecular Properties in xml file.
    """
    def __init__(self, case):
        """FluidCharacteristicsModel Constuctor."""
        self.case = case
        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_prop   = self.case.xmlGetNode('physical_properties')
        self.node_fluid  = self.node_prop.xmlInitNode('fluid_properties')
        self.node_comp   = self.node_models.xmlInitNode('compressible_model', 'model')
        self.node_gas    = self.node_models.xmlInitNode('gas_combustion',     'model')
        self.node_coal   = self.node_models.xmlInitNode('solid_fuels',        'model')
        self.node_hgn    = self.node_models.xmlInitNode('hgn_model',          'model')

        # Info on available libraries

        self.mask_builtin   = 1 << 0
        self.mask_CoolProp  = 1 << 1
        self.mask_EOS       = 1 << 2
        self.mask_freesteam = 1 << 3

        self.tables = 0

        # Properties id cards
        self.properties_info = propertiesIdInfo()
        self.properties_info.add_property(name="density",
                                          output="Density",
                                          ref_tag="rho0",
                                          unit="kg/m^3")
        self.properties_info.add_property(name="molecular_viscosity",
                                          output="Molecular viscosity",
                                          ref_tag="mu0",
                                          unit="Pa.s")
        self.properties_info.add_property(name="surface_tension",
                                          output="Surface tension",
                                          ref_tag="sigma0",
                                          unit="N/m")
        self.properties_info.add_property(name="specific_heat",
                                          output="Sepcific heat",
                                          ref_tag="cp0",
                                          unit="J/kg/K")
        self.properties_info.add_property(name="thermal_conductivity",
                                          output="Thermal conductivity",
                                          ref_tag="lambda0",
                                          unit="W/m/K")
        self.properties_info.add_property(name="volume_viscosity",
                                          output="Volume viscosity",
                                          ref_tag="viscv0",
                                          unit="kg/(m.s)")

        self.lib_properties = {}
        self.lib_properties['user_material'] = self.mask_builtin

        self.eos = eosWrapper()
        if self.eos.isActive():
            self.tables += self.mask_EOS
            # suppress perfect gas
            fls = self.eos.getListOfFluids()
            for fli in fls:
                if fli in eos_excl:
                    continue
                if fli not in self.lib_properties.keys():
                    self.lib_properties[fli] = self.mask_EOS
                else:
                    self.lib_properties[fli] += self.mask_EOS

        try:
            cfg = case.case['package'].config

            if cfg.libs['freesteam'].have:
                self.tables += self.mask_freesteam
                fli = 'Water'
                if fli not in self.lib_properties.keys():
                    self.lib_properties[fli] = self.mask_freesteam
                else:
                    self.lib_properties[fli] += self.mask_freesteam

            if cfg.libs['coolprop']:
                if not coolprop_fluids:
                    load_coolprop_fluids(cfg)
                if coolprop_fluids:
                    self.tables += self.mask_CoolProp

        except Exception:  # if case/package not available (should not happen)
            print("Warning: package configuration not available")
            pass

        for fli in coolprop_fluids:
            if fli not in self.lib_properties.keys():
                self.lib_properties[fli] = self.mask_CoolProp
            else:
                self.lib_properties[fli] += self.mask_CoolProp

        # Base model needs density and molecular viscosity

        self.lst = [('density', 'Rho'),('molecular_viscosity', 'Mu')]
        self.node_density   = self.setNewFluidProperty(self.node_fluid, \
                                                       'density')
        self.node_viscosity = self.setNewFluidProperty(self.node_fluid, \
                                                       'molecular_viscosity')
        self.node_lst = [self.node_density, self.node_viscosity]

        self.node_heat = None
        self.node_cond = None
        self.node_vol_visc = None
        self.node_dyn = None

        # Get surface tension for Vof model

        if self.node_hgn['model'] not in [None, "off"]:
            self.lst.append(('surface_tension', 'Sigma'))
            self.node_surface_tension   = self.setNewFluidProperty(self.node_fluid, \
                                                                   'surface_tension')
            self.node_lst.append(self.node_surface_tension)

        # Get thermal scalar and model

        thm = ThermalScalarModel(self.case)
        tsn = thm.getThermalScalarName()
        self.tsm = thm.getThermalScalarModel()

        # If thermal model enabled, add thermal conductivity and specific heat

        if self.tsm != "off":
            self.lst.extend([('specific_heat', 'Cp'), \
                             ('thermal_conductivity', 'Al')])
            self.node_heat = self.setNewFluidProperty(self.node_fluid, \
                                                      'specific_heat')
            self.node_cond = self.setNewFluidProperty(self.node_fluid, \
                                                      'thermal_conductivity')
            self.node_lst.extend([self.node_heat, self.node_cond])

        # Define volume viscosity for compressible model

        if self.node_comp['model'] not in [None, "off"]:
            self.lst.append(('volume_viscosity', 'Viscv0'))
            self.node_vol_visc  = self.setNewFluidProperty(self.node_fluid, \
                                                           'volume_viscosity')
            self.node_lst.append(self.node_vol_visc)

        # Define dynamic diffusion for reactive flow
        elif self.node_coal['model'] not in [None, "off"] \
             or self.node_gas['model'] not in [None, "off"]:
            self.lst.append(('dynamic_diffusion', 'Diftl0'))
            self.node_dyn = self.setNewFluidProperty(self.node_fluid, \
                                                     'dynamic_diffusion')
            self.node_lst.append(self.node_dyn)

        # Build scalars list

        self.list_scalars = []

        if self.tsm == "temperature_celsius":
            self.list_scalars.append((tsn, self.tr("Thermal scalar: temperature (\xB0 C)")))
        elif self.tsm == "temperature_kelvin":
            self.list_scalars.append((tsn, self.tr("Thermal scalar: temperature (K)")))
        elif self.tsm != "off":
            self.list_scalars.append((tsn, self.tr("Thermal scalar")))

        self.m_sca = DefineUserScalarsModel(self.case)
        for s in self.m_sca.getUserScalarNameList():
            self.list_scalars.append((s, self.tr("Additional scalar")))

        # Notebook

        self.notebook = NotebookModel(self.case)


    def __nodeFromTag(self, name):
        """
        Private method : return node with attibute name 'name'
        """
        for node in self.node_lst:
            if node['name'] == name:
                return node

    def getListOfRefProperties(self, first_elt = None):
        """
        Returns a list of physical properties which need definition.
        An optional argument 'first_elt' allows setting a given property
        as the first member of the list.
        """

        l = ['density', 'molecular_viscosity']

        # Vof properties
        if self.node_hgn['model'] not in [None, "off"]:
            l.append('surface_tension')

        # Thermal properties
        if self.tsm != 'off':
            l.extend(['specific_heat', 'thermal_conductivity'])

        # Compressible properties
        if self.node_comp['model'] not in [None, "off"]:
            l.append('volume_viscosity')

        if first_elt in l:
            l.insert(0, l.pop(l.index(first_elt)))

        return l

    def getLibPropertiesDict(self):
       """
       Return dictionnary with available properties and matching flags
       """
       return self.lib_properties


    def getLibPropertyMethods(self, material):
        """
        Return list of methods for a given library and material.
        Values are returned as a tuple: (name, available)
        """
        methods = []

        if material not in self.lib_properties.keys():
            methods = [('unknown', False)]
            return methods

        material_flags = self.lib_properties[material]

        if material_flags & self.mask_builtin:
            methods.append(('user_properties', True))

        if material_flags & self.mask_EOS:
            if material not in eos_excl:
                fls = self.eos.getFluidMethods(material)
                for fli in fls:
                    if fli != 'PerfectGas':
                        methods.append((fli, True))

        if material_flags & self.mask_freesteam:
            avail = self.tables & self.mask_freesteam != 0
            methods.append(("freesteam", avail))

        if material_flags & self.mask_CoolProp:
            avail = self.tables & self.mask_CoolProp != 0
            avail = False
            methods.append(("CoolProp", avail))

        return methods


    def defaultFluidCharacteristicsValues(self):
        """
        Return in a dictionnary which contains default values.
        (also called by CoalCombustionModel)
        """
        default = {}

        # Particular default values init. taking into account chosen model
        mdl_atmo, mdl_joule, mdl_thermal, mdl_gas, mdl_coal, mdl_comp, mdl_hgn=\
            self.getThermoPhysicalModel()

        default['reference_pressure']    = 1.01325e+5

        if self.tsm == "temperature_celsius":
           default['reference_temperature'] = 20.
        else:
           default['reference_temperature'] = 293.15

        default['reference_fuel_temperature']      = 436.
        default['reference_oxydant_temperature']   = 353.
        # molar mass for dry air
        default['reference_molar_mass'] = 28.966e-3

        # Initial values for properties: 20 Celsius degrees air at atmospheric
        # pressure except for VoF

        default['density']                = 1.17862
        default['density_0']              = 1000.
        default['density_1']              = 1.
        default['molecular_viscosity']    = 1.83e-05
        default['molecular_viscosity_0']  = 1.e-03
        default['molecular_viscosity_1']  = 1.e-05
        default['surface_tension']        = 0.
        if mdl_hgn != "off":
           default['density']             = default['density_1']
           default['molecular_viscosity'] = default['molecular_viscosity_1']
        default['specific_heat']          = 1017.24
        default['thermal_conductivity']   = 0.02495
        default['dynamic_diffusion']      = 0.01
        default['volume_viscosity']       = 0.
        default['material']               = "user_material"
        default['method']                 = "user_properties"
        default['reference']              = None

        return default

    @Variables.undoLocal
    def setPressure(self, value):
        """
        Set value of reference pressure into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_fluid.xmlSetData('reference_pressure', value)


    @Variables.noUndo
    def getPressure(self):
        """
        Return the value of reference pressure.
        """
        value = self.node_fluid.xmlGetDouble('reference_pressure')
        if value is None:
           p_str = 'reference_pressure'
           value = self.defaultFluidCharacteristicsValues()[p_str]
           self.setPressure(value)

        return value


    @Variables.undoLocal
    def setTemperature(self, value):
        """
        Set reference temperature.
        """
        self.isGreater(value, 0.0)
        self.node_fluid.xmlSetData('reference_temperature', value)


    @Variables.noUndo
    def getTemperature(self):
        """
        Get reference temperature.
        """
        value = self.node_fluid.xmlGetDouble('reference_temperature')
        if value is None:
           t_str = 'reference_temperature'
           value = self.defaultFluidCharacteristicsValues()[t_str]
           self.setTemperature(value)
        return value


    @Variables.undoLocal
    def setTempOxydant(self, value):
        """
        Set reference temperature for Oxydant.
        """
        self.isGreater(value, 0.0)
        self.node_fluid.xmlSetData('reference_oxydant_temperature', value)


    @Variables.noUndo
    def getTempOxydant(self):
        """
        Get reference temperaturefor Oxydant.
        """
        value = self.node_fluid.xmlGetDouble('reference_oxydant_temperature')
        if value is None:
           ot_str = 'reference_oxydant_temperature'
           value = self.defaultFluidCharacteristicsValues()[ot_str]
           self.setTempOxydant(value)
        return value


    @Variables.undoLocal
    def setTempFuel(self, value):
        """
        Set reference temperature.
        """
        self.isGreater(value, 0.0)
        self.node_fluid.xmlSetData('reference_fuel_temperature', value)


    @Variables.noUndo
    def getTempFuel(self):
        """
        Get reference temperature.
        """
        value = self.node_fluid.xmlGetDouble('reference_fuel_temperature')
        if value is None:
           ft_str = 'reference_fuel_temperature'
           value = self.defaultFluidCharacteristicsValues()[ft_str]
           self.setTempFuel(value)
        return value


    @Variables.undoLocal
    def setMassemol(self, value):
        """
        Set reference molar mass.
        """
        self.isGreater(value, 0.0)
        self.node_fluid.xmlSetData('reference_molar_mass', value)


    @Variables.noUndo
    def getMassemol(self):
        """
        Get reference molar mass.
        """
        value = self.node_fluid.xmlGetDouble('reference_molar_mass')
        if value is None:
           mm_str = 'reference_molar_mass'
           value = self.defaultFluidCharacteristicsValues()[mm_str]
           self.setMassemol(value)
        return value


    @Variables.noUndo
    def getThermoPhysicalModel(self):
        """
        Return values of attribute "model" of all thermophysical model nodes.
        (also called by NumericalParamGlobalView and TimeStepView)
        """
        d = {}

        for mdl in ['atmospheric_flows', 'compressible_model',
                    'gas_combustion', 'joule_effect', 'solid_fuels',
                    'thermal_scalar', 'hgn_model']:
            d[mdl] = 'off'

            node = self.node_models.xmlGetNode(mdl, 'model')

            if node:
                if node['model'] == "":
                    node['model'] = "off"
                if node['model'] != 'off':
                    d[mdl] = node['model']

        return d['atmospheric_flows'], \
               d['joule_effect'],      \
               d['thermal_scalar'],    \
               d['gas_combustion'],    \
               d['solid_fuels'],       \
               d['compressible_model'],\
               d['hgn_model']


    @Variables.noUndo
    def getMaterials(self):
        """
        get the nature of materials
        """
        nodem = self.node_fluid.xmlGetNode('material')
        if nodem is None:
            material = self.defaultFluidCharacteristicsValues()['material']
            self.setMaterials(material)
            nodem = self.node_fluid.xmlGetNode('material')
        material = nodem['choice']
        return material


    @Variables.undoLocal
    def setMaterials(self, material):
        """
        set the nature of materials
        """
        childNode = self.node_fluid.xmlInitChildNode('material')
        m = childNode.xmlGetNode('material')
        oldMaterial = None
        if m != None:
            oldMaterial = m['choice']
        childNode.xmlSetAttribute(choice = material)
        self.updateMethod(oldMaterial)

        if material == "user_material":
            for node in self.node_fluid.xmlGetChildNodeList('property'):
                if node['choice'] == "thermal_law":
                    node['choice'] = "user_law"


    @Variables.noUndo
    def getMethod(self):
        """
        get the method used to compute properties
        """
        nodem = self.node_fluid.xmlGetNode('method')
        if nodem is None:
            method = self.updateMethod("")
            nodem = self.node_fluid.xmlGetNode('method')
        method = nodem['choice']
        return method


    @Variables.undoLocal
    def setMethod(self, method):
        """
        update reference value for EOS
        """
        childNode = self.node_fluid.xmlInitChildNode('method')
        childNode.xmlSetAttribute(choice = method)
        self.getReference()

        # suppress phase choice (not used)
        nodem = self.node_fluid.xmlGetNode('phas')
        if nodem:
            nodem.xmlRemoveNode()

        # suppress reference choice if not EOS
        if method in ("user_properties", "freesteam", "CoolProp"):
            nodem = childNode.xmlGetNode('reference')
            if nodem:
                nodem.xmlRemoveNode()


    @Variables.undoGlobal
    def updateMethod(self, oldMaterial):
        """
        update reference value for EOS
        """
        material = self.getMaterials()
        if oldMaterial == material:
            return

        old_method = self.getMethod

        # If fluid not known, do no change anything (error may be caught at
        # computation time, or the GUI may simply be incomplete)
        if material not in self.lib_properties.keys():
            return

        material_flags = self.lib_properties[material]
        methods = []

        if material_flags & self.mask_builtin:
            methods.append('user_properties')

        if material_flags & self.mask_EOS:
            fls = self.eos.getFluidMethods(material)
            for fli in fls:
                methods.append(fli)

        if material_flags & self.mask_CoolProp:
            methods.append("CoolProp")

        if material_flags & self.mask_freesteam:
            methods.append("freesteam")

        if old_method not in methods and methods:
            methods.append("unknown")
            self.setMethod(methods[0])

        reference = self.getReference() # to force update


    @Variables.noUndo
    def getAvailReferences(self, material, method):
        """
        return available reference value for EOS
        """
        if method in ("user_properties", "freesteam", "CoolProp"):
            return None

        references = self.eos.getFluidReferences(material, self.getMethod())

        return references


    @Variables.noUndo
    def getReference(self):
        """
        return reference value for EOS
        """
        reference = ""
        material = self.getMaterials()
        method = self.getMethod()
        if method in ("user_properties", "freesteam", "CoolProp"):
            return None

        nodem = self.node_fluid.xmlGetChildNode('method')
        reference = nodem.xmlGetString('reference')

        if not reference and self.lib_properties[material] & self.mask_EOS and self.eos.isActive():
            ref = self.eos.getFluidReferences(material, self.getMethod())
            if ref:
                reference = ref[0]

            # update XML
            self.setReference('reference')

        return reference


    @Variables.noUndo
    def setReference(self, reference):
        """
        set reference value for EOS
        """
        nodem = self.node_fluid.xmlInitChildNode('method')
        if not reference:
            if nodem.xmlGetNode('reference'):
                node.xmlRemoveChild('reference')
        else:
            nodem.xmlSetData('reference', reference)


    @Variables.noUndo
    def getInitialValue(self, tag):
        """
        Return initial value of the markup tag : 'density', or
        'molecular_viscosity', or 'surface_tension', or 'specific_heat' or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity', 'surface_tension',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        pp = node.xmlGetDouble('initial_value')
        if pp is None:
            pp = self.defaultFluidCharacteristicsValues()[tag]
            self.setInitialValue(tag, pp)
        return pp


    @Variables.noUndo
    def getValue(self, f_id, tag):
        """
        Return value of the markup tag : 'density', or
        'molecular_viscosity' for fluid of id f_id
        """
        self.isInList(tag, ('density', 'molecular_viscosity'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        pp = node.xmlGetDouble('value' + '_' + str(f_id))
        if pp is None:
            pp = self.defaultFluidCharacteristicsValues()[tag+'_'+str(f_id)]
            self.setValue(f_id, tag, pp)
        return pp


    @Variables.undoLocal
    def setInitialValue(self, tag, val):
        """
        Set initial value for the markup tag : 'density', or
        'molecular_viscosity', or 'surface_tension', or 'specific_heat'or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity', 'surface_tension',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
        if tag != 'volume_viscosity' and tag != 'surface_tension':
            self.isGreater(val, 0.)
        else:
            self.isPositiveFloat(val)
        node = self.node_fluid.xmlGetNode('property', name=tag)
        node.xmlSetData('initial_value', val)


    @Variables.undoLocal
    def setValue(self, f_id, tag, val):
        """
        Set initial value for the markup tag : 'density', or
        'molecular_viscosity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        node.xmlSetData('value'+'_'+str(f_id), val)


    @Variables.noUndo
    def getInitialValueDensity(self):
        """Return initial value of density"""
        return self.getInitialValue('density')


    @Variables.noUndo
    def getVofValueDensity(self, f_id):
        """Return value of density of fluid of id f_id"""
        return self.getValue(f_id, 'density')


    @Variables.undoLocal
    def setInitialValueDensity(self, val):
        """Put initial value for density"""
        self.setInitialValue('density', val)


    @Variables.undoLocal
    def setVofValueDensity(self, f_id, val):
        """Set value for density of fluid of id f_id"""
        self.setValue(f_id, 'density', val)


    @Variables.noUndo
    def getInitialValueViscosity(self):
        """Return initial value of viscosity"""
        return self.getInitialValue('molecular_viscosity')


    @Variables.noUndo
    def getVofValueViscosity(self, f_id):
        """Return value of viscosity of fluid of id f_id"""
        return self.getValue(f_id, 'molecular_viscosity')


    @Variables.undoLocal
    def setInitialValueViscosity(self, val):
        """Put initial value for viscosity"""
        self.setInitialValue('molecular_viscosity', val)


    @Variables.undoLocal
    def setVofValueViscosity(self, f_id, val):
        """Set value for viscosity of fluid of id f_id"""
        self.setValue(f_id, 'molecular_viscosity', val)


    @Variables.noUndo
    def getInitialValueVolumeViscosity(self):
        """Return initial value of volume viscosity"""
        return self.getInitialValue('volume_viscosity')


    @Variables.undoLocal
    def setInitialValueVolumeViscosity(self, val):
        """Put initial value for volume viscosity"""
        self.setInitialValue('volume_viscosity', val)


    @Variables.noUndo
    def getInitialValueSurfaceTension(self):
        """Return value of surface tension"""
        return self.getInitialValue('surface_tension')


    @Variables.undoLocal
    def setInitialValueSurfaceTension(self, val):
        """Set value for surface tension"""
        self.setInitialValue('surface_tension', val)


    @Variables.noUndo
    def getInitialValueHeat(self):
        """Return initial value of specific heat"""
        return self.getInitialValue('specific_heat')


    @Variables.undoLocal
    def setInitialValueHeat(self, val):
        """Put initial value for specific heat"""
        self.setInitialValue('specific_heat', val)


    @Variables.noUndo
    def getInitialValueCond(self):
        """Return initial value of conductivity"""
        return self.getInitialValue('thermal_conductivity')


    @Variables.undoLocal
    def setInitialValueCond(self, val):
        """Put initial value for conductivity"""
        self.setInitialValue('thermal_conductivity', val)


    @Variables.noUndo
    def getInitialValueDyn(self):
        """Return initial value of conductivity"""
        return self.getInitialValue('dynamic_diffusion')


    @Variables.undoLocal
    def setInitialValueDyn(self, val):
        """Put initial value for conductivity"""
        self.setInitialValue('dynamic_diffusion', val)


    @Variables.noUndo
    def getFormulaTry(self, tag, zone="1"):
        """
        Return a formula for I{tag} 'density', 'molecular_viscosity',
        'specific_heat' or 'thermal_conductivity'
        only if actually defined, None otherwise.
        """
        self.isInList(tag, ('density', 'dynamic_diffusion',
                            'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        if zone != "1":
            if node.xmlGetChildNode("zone", zone_id=zone):
                node = node.xmlGetChildNode("zone", zone_id=zone)
            else:
                node = node.xmlInitChildNode("zone", zone_id=zone)
        formula = node.xmlGetChildString('formula')
        return formula


    @Variables.noUndo
    def getFormula(self, tag, zone="1"):
        """
        Return a formula for I{tag} 'density', 'molecular_viscosity',
        'specific_heat' or 'thermal_conductivity'
        """
        formula = self.getFormulaTry(tag, zone)
        if not formula:
            formula = self.getDefaultFormula(tag)
            self.setFormula(tag, formula, zone)
        return formula


    @Variables.noUndo
    def getDefaultFormula(self, tag):
        """
        Return default formula
        """
        self.isInList(tag, self.properties_info.get_property_list())

        prop_info = self.properties_info.get_property(name=tag)
        formula = prop_info['name'] + " = " + prop_info['ref_tag'] + ";"

        return formula


    @Variables.undoLocal
    def setFormula(self, tag, formula, zone="1"):
        """
        Gives a formula for 'density', 'molecular_viscosity',
        'specific_heat' or 'thermal_conductivity'
        """
        self.isInList(tag, self.properties_info.get_property_list())

        node = self.node_fluid.xmlGetNode('property', name=tag)
        if zone != "1":
            if node.xmlGetChildNode("zone", zone_id=zone):
                node = node.xmlGetChildNode("zone", zone_id=zone)
            else:
                node = node.xmlInitChildNode("zone", zone_id=zone)
        node.xmlSetData('formula', formula)


    @Variables.noUndo
    def getPropertyMode(self, tag):
        """Return choice of node I{tag}. Choice is constant or variable"""
        self.isInList(tag, self.properties_info.get_property_list()+['dynamic_diffusion'])
        node = self.__nodeFromTag(tag)
        c = node['choice']
        self.isInList(c, ('constant', 'thermal_law', 'user_law', 'predefined_law'))
        return c


    @Variables.undoGlobal
    def setPropertyMode(self, tag, choice):
        """Put choice in xml file's node I{tag}"""
        self.isInList(tag, self.properties_info.get_property_list()+['dynamic_diffusion'])
        self.isInList(choice, ('constant', 'thermal_law',
                               'user_law', 'predefined_law'))

        node = self.__nodeFromTag(tag)
        node['choice'] = choice

        if choice == 'constant':

            if tag == 'density':
                from code_saturne.model.TimeStepModel import TimeStepModel
                TimeStepModel(self.case).RemoveThermalTimeStepNode()
                del TimeStepModel
        else:
            if node.xmlGetNode('listing_printing'):
                node.xmlRemoveChild('listing_printing')
            if node.xmlGetNode('postprocessing_recording'):
                node.xmlRemoveChild('postprocessing_recording')

    # MEG Generation related functions
    def getFormulaComponents(self, tag, scalar=None, zone="1"):
        """
        Get the formula components for a given tag
        """

        if tag in self.getListOfRefProperties():
            return self.getFormulaPropertyComponents(tag, zone)

        elif tag == 'scalar_diffusivity' and scalar != None:
            return self.getFormulaDiffComponents(scalar, zone)

        else:
            msg = 'Formula is not available for field %s in MEG' % tag
            raise Exception(msg)


    def getFormulaPropertyComponents(self, prop_name, zone="1"):
        """
        User formula for a given physical property other than diffusion.
        """

        prop_info = self.properties_info.get_property(name=prop_name)

        exp = self.getFormula(prop_name, zone)
        req = [(prop_info['name'], prop_info['output'])]

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for pn in self.getListOfRefProperties(first_elt=prop_name):
            p = self.properties_info.get_property(name=pn)

            vstr  = str(self.getInitialValue(tag=p['name']))
            tstr  = p['ref_tag']
            ostr  = p['output']
            ustr  = p['unit']
            symbols.append((tstr, '{} (reference value, {}) = {}'.format(ostr, ustr, vstr)))

        ref_pressure = self.getPressure()
        symbols.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        symbols.append(('volume', 'Zone volume'))
        symbols.append(('fluid_volume', 'Zone fluid volume (different if porous)'))

        for _c in ['x', 'y', 'z']:
            symbols.append(('g%s' % _c,
                            'Gravity component along the %s axis' % _c))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaDiffComponents(self, scalar, zone="1"):
        """
        User formula for the diffusion coefficient
        """
        name = self.m_sca.getScalarDiffusivityName(scalar)
        exp  = self.m_sca.getDiffFormula(scalar, zone)
        req = [(str(name), str(scalar)+' diffusion coefficient')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('volume', 'Zone volume'),
               ('fluid_volume', 'Zone fluid volume (different if porous)')]
        sym.append((str(scalar),str(scalar)))
        diff0_value = self.m_sca.getScalarDiffusivityInitialValue(scalar)
        sym.append((str(name)+'_ref', str(scalar)+' diffusion coefficient (reference value, m^2/s) = '+str(diff0_value)))

        sym.append(('volume', 'Zone volume'))
        sym.append(('fluid_volume', 'Zone fluid volume (different if porous)'))

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, sym


##    def RemoveThermoConductNode(self):
##        """Remove balise property for thermal_conductivity"""
##        self.node_fluid.xmlRemoveChild('property', name='thermal_conductivity')


    def tr(self, text):
        """
        translation
        """
        return text


#-------------------------------------------------------------------------------
# FluidCharacteristicsModel test case
#-------------------------------------------------------------------------------


class FluidCharacteristicsModelTestCase(ModelTest):
    """
    """
    def checkFluidCharacteristicsInstantiation(self):
        """Check whether the FluidCaracteristicsModel class could be instantiated"""
        model = None
        model = FluidCharacteristicsModel(self.case)
        assert model != None, 'Could not instantiate FluidCaracteristicsModel'

    def checkGetThermoPhysicalModel(self):
        """Check whether thermal physical models could be get"""
        mdl = FluidCharacteristicsModel(self.case)
        from code_saturne.model.ThermalScalarModel import ThermalScalarModel
        ThermalScalarModel(self.case).setThermalModel('temperature_celsius')
        del ThermalScalarModel
        assert mdl.getThermoPhysicalModel() == ('off', 'off', 'temperature_celsius', 'off', 'off', 'off', 'off'),\
        'Could not get thermophysical models in FluidCaracteristicsModel'

    def checkSetandGetInitialValue(self):
        """Check whether the initial value for the properties could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValue('density', 123.0)
        mdl.setInitialValue('molecular_viscosity', 1.5e-5)
        mdl.setInitialValue('surface_tension', 0.)
        mdl.setInitialValue('specific_heat', 1212)
        mdl.setInitialValue('thermal_conductivity', 0.04)
        doc = '''<fluid_properties>
                    <property choice="constant" label="Density" name="density">
                        <initial_value>123</initial_value>
                    </property>
                    <property choice="constant" label="LamVisc" name="molecular_viscosity">
                        <initial_value>1.5e-05</initial_value>
                    </property>
                    <property choice="constant" label="SurfTen" name="surface_tension">
                        <initial_value>0.</initial_value>
                    </property>
                    <property choice="constant" label="Sp. heat" name="specific_heat">
                        <initial_value>1212</initial_value>
                    </property>
                    <property choice="constant" label="Th. cond" name="thermal_conductivity">
                        <initial_value>0.04</initial_value>
                    </property>
                 </fluid_properties>'''
        assert mdl.node_fluid == self.xmlNodeFromString(doc),\
        'Could not set initial value for properties in FluidCharacteristicsModel model'
        assert mdl.getInitialValue('specific_heat') == 1212.,\
         'Could not get initial value for properties in FluidCharacteristicsModel model'

    def checkSetandGetInitialValueDensity(self):
        """Check whether the initial value for density could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValueDensity(89.98)
        doc = '''<property choice="constant" label="Density" name="density">
                    <initial_value>89.98</initial_value>
                 </property>'''
        assert mdl.node_density == self.xmlNodeFromString(doc),\
        'Could not set initial value of density'
        assert mdl.getInitialValueDensity() == 89.98,\
        'Could not get initial value of density'

    def checkSetandGetInitialValueViscosity(self):
        """Check whether the initial value for molecular_viscosity could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValueViscosity(1.2e-4)
        doc = '''<property choice="constant" label="LamVisc" name="molecular_viscosity">
                    <initial_value>1.2e-4</initial_value>
                 </property>'''
        assert mdl.node_viscosity == self.xmlNodeFromString(doc),\
        'Could not set initial value of molecular_viscosity'
        assert mdl.getInitialValueViscosity() == 1.2e-4,\
        'Could not get initial value of molecular_viscosity'

    def checkSetandGetInitialValueSurfaceTension(self):
        """Check whether the initial value for surface_tension could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValueSurfaceTension(0.)
        doc = '''<property choice="constant" label="SurfTen" name="surface_tension">
                    <initial_value>0.</initial_value>
                 </property>'''
        assert mdl.node_surface_tension == self.xmlNodeFromString(doc),\
        'Could not set initial value of surface_tension'
        assert mdl.getInitialValueSurfaceTension() == 0.,\
        'Could not get initial value of surface_tension'

    def checkSetandGetInitialValueHeat(self):
        """Check whether the initial value for specific_heat could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValueHeat(456.9)
        doc = '''<property choice="constant" label="Sp. heat" name="specific_heat">
                    <initial_value>456.9</initial_value>
                 </property>'''
        assert mdl.node_heat == self.xmlNodeFromString(doc),\
        'Could not set initial value of specific_heat'
        assert mdl.getInitialValueHeat() == 456.9,\
        'Could not get initial value of specific_heat'

    def checkSetandGetInitialValueCond(self):
        """Check whether the initial value for thermal_conductivity could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValueCond(456.9)
        doc = '''<property choice="constant" label="Th. cond" name="thermal_conductivity">
                    <initial_value>456.9</initial_value>
                 </property>'''
        assert mdl.node_cond == self.xmlNodeFromString(doc),\
        'Could not set initial value of thermal_conductivity'
        assert mdl.getInitialValueCond() == 456.9,\
        'Could not get initial value of thermal_conductivity'

    def checkSetandGetPropertyMode(self):
        """Check whether choice constant or variable could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setPropertyMode('density', 'user_law')
        doc = '''<property choice="user_law" label="Density" name="density">
                    <initial_value>1.17862</initial_value>
                 </property>'''
        assert mdl.node_density == self.xmlNodeFromString(doc),\
        'Could not set choice constant or variable for property'
        assert mdl.getPropertyMode('density') == 'user_law',\
        'Could not get choice constant or variable for property'

        mdl.setPropertyMode('density', 'constant')
        doc2 = '''<property choice="constant" label="Density" name="density">
                    <initial_value>1.17862</initial_value>
                 </property>'''
        assert mdl.node_density == self.xmlNodeFromString(doc2),\
        'Could not set listing and recording status for property and remove this nodes'

##    def checkRemoveThermoConductNode(self):
##        """Check whether node thermal conductivity could be removed"""
##        mdl = FluidCharacteristicsModel(self.case)
##        mdl.setPropertyMode('density', 'variable')
##        mdl.setPropertyMode('molecular_viscosity', 'variable')
##        mdl.setPropertyMode('specific_heat', 'variable')
##        mdl.RemoveThermoConductNode()
##        doc = '''<fluid_properties>
##                    <property choice="variable" label="Density" name="density">
##                        <initial_value>1.17862</initial_value>
##                    </property>
##                    <property choice="variable" label="LamVisc" name="molecular_viscosity">
##                        <initial_value>1.83e-05</initial_value>
##                    </property>
##                    <property choice="variable" label="Sp. heat" name="specific_heat">
##                        <initial_value>1017.24</initial_value>
##                    </property>
##                 </fluid_properties>'''
##        assert mdl.node_fluid == self.xmlNodeFromString(doc),\
##        'Could not remove thermal_conductivity node'


def suite():
    testSuite = unittest.makeSuite(FluidCharacteristicsModelTestCase, "check")
    return testSuite


def runTest():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
