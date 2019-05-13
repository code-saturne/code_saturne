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

EOS = 1
try:
   import eosAva
except:
   EOS = 0
else :
   import eosAva

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

        # Look for thermo tables

        import code_saturne.cs_config as cs_config
        cfg = cs_config.config()
        self.freesteam = 0
        if cfg.libs['freesteam'].have != "no":
            self.freesteam = 1

        self.coolprop = 0
        if cfg.libs['coolprop'].have != "no":
            self.coolprop = 1

        # Notebook

        self.notebook = NotebookModel(self.case)


    def __nodeFromTag(self, name):
        """
        Private method : return node with attibute name 'name'
        """
        for node in self.node_lst:
            if node['name'] == name:
                return node


    def defaultFluidCharacteristicsValues(self):
        """
        Return in a dictionnary which contains default values.
        (also called by CoalCombustionModel)
        """
        default = {}

        default['reference_pressure']    = 1.01325e+5

        if self.tsm == "temperature_celsius":
           default['reference_temperature'] = 20.
        else:
           default['reference_temperature'] = 293.15

        default['fuel_temperature']      = 436.
        default['oxydant_temperature']   = 353.
        # molar mass for dry air
        default['reference_molar_mass'] = 28.966e-3

        #Initial values for properties: 20 Celsius degrees air at atmospheric
        #pressure

        default['density']              = 1.17862
        default['molecular_viscosity']  = 1.83e-05
        default['specific_heat']        = 1017.24
        default['thermal_conductivity'] = 0.02495
        default['dynamic_diffusion']    = 0.01
        default['volume_viscosity']     = 0.
        default['material']             = "user_material"
        default['method']               = "user_properties"
        default['phas']                 = "liquid"

        return default

    @Variables.undoLocal
    def setPressure(self, value):
        """
        Set value of reference pressure into xml file.
        """
        self.isGreaterOrEqual(value, 0.0)
        self.node_fluid.xmlSetData('reference_pressure',value)


    @Variables.noUndo
    def getPressure(self):
        """
        Return the value of reference pressure.
        """
        value = self.node_fluid.xmlGetDouble('reference_pressure')
        if value == None:
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
        if not value :
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
        if not value :
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
        if not value :
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
        if not value :
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
                    'thermal_scalar']:
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
               d['compressible_model']


    @Variables.noUndo
    def getMaterials(self):
        """
        get the nature of materials
        """
        nodem = self.node_fluid.xmlGetNode('material')
        if nodem == None:
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
        get the nature of materials
        """
        nodem = self.node_fluid.xmlGetNode('method')
        if nodem == None:
            method = self.defaultFluidCharacteristicsValues()['method']
            self.setMethod(method)
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

        # suppress phas choice if not EOS
        if method == self.defaultFluidCharacteristicsValues()['method'] or \
           method == "freesteam" or method == "CoolProp":
            nodem = self.node_fluid.xmlGetNode('phas')
            if nodem:
                nodem.xmlRemoveNode()


    @Variables.undoGlobal
    def updateMethod(self, oldMaterial):
        """
        update reference value for EOS
        """
        material = self.getMaterials()
        if oldMaterial != material:
            if material == self.defaultFluidCharacteristicsValues()['material'] :
               self.setMethod(self.defaultFluidCharacteristicsValues()['method'])
            elif EOS == 1:
               self.ava = eosAva.EosAvailable()
               self.ava.setMethods(material)
               self.setMethod(self.ava.whichMethods()[0])
            elif self.coolprop == 1:
                self.setMethod("CoolProp")
            else:
                self.setMethod("freesteam")
            # suppress phas choice if not EOS
            if self.getMethod() == self.defaultFluidCharacteristicsValues()['method'] or \
               self.getMethod() == "freesteam" or self.getMethod() == "CoolProp":
                nodem = self.node_fluid.xmlGetNode('phas')
                if nodem:
                   nodem.xmlRemoveNode()


    @Variables.noUndo
    def getReference(self):
        """
        return reference value for EOS
        """
        reference = ""
        material = self.getMaterials()
        method = self.getMethod()
        if material == "user_material":
            reference = material
        else:
            if self.freesteam == 1:
                if self.getMethod() == "freesteam":
                    reference = "freesteam"
            if self.coolprop == 1:
                if self.getMethod() == "CoolProp":
                    reference = "CoolProp"
            if EOS == 1:
                if self.getMethod() != "freesteam" and self.getMethod() != "CoolProp":
                    phas = self.getFieldNature()
                    self.ava = eosAva.EosAvailable()
                    self.ava.setMethods(material)
                    self.ava.setReferences(material, self.getMethod())
                    ref = self.ava.whichReferences()

                    if phas == "liquid":
                        reference = ref[0]
                    elif phas == "gas":
                        reference = ref[1]
        # update XML
        childNode = self.node_fluid.xmlInitChildNode('reference')
        childNode.xmlSetAttribute(choice = reference)

        return reference


    @Variables.noUndo
    def getFieldNature(self):
        """
        get field nature for EOS
        """
        nodem = self.node_fluid.xmlGetNode('phas')
        if nodem == None:
            nature = self.defaultFluidCharacteristicsValues()['phas']
            self.setFieldNature(nature)
            nodem = self.node_fluid.xmlGetNode('phas')
        nature = nodem['choice']
        return nature


    @Variables.noUndo
    def setFieldNature(self, nature):
        """
        set field nature for EOS
        """
        self.isInList(nature, ('liquid', 'gas'))
        childNode = self.node_fluid.xmlInitChildNode('phas')
        childNode.xmlSetAttribute(choice = nature)


    @Variables.noUndo
    def getInitialValue(self, tag):
        """
        Return initial value of the markup tag : 'density', or
        'molecular_viscosity', or 'specific_heat'or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        pp = node.xmlGetDouble('initial_value')
        if pp == None:
            pp = self.defaultFluidCharacteristicsValues()[tag]
            self.setInitialValue(tag, pp)
        return pp


    @Variables.undoLocal
    def setInitialValue(self, tag, val):
        """
        Put initial value for the markup tag : 'density', or
        'molecular_viscosity', or 'specific_heat'or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
        if tag != 'volume_viscosity':
            self.isGreater(val, 0.)
        else:
            self.isPositiveFloat(val)
        node = self.node_fluid.xmlGetNode('property', name=tag)
        node.xmlSetData('initial_value', val)


    @Variables.noUndo
    def getInitialValueDensity(self):
        """Return initial value of density"""
        return self.getInitialValue('density')


    @Variables.undoLocal
    def setInitialValueDensity(self, val):
        """Put initial value for density"""
        self.setInitialValue('density', val)


    @Variables.noUndo
    def getInitialValueViscosity(self):
        """Return initial value of viscosity"""
        return self.getInitialValue('molecular_viscosity')


    @Variables.undoLocal
    def setInitialValueViscosity(self, val):
        """Put initial value for viscosity"""
        self.setInitialValue('molecular_viscosity', val)


    @Variables.noUndo
    def getInitialValueVolumeViscosity(self):
        """Return initial value of volume viscosity"""
        return self.getInitialValue('volume_viscosity')


    @Variables.undoLocal
    def setInitialValueVolumeViscosity(self, val):
        """Put initial value for volume viscosity"""
        self.setInitialValue('volume_viscosity', val)


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
    def getFormula(self, tag):
        """
        Return a formula for I{tag} 'density', 'molecular_viscosity',
        'specific_heat' or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        formula = node.xmlGetString('formula')
        if not formula:
            formula = self.getDefaultFormula(tag)
            self.setFormula(tag, formula)
        return formula


    @Variables.noUndo
    def getDefaultFormula(self, tag):
        """
        Return default formula
        """
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity'))

        formula = tag + " = -1.;"

        return formula


    @Variables.undoLocal
    def setFormula(self, tag, str):
        """
        Gives a formula for 'density', 'molecular_viscosity',
        'specific_heat'or 'thermal_conductivity'
        """
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity'))
        node = self.node_fluid.xmlGetNode('property', name=tag)
        node.xmlSetData('formula', str)


    @Variables.noUndo
    def getPropertyMode(self, tag):
        """Return choice of node I{tag}. Choice is constant or variable"""
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
        node = self.__nodeFromTag(tag)
        c = node['choice']
        self.isInList(c, ('constant', 'thermal_law', 'user_law', 'predefined_law'))
        return c


    @Variables.undoGlobal
    def setPropertyMode(self, tag, choice):
        """Put choice in xml file's node I{tag}"""
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity',
                            'volume_viscosity', 'dynamic_diffusion'))
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
    def getFormulaComponents(self, tag, scalar=None):
        """
        Get the formula components for a given tag
        """

        if tag == 'density':
            return self.getFormulaRhoComponents()

        elif tag == 'molecular_viscosity':
            return self.getFormulaMuComponents()

        elif tag == 'specific_heat':
            return self.getFormulaCpComponents()

        elif tag == 'thermal_conductivity':
            return self.getFormulaAlComponents()

        elif tag == 'volume_viscosity':
            return self.getFormulaViscv0Components()

        elif tag == 'scalar_diffusivity' and scalar != None:
            return self.getFormulaDiffComponents(scalar)

        else:
            msg = 'Formula is not available for field %s in MEG' % tag
            raise Exception(msg)


    def getFormulaRhoComponents(self):
        """
        User formula for density
        """

        exp = self.getFormula('density')
        req = [('density', 'Density')]

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        rho0_value = self.getInitialValueDensity()
        symbols.append(('rho0', 'Density (reference value) = ' + str(rho0_value)))

        ref_pressure = self.getPressure()
        symbols.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaMuComponents(self):
        """
        User formula for molecular viscosity
        """

        exp = self.getFormula('molecular_viscosity')
        req = [('molecular_viscosity', 'Molecular Viscosity')]

        symbols = []

        for s in self.list_scalars:
           symbols.append(s)
        mu0_value = self.getInitialValueViscosity()
        symbols.append(('mu0', 'Viscosity (reference value) = ' + str(mu0_value)))

        rho0_value = self.getInitialValueDensity()
        symbols.append(('rho0', 'Density (reference value) = ' + str(rho0_value)))

        ref_pressure = self.getPressure()
        symbols.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        symbols.append(('rho', 'Density'))

        from code_saturne.model.CompressibleModel import CompressibleModel
        if CompressibleModel(self.case).getCompressibleModel() == 'on':
            symbols.append(('T', 'Temperature'))
            ref_temperature = self.getTemperature()
            symbols.append(('t0', 'Reference temperature = '+str(ref_temperature)+' K'))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaCpComponents(self):
        """
        User formula for specific heat
        """
        exp = self.getFormula('specific_heat')
        req = [('specific_heat', 'Specific heat')]

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        cp0_value = self.getInitialValueHeat()
        symbols.append(('cp0', 'Specific heat (reference value) = ' + str(cp0_value)))

        ref_pressure = self.getPressure()
        symbols.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaViscv0Components(self):
        """
        User formula for volumic viscosity
        """
        exp = self.getFormula('volume_viscosity')
        req = [('volume_viscosity', 'Volume viscosity')]

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        viscv0_value = self.getInitialValueVolumeViscosity()
        ref_pressure = self.getPressure()
        ref_temperature = self.getTemperature()
        symbols.append(('viscv0', 'Volume viscosity (J/kg/K) = '+str(viscv0_value)))
        symbols.append(('p0', 'Reference pressure (Pa) = '+str(ref_pressure)))
        symbols.append(('t0', 'Reference temperature (K) = '+str(ref_temperature)))
        symbols.append(('T', 'Temperature'))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaAlComponents(self):
        """
        User formula for thermal conductivity
        """
        exp = self.getFormula('thermal_conductivity')
        req = [('thermal_conductivity', 'Thermal conductivity')]

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        lambda0_value = self.getInitialValueCond()
        ref_pressure = self.getPressure()
        symbols.append(('lambda0', 'Thermal conductivity (reference value) = ' + str(lambda0_value)))
        symbols.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols;


    def getFormulaDiffComponents(self, scalar):
        """
        User formula for the diffusion coefficient
        """
        name = self.m_sca.getScalarDiffusivityName(scalar)
        exp  = self.m_sca.getDiffFormula(scalar)
        req = [(str(name), str(scalar)+' diffusion coefficient')]
        sym = [('x','cell center coordinate'),
               ('y','cell center coordinate'),
               ('z','cell center coordinate'),]
        sym.append((str(scalar),str(scalar)))
        diff0_value = self.m_sca.getScalarDiffusivityInitialValue(scalar)
        sym.append((str(name)+'_ref', str(scalar)+' diffusion coefficient (reference value, m^2/s) = '+str(diff0_value)))

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
        assert mdl.getThermoPhysicalModel() == ('off', 'off', 'temperature_celsius', 'off', 'off', 'off'),\
        'Could not get thermophysical models in FluidCaracteristicsModel'

    def checkSetandGetInitialValue(self):
        """Check whether the initial value for the properties could be set and get"""
        mdl = FluidCharacteristicsModel(self.case)
        mdl.setInitialValue('density', 123.0)
        mdl.setInitialValue('molecular_viscosity', 1.5e-5)
        mdl.setInitialValue('specific_heat', 1212)
        mdl.setInitialValue('thermal_conductivity', 0.04)
        doc = '''<fluid_properties>
                    <property choice="constant" label="Density" name="density">
                        <initial_value>123</initial_value>
                    </property>
                    <property choice="constant" label="LamVisc" name="molecular_viscosity">
                        <initial_value>1.5e-05</initial_value>
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
