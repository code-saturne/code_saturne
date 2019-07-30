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

import sys, unittest
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.OutputFieldsModel import OutputFieldsModel
from code_saturne.model.SpeciesModel import SpeciesModel

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
# Constructor
#-------------------------------------------------------------------------------


class ThermodynamicsModel(MainFieldsModel, Variables, Model):
    """
    This class manages the Field objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case            = case
        self.XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
        self.__XMLNodefields = self.XMLNodethermo.xmlInitNode('fields')
        self.__XMLThermo     = self.XMLNodethermo.xmlInitNode('thermodynamics')
        self.XMLNodeproperty = self.XMLNodethermo.xmlInitNode('properties')

        self.m_spe = SpeciesModel(self.case)
        self.notebook = NotebookModel(self.case)
        self.list_scalars = []

        self.m_out = OutputFieldsModel(self.case)
        label = self.m_out.getVariableLabel("none", "pressure")
        self.list_scalars.append(('pressure', label))


    def defaultValues(self):
        default = {}
        default['density']              = 1.8
        default['molecular_viscosity']  = 0.0000456
        default['specific_heat']        = 4000
        default['thermal_conductivity'] = 1.e-5
        default['surface_tension']      = 0.075
        default['emissivity']           = 0.
        default['elasticity']           = 0.9
        default['radiative']            = "off"
        default['material']             = "user_material"
        default['method']               = "user_properties"
        default['user_property']        = "constant"
        default['thetisClipping']       = "on"
        default['cathareClipping']      = "on"
        default['propertyChoice']       = "constant"

        return default


    def propertiesFormulaList(self):
        lst = ('density', 'molecular_viscosity',
                'specific_heat', 'thermal_conductivity',
                'surface_tension',
                'temperature','d_rho_d_P', 'd_rho_d_h',
                'SaturationEnthalpyLiquid', 'SaturationEnthalpyGas',
                'd_Hsat_d_P_Liquid', 'd_Hsat_d_P_Gas',
                'SaturationTemperature','d_Tsat_d_P', 'LatentHeat')
        return lst


    @Variables.undoLocal
    def setMaterials(self, fieldId, material):
        """
        set the nature of materials
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('material')
        m = childNode.xmlGetNode('material')
        oldMaterial = None
        if m != None:
            oldMaterial = m['choice']
        childNode.xmlSetAttribute(choice = material)
        if material != "user_material":
            for prop in self.propertiesFormulaList():
                node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=prop)
                if node:
                    node.xmlRemoveChild('formula')
            # add node enthalpy if needed
            XMLNodeVariable = self.XMLNodethermo.xmlGetNode('variables')
            node = XMLNodeVariable.xmlGetNode('variable', field_id=fieldId, name="enthalpy")
            if not node:
                Variables(self.case).setNewVariableProperty("variable", "", XMLNodeVariable, fieldId, "enthalpy", "enthalpy_"+field_name)

            if self.getPredefinedFlow() == 'None' or self.getPredefinedFlow() == 'particles_flow':
                if oldMaterial == "user_material" and (str(fieldId) == '1' or str(fieldId) == '2'):
                    if str(fieldId) == '1':
                        field2 = '2'
                    else:
                        field2 = '1'
                    node2 = self.__XMLNodefields.xmlGetNode('field', field_id = field2)
                    childNode2 = node2.xmlInitChildNode('material')
                    m2 = childNode2.xmlGetNode('material')
                    oldMaterial2 = None
                    if m2 != None:
                        oldMaterial2 = m2['choice']
                    childNode2.xmlSetAttribute(choice = material)
                    # update method
                    self.updateMethod(field2, oldMaterial2)
        else:
            # if predefine flow same material for 2 phases
            if self.getPredefinedFlow() != 'None' and self.getPredefinedFlow() != 'particles_flow':
                if str(fieldId) == '1' or str(fieldId) == '2':
                    if str(fieldId) == '1':
                        field2 = '2'
                    else:
                        field2 = '1'
                    node2 = self.__XMLNodefields.xmlGetNode('field', field_id = field2)
                    childNode2 = node2.xmlInitChildNode('material')
                    m2 = childNode2.xmlGetNode('material')
                    oldMaterial2 = None
                    if m2 != None:
                        oldMaterial2 = m2['choice']
                    childNode2.xmlSetAttribute(choice = material)
                    # update method
                    self.updateMethod(field2, oldMaterial2)


    @Variables.noUndo
    def getMaterials(self, fieldId):
        """
        get the nature of materials
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodem = node.xmlGetNode('material')
        if nodem == None :
            material = self.defaultValues()['material']
            self.setMaterials(fieldId, material)
        material = node.xmlGetNode('material')['choice']
        return material


    @Variables.undoLocal
    def setMethod(self, fieldId, method):
        """
        set the nature of materials
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('method')
        childNode.xmlSetAttribute(choice = method)
        self.updateReference(fieldId)


    @Variables.noUndo
    def getMethod(self, fieldId):
        """
        get the nature of materials
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodem = node.xmlGetNode('method')
        if nodem == None :
            method = self.defaultValues()['method']
            self.setMethod(fieldId, method)
        method = node.xmlGetNode('method')['choice']
        return method


    @Variables.undoGlobal
    def updateMethod(self, fieldId, oldMaterial):
        """
        update reference value for EOS
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        material = self.getMaterials(fieldId)
        if oldMaterial != material :
            if material == self.defaultValues()['material'] :
               self.setMethod(fieldId, self.defaultValues()['method'])
            elif EOS == 1 :
                self.ava = eosAva.EosAvailable()
                self.ava.setMethods(material)
                fls = self.ava.whichMethods()
                if "Cathare" in fls:
                    self.setMethod(fieldId, "Cathare")
                else:
                    for fli in fls:
                        if fli != "Ovap" and fli != "Flica4" and fli != "StiffenedGas":
                            self.setMethod(fieldId, fli)
                            break

            if material == 'user_material' :
                choice = 'constant'
            else :
                choice = 'user_law'
            for tag in ['density', 'molecular_viscosity', 'specific_heat', 'thermal_conductivity'] :
                self.setPropertyMode(fieldId, tag, choice)

            if self.getFieldNature(fieldId) == "gas":
                tag = 'surface_tension'
                fieldId = "none"
                self.setPropertyMode(fieldId, tag, choice)


    @Variables.undoGlobal
    def updateReference(self, fieldId):
        """
        return reference value for EOS
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        reference = ""
        material = self.getMaterials(fieldId)
        method = self.getMethod(fieldId)
        if material == "user_material" :
            reference = material
        else :
            phase = self.getFieldNature(fieldId)
            if EOS == 1 :
                self.ava = eosAva.EosAvailable()
                self.ava.setMethods(material)
                self.ava.setReferences(material, self.getMethod(fieldId))
                ref = self.ava.whichReferences()

                # Reference for methods
                liqidx = 0
                gasidx = 1
                # if Refprop, the list has an additional element at position 0
                if self.getMethod(fieldId) == "Refprop":
                    liqidx += 1
                    gasidx += 1

                if len(ref) == 1 :
                   # cas des gaz par exemple
                   reference = ref[0]
                else :
                   if phase == "liquid" :
                      reference = ref[liqidx]
                   elif phase == "gas" :
                      reference = ref[gasidx]
            else :
                reference = material + phase
        # update XML
        childNode = node.xmlInitChildNode('reference')
        childNode.xmlSetAttribute(choice = reference)

        return reference


    @Variables.noUndo
    def getInitialValue(self, fieldId, tag):
        """
        Return initial value of the markup tag : 'density', or
        'molecular_viscosity', 'specific_heat', 'thermal_conductivity',
        'surface_tension', 'emissivity' or 'elasticity'
        """
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId),fieldIdList)

        lst = ('density', 'molecular_viscosity',
                'specific_heat', 'thermal_conductivity',
                'surface_tension', 'emissivity','elasticity')
        self.isInList(tag, lst)
        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)
        pp = node.xmlGetDouble('initial_value')
        if pp == None:
            pp = self.defaultValues()[tag]
            self.setInitialValue(fieldId, tag, pp)
        return pp


    @Variables.undoLocal
    def setInitialValue(self, fieldId, tag, val):
        """
        Put initial value for the markup tag : 'density', or
        'molecular_viscosity', 'specific_heat', 'thermal_conductivity',
        'surface_tension', 'emissivity' or 'elasticity'
        """
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId),fieldIdList)
        lst = ('density', 'molecular_viscosity',
                'specific_heat', 'thermal_conductivity',
                'surface_tension', 'emissivity','elasticity')
        self.isInList(tag, lst)
        self.isFloat(val)
        if tag != 'emissivity' and tag != 'elasticity':
            self.isGreater(val, 0.)
        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)
        node.xmlSetData('initial_value', val)


    @Variables.noUndo
    def getInitialValueDensity(self, fieldId):
        """Return initial value of density"""
        return self.getInitialValue(fieldId, 'density')


    @Variables.undoLocal
    def setInitialValueDensity(self, fieldId, val):
        """Put initial value for density"""
        self.setInitialValue(fieldId, 'density', val)


    @Variables.noUndo
    def getInitialValueViscosity(self, fieldId):
        """Return initial value of viscosity"""
        return self.getInitialValue(fieldId, 'molecular_viscosity')


    @Variables.undoLocal
    def setInitialValueViscosity(self, fieldId, val):
        """Put initial value for viscosity"""
        self.setInitialValue(fieldId, 'molecular_viscosity', val)


    @Variables.noUndo
    def getInitialValueHeat(self, fieldId):
        """Return initial value of specific heat"""
        return self.getInitialValue(fieldId, 'specific_heat')


    @Variables.undoLocal
    def setInitialValueHeat(self, fieldId, val):
        """Put initial value for specific heat"""
        self.setInitialValue(fieldId, 'specific_heat', val)


    @Variables.noUndo
    def getInitialValueCond(self, fieldId):
        """Return initial value of conductivity"""
        return self.getInitialValue(fieldId, 'thermal_conductivity')


    @Variables.undoLocal
    def setInitialValueCond(self, fieldId, val):
        """Put initial value for conductivity"""
        self.setInitialValue(fieldId, 'thermal_conductivity', val)


    @Variables.noUndo
    def getInitialValueTens(self):
        """Return initial value of surface tension"""
        fieldId = 'none'
        return self.getInitialValue(fieldId, 'surface_tension')


    @Variables.undoLocal
    def setInitialValueTens(self, val):
        """Put initial value for surface tension"""
        fieldId = 'none'
        self.setInitialValue(fieldId, 'surface_tension', val)


    @Variables.noUndo
    def getInitialValueEmissivity(self, fieldId):
        """Return initial value of emissivity"""
        return self.getInitialValue(fieldId, 'emissivity')


    @Variables.undoLocal
    def setInitialValueEmissivity(self, fieldId, val):
        """Put initial value for emissivity"""
        self.setInitialValue(fieldId, 'emissivity', val)


    @Variables.noUndo
    def getInitialValueElastCoef(self, fieldId):
        """Return initial value of elasticity coefficient"""
        return self.getInitialValue(fieldId, 'elasticity')


    @Variables.undoLocal
    def setInitialValueElastCoef(self, fieldId, val):
        """Put initial value for elasticity coefficient"""
        self.setInitialValue(fieldId, 'elasticity', val)


    @Variables.noUndo
    def getFormula(self, fieldId, tag):
        """
        Return a formula for properties
        """
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId), fieldIdList)
        self.isInList(tag, self.propertiesFormulaList())
        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)
        if node:
            return node.xmlGetString('formula')
        else:
            return None


    @Variables.undoLocal
    def setFormula(self, fieldId, tag, strg):
        """
        Gives a formula for properties
        """
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId), fieldIdList)
        self.isInList(tag, self.propertiesFormulaList())
        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)
        node.xmlSetData('formula', strg)


    @Variables.undoLocal
    def setRadiativeTransferStatus(self, fieldId, status):
        """
        set status for radiative resolution transfer
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isOnOff(status)

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        childNode = node.xmlInitChildNode('particles_radiative_transfer')
        childNode.xmlSetAttribute(status = status)


    @Variables.noUndo
    def getRadiativeTransferStatus(self, fieldId):
        """
        return status for radiative resolution transfer
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.__XMLNodefields.xmlGetNode('field', field_id = fieldId)
        nodeh= node.xmlGetNode('particles_radiative_transfer')
        if nodeh == None :
            rad = self.defaultValues()['radiative']
            self.setRadiativeTransferStatus(fieldId,rad)
        rad = node.xmlGetNode('particles_radiative_transfer')['status']
        return rad


    @Variables.noUndo
    def getPropertyMode(self, fieldId, tag):
        """Return choice of node I{tag}. Choice is constant or variable"""
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId), fieldIdList)
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity', 'surface_tension'))
        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)

        c = None
        if node:
            if node['choice'] != "" and node['choice'] != None:
                c = node['choice']
                self.isInList(c, ('constant', 'user_law', 'table_law'))

        if c == None:
            c = self.defaultValues()['propertyChoice']
            if node:
                self.setPropertyMode(fieldId, tag, c)

        return c


    @Variables.undoLocal
    def setPropertyMode(self, fieldId, tag, choice):
        """Put choice in xml file's node I{tag}"""
        fieldIdList = self.getFieldIdList()
        fieldIdList.append('none')
        self.isInList(str(fieldId), fieldIdList)
        self.isInList(tag, ('density', 'molecular_viscosity',
                            'specific_heat', 'thermal_conductivity', 'surface_tension'))
        self.isInList(choice, ('constant', 'user_law', 'table_law'))

        node = self.XMLNodeproperty.xmlGetNode('property', field_id = fieldId, name=tag)
        node['choice'] = choice

        if choice == 'constant':
            node.xmlRemoveChild('formula')


    def getXMLNodefieldsNode(self):
        return self.__XMLNodefields


    def getXMLThermo(self):
        return self.__XMLThermo


    # MEG Generation related functions
    def getFormulaComponents(self, fieldId, tag):
        """
        Get the formula components for a given tag
        """

        if tag == 'density':
            return self.getFormulaRhoComponents(fieldId)

        elif tag == 'molecular_viscosity':
            return self.getFormulaMuComponents(fieldId)

        elif tag == 'specific_heat':
            return self.getFormulaCpComponents(fieldId)

        elif tag == 'thermal_conductivity':
            return self.getFormulaAlComponents(fieldId)

        elif tag == 'd_rho_d_P':
            return self.getFormuladrodpComponents(fieldId)

        elif tag == 'd_rho_d_h':
            return self.getFormuladrodhComponents(fieldId)

        elif tag == 'temperature':
            return self.getFormulaTemperatureComponents(fieldId)

        elif tag == 'surface_tension':
            return self.getFormulaStComponents()

        elif tag == 'd_Tsat_d_P':
            return self.getFormuladTsatdpComponents()

        elif tag == 'LatentHeat':
            return self.getFormulaHlatComponents()

        elif tag == 'SaturationTemperature':
            return self.getFormulaTsatComponents()

        elif 'd_Hsat_d_P_' in tag:
            return self.getFormuladHsatdpComponents(tag)

        elif 'SaturationEnthalpy' in tag:
            return self.getFormulaHsatComponents(tag)

        else:
            msg = 'Formula is not available for field %s_%s in MEG' % (tag,str(fieldId))
            raise Exception(msg)


    def getFormulaRhoComponents(self, fieldId):
        """
        User formula for density
        """
        exp = self.getFormula(fieldId, 'density')
        if not exp:
            exp = "rho = 1.8;"
        req = [('rho', 'Density')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, "enthalpy_"+str(fieldId)))
            known_fields.append((label, "enthalpy_"+str(fieldId)))

        rho0_value = self.getInitialValue(fieldId, 'density')
        symbols.append(('rho0', 'Density (reference value) = '+str(rho0_value)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormulaMuComponents(self, fieldId):
        """
        User formula for molecular viscosity
        """
        exp = self.getFormula(fieldId, 'molecular_viscosity')
        if not exp:
            exp = "mu = 4.56e-05;"
        req = [('mu', 'Molecular Viscosity')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, 'enthalpy_'+str(fieldId)))
            known_fields.append((label, 'enthalpy_'+str(fieldId)))

        mu0_val = self.getInitialValue(fieldId, 'molecular_viscosity')
        symbols.append(('mu0', 'Viscosity (reference value) = '+str(mu0_val)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormulaCpComponents(self, fieldId):
        """
        User formula for specific heat
        """
        exp = self.getFormula(fieldId, 'specific_heat')

        if not exp:
            exp = "cp = 4000.;"
        req = [('cp', 'Specific heat')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, "enthalpy_"+str(fieldId)))
            known_fields.append((label, "enthalpy_"+str(fieldId)))

        cp0_val = self.getInitialValue(fieldId, "specific_heat")
        symbols.append(('cp0', 'Specific heat (reference value) = '+str(cp0_val)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormulaAlComponents(self, fieldId):
        """
        User formula for thermal conductivity
        """
        exp = self.getFormula(fieldId, 'thermal_conductivity')
        if not exp:
            exp = "lambda = 1.e-5;"
        req = [('lambda', 'Thermal conductivity')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, 'enthalpy_'+str(fieldId)))
            known_fields.append((label, 'enthalpy_'+str(fieldId)))

        l0_val = self.getInitialValue(fieldId, 'thermal_conductivity')
        symbols.append(('lambda0', 'Thermal conductivity (reference value) = '+str(l0_val)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormulaStComponents(self):
        """
        User formula for surface tension
        """
        exp = self.getFormula('none', 'surface_tension')
        if not exp:
            exp = "sigma = 0.075;"
        req = [('sigma', 'Surface Tension')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        for fieldId in self.getFieldIdList():
            if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
                label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
                symbols.append((label, 'enthalpy_'+str(fieldId)))
                known_fields.append((label, 'enthalpy_'+str(fieldId)))

        s0_val = self.getInitialValue('none', 'surface_tension')
        symbols.append(('sigma0', 'Surface tension (reference value) = '+str(s0_val)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
              symbols.append((s, s))
              known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormulaTemperatureComponents(self, fieldId):
        """
        User formula for temperature as a function of enthalpy
        """
        label = self.m_out.getVariableLabel(str(fieldId), 'temperature')
        exp = self.getFormula(fieldId, 'temperature')
        if not exp:
            exp = label + " = 273.15;"
        req = [(label, 'temperature')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, 'enthalpy_'+str(fieldId)))
            known_fields.append((label, 'enthalpy_'+str(fieldId)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormuladrodpComponents(self, fieldId):
        """
        User formula for d(ro) / dp (compressible flow)
        """
        exp = self.getFormula(fieldId, 'd_rho_d_P')
        if not exp:
            exp = "d_rho_d_P = 0.;"
        req = [('d_rho_d_P', 'Partial derivative of density with respect to pressure')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, 'enthalpy_'+str(fieldId)))
            known_fields.append((label, 'enthalpy_'+str(fieldId)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return  exp, req, known_fields, symbols


    def getFormuladrodhComponents(self, fieldId):
        """
        User formula for d(ro) / dh (compressible flow)
        """
        exp = self.getFormula(fieldId, 'd_rho_d_h')
        if not exp:
            exp = "d_rho_d_h = 0.;"
        req = [('d_rho_d_h', 'Partial derivative of density with respect to enthalpy')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
           symbols.append(s)
           known_fields.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, 'enthalpy_'+str(fieldId)))
            known_fields.append((label, 'enthalpy_'+str(fieldId)))

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, s))
            known_fields.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, known_fields, symbols


    def getFormuladTsatdpComponents(self):

        exp = self.getFormula('none', 'd_Tsat_d_P')
        label = self.m_out.getVariableLabel('none', 'd_Tsat_d_P')
        req = [(label, 'Partial derivative of Saturation temperature with respect to pressure')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
            symbols.append(s)
            known_fields.append(s)

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))
            known_fields.append((s, self.tr("Additional species")))

        return exp, req, known_fields, symbols


    def getFormulaHlatComponents(self):

        exp = self.getFormula('none', 'LatentHeat')
        label = self.m_out.getVariableLabel('none', 'LatentHeat')
        req = [(label, 'latent heat')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
            symbols.append(s)
            known_fields.append(s)

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr('Additional species')))
            known_fields.append((s, self.tr('Additional species')))

        return exp, req, known_fields, symbols


    def getFormulaTsatComponents(self):

        exp = self.getFormula('none', 'SaturationTemperature')
        label = self.m_out.getVariableLabel('none', 'SaturationTemperature')
        req = [(label, 'SaturationTemperature')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
            symbols.append(s)
            known_fields.append(s)

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr('additional species')))
            known_fields.append((s, self.tr('additional species')))

        return exp, req, known_fields, symbols

    def getFormuladHsatdpComponents(self, tag):

        exp = self.getFormula('none', tag)

        label = self.m_out.getVariableLabel('none', tag)
        req  = [(label, 'Partial derivative of enthalpy of saturation with respect to pressure')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
            symbols.append(s)
            known_fields.append(s)

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr('Additional species')))
            known_fields.append((s, self.tr('Additional species')))

        return exp, req, known_fields, symbols


    def getFormulaHsatComponents(self, tag):

        label = self.m_out.getVariableLabel('none', tag)
        exp = self.getFormula('none', tag)

        req = [(label, 'enthalpy of saturation')]

        # Predefined Symbols
        symbols = []

        # Known fields (for equation->C translation) within the code
        known_fields = []

        for s in self.list_scalars:
            symbols.append(s)
            known_fields.append(s)

        symbols.append(('volume', 'Zone volume'))

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr('additional species')))
            known_fields.append((s, self.tr('additional species')))

        return exp, req, known_fields, symbols


    def tr(self, text):
        """
        translation
        """
        return text



#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class ThermodynamicsTestCase(ModelTest):
    """
    """
    def checkThermodynamicsInstantiation(self):
        """Check whether the ThermodynamicsModel class could be instantiated"""
        model = None
        model = ThermodynamicsModel(self.case)
        assert model != None, 'Could not instantiate ThermodynamicsModel'


    def checkGetandSetMaterials(self):
        """Check whether the ThermodynamicsModel class could set and get Materials"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setMaterials('1','user_material')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                                 <material choice="user_material"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefieldsNode() == self.xmlNodeFromString(doc),\
            'Could not set Materials'
        assert mdl.getMaterials('1') == 'user_material',\
            'Could not get Materials'


    def checkGetandSetMethod(self):
        """Check whether the ThermodynamicsModel class could set and get Method"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setMaterials('1','user_material')
        mdl.setMethod('1','user_properties')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                                 <method choice="user_properties"/>
                                 <material choice="user_material"/>
                                 <reference choice="user_material"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefieldsNode() == self.xmlNodeFromString(doc),\
            'Could not set Method'
        assert mdl.getMethod('1') == 'user_properties',\
            'Could not get Method'


    def checkGetandSetInitialValue(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValue"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValue('2','density',1.24)
        mdl.setInitialValue('2','molecular_viscosity',4.21)
        mdl.setInitialValue('2','specific_heat',6.23)
        mdl.setInitialValue('2','thermal_conductivity',885.21)
        mdl.setInitialValue('none','surface_tension',0.075)
        mdl.setInitialValue('2','emissivity',15446.2)
        mdl.setInitialValue('2','elasticity',22.2)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="Diam2" name="Diameter">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="emissivity2" name="emissivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         15446.2
                                 </initial_value>
                         </property>
                         <property choice="" field_id="2" label="elasticity2" name="elasticity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         22.2
                                 </initial_value>
                         </property>
                         <property choice="" field_id="2" label="temp2" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="density2" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                                  <initial_value>
                                         1.24
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="2" label="molecular_viscosity_2" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         4.21
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="2" label="specific_heat_2" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         6.23
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="2" label="thermal_conductivity_2" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         885.21
                                 </initial_value>
                         </property>
                         <property choice="" field_id="2" label="mass_trans2" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="none" label="Surf_tens" name="surface_tension">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         0.075
                                 </initial_value>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValue'
        assert mdl.getInitialValue('2','density') == 1.24,\
            'Could not get InitialValue density'
        assert mdl.getInitialValue('2','molecular_viscosity') == 4.21,\
            'Could not get InitialValue molecular_viscosity'
        assert mdl.getInitialValue('2','specific_heat') == 6.23,\
            'Could not get InitialValue specific_heat'
        assert mdl.getInitialValue('2','thermal_conductivity') == 885.21,\
            'Could not get InitialValue thermal_conductivity'
        assert mdl.getInitialValue('2','emissivity') == 15446.2,\
            'Could not get InitialValue emissivity'
        assert mdl.getInitialValue('2','elasticity') == 22.2,\
            'Could not get InitialValue elasticity'
        assert mdl.getInitialValue('none','surface_tension') == 0.075,\
            'Could not get InitialValue surface_tension'


    def checkGetandSetInitialValueDensity(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueDensity"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueDensity('1',8.8)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         8.8
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueDensity'
        assert mdl.getInitialValueDensity('1') == 8.8,\
            'Could not get InitialValueDensity'


    def checkGetandSetInitialValueViscosity(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueViscosity"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueViscosity('1',7.7)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         7.7
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueViscosity'
        assert mdl.getInitialValueViscosity('1') == 7.7,\
            'Could not get InitialValueViscosity'


    def checkGetandSetInitialValueHeat(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueHeat"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueHeat('1',118.712)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         118.712
                                 </initial_value>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueHeat'
        assert mdl.getInitialValueHeat('1') == 118.712,\
            'Could not get InitialValueHeat'


    def checkGetandSetInitialValueCond(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueCond"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueCond('1',456.1)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         456.1
                                 </initial_value>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueCond'
        assert mdl.getInitialValueCond('1') == 456.1,\
            'Could not get InitialValueCond'


    def checkGetandSetInitialValueTens(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueTens"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueTens('none',0.075)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="none" label="Surf_tens1" name="surface_tension">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         0.075
                                 </initial_value>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueTens'
        assert mdl.getInitialValueTens('none') == 0.075,\
            'Could not get InitialValueTens'


    def checkGetandSetInitialValueEmissivity(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueEmissivity"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueEmissivity('2',0.008)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="Diam2" name="Diameter">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="emissivity2" name="emissivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         0.008
                                 </initial_value>
                         </property>
                         <property choice="" field_id="2" label="elasticity2" name="elasticity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="temp2" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="density2" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="molecular_viscosity_2" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="specific_heat_2" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="thermal_conductivity_2" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="mass_trans2" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueEmissivity'
        assert mdl.getInitialValueEmissivity('2') == 0.008,\
            'Could not get InitialValueEmissivity'


    def checkGetandSetInitialValueElastCoef(self):
        """Check whether the ThermodynamicsModel class could set and get InitialValueElastCoef"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'solid', 'on', 'on', 'off', 2)
        mdl = ThermodynamicsModel(self.case)
        mdl.setInitialValueElastCoef('2',0.42)
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="Diam2" name="Diameter">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="emissivity2" name="emissivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="elasticity2" name="elasticity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <initial_value>
                                         0.42
                                 </initial_value>
                         </property>
                         <property choice="" field_id="2" label="temp2" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="density2" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="molecular_viscosity_2" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="specific_heat_2" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="2" label="thermal_conductivity_2" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="2" label="mass_trans2" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set InitialValueElastCoef'
        assert mdl.getInitialValueElastCoef('2') == 0.42,\
            'Could not get InitialElastValueCoef'


    def checkGetandSetFormula(self):
        """Check whether the ThermodynamicsModel class could set and get Formula"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setFormula('1','density','2.123')
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <formula>
                                         2.123
                                 </formula>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set Formula'
        assert mdl.getFormula('1', 'density') == '2.123',\
            'Could not get Formula'


    def checkGetandSetRadiativeTransferStatus(self):
        """Check whether the ThermodynamicsModel class could set and get RadiativeTransferStatus"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setRadiativeTransferStatus('1','on')
        doc = '''<fields>
                         <field field_id="1" label="Field1">
                                 <type choice="continuous"/>
                                 <carrier_field field_id="off"/>
                                 <phase choice="liquid"/>
                                 <hresolution status="on"/>
                                 <compressible status="off"/>
                                 <particles_radiative_transfer status="on"/>
                         </field>
                 </fields>'''
        assert mdl.getXMLNodefieldsNode() == self.xmlNodeFromString(doc),\
            'Could not set RadiativeTransferStatus'
        assert mdl.getRadiativeTransferStatus('1') == 'on',\
            'Could not get RadiativeTransferStatus'


    def checkGetandSetPropertyMode(self):
        """Check whether the ThermodynamicsModel class could set and get PropertyMode"""
        MainFieldsModel(self.case).addField()
        mdl = ThermodynamicsModel(self.case)
        mdl.setPropertyMode('1','density','user_law')
        doc = '''<properties>
                         <property choice="" field_id="1" label="temp1" name="temperature">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="user_law" field_id="1" label="density1" name="density">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="molecular_viscosity_1" name="molecular_viscosity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="specific_heat_1" name="specific_heat">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="constant" field_id="1" label="thermal_conductivity_1" name="thermal_conductivity">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                         <property choice="" field_id="1" label="mass_trans1" name="mass_trans">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </property>
                 </properties>'''
        assert mdl.XMLNodeproperty == self.xmlNodeFromString(doc),\
            'Could not set PropertyMode'
        assert mdl.getPropertyMode('1','density') == 'user_law',\
            'Could not get PropertyMode'


def suite():
    testSuite = unittest.makeSuite(ThermodynamicsTestCase, "check")
    return testSuite


def runTest():
    print("ThermodynamicsTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())


