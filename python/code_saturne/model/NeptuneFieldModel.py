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

import traceback
import sys, unittest, copy

from code_saturne.model.Common import LABEL_LENGTH_MAX
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *

class NeptuneField(Variables, Model):
    """
    This class defines the properties of one single neptune_cfd field
    """
    default_carrier_id = "off"

    def __init__(self, case, f_id):

        self.case = case # Could be a global variable (is actually a singleton)
        self._f_id = str(f_id)

        self._label = ""
        self._phase = ""
        self._flow_type = ""
        self._carrier_id =  ""
        self._enthalpy_model = ""
        self._compressible = ""

        if (self._f_id not in ["none", "off", "all"]):
            thermo_node = self.case.xmlGetNode('thermophysical_models')
            fields_node = thermo_node.xmlInitNode('fields')
            self._xml_node = fields_node.xmlInitNode('field', field_id = self.f_id)
            self._xml_variable_node = thermo_node.xmlInitNode('variables')
            self._xml_property_node = thermo_node.xmlInitNode('properties')
        else:
            self._label = self._f_id
            if self._f_id == "all":
                self._flow_type = "continuous"


    def __str__(self):
        return str(self._label)

    def __repr__(self):
        return ("Field(f_id={0}, label={1}, phase={2}, flow_type={3},"
               "carrier_id={4}, enthalpy_model={5}, compressible={6})").format(
                       self.f_id, self._label, self._phase, self._flow_type, self._carrier_id,
                       self._enthalpy_model, self._compressible)

    def load_from_xml(self):
        try:
            self._label = self.label
            self._phase = self.phase
            self._flow_type = self.flow_type
            self._carrier_id = self.carrier_id
            self._enthalpy_model = self.enthalpy_model
            self._compressible = self.compressible
        except Exception as e:
            raise Exception('problem during load_from_xml') from e

    @property
    def f_id(self):
        return str(self._f_id)

    @f_id.setter
    def f_id(self, value):
        try:
            self.isPositiveInt(value)
        except ValueError:
            assert(type(value) == str)
            assert(value.isdigit())
        for node in self.case.xmlGetNodeWithAttrList("field_id", field_id=self._f_id):
            node["field_id"] = str(value)
        self._f_id = str(value)

    @property
    @Variables.noUndo
    def label(self):
        if self.f_id in ["none", "off", "all"]:
            return self.f_id
        return self._xml_node.xmlGetAttribute("label", default=self._label)

    @label.setter
    @Variables.undoLocal
    def label(self, value):
        #TODO add rules from MainFieldsModel
        if self.isStr(value):
            old_label = self._label
            self._label = value[:LABEL_LENGTH_MAX] # Why ?
            # Find a way to check that value is not in existing labels
            self._xml_node.xmlSetAttribute(label=value)
            for node in self.case.xmlGetNodeList('variable') \
                    + self.case.xmlGetNodeList('property'):
                if node['field_id'] == str(self.f_id):
                    li = node['label'].rsplit(old_label, 1)
                    node['label'] = value.join(li)

    @property
    def flow_type(self):
        if self._label == "all":
            return "continuous" # ugly hack...
        node = self._xml_node.xmlGetNode("type")
        if node != None:
            return node.xmlGetAttribute("choice", default=self._flow_type)
        else:
            return self._flow_type

    @flow_type.setter
    def flow_type(self, value):
        #TODO add rules from MainFieldsModel
        self.isInList(value, ("continuous", "dispersed", "auto"))
        old_value = self._flow_type
        self._flow_type = value
        child_node = self._xml_node.xmlInitChildNode("type")
        child_node.xmlSetAttribute(choice = value)

        if (value != "continuous" and self.carrier_id == "off"):
            self.carrier_id = NeptuneField.default_carrier_id

    @property
    def phase(self):
        node = self._xml_node.xmlGetNode("phase")
        if node != None:
            return self._xml_node.xmlGetNode("phase").xmlGetAttribute("choice", default=self._phase)
        else:
            return self._phase

    @phase.setter
    def phase(self, value):
        #TODO add rules from MainFieldsModel
        self.isInList(value, ("liquid", "solid", "gas"))
        self._phase = value
        child_node = self._xml_node.xmlInitChildNode("phase")
        child_node.xmlSetAttribute(choice = value)
        if value == "solid":
            self.flow_type = "dispersed"

    @property
    def enthalpy_model(self):
        node = self._xml_node.xmlGetNode("hresolution")
        if node != None:
            return node.xmlGetAttribute("model", default=self._enthalpy_model)
        else:
            return self._enthalpy_model

    @enthalpy_model.setter
    def enthalpy_model(self, value):
        self.isInList(value, ("off", "total_enthalpy", "specific_enthalpy"))

        self._enthalpy_model = value
        child_node = self._xml_node.xmlInitChildNode("hresolution")
        child_node.xmlSetAttribute(model = value)

        if value == "off":
            child_node.xmlSetAttribute(status = "off")
            self.removeVariableProperty("variable", self._xml_variable_node, self.f_id, "enthalpy")
            self.removeVariableProperty("property", self._xml_property_node, self.f_id, "temperature")
        else:
            child_node.xmlSetAttribute(status = "on")
            self.setNewVariableProperty("variable", "", self._xml_variable_node, self.f_id, "enthalpy", "enthalpy_"+self.label, post = True)
            self.setNewVariableProperty("property", "", self._xml_property_node, self.f_id, "temperature", "temp_"+self.label, post = True)

    @property
    def carrier_id(self):
        node = self._xml_node.xmlGetNode("carrier_field")
        if node != None:
            return node.xmlGetAttribute("field_id", default=self._carrier_id)
        else:
            return self._carrier_id # or return "off" directly ?

    @carrier_id.setter
    def carrier_id(self, field_id):
        #TODO add check on existence of carrier field
        child_node = self._xml_node.xmlInitChildNode("carrier_field")
        self._carrier_id = field_id
        child_node.xmlSetAttribute(field_id = str(field_id))

    @property
    def compressible(self):
        node = self._xml_node.xmlGetNode("compressible")
        if node != None:
            return self._xml_node.xmlGetNode("compressible").xmlGetAttribute("status", default=self._compressible)
        else:
            return self._compressible

    @compressible.setter
    def compressible(self, value):
        #TODO add rules from MainFieldsModel
        self.isOnOff(value)
        self._compressible = value
        child_node = self._xml_node.xmlInitChildNode("compressible")
        child_node.xmlSetAttribute(status = value)
        if value == "on":
            self.setNewVariableProperty("property", "", self._xml_property_node, str(self.f_id), "d_rho_d_P", "drho_dP_"+self.label)
            self.setNewVariableProperty("property", "", self._xml_property_node, str(self.f_id), "d_rho_d_h", "drho_dh_"+self.label)
        else :
            self.removeVariableProperty("property", self._xml_property_node, str(self.f_id), "d_rho_d_P")
            self.removeVariableProperty("property", self._xml_property_node, str(self.f_id), "d_rho_d_h")

    @Variables.undoLocal
    def createVariableProperties(self, multiregime=False):
        """
        add XML variable and properties
        """

        field_name = self.label
        field_number = self.f_id

        self.setNewVariableProperty("variable", "", self._xml_variable_node, field_number,
                                                    "volume_fraction", "vol_f_" + field_name, post=True)
        self.setNewVariableProperty("variable", "", self._xml_variable_node, field_number, "velocity",
                                                    "U_" + field_name, dim='3', post=True)
        if self.enthalpy_model != "off":
            self.setNewVariableProperty("variable", "", self._xml_variable_node, field_number, "enthalpy",
                                                        "enthalpy_" + field_name, post=True)

        # Physical properties are set by default to "constant" to avoid uninitialized states with the GUI
        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "density", "density_"+field_name)
        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "molecular_viscosity", "molecular_viscosity_"+field_name)
        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "specific_heat", "specific_heat_"+field_name)
        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "thermal_conductivity", "thermal_conductivity_"+field_name)

        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "mass_trans", "mass_trans_"+field_name)
        self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "wall_distance", "y_plus_"+field_name, support = "boundary")
        if self.compressible == "on":
           self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "drodp", "drodp_"+field_name)
           self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "drodh", "drodh_"+field_name)
        if self.enthalpy_model != "off":
           self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "temperature", "temp_"+field_name, post = True)
        if self.flow_type == "dispersed" or multiregime:
           self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "diameter", "diam_"+field_name)
           self.setNewVariableProperty("property", "", self._xml_property_node, field_number, "drift_component", "drift_component_"+field_name, dim='3')

#
#    @property
#    def thermodynamic_state(self):
#        return self.thermodynamic_state
#
#    @thermodynamic_state.setter
#    def thermodynamic_state(self, phase):
#        self.isInList(phase, ("liquid", "solid", "gas"))
#        if self.has_xml:
#            childNode = self._xml_node.xmlInitChildNode('phase')
#            if childNode != None:
#                oldstatus = childNode['choice']
#                if phase != oldstatus:
#                    if phase == "solid": #TODO : find a way to gather clearly this type of rules
#                      self.setCriterion(self.f_id, "dispersed")
#
#            childNode.xmlSetAttribute(choice = phase)
