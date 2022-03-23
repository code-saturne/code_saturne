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

import sys, unittest
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import MainFieldsModel

from code_saturne.model.ProfilesModel import ProfilesModel
from code_saturne.model.TimeAveragesModel import TimeAveragesModel

#-------------------------------------------------------------------------------
# Constructor
#-------------------------------------------------------------------------------


class NonCondensableModel(MainFieldsModel, Variables, Model):
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
        self.XMLNodeNonCondensable = self.XMLNodethermo.xmlInitNode('non_condensable_list')
        self.node_anal       = self.case.xmlInitNode('analysis_control')
        self.node_average    = self.node_anal.xmlInitNode('time_averages')
        self.node_profile    = self.node_anal.xmlInitNode('profiles')


    def defaultValues(self):
        default = {}

        default['typeNonCond'] = "Air"
        default['MolarMass']   = 0.02896
        default['Cobin1']      = 2.2e-5
        default['Cobin2']      = 1.5
        return default


    def HydrogenValues(self):
        """
        Return in a dictionnary which contains default values
        """
        hydrogen = {}
        hydrogen['MolarMass'] = 0.002016
        hydrogen['Cobin1'] = 0.78e-4
        hydrogen['Cobin2'] = 1.75
        return hydrogen


    def AirValues(self):
        """
        Return in a dictionnary which contains default values
        """
        air = {}
        air['MolarMass']   = 0.02896
        air['Cobin1']      = 2.2e-5
        air['Cobin2']      = 1.5
        return air


    def NitrogenValues(self):
        """
        Return in a dictionnary which contains default values
        """
        nitrogen = {}
        nitrogen['MolarMass'] = 28.0134e-3
        nitrogen['Cobin1'] = 0.227e-4
        nitrogen['Cobin2'] = 1.75
        return nitrogen


    def HeliumValues(self):
        """
        Return in a dictionnary which contains default values
        """
        helium = {}
        helium['MolarMass'] = 4.003e-3
        helium['Cobin1'] = 0.73e-4
        helium['Cobin2'] = 1.75
        return helium

    def OxygenValues(self):
        """
        Return in a dictionnary which contains default values
        """
        oxygen = {}
        oxygen['MolarMass'] = 32.e-3
        oxygen['Cobin1'] = 0.24e-4
        oxygen['Cobin2'] = 1.71
        return oxygen

    def checkNonCondensableRequirements(self):
        """
        Check if the requirements for non condensable gas addition are met
        - Energy equation must be activated for all gases
        """
        if self.getGasPhaseList() in [[], None]:
            return False
        for field_id in self.getGasPhaseList():
            if self.getEnergyModel(field_id) != "off":
                return True
        return False

    def getNonCondensableLabelList(self):
        """
        Return the non condensable list
        """
        list = []
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            list.append(node['label'])
        return list


    def getNonCondensableNameList(self):
        """
        Return the non condensable list
        """
        list = []
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            list.append(node['name'])
        return list


    def getNonCondensableByFieldId(self, FieldId):
        """
        Return the non condensable name list for a fieldId
        """
        self.isInList(str(FieldId),self.getFieldIdList())
        list = []
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if self.getNonCondFieldId(node['name']) == str(FieldId) :
                list.append(node['name'])
        return list


    @Variables.undoLocal
    def addNonCondensable(self):
        """
        Add a non condensable
        """
        fieldId = self.getFirstGasField()
        field_name = self.getFieldLabelsList()[int(fieldId)-1]
        type = self.defaultValues()['typeNonCond']
        label   = type + "_1" + "_" + field_name
        if label in self.getNonCondensableLabelList() :
           labelNumber = 1
           label = type + "_" + str(labelNumber) + "_" + field_name
           while label in self.getNonCondensableLabelList() :
               labelNumber += 1
               label = type + "_" + str(labelNumber)+ "_" + field_name

        name = "mass_fraction_non_condensable_gas_" + str(len(self.getNonCondensableNameList())+1)

        Variables(self.case).setNewVariableProperty("variable", "", self.XMLNodeNonCondensable, fieldId, name, label)

        # for non condensable we need use cathare2 or cathare tables
        from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
        ref_material = "Water"
        for m in ("Cathare2", "Cathare"):
            if m in self.eos.getFluidMethods(ref_material):
                for i in (1,2):
                    ThermodynamicsModel(self.case).setMaterials(i, ref_material)
                    ThermodynamicsModel(self.case).setMethod(i, m)

                break
        del ThermodynamicsModel

        return name


    @Variables.undoLocal
    def setNonCondLabel(self, name, label):
        """
        set label for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               node['label'] = label


    @Variables.noUndo
    def getNonCondLabel(self, name):
        """
        return label for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               return node['label']


    @Variables.undoLocal
    def setNonCondFieldId(self, name, field):
        """
        set field Id for non condensable
        """
        fieldId = -1
        for id in self.getGasPhaseList() :
            if field == self.getLabel(id) :
               fieldId = id
        self.isInList(fieldId, self.getGasPhaseList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               node['field_id'] = fieldId


    @Variables.noUndo
    def getNonCondFieldId(self, name):
        """
        return field Id for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               return node['field_id']


    @Variables.undoLocal
    def setNonCondType(self, name, type):
        """
        set type for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())

        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               childNode = node.xmlInitChildNode('Type')
               childNode.xmlSetAttribute(choice = type)
               self.setIncTyp(name,type)


    @Variables.noUndo
    def getNonCondType(self, name):
        """
        return type for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())

        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               childNode = node.xmlGetNode('Type')
               if childNode is None :
                   type = self.defaultValues()['typeNonCond']
                   self.setNonCondType(name, type)
               type = node.xmlGetNode('Type')['choice']
               return type


    @Variables.undoGlobal
    def setIncTyp(self, name, Inctyp):
        masmol = 0.
        cobin1 = 0.
        cobin2 = 0.
        if Inctyp == 'H2' :
            masmol = self.HydrogenValues()['MolarMass']
            cobin1 = self.HydrogenValues()['Cobin1']
            cobin2 = self.HydrogenValues()['Cobin2']
        elif Inctyp == 'N2' :
            masmol = self.NitrogenValues()['MolarMass']
            cobin1 = self.NitrogenValues()['Cobin1']
            cobin2 = self.NitrogenValues()['Cobin2']
        elif Inctyp == 'HE' :
            masmol = self.HeliumValues()['MolarMass']
            cobin1 = self.HeliumValues()['Cobin1']
            cobin2 = self.HeliumValues()['Cobin2']
        elif Inctyp == 'O2' :
            masmol = self.OxygenValues()['MolarMass']
            cobin1 = self.OxygenValues()['Cobin1']
            cobin2 = self.OxygenValues()['Cobin2']
        elif Inctyp == 'Air' :
            masmol = self.AirValues()['MolarMass']
            cobin1 = self.AirValues()['Cobin1']
            cobin2 = self.AirValues()['Cobin2']

        self.setNonCondMassMol(name, masmol)
        self.setNonCondCobin1( name, cobin1)
        self.setNonCondCobin2( name, cobin2)


    @Variables.undoLocal
    def setNonCondMassMol(self, name, MassMol):
        """
        set molar mass for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        self.isPositiveFloat(MassMol)

        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               node.xmlSetData('MolarMass', MassMol)


    @Variables.noUndo
    def getNonCondMassMol(self, name):
        """
        return mass molar for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               massmol = node.xmlGetDouble('MolarMass')
               if massmol is None :
                   massmol = self.defaultValues()['MolarMass']
                   self.setNonCondMassMol(name, massmol)
               return massmol


    @Variables.undoLocal
    def setNonCondCobin1(self, name, cobin1):
        """
        set cobin 1 for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        self.isPositiveFloat(cobin1)

        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               node.xmlSetData('Cobin1', cobin1)


    @Variables.noUndo
    def getNonCondCobin1(self, name):
        """
        return cobin 1 for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               cobin1 = node.xmlGetDouble('Cobin1')
               if cobin1 is None :
                   cobin1 = self.defaultValues()['Cobin1']
                   self.setNonCondCobin1(name, cobin1)
               return cobin1


    @Variables.undoLocal
    def setNonCondCobin2(self, name, cobin2):
        """
        set cobin 2 for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        self.isPositiveFloat(cobin2)

        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               node.xmlSetData('Cobin2', cobin2)


    @Variables.noUndo
    def getNonCondCobin2(self, name):
        """
        return cobin 2 for non condensable
        """
        self.isInList(name,self.getNonCondensableNameList())
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            if node['name'] == name :
               cobin2 = node.xmlGetDouble('Cobin2')
               if cobin2 is None :
                   cobin2 = self.defaultValues()['Cobin2']
                   self.setNonCondCobin2(name, cobin2)
               return cobin2


    @Variables.undoGlobal
    def deleteNonCondensable(self, number):
        """
        Suppress non condensable from list
        """
        name = "mass_fraction_non_condensable_gas_" + str(number)
        self.isInList(name,self.getNonCondensableNameList())

        # delete non condensable
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            try :
               if node['name'] == name :
                  node.xmlRemoveNode()
            except :
               pass

        #suppress profile
        for node in reversed(self.node_profile.xmlGetNodeList('profile')) :
            suppress = 0
            for child in node.xmlGetNodeList('var_prop') :
                if (child['name'] == name) :
                    suppress = 1
            if suppress == 1 :
                name = node['name']
                ProfilesModel(self.case).deleteProfile(name)

        #suppress average
        for node in reversed(self.node_average.xmlGetNodeList('time_average')) :
            suppress = 0
            for child in node.xmlGetNodeList('var_prop') :
                if (child['name'] == name) :
                    suppress = 1
                    break
            if suppress == 1 :
                name = node['name']
                TimeAveragesModel(self.case).deleteTimeAverage(name)

        # update name for other non condensable in XML file
        index = 1
        for node in self.XMLNodeNonCondensable.xmlGetNodeList('variable'):
            try :
               if index >= number :
                  oldname = "mass_fraction_non_condensable_gas_" + str(index+1)
                  name = "mass_fraction_non_condensable_gas_" + str(index)
                  node['name'] = name
                  for n in self.case.xmlGetNodeList('var_prop', name=oldname):
                      n['name'] = name
            except :
               pass
            index += 1

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class NonCondensableTestCase(ModelTest):
    """
    """
    def checkNonCondensableInstantiation(self):
        """Check whether the NonCondensableModel class could be instantiated"""
        model = None
        model = NonCondensableModel(self.case)
        assert model != None, 'Could not instantiate NonCondensableModel'


    def checkGetNonCondensableLabelList(self):
        """Check whether the  NonCondensableModel class could get the NonCondensableLabelList"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        assert mdl.getNonCondensableLabelList() == ['Air_1'],\
            'Could not get NonCondensableLabelList'


    def checkGetNonCondensableNameList(self):
        """Check whether the  NonCondensableModel class could get the NonCondensableNameList"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        assert mdl.getNonCondensableNameList() == ['mass_fraction_non_condensable_gas_1'],\
            'Could not get NonCondensableNameList'


    def checkGetNonCondensableByFieldId(self):
        """Check whether the  NonCondensableModel class could get the NonCondensableByFieldId"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        assert mdl.getNonCondensableByFieldId('1') == ['mass_fraction_non_condensable_gas_1'],\
            'Could not get NonCondensableByFieldId'


    def checkaddNonCondensable(self):
        """Check whether the  NonCondensableModel class could addNonCondensable"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not addNonCondensable'


    def checkGetandSetNonCondLabel(self):
        """Check whether the NonCondensableModel class could set and get NonCondLabel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setNonCondLabel('mass_fraction_non_condensable_gas_1','example_label')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="example_label" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not set NonCondLabel'
        assert mdl.getNonCondLabel('mass_fraction_non_condensable_gas_1') == 'example_label',\
            'Could not get NonCondLabel'


    def checkGetandSetNonCondFieldId(self):
        """Check whether the NonCondensableModel class could set and get NonCondFieldId"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setNonCondFieldId('mass_fraction_non_condensable_gas_1','field1')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not set NonCondFieldId'
        assert mdl.getNonCondFieldId('mass_fraction_non_condensable_gas_1') == '1',\
            'Could not get NonCondFieldId'


    def checkGetandSetNonCondType(self):
        """Check whether the NonCondensableModel class could set and get NonCondType"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setNonCondType('mass_fraction_non_condensable_gas_1','N2')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <Type choice="N2"/>
                                 <MolarMass>
                                         0.028
                                 </MolarMass>
                                 <Cobin1>
                                         2.27e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.75
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not set NonCondType'
        assert mdl.getNonCondType('mass_fraction_non_condensable_gas_1') == 'N2',\
            'Could not get NonCondType'


    def checkGetandSetNonCondMassMol(self):
        """Check whether the NonCondensableModel class could set and get NonCondMassMoll"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setNonCondMassMol('mass_fraction_non_condensable_gas_1',8510.1)
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         8510.1
                                 </MolarMass>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not set NonCondMassMol'
        assert mdl.getNonCondMassMol('mass_fraction_non_condensable_gas_1') == 8510.1,\
            'Could not get NonCondMassMol'


    def checkGetandSetNonCondCobin1and2(self):
        """Check whether the NonCondensableModel class could set and get NonCondCobin1and2"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setNonCondCobin1('mass_fraction_non_condensable_gas_1',5.5)
        mdl.setNonCondCobin2('mass_fraction_non_condensable_gas_1',6.6)
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <Cobin1>
                                         5.5
                                 </Cobin1>
                                 <Cobin2>
                                         6.6
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not set NonCondCobin1and2'
        assert mdl.getNonCondCobin1('mass_fraction_non_condensable_gas_1') == (5.5),\
            'Could not get NonCondCobin1'
        assert mdl.getNonCondCobin2('mass_fraction_non_condensable_gas_1') == (6.6),\
            'Could not get NonCondCobin2'


    def checkdeleteNonCondensable(self):
        """Check whether the NonCondensableModel class could deleteNonCondensable"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.addNonCondensable()
        mdl.deleteNonCondensable(1)
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_2" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not delete noncondensable'


    def checksetIncTyp(self):
        """Check whether the NonCondensableModel class could deleteNonCondensable"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'gas', 'on', 'on', 'off', 1)
        mdl = NonCondensableModel(self.case)
        mdl.addNonCondensable()
        mdl.setIncTyp('mass_fraction_non_condensable_gas_1','N2')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         0.028
                                 </MolarMass>
                                 <Cobin1>
                                         2.27e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.75
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not setIncTyp for N2'

        mdl.setIncTyp('mass_fraction_non_condensable_gas_1','H2')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         0.002016
                                 </MolarMass>
                                 <Cobin1>
                                         7.8e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.75
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not setIncTyp for H2'

        mdl.setIncTyp('mass_fraction_non_condensable_gas_1','HE')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         0.004003
                                 </MolarMass>
                                 <Cobin1>
                                         7.3e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.75
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not setIncTyp for HE'

        mdl.setIncTyp('mass_fraction_non_condensable_gas_1','O2')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         0.032
                                 </MolarMass>
                                 <Cobin1>
                                         2.4e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.71
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not setIncTyp for O2'

        mdl.setIncTyp('mass_fraction_non_condensable_gas_1','Air')
        doc = '''<non_condensable_list>
                         <variable field_id="1" label="Air_1" name="mass_fraction_non_condensable_gas_1">
                                 <listing_printing status="on"/>
                                 <postprocessing_recording status="on"/>
                                 <MolarMass>
                                         0.02896
                                 </MolarMass>
                                 <Cobin1>
                                         2.2e-05
                                 </Cobin1>
                                 <Cobin2>
                                         1.5
                                 </Cobin2>
                         </variable>
                 </non_condensable_list>'''
        assert mdl.XMLNodeNonCondensable == self.xmlNodeFromString(doc),\
            'Could not setIncTyp for Air'


def suite():
    testSuite = unittest.makeSuite(NonCondensableTestCase, "check")
    return testSuite


def runTest():
    print("NonCondensableTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
