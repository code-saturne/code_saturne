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
from code_saturne.Base.XMLvariables import Model
from code_saturne.Base.XMLengine import *
from code_saturne.Pages.MainFieldsModel import *       # TODO change
import copy

#-------------------------------------------------------------------------------
# Turbulence model description
#-------------------------------------------------------------------------------

class TurbulenceModelsDescription:
    """
    """
    continuousTurbulenceModels = ['none', 'mixing_length',
                                  'k-epsilon', 'k-epsilon_linear_production',
                                  'rij-epsilon_ssg', 'rij-epsilon_ebrsm',
                                  'les_smagorinsky', 'les_wale']

    dispersedTurbulenceModels  = ['none','tchen','q2-q12', 'r2-q12', 'r2-r12-tchen']

    continuousCouplingModels = ['none','separate_phase','separate_phase_cond']
    dispersedCouplingModels  = ['none','small_inclusions','large_inclusions']

    ThermalTurbFluxModels = ['sgdh', 'ggdh']

    turbulenceVariables = {}
    turbulenceVariables['none'] = []
    turbulenceVariables['mixing_length'] = []
    turbulenceVariables['k-epsilon'] = ['TurbKineEner_k', 'TurbDissip']
    turbulenceVariables['k-epsilon_linear_production'] = ['TurbKineEner_k', 'TurbDissip']
    turbulenceVariables['rij-epsilon_ssg'] = ['ReynoldsStress', 'TurbDissip']
    turbulenceVariables['rij-epsilon_ebrsm'] = ['ReynoldsStress', 'TurbDissip']
    turbulenceVariables['les_smagorinsky'] = []
    turbulenceVariables['les_wale'] = []
    turbulenceVariables['tchen'] = []
    turbulenceVariables['q2-q12'] = ['TurbKineEner_q2', 'Covariance_q12']
    turbulenceVariables['r2-q12'] = ['ReynoldsStress','Covariance_q12']
    turbulenceVariables['r2-r12-tchen'] = ['ReynoldsStress',
                                           'R12XX','R12XY','R12XZ','R12YY','R12YZ','R12ZZ']

    turbulenceVariables['all'] = turbulenceVariables['k-epsilon'] \
                               + turbulenceVariables['k-epsilon_linear_production'] \
                               + turbulenceVariables['rij-epsilon_ssg'] \
                               + turbulenceVariables['rij-epsilon_ebrsm'] \
                               + turbulenceVariables['les_smagorinsky'] \
                               + turbulenceVariables['les_wale'] \
                               + turbulenceVariables['q2-q12'] \
                               + turbulenceVariables['r2-q12'] \
                               + turbulenceVariables['r2-r12-tchen']

    turbulenceProperties = {}
    turbulenceProperties['none'] = []
    turbulenceProperties['mixing_length'] = ["turb_viscosity"]
    turbulenceProperties['k-epsilon'] = ["turb_viscosity"]
    turbulenceProperties['k-epsilon_linear_production'] = ["turb_viscosity"]
    turbulenceProperties['rij-epsilon_ssg'] = ["turb_viscosity"]
    turbulenceProperties['rij-epsilon_ebrsm'] = ["turb_viscosity"]
    turbulenceProperties['les_smagorinsky'] = ["turb_viscosity"]
    turbulenceProperties['les_wale'] = ["turb_viscosity"]
    turbulenceProperties['tchen'] = ["TurbKineEner_q2", "Covariance_q12", "turb_viscosity"]
    turbulenceProperties['q2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-r12-tchen'] = ["turb_viscosity"]


#-------------------------------------------------------------------------------
# Main turbulence model class
#-------------------------------------------------------------------------------

class TurbulenceModel(MainFieldsModel):
    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constuctor.
        """
        #
        # XML file parameters
        MainFieldsModel.__init__(self, case)
        self.case = case
        self.XMLClosure      = self.case.xmlGetNode('closure_modeling')
        self.XMLturbulence   = self.XMLClosure.xmlInitNode('turbulence')
        self.XMLNodeVariable = self.XMLturbulence.xmlInitNode('variables')
        self.XMLNodeproperty = self.XMLturbulence.xmlInitNode('properties')


    def defaultValues(self):
        default = {}

        default['length_scale']     = 1.0
        default['two_way_coupling'] = "none"
        default['model']            = TurbulenceModelsDescription.continuousTurbulenceModels[0]
        default['turb_flux']        = TurbulenceModelsDescription.ThermalTurbFluxModels[0]
        default['length']           = 10.0
        return default


    @Variables.undoGlobal
    def setTurbulenceModel(self, fieldId, model):
        """
        put turbulence model for fieldId
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        field_name = self.getFieldLabelsList()[int(fieldId)-1]

        criterion = self.getCriterion(fieldId)
        if criterion == "continuous":
           self.isInList(model, TurbulenceModelsDescription.continuousTurbulenceModels)
        else:
           self.isInList(model, TurbulenceModelsDescription.dispersedTurbulenceModels)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        oldmodel = ""
        if node != None:
            oldmodel = self.defaultValues()['model']

        if node == None:
            if criterion == "continuous":
                self.XMLturbulence.xmlInitChildNode('field', field_id = fieldId,
                                                    model = model,
                                                    turb_flux = self.defaultValues()['turb_flux'],
                                                    two_way_coupling = TurbulenceModelsDescription.continuousCouplingModels[0])
            else:
                self.XMLturbulence.xmlInitChildNode('field', field_id = fieldId,
                                                    model = model,
                                                    turb_flux = self.defaultValues()['turb_flux'],
                                                    two_way_coupling = TurbulenceModelsDescription.dispersedCouplingModels[0])
        else:
            oldmodel = node['model']
            node['model'] = model

        if oldmodel != model:
           if oldmodel != "":
               # erase old variables and properties from XML
               for var in TurbulenceModelsDescription.turbulenceVariables[oldmodel]:
                   Variables(self.case).removeVariableProperty("variable", self.XMLNodeVariable, fieldId, var)

               for var in TurbulenceModelsDescription.turbulenceProperties[oldmodel]:
                   Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldId, var)

           # add new variables and properties from XML
           for var in TurbulenceModelsDescription.turbulenceVariables[model]:
             if var == "ReynoldsStress":
                 Variables(self.case).setNewVariableProperty("variable", "",
                                                             self.XMLNodeVariable,
                                                             fieldId,
                                                             var,
                                                             var+"_"+field_name,
                                                             dim=6)
             else:
                 Variables(self.case).setNewVariableProperty("variable", "",
                                                             self.XMLNodeVariable,
                                                             fieldId,
                                                             var,
                                                             var+"_"+field_name)

           for var in TurbulenceModelsDescription.turbulenceProperties[model]:
               Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, var, var+"_"+field_name)

           if oldmodel == "mixing_length":
               node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
               noder = node.xmlGetNode('length_scale')
               noder.xmlRemoveNode()

           # update other field if continuous and set to none
           if model == "none" and criterion == "continuous":
               for id in self.getFieldIdList():
                   if self.getCriterion(id) == "dispersed":
                       if self.getCarrierField(id) == str(fieldId):
                           self.setTurbulenceModel(id, TurbulenceModelsDescription.dispersedTurbulenceModels[0])
                           if criterion == "continuous":
                               self.setTwoWayCouplingModel(id, TurbulenceModelsDescription.continuousCouplingModels[0])
                           else:
                               self.setTwoWayCouplingModel(id, TurbulenceModelsDescription.continuousCouplingModels[0])


    @Variables.noUndo
    def getTurbulenceModel(self, fieldId):
        """
        get turbulence model for fieldId
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        criterion = self.getCriterion(fieldId)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node == None:
            model = ""
            if criterion == "continuous":
               model = TurbulenceModelsDescription.continuousTurbulenceModels[0]
            else:
               model = TurbulenceModelsDescription.dispersedTurbulenceModels[0]
            self.setTurbulenceModel(fieldId, model)
            node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        model = node['model']

        return model

    @Variables.undoLocal
    def setThermalTurbulentFlux(self, fieldId, model):

        self.isInList(str(fieldId), self.getFieldIdList())

        self.isInList(model,TurbulenceModelsDescription.ThermalTurbFluxModels)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)

        critrerion =  self.getCriterion(fieldId)
        if node == None:
            if criterion == "continuous":
                self.XMLturbulence.xmlInitChildNode('field',
                                                    field_id = fieldId,
                                                    model = self.defaultValues()['model'],
                                                    turb_flux = model,
                                                    two_way_coupling = TurbulenceModelsDescription.continuousCouplingModels[0])
            else:
                self.XMLturbulence.xmlInitChildNode('field',
                                                    field_id = fieldId,
                                                    model = self.defaultValues()['model'],
                                                    turb_flux = model,
                                                    two_way_coupling = TurbulenceModelsDescription.dispersedCouplingModels[0])

        node['turb_flux'] = model


    @Variables.noUndo
    def getThermalTurbulentFlux(self, fieldId):

        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node == None:
            if self.getEnergyResolution(fieldId) == 'on':
                self.setThermalTurbulentFlux(fieldId,
                                             self.defaultValues()['turb_flux'])
            else:
                self.setThermalTurbulentFlux(fieldId, 'none')
        else:
            if self.getEnergyResolution(fieldId) != 'on':
                node['turb_flux'] = 'none'
            elif node['turb_flux'] == 'none':
               node['turb_flux'] = 'sgdh'


        model = node['turb_flux']
        return model


    @Variables.undoLocal
    def setTwoWayCouplingModel(self, fieldId, model):
        """
        put two way coupling model for fieldId dispersed field
        """
        self.isInList(str(fieldId),self.getDispersedFieldList())
        self.isInList(model, TurbulenceModelsDescription.dispersedCouplingModels)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node == None:
            self.XMLturbulence.xmlInitChildNode('field', field_id = fieldId,
                                                               model = self.defaultValues()['model'],
                                                               two_way_coupling = model)
        node['two_way_coupling'] = model


    @Variables.noUndo
    def getTwoWayCouplingModel(self, fieldId):
        """
        get two way coupling for fieldId dispersed field
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node == None:
            self.setTwoWayCouplingModel(fieldId, self.defaultValues()['two_way_coupling'])
        model = node['two_way_coupling']

        return model


    @Variables.undoLocal
    def setContinuousCouplingModel(self, model):
        """
        put two way coupling model for continuous fields
        """
        self.isInList(model, TurbulenceModelsDescription.continuousCouplingModels)

        node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        if node == None:
            node = self.XMLturbulence.xmlInitChildNode('continuous_field_coupling')
        node['model'] = model


    @Variables.noUndo
    def getContinuousCouplingModel(self):
        """
        get two way coupling for continuous fields
        """
        node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        if node == None:
            self.setContinuousCouplingModel(self.defaultValues()['two_way_coupling'])
            node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        model = node['model']

        return model


    @Variables.undoLocal
    def setMixingLength(self, fieldId, value):
        """
        put value for mixing length
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        self.isFloat(value)

        fieldNode = self.XMLturbulence.xmlGetNode('field', field_id = str(fieldId))
        fieldNode.xmlSetData('length_scale', value)


    @Variables.noUndo
    def getMixingLength(self, fieldId):
        """
        get value for mixing length
        """
        self.isInList(str(fieldId),self.getFieldIdList())

        fieldNode = self.XMLturbulence.xmlGetNode('field', field_id = str(fieldId))
        lengthNode = fieldNode.xmlGetNode('length_scale')

        if lengthNode == None:
            value = self.defaultValues()['length_scale']
            self.setMixingLength(fieldId, value)
            lengthNode = fieldNode.xmlGetNode('length_scale')

        value = fieldNode.xmlGetDouble('length_scale')

        return value


    def isSecondOrderTurbulenceModel(self, fieldId):
        """
        return 1 if turbulent model of field is k-eps or Rij
        """
        self.isInList(str(fieldId),self.getFieldIdList())
        flag = 0
        if (self.getTurbulenceModel(fieldId) == "k-epsilon" \
         or self.getTurbulenceModel(fieldId) == "k-epsilon_linear_production" \
         or self.getTurbulenceModel(fieldId) == "rij-epsilon_ssg" \
         or self.getTurbulenceModel(fieldId) == "rij-epsilon_ebrsm" \
         or self.getTurbulenceModel(fieldId) == "les_smagorinsky" \
         or self.getTurbulenceModel(fieldId) == "les_wale" ):
            flag = 1
        return flag

    def useAdvancedThermalFluxes(self, fieldId):
        flag = False

        if self.getCriterion(fieldId) == 'continuous' and \
           'rij-epsilon' in self.getTurbulenceModel(fieldId):

            flag = True

        return flag

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------

class TurbulenceTestCase(ModelTest):
    """
    """
    def checkTurbulenceInstantiation(self):
        """Check whether the TurbulenceModel class could be instantiated"""
        model = None
        model = TurbulenceModel(self.case)
        assert model != None, 'Could not instantiate TurbulenceModel'


    def checkGetandSetTurbulenceModel(self):
        """Check whether the TurbulenceModel class could set and get TurbulenceModel"""
        MainFieldsModel(self.case).addField()
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('1','mixing_length')
        doc = '''<turbulence>
                         <variables/>
                         <properties/>
                         <field field_id="1" model="mixing_length" two_way_coupling="none"/>
                 </turbulence>'''
        assert mdl.XMLturbulence == self.xmlNodeFromString(doc),\
            'Could not set TurbulenceModel'
        assert mdl.getTurbulenceModel('1') == 'mixing_length',\
            'Could not get TurbulenceModel'


    def checkGetandSetTwoWayCouplingModel(self):
        """Check whether the TurbulenceModel class could set and get TwoWayCouplingModel"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField("2", "field2", 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = TurbulenceModel(self.case)

        mdl.setTurbulenceModel('2','tchen')
        mdl.setTwoWayCouplingModel('2','small_inclusions')
        doc = '''<turbulence>
                         <variables/>
                         <properties>
                                 <property choice="" field_id="2" label="TurbKineEner_q22" name="TurbKineEner_q2">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                                 <property choice="" field_id="2" label="Covariance_q122" name="Covariance_q12">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                                 <property choice="" field_id="2" label="turb_viscosity2" name="turb_viscosity">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                         </properties>
                         <field field_id="2" model="tchen" two_way_coupling="small_inclusions"/>
                 </turbulence>'''
        assert mdl.XMLturbulence == self.xmlNodeFromString(doc),\
            'Could not set TwoWayCouplingModel'
        assert mdl.getTwoWayCouplingModel('2') == 'small_inclusions',\
            'Could not get TwoWayCouplingModel'


    def checkGetandSetContinuousCouplingModel(self):
        """Check whether the TurbulenceModel class could set and get ContinuousCouplingModel"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField("2", "field2", 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = TurbulenceModel(self.case)
        mdl.setContinuousCouplingModel('separate_phase')
        doc = '''<turbulence>
                         <variables/>
                         <properties/>
                         <continuous_field_coupling model="separate_phase"/>
                 </turbulence>'''
        assert mdl.XMLturbulence == self.xmlNodeFromString(doc),\
            'Could not set ContinuousCouplingModel'
        assert mdl.getContinuousCouplingModel() == 'separate_phase',\
            'Could not get ContinuousCouplingModel'


    def checkGetandSetMixingLength(self):
        """Check whether the TurbulenceModel class could set and get MixingLenght"""
        MainFieldsModel(self.case).addField()
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('1','mixing_length')
        mdl.setMixingLength('1',8.4)
        doc = '''<turbulence>
                         <variables/>
                         <properties/>
                         <field field_id="1" model="mixing_length" two_way_coupling="none">
                                 <length_scale>
                                         8.4
                                 </length_scale>
                         </field>
                 </turbulence>'''
        assert mdl.XMLturbulence == self.xmlNodeFromString(doc),\
            'Could not set MixingLenght'
        assert mdl.getMixingLength('1') == 8.4,\
            'Could not get MixingLenght'


    def checkisSecondOrderTurbulenceModel(self):
        """Check whether the TurbulenceModel class could get isSecondOrderTurbulenceModel"""
        MainFieldsModel(self.case).addField()
        mdl = TurbulenceModel(self.case)
        mdl.setTurbulenceModel('1','k-epsilon')
        assert mdl.isSecondOrderTurbulenceModel('1') == 1,\
            'Could not get isSecondOrderTurbulenceModel'


    def checkGetandSetDomainLength(self):
        """Check whether the GeneralitiesModel class could set and get Domain Lenght"""
        mdl = TurbulenceModel(self.case)
        mdl.setDomainLength(1)
        doc = '''<turbulence>
                         <variables/>
                         <properties/>
                         <eddy_length_scale>
                                 1
                         </eddy_length_scale>
                         </field>
                 </turbulence>'''
        assert mdl.XMLturbulence == self.xmlNodeFromString(doc),\
            'Could not set Domain Length'
        assert mdl.getDomainLength() == (1),\
            'Could not get Domain Length'


def suite():
    testSuite = unittest.makeSuite(TurbulenceTestCase, "check")
    return testSuite


def runTest():
    print("TurbulenceTestCase")

    runner = unittest.TextTestRunner()
    runner.run(suite())
