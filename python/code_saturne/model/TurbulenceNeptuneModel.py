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
from code_saturne.model.MainFieldsModel import *       # TODO change
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

    reverseCouplingModels = ['k-epsilon', 'k-epsilon_linear_production',
                             'rij-epsilon_ssg', 'rij-epsilon_ebrsm']

    dispersedTurbulenceModels = ['none', 'q2-q12-tchen', 'q2-q12', 'r2-q12', 'r2-r12-tchen']

    bubblyFlowsTurbulenceModels = ["none", "q2-q12-tchen", "r2-r12-tchen"]
    dropletFlowsTurbulenceModels = ["none", "q2-q12", "r2-q12"]

    continuousCouplingModels = ['none', 'separate_phase', 'separate_phase_cond']
    dispersedCouplingModels = ['none', 'small_inclusions', 'large_inclusions']

    ThermalTurbFluxModels = ['sgdh', 'ggdh']

    turbulenceVariables = {}
    turbulenceVariables['none'] = []
    turbulenceVariables['mixing_length'] = []
    turbulenceVariables['k-epsilon'] = ['k', 'epsilon']
    turbulenceVariables['k-epsilon_linear_production'] = ['k', 'epsilon']
    turbulenceVariables['rij-epsilon_ssg'] = ['reynolds_stress', 'epsilon']
    turbulenceVariables['rij-epsilon_ebrsm'] = ['reynolds_stress', 'epsilon', 'alpha']
    turbulenceVariables['les_smagorinsky'] = []
    turbulenceVariables['les_wale'] = []
    turbulenceVariables['q2-q12-tchen'] = []
    turbulenceVariables['q2-q12'] = ['TurbKineEner_qp', 'covariance_qfp']
    turbulenceVariables['r2-q12'] = ['reynolds_stress','covariance_qfp']
    turbulenceVariables['r2-r12-tchen'] = ['reynolds_stress',
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
    turbulenceProperties['q2-q12-tchen'] = ["TurbKineEner_qp", "covariance_qfp", "turb_viscosity"]
    turbulenceProperties['q2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-q12'] = ["turb_viscosity"]
    turbulenceProperties['r2-r12-tchen'] = ["turb_viscosity"]


#-------------------------------------------------------------------------------
# Main turbulence model class
#-------------------------------------------------------------------------------

class TurbulenceModel(Variables, Model):
    """
    This class manages the turbulence objects in the XML file
    """

    def __init__(self, case):
        """
        Constructor.
        """
        #
        # XML file parameters
        self.mainFieldsModel = MainFieldsModel(case)
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
        self.mainFieldsModel.isFieldIdValid(fieldId)
        field = self.mainFieldsModel.getFieldFromId(fieldId)

        field_name = field.label
        criterion = field.flow_type
        if criterion == "continuous":
           self.isInList(model, TurbulenceModelsDescription.continuousTurbulenceModels)
        else:
           self.isInList(model, TurbulenceModelsDescription.dispersedTurbulenceModels)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        oldmodel = ""
        if node != None:
            oldmodel = self.defaultValues()['model']

        if node is None:
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
                   self.removeVariableProperty("variable", self.XMLNodeVariable, fieldId, var)

               for var in TurbulenceModelsDescription.turbulenceProperties[oldmodel]:
                   self.removeVariableProperty("property", self.XMLNodeproperty, fieldId, var)

           # add new variables and properties from XML
           for var in TurbulenceModelsDescription.turbulenceVariables[model]:
             if var == "reynolds_stress":
                 self.setNewVariableProperty("variable", "",
                                                             self.XMLNodeVariable,
                                                             fieldId,
                                                             var,
                                                             var+"_"+field_name,
                                                             dim=6)
             elif var == 'covariance_qfp':
                 if field.carrier_id == 'all':
                     _carr_lst = self.mainFieldsModel.getContinuousFieldList()
                 else:
                     _carr_lst = [self.mainFieldsModel.getFieldFromId(field.carrier_id)]

                 for _f in _carr_lst:
                     _vname = "_".join([var, _f.f_id])
                     _vlabel = "_".join([var, _f.label, field_name])
                     self.setNewVariableProperty("variable", "",
                                                 self.XMLNodeVariable,
                                                 fieldId,
                                                 _vname,
                                                 _vlabel)
             else:
                 if var == 'covariance_qfp':
                     var += "_" + self.mainFieldsModel.getFieldLabelsList()[0]
                 self.setNewVariableProperty("variable", "",
                                                             self.XMLNodeVariable,
                                                             fieldId,
                                                             var,
                                                             var+"_"+field_name)

           for var in TurbulenceModelsDescription.turbulenceProperties[model]:
               self.setNewVariableProperty("property", "", self.XMLNodeproperty, fieldId, var, var+"_"+field_name)

           if oldmodel == "mixing_length":
               node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
               noder = node.xmlGetNode('length_scale')
               noder.xmlRemoveNode()

           # update other field if continuous and set to none
           if model == "none" and criterion == "continuous":
               for fld in self.mainFieldsModel.list_of_fields:
                   if fld.flow_type == "dispersed":
                       if fld.carrier_id == str(fieldId):
                           self.setTurbulenceModel(fld.f_id, TurbulenceModelsDescription.dispersedTurbulenceModels[0])
                           if criterion == "continuous":
                               self.setTwoWayCouplingModel(fld.f_id, TurbulenceModelsDescription.continuousCouplingModels[0])
                           else:
                               self.setTwoWayCouplingModel(fld.f_id, TurbulenceModelsDescription.continuousCouplingModels[0])


    @Variables.noUndo
    def getTurbulenceModel(self, fieldId):
        """
        get turbulence model for fieldId
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        field = self.mainFieldsModel.getFieldFromId(fieldId)

        criterion = field.flow_type

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node is None:
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

        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isInList(model,TurbulenceModelsDescription.ThermalTurbFluxModels)
        field = self.mainFieldsModel.getFieldFromId(fieldId)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)

        critrerion =  field.flow_type
        if node is None:
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

        self.mainFieldsModel.isFieldIdValid(fieldId)
        field = self.mainFieldsModel.getFieldFromId(fieldId)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node is None:
            if field.enthalpy_model != 'off':
                self.setThermalTurbulentFlux(fieldId,
                                             self.defaultValues()['turb_flux'])
            else:
                self.setThermalTurbulentFlux(fieldId, 'none')
        else:
            if field.enthalpy_model == 'off':
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
        field = self.mainFieldsModel.getFieldFromId(fieldId)
        assert(field.flow_type == "dispersed")
        self.isInList(model, TurbulenceModelsDescription.dispersedCouplingModels)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node is None:
            node = self.XMLturbulence.xmlInitChildNode('field', field_id=fieldId,
                                                       model=self.defaultValues()['model'],
                                                       two_way_coupling=model)
        node['two_way_coupling'] = model


    @Variables.noUndo
    def getTwoWayCouplingModel(self, fieldId):
        """
        get two way coupling for fieldId dispersed field
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        node = self.XMLturbulence.xmlGetNode('field', field_id = fieldId)
        if node is None:
            self.setTwoWayCouplingModel(fieldId, self.defaultValues()['two_way_coupling'])
            node = self.XMLturbulence.xmlGetNode('field', field_id=fieldId)
        model = node['two_way_coupling']

        return model


    @Variables.undoLocal
    def setContinuousCouplingModel(self, model):
        """
        put two way coupling model for continuous fields
        """
        self.isInList(model, TurbulenceModelsDescription.continuousCouplingModels)

        node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        if node is None:
            node = self.XMLturbulence.xmlInitChildNode('continuous_field_coupling')
        node['model'] = model


    @Variables.noUndo
    def getContinuousCouplingModel(self):
        """
        get two way coupling for continuous fields
        """
        node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        if node is None:
            self.setContinuousCouplingModel(self.defaultValues()['two_way_coupling'])
            node = self.XMLturbulence.xmlGetNode('continuous_field_coupling')
        model = node['model']

        return model


    @Variables.undoLocal
    def setMixingLength(self, fieldId, value):
        """
        put value for mixing length
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)
        self.isFloat(value)

        fieldNode = self.XMLturbulence.xmlGetNode('field', field_id = str(fieldId))
        fieldNode.xmlSetData('length_scale', value)


    @Variables.noUndo
    def getMixingLength(self, fieldId):
        """
        get value for mixing length
        """
        self.mainFieldsModel.isFieldIdValid(fieldId)

        fieldNode = self.XMLturbulence.xmlGetNode('field', field_id = str(fieldId))
        lengthNode = fieldNode.xmlGetNode('length_scale')

        if lengthNode is None:
            value = self.defaultValues()['length_scale']
            self.setMixingLength(fieldId, value)
            lengthNode = fieldNode.xmlGetNode('length_scale')

        value = fieldNode.xmlGetDouble('length_scale')

        return value


    def modelLevelIsAboveTwoEquations(self, fieldId):
        """
        return 1 if turbulent model of field is k-eps or Rij
        """
        if fieldId == "all":
            flag = 1
            for continuous_field in self.mainFieldsModel.getContinuousFieldList():
                flag = flag and self.modelLevelIsAboveTwoEquations(continuous_field.f_id)
        else:
            flag = 0
            self.mainFieldsModel.isFieldIdValid(fieldId)
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

        field = self.mainFieldsModel.getFieldFromId(fieldId)
        if field.flow_type == 'continuous' and \
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

        mdl.setTurbulenceModel('2','q2-q12-tchen')
        mdl.setTwoWayCouplingModel('2','small_inclusions')
        doc = '''<turbulence>
                         <variables/>
                         <properties>
                                 <property choice="" field_id="2" label="TurbKineEner_q22" name="TurbKineEner_qp">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                                 <property choice="" field_id="2" label="Covariance_q122" name="covariance_qfp">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                                 <property choice="" field_id="2" label="turb_viscosity2" name="turb_viscosity">
                                         <listing_printing status="on"/>
                                         <postprocessing_recording status="on"/>
                                 </property>
                         </properties>
                         <field field_id="2" model="q2-q12-tchen" two_way_coupling="small_inclusions"/>
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
