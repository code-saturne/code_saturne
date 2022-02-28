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

# -------------------------------------------------------------------------------

import unittest
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
from code_saturne.model.NotebookModel import NotebookModel


class InterfacialForcesModel(MainFieldsModel, Variables, Model):
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
        self.turb_m   = TurbulenceModel(case)
        self.notebook = NotebookModel(case)

        self.case = case
        self.XMLClosure = self.case.xmlGetNode('closure_modeling')
        self.XMLInterForce = self.XMLClosure.xmlInitNode('interfacial_forces')
        self.XMLNodethermo = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeproperty = self.XMLNodethermo.xmlInitNode('properties')

        self.__availableturbulentedispersionModelsList = ["none", "LLB_model", "GTD_model"]
        self.__availablewallforcesModelList = ["none", "antal", "tomiyama"]
        self.__availableContinuousDragModelList = ["Large_Interface_Model", "G_Large_Interface_Model",
                                                   "Large_Bubble_Model"]
        self.__availableGasDispersedDragModelList = ["ishii", "Wen_Yu"]
        self.__availableSolidLiquidDispersedDragModelList = ["Gobin"]

        self.__availableAddedMassModelsLists = ["none", "standard", "zuber"]
        self.__availableLiftModelsLists = ["none", "coef_cst", "Tomiyama_SMD", "Zeng_Baalbaki"]

        # Init freeCouples for forces : criterion checking !
        self.__allCouples = []
        for i, fieldaId in enumerate(self.getContinuousFieldList()):
            for fieldbId in self.getContinuousFieldList()[i + 1:]:
                self.__allCouples.append([self.getLabel(fieldaId), self.getLabel(fieldbId), "continuous"])
        for fieldbId in self.getDispersedFieldList():
            fieldaId = self.getCarrierField(fieldbId)
            self.__allCouples.append([self.getLabel(fieldaId), self.getLabel(fieldbId), "dispersed"])

    def getAllCouples(self):
        return self.__allCouples

    # TODO : to remove once the "auto" type of interaction is enabled ?
    def getGLIMfields(self):
        fields = []
        for i, fieldaId in enumerate(self.getContinuousFieldList()):
            for fieldbId in self.getContinuousFieldList()[i + 1:]:
                if self.getContinuousCouplingModel(fieldaId, fieldbId) == "G_Large_Interface_Model":
                    fields += [fieldaId, fieldbId]
        return fields

    def getAvailableContinuousDragModelList(self):
        """
        Get available turbulente dispersion list model
        """
        return self.__availableContinuousDragModelList

    def getAvailableTurbulenteDispersionModelList(self, fieldaId, fieldbId):
        """
        Get available turbulente dispersion list model
        """
        # Only if : fieldaId is continuous and k-eps or Rij-eps
        #           fieldbId is dispersed without turbulence
        lst = []
        if (fieldaId in self.getContinuousFieldList() and fieldbId in self.getDispersedFieldList()) :
            if (self.turb_m.isSecondOrderTurbulenceModel(fieldaId) and self.turb_m.getTurbulenceModel(fieldbId) == "none") :
                lst = self.__availableturbulentedispersionModelsList
            else :
                lst = ["none"]
        else :
            lst = ["none"]

        return lst


    def getAvailableWallForcesModelList(self, fieldaId, fieldbId):
        """
        Get available wall forces list model
        """
        # only for bubbles in water
        lst = []
        if ((fieldbId in self.getDispersedFieldList()) and self.getFieldNature(fieldbId) == "gas") :
            lst = self.__availablewallforcesModelList
        else :
            lst = ["none"]
        return lst


    def getAvailableDragModels(self, fieldaId, fieldbId) :
        """
        return Model list according to fields criterion
        """
        # field A is continuous and field B is dispersed
        self.isInList(fieldaId,self.getContinuousFieldList())
        self.isInList(fieldbId,self.getDispersedFieldList())

        if (self.getFieldNature(fieldaId) == "liquid") and (self.getFieldNature(fieldbId) == "solid") :
            return self.__availableSolidLiquidDispersedDragModelList
        else :
            predefined_flow = self.getPredefinedFlow()
            if predefined_flow == "boiling_flow":
                return ["ishii"]
            elif predefined_flow == "droplet_flow":
                return ["Wen_Yu"]
            else:
                return self.__availableGasDispersedDragModelList


    def getAvailableAddedMassModels(self) :
        """
        return Model list according to fields criterion
        """
        return self.__availableAddedMassModelsLists


    def getAvailableLiftModels(self) :
        """
        return Model list according to fields criterion
        """
        return self.__availableLiftModelsLists


    def defaultValues(self):
        default = {}
        predefined_flow = self.getPredefinedFlow()

        default['gasdisperseddragmodel'] = 'ishii'
        default['liquidsoliddisperseddragmodel'] = 'Gobin'
        default['addedmassmodel'] = 'zuber'
        default['liftmodel'] = 'Tomiyama_SMD'
        default['turbulent_dispersion_model'] = "none"
        default['wallforcemodel'] = 'tomiyama'
        default['nowallforcemodel'] = 'none'

        if predefined_flow == "boiling_flow":
            GTD_condition_1 = (self.turb_m.getTurbulenceModel("1") in
                               ["k-epsilon",
                                "k-epsilon_linear_production",
                                "rij-epsilon_ssg'",
                                "rij-epsilon_ebrsm"])

            GTD_condition_2 = (self.turb_m.getTurbulenceModel("2") == "none")
            if GTD_condition_1 and GTD_condition_2:
                default['turbulent_dispersion_model'] = "GTD_model"
        elif predefined_flow == "droplet_flow":
            default['gasdisperseddragmodel'] = "Wen_Yu"
            default['liftmodel'] = "Zeng_Baalbaki"
            default['addedmassmodel'] = "none"
            default['wallforcemodel'] = "none"
        elif predefined_flow == "particles_flow":
            default['gasdisperseddragmodel'] = 'ishii'
            default['liftmodel'] = "Zeng_Baalbaki"
            default['wallforcemodel'] = 'none'
        else:
            pass

        return default


    def defaultValuesContinuous(self):
        default = {}
        predefined_flow = self.getPredefinedFlow()

        default['continuousdragmodel']           = 'Large_Interface_Model'
        default['BubblesForLIM']                 = 'off'
        default['InterfaceSharpening'] = 'none'
        default["SurfaceTension"] = "none"

        if predefined_flow == "free_surface":
            pass
        elif predefined_flow == "multiregime":
            default['continuousdragmodel'] = "G_Large_Interface_Model"
            default['BubblesForLIM'] = "on"
        else:
            pass

        return default

    def setDefaultParameters(self, field_id_a, field_id_b):
        predefined_flow = self.getPredefinedFlow()

        # Dispersed models
        if predefined_flow in ["boiling_flow", "droplet_flow", "particles_flow"]:
            default = self.defaultValues()
            if (self.getFieldNature(field_id_a) == "liquid") and (self.getFieldNature(field_id_b) == "solid"):
                self.setDragModel(field_id_a, field_id_b, default["liquidsoliddisperseddragmodel"])
            else:
                self.setDragModel(field_id_a, field_id_b, default["gasdisperseddragmodel"])
            self.setLiftModel(field_id_a, field_id_b, default["liftmodel"])
            self.setAddMassModel(field_id_a, field_id_b, default["addedmassmodel"])
            self.setTurbDispModel(field_id_a, field_id_b, default["turbulent_dispersion_model"])
            self.setWallForceModel(field_id_a, field_id_b, default["wallforcemodel"])

        # Continuous models
        elif predefined_flow in ["free_surface", "multiregime"]:
            default = self.defaultValuesContinuous()
            self.setContinuousCouplingModel(field_id_a, field_id_b, default["continuousdragmodel"])
            self.setBubblesForLIMStatus(field_id_a, field_id_b, default["BubblesForLIM"])
            self.setInterfaceSharpeningModel(field_id_a, field_id_b, default["InterfaceSharpening"])
            self.setSurfaceTensionModel(field_id_a, field_id_b, default["SurfaceTension"])

        # User defined flow : do nothing for now
        else:
            pass


    @Variables.noUndo
    def getForceList(self) :
        """
        return list of field couple for forces
        """
        lst = []
        for node in self.XMLInterForce.xmlGetNodeList('force') :
            couple=[node['field_id_a'], node['field_id_b']]
            lst.append(couple)
        return lst

    def __updatedrag(self, fieldaId, fieldbId) :
        """
        update drag model after fieldIdChange
        """
        model = ""
        model = self.defaultValues()['gasdisperseddragmodel']
        if (self.getFieldNature(fieldaId) == "liquid") and (self.getFieldNature(fieldbId == "solid")):
            model = self.defaultValues()['liquidsoliddisperseddragmodel']
        self.setDragModel(fieldaId, fieldbId, model)
        if model != "none" :
            Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldbId, "drag_coefficient", "drag_coef"+str(fieldbId))
        else :
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldbId, "drag_coefficient")


    def __updateWallForce(self, fieldaId, fieldbId) :
        """
        update wall force after fieldIdChange
        """
        model = self.getWallForceModel(fieldaId, fieldbId)
        if model not in self.getAvailableWallForcesModelList(fieldaId, fieldbId) :
            model = self.getAvailableWallForcesModelList(fieldaId, fieldbId)[0]
            self.setWallForceModel(fieldaId, fieldbId, model)


    def __updateTurbulentDispersion(self, fieldaId, fieldbId) :
        """
        update turbulent dispersion after fieldIdChange
        """
        model = self.getTurbDispModel(fieldaId, fieldbId)
        if model not in self.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId) :
            model = self.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId)[0]
            self.setTurbDispModel(fieldaId, fieldbId, model)


    @Variables.undoGlobal
    def setDragModel(self, fieldaId, fieldbId, model) :
        """
        set drag model for a couple fieldId
        """
        self.isInList(model, self.getAvailableDragModels(fieldaId, fieldbId))
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        ChildNode = node.xmlInitChildNode('drag_model')
        ChildNode['model'] = model
        if model != "none":
            Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldbId,
                                                        "drag_coefficient", "drag_coef" + str(fieldbId))
        else:
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldbId, "drag_coefficient")


    @Variables.noUndo
    def getDragModel(self, fieldaId, fieldbId) :
        """
        get drag model for a couple fieldId
        """
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        childNode = node.xmlGetNode('drag_model')
        if childNode is None :
            model = self.defaultValues()['gasdisperseddragmodel']
            if (self.getFieldNature(fieldaId) == "liquid") and (self.getFieldNature(fieldbId) == "solid"):
                model = self.defaultValues()['liquidsoliddisperseddragmodel']
            self.setDragModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('drag_model')['model']
        return model


    @Variables.undoLocal
    def setAddMassModel(self, fieldaId, fieldbId, model) :
        """
        set added mass model for a couple fieldId
        """
        self.isInList(model, self.getAvailableAddedMassModels())
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        ChildNode = node.xmlInitChildNode('added_mass')
        ChildNode['model'] = model


    @Variables.noUndo
    def getAddMassModel(self, fieldaId, fieldbId) :
        """
        get added mass model for a couple fieldId
        """
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        childNode = node.xmlGetNode('added_mass')
        if childNode is None :
            model = self.defaultValues()['addedmassmodel']
            self.setAddMassModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('added_mass')['model']
        return model


    @Variables.undoLocal
    def setLiftModel(self, fieldaId, fieldbId, model) :
        """
        set lift model for a couple fieldId
        """
        self.isInList(model, self.getAvailableLiftModels())
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        ChildNode = node.xmlInitChildNode('lift_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getLiftModel(self, fieldaId, fieldbId) :
        """
        get lift model for a couple fieldId
        """
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        childNode = node.xmlGetNode('lift_model')
        if childNode is None :
            model = self.defaultValues()['liftmodel']
            self.setLiftModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('lift_model')['model']
        return model


    @Variables.undoLocal
    def setWallForceModel(self, fieldaId, fieldbId, model) :
        """
        set wall force model for a couple fieldId
        """
        self.isInList(model, self.getAvailableWallForcesModelList(fieldaId, fieldbId))
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        ChildNode = node.xmlInitChildNode('wall_force_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getWallForceModel(self, fieldaId, fieldbId) :
        """
        get wall force model for a couple fieldId
        """
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        childNode = node.xmlGetNode('wall_force_model')
        if childNode is None :
            if len(self.getAvailableWallForcesModelList(fieldaId, fieldbId)) > 1:
                model = self.defaultValues()['wallforcemodel']
            else:
                model = self.defaultValues()['nowallforcemodel']
            self.setWallForceModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('wall_force_model')['model']
        return model


    @Variables.undoLocal
    def setTurbDispModel(self, fieldaId, fieldbId, model) :
        """
        set turbulent dispersion model for a couple fieldId
        """
        self.isInList(model, self.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId))
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        ChildNode = node.xmlInitChildNode('turbulent_dispersion_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getTurbDispModel(self, fieldaId, fieldbId) :
        """
        get turbulent dispersion model for a couple fieldId
        """
        node = self.XMLInterForce.xmlInitChildNode('force', field_id_a=fieldaId, field_id_b=fieldbId)
        childNode = node.xmlGetNode('turbulent_dispersion_model')
        if childNode is None :
            model = self.defaultValues()['turbulent_dispersion_model']
            self.setTurbDispModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('turbulent_dispersion_model')['model']
        return model


    @Variables.undoGlobal
    def deleteForce(self, fieldaId, fieldbId) :
        """
        suppress force
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())
        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        node.xmlRemoveNode()

    @Variables.undoGlobal
    def setContinuousCouplingModel(self, field_id_a, field_id_b, model):
        """
        set turbulent dispersion model for a couple fieldId
        """
        self.isInList(model, self.getAvailableContinuousDragModelList())
        ChildNode = self.XMLInterForce.xmlInitChildNode('continuous_field_momentum_transfer',
                                                        field_id_a=field_id_a,
                                                        field_id_b=field_id_b)
        ChildNode['model'] = model
        if model != "separate_phases":
            ChildNode2 = ChildNode.xmlGetChildNode('gradP_correction')
            if ChildNode2 != None:
                ChildNode2.xmlRemoveNode()
            ChildNode2 = ChildNode.xmlGetChildNode('gradP_correction_model')
            if ChildNode2 != None:
                ChildNode2.xmlRemoveNode()

        if model not in ["Large_Interface_Model", "G_Large_Interface_Model"]:
            ChildNode2 = ChildNode.xmlGetChildNode('BubblesForLIM')
            if ChildNode2 != None:
                ChildNode2.xmlRemoveNode()

        if model == "none":
            ChildNode2 = ChildNode.xmlGetChildNode('InterfaceSharpening')
            if ChildNode2 != None:
                ChildNode2.xmlRemoveNode()

        if model == "G_Large_Interface_Model":
            self.setBubblesForLIMStatus(field_id_a, field_id_b, "on")

    @Variables.noUndo
    def getContinuousCouplingModel(self, field_id_a, field_id_b):
        """
        get turbulent dispersion model for a couple fieldId
        """
        ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer', field_id_a=field_id_a,
                                                       field_id_b=field_id_b)

        # separate_phases model has been removed from the GUI!
        # Hence if the test leads to separate phases we set the model to none...
        if ChildNode is None:
            model = self.defaultValuesContinuous()['continuousdragmodel']
            self.setContinuousCouplingModel(field_id_a, field_id_b, model)
            ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer', field_id_a=field_id_a,
                                                           field_id_b=field_id_b)
        elif ChildNode['model'] == 'separate_phases':
            self.setContinuousCouplingModel(field_id_a, field_id_b, 'none')
            ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer', field_id_a=field_id_a,
                                                           field_id_b=field_id_b)

        model = ChildNode['model']
        if model is None:
            model = "none"
        return model

    @Variables.undoLocal
    def setBubblesForLIMStatus(self, field_id_a, field_id_b, status):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        old_status = node.xmlGetString('BubblesForLIM')
        node.xmlSetData('BubblesForLIM', status)

        mfm = MainFieldsModel(self.case)
        if mfm.getPredefinedFlow() == 'free_surface':
            fieldId = self.getContinuousFieldList()[1]

            if status == "off" and old_status:
                node = self.XMLClosure.xmlGetChildNode('interfacial_area_diameter')
                if node:
                    node.xmlRemoveNode()

                Variables(self.case).removeVariableProperty("property",
                                                            self.XMLNodeproperty,
                                                            fieldId,
                                                            "diameter")
                Variables(self.case).removeVariableProperty("property",
                                                            self.XMLNodeproperty,
                                                            fieldId,
                                                            "drift_component")

            elif status == "on":
                field_name = MainFieldsModel(self.case).getFieldLabelsList()[int(fieldId)-1]
                Variables(self.case).setNewVariableProperty('property', '',
                                                            self.XMLNodeproperty,
                                                            fieldId,
                                                            'diameter',
                                                            'diam_'+field_name)
                Variables(self.case).setNewVariableProperty('property', '',
                                                            self.XMLNodeproperty,
                                                            fieldId,
                                                            'drift_component',
                                                            'drift_component_'+field_name)

    @Variables.noUndo
    def getBubblesForLIMStatus(self, field_id_a, field_id_b):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        if node is None:
            node = self.XMLInterForce.xmlInitChildNode('continuous_field_momentum_transfer',
                                                       field_id_a=field_id_a,
                                                       field_id_b=field_id_b)
        ChildNode = node.xmlGetChildNode('BubblesForLIM')
        if ChildNode is None:
            status = self.defaultValuesContinuous()['BubblesForLIM']
            self.setBubblesForLIMStatus(field_id_a, field_id_b, status)
        status = node.xmlGetString('BubblesForLIM')
        return status

    @Variables.undoLocal
    def setInterfaceSharpeningModel(self, field_id_a, field_id_b, status):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        node.xmlSetData('InterfaceSharpening', status)

    @Variables.noUndo
    def getInterfaceSharpeningModel(self, field_id_a, field_id_b):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        ChildNode = node.xmlGetChildNode('InterfaceSharpening')
        if ChildNode is None:
            status = self.defaultValuesContinuous()['InterfaceSharpening']
            self.setInterfaceSharpeningModel(field_id_a, field_id_b, status)
        status = node.xmlGetString('InterfaceSharpening')
        return status

    @Variables.undoLocal
    def setSurfaceTensionModel(self, field_id_a, field_id_b, model):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        st_node = node.xmlInitChildNode('SurfaceTension')
        st_node.xmlSetData('model', model)

    @Variables.noUndo
    def getSurfaceTensionModel(self, field_id_a, field_id_b):
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer',
                                                  field_id_a=field_id_a,
                                                  field_id_b=field_id_b)
        st_node = node.xmlGetChildNode('SurfaceTension')
        if st_node is None:
            st_node = node.xmlInitChildNode('SurfaceTension')
            model = self.defaultValuesContinuous()['SurfaceTension']
            self.setSurfaceTensionModel(field_id_a, field_id_b, model)
        model = st_node.xmlGetString('model')
        return model

#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class InterfacialForcesTestCase(ModelTest):
    """
    """
    def checkInterfacialForcesInstantiation(self):
        """Check whether the InterfaciallForcesModel class could be instantiated"""
        model = None
        model = InterfacialForcesModel(self.case)
        assert model != None, 'Could not instantiate InterfacialForcesModel'


    def checkGetAvailableContinuousDragModelList(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableContinuousDragModelList"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableContinuousDragModelList() == ["none", "Large_Interface_Model", "G_Large_Interface_Model",
                                                             "Large_Bubble_Model"], \
            'Could not get AvailableContinuousDragModelList'


    def checkGetAvailableTurbulenteDispersionModelList(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableTurbulenteDispersionModelList"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        TurbulenceModel(self.case).setTurbulenceModel('1','k-epsilon')
        TurbulenceModel(self.case).setTurbulenceModel('2','none')
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableTurbulenteDispersionModelList('1','2') == ["none","LLB_model","GTD_model"],\
            'Could not get AvailableTurbulenteDispersionModelList'


    def checkGetAvailableWallForcesModelList(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableWallForcesModelList"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        ThermodynamicsModel(self.case).setMaterials('2','Water')
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableWallForcesModelList('1','2') == ["none", "antal", "tomiyama"],\
            'Could not get AvailableWallForcesModelList'


    def checkGetAvailableDragModels(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableDragModels"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableDragModels('1','2') == ["none", "ishii", "Wen_Yu"],\
            'Could not get AvailableDragModels'
        MainFieldsModel(self.case).setFieldNature('2','solid')
        assert mdl.getAvailableDragModels('1','2') == ["none", "Gobin"],\
            'Could not get AvailableDragModels'


    def checkGetAvailableAddedMassModels(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableAddedMassModels"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableAddedMassModels() == ["none", "standard", "zuber"],\
            'Could not get AvailableAddedMassModels'


    def checkGetAvailableLiftModels(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableLiftModels"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableLiftModels() == ["none", "coef_cst", "Tomiyama_SMD", "Zeng_Baalbaki"], \
            'Could not get AvailableLiftModels'


    def checkGetForceList(self):
        """Check whether the InterfacialEnthalpyModel class could get force list"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getForceList() == [['1', '2']],\
            'Could not get force list'

    def checkGetandSetDragModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get DragModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setDragModel('1','2','ishii')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set DragModel'
        assert mdl.getDragModel('1','2') == 'ishii',\
            'Could not get DragModel'


    def checkGetandSetAddMassModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get AddMassModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setAddMassModel('1','2','standard')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="standard"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set AddMassModel'
        assert mdl.getAddMassModel('1','2') == 'standard',\
            'Could not set AddMassModel'


    def checkGetandSetLiftModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get LiftModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setLiftModel('1','2','coef_cst')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="coef_cst"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set LiftModel'
        assert mdl.getLiftModel('1','2') == 'coef_cst',\
            'Could not set LiftModel'


    def checkGetandSetWallForceModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get WallForceModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        ThermodynamicsModel(self.case).setMaterials('2','Water')
        mdl = InterfacialForcesModel(self.case)
        mdl.setWallForceModel('1','2','antal')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="antal"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set WallForceModel'
        assert mdl.getWallForceModel('1','2') == 'antal',\
            'Could not set WallForceModel'


    def checkGetandSetTurbDispModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get TurbDispModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        TurbulenceModel(self.case).setTurbulenceModel('1','k-epsilon')
        TurbulenceModel(self.case).setTurbulenceModel('2','none')
        mdl = InterfacialForcesModel(self.case)
        mdl.setTurbDispModel('1','2','LLB_model')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="LLB_model"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set TurbDispModel'
        assert mdl.getTurbDispModel('1','2') == 'LLB_model',\
            'Could not set TurbDispModel'


    def checkdeleteForce(self):
        """Check whether the InterfacialEnthalpyModel class could deleteForce"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.deleteForce('1','2')
        doc = '''<interfacial_forces/>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not delete Force'


    def checkGetandSetContinuousCouplingModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get ContinuousCouplingModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('1', '2', 'Large_Interface_Model')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="Large_Interface_Model"/>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc), \
            'Could not set ContinuousCouplingModel'
        assert mdl.getContinuousCouplingModel('1', '2') == 'Large_Interface_Model', \
            'Could not get ContinuousCouplingModel'

    def checkGetandSetTurbulenceEffectStatus(self):
        """Check whether the InterfacialEnthalpyModel class could set and get TurbulenceEffectStatus"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('1', '2', 'Large_Interface_Model')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="Large_Interface_Model">
                                 <TurbulenceEffect>
                                         on
                                 </TurbulenceEffect>
                         </continuous_field_momentum_transfer>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set TurbulenceEffectStatus'
        assert mdl.getTurbulenceEffectStatus() == 'on',\
            'Could not set TurbulenceEffectStatus'


    def checkGetandSetFrictionStatus(self):
        """Check whether the InterfacialEnthalpyModel class could set and get FrictionStatus"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('1', '2', 'Large_Interface_Model')
        mdl.setFrictionStatus('on')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="Large_Interface_Model">
                                 <WaweRoughness>
                                         on
                                 </WaweRoughness>
                         </continuous_field_momentum_transfer>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set FrictionStatus'
        assert mdl.getFrictionStatus() == 'on',\
            'Could not set FrictionStatus'


def suite():
    testSuite = unittest.makeSuite(InterfacialForcesTestCase, "check")
    return testSuite


def runTest():
    print("InterfacialForcesTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite())
