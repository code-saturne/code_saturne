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
from code_saturne.Base.XMLmodel import *
from code_saturne.Pages.MainFieldsModel import MainFieldsModel
from code_saturne.Pages.TurbulenceNeptuneModel import TurbulenceModel
from code_saturne.Pages.ThermodynamicsModel import ThermodynamicsModel
import copy

class InterfacialForcesModel(TurbulenceModel,
                             ThermodynamicsModel):

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
        TurbulenceModel.__init__(self, case)
        ThermodynamicsModel.__init__(self, case)
        self.case = case
        self.XMLClosure      = self.case.xmlGetNode('closure_modeling')
        self.XMLInterForce   = self.XMLClosure.xmlInitNode('interfacial_forces')
        self.XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
        self.XMLNodeproperty = self.XMLNodethermo.xmlInitNode('properties')

        self.__availableturbulentedispersionModelsList    = ["none","LLB_model","GTD_model"]
        self.__availablewallforcesModelList               = ["none", "antal", "tomiyama"]
        self.__availableContinuousDragModelList           = ["none", "separate_phases" ,"Large_Interface_Model", "Large_Bubble_Model"]
        self.__availableGasDispersedDragModelList         = ["none", "ishii"]
        self.__availableSolidLiquidDispersedDragModelList = ["none", "inclusions", "Wen_Yu"]
        self.__availableAddedMassModelsLists              = ["none", "standard", "zuber"]
        self.__availableLiftModelsLists                   = ["none", "coef_cst", "Tomiyama_SMD"]

        # Init freeCouples for forces : criterion checking !
        self.__freeCouples = []

        for fieldaId in self.getContinuousFieldList() :
            for fieldbId in self.getDispersedFieldList() :
                self.__freeCouples.append((fieldaId, fieldbId))

        XMLNodesForces = self.XMLInterForce.xmlGetNodeList('force', 'field_id_a', 'field_id_b')
        for node in XMLNodesForces :
            # Get fields
            fieldaId = node['field_id_a']
            fieldbId = node['field_id_b']
            #
            # Update free couples
            try :
                self.__freeCouples.remove((fieldaId, fieldbId))
            except:
                pass


    def getFreeCouples (self) :
         """
         return list of free couples
         """
         return self.__freeCouples


    def getFieldIdaList (self, fieldaId) :
         """
         return list of free couples
         """
         list = []
         list.append(fieldaId)
         for fielda, fieldb in self.__freeCouples :
             if fielda not in list :
                 list.append(fielda)
         return list


    def getFieldIdbList (self, fieldaId) :
         """
         return list of free couples
         """
         list = []
         for fielda, fieldb in self.__freeCouples :
             if str(fielda) == str(fieldaId) :
                 list.append(fieldb)
         return list


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
        list = []
        if (fieldaId in self.getContinuousFieldList() and fieldbId in self.getDispersedFieldList()) :
            if (self.isSecondOrderTurbulenceModel(fieldaId) and self.getTurbulenceModel(fieldbId) == "none") :
                list = self.__availableturbulentedispersionModelsList
            else :
                list = ["none"]
        else :
            list = ["none"]
        return list


    def getAvailableWallForcesModelList(self, fieldaId, fieldbId):
        """
        Get available wall forces list model
        """
        # only for bubbles in water
        list = []
        if ((fieldbId in self.getDispersedFieldList()) and self.getFieldNature(fieldbId) == "gas") :
            list = self.__availablewallforcesModelList
        else :
            list = ["none"]
        return list


    def getAvailableDragModels(self, fieldaId, fieldbId) :
        """
        return Model list according to fields criterion
        """
        # field A is continuous and field B is dispersed
        self.isInList(fieldaId,self.getContinuousFieldList())
        self.isInList(fieldbId,self.getDispersedFieldList())

        if self.getFieldNature(fieldbId) == "gas" :
            return self.__availableGasDispersedDragModelList
        else :
            return self.__availableSolidLiquidDispersedDragModelList


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


    def defaultValues(self, fieldaId, fieldbId):
        default = ThermodynamicsModel.defaultValues(self)

        # Drag model
        default['gasdisperseddragmodel']         = 'ishii'
        default['liquidsoliddisperseddragmodel'] = 'inclusions'

        # Added mass model
        default['addedmassmodel']                = 'zuber'

        # lift
        default['liftmodel']                     = 'Tomiyama_SMD'

        # turbulent dispersion
        default['turbulent_dispersion_model']    = self.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId)[0]

        # wall forces
        default['wallforcemodel']                = 'tomiyama'
        default['nowallforcemodel']              = 'none'

        return default


    def defaultValuesContinuous(self):
        default = ThermodynamicsModel.defaultValues(self)
        # Drag model
        default['continuousdragmodel']           = 'Large_Interface_Model'

        # gradP correction
        default['gradP_correction_status']       = 'off'
        default['gradP_correction_model']        = 'refined_gradient'
        default['TurbulenceEffect']              = 'off'
        default['FrictionEffect']                = 'off'

        # buubbles forces for LIM
        default['BubblesForLIM']                 = 'off'

        # interface sharpening options
        default['InterfaceSharpening']           = 'off'
        default['UnsharpenedCells']              = 'off'
        default['SurfaceTension']                = 'off'

        return default


    @Variables.noUndo
    def getForceList(self) :
        """
        return list of field couple for forces
        """
        list = []
        for node in self.XMLInterForce.xmlGetNodeList('force') :
            couple=[node['field_id_a'], node['field_id_b']]
            list.append(couple)
        return list


    @Variables.undoGlobal
    def addForce(self) :
        """
        add a new force
        """
        couple = []
        if len(self.getFreeCouples()) > 0 :
            couple = self.getFreeCouples()[0]
            node = self.XMLInterForce.xmlInitChildNode('force', field_id_a = couple[0], field_id_b = couple[1])

            model = ""
            if self.getFieldNature(couple[1]) == "gas" :
                model = self.defaultValues(couple[0], couple[1])['gasdisperseddragmodel']
            else :
                model = self.defaultValues(couple[0], couple[1])['liquidsoliddisperseddragmodel']
            node.xmlInitChildNode('drag_model', model = model)
            if node != "none" :
                Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, couple[1], "drag_coefficient", "drag_coef"+str(couple[1]))
            model = self.defaultValues(couple[0], couple[1])['addedmassmodel']
            node.xmlInitChildNode('added_mass', model = model)

            model = self.defaultValues(couple[0], couple[1])['liftmodel']
            node.xmlInitChildNode('lift_model', model = model)

            if len(self.getAvailableWallForcesModelList(couple[0], couple[1])) > 1:
                model = self.defaultValues(couple[0], couple[1])['wallforcemodel']
            else:
                model = self.defaultValues(couple[0], couple[1])['nowallforcemodel']
            node.xmlInitChildNode('wall_force_model', model = model)

            model = self.defaultValues(couple[0], couple[1])['turbulent_dispersion_model']
            node.xmlInitChildNode('turbulent_dispersion_model', model = model)

            try :
                self.__freeCouples.remove((couple[0], couple[1]))
            except:
                pass
        return couple


    @Variables.undoGlobal
    def setFielda(self, oldfieldaId, fieldbId, fieldaId):
        """
        put field id for field a
        """
        self.isInList([oldfieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = oldfieldaId, field_id_b = fieldbId)
        node['field_id_a'] = fieldaId

        if oldfieldaId != fieldaId :
            if (fieldaId, fieldbId) not in self.__freeCouples :
                new_fieldb = self.getFieldIdbList(fieldaId)[0]
                node['field_id_b'] = new_fieldb
                if self.getFieldNature(fieldbId) != self.getFieldNature(new_fieldb) :
                    self.__updatedrag(fieldaId, new_fieldb)
                self.__freeCouples.remove((fieldaId, new_fieldb))
            else :
                self.__freeCouples.remove((fieldaId, fieldbId))

            self.__freeCouples.append((oldfieldaId, fieldbId))
            self.__updatedrag(fieldaId, fieldbId)
            self.__updateWallForce(fieldaId, fieldbId)


    @Variables.undoGlobal
    def setFieldb(self, fieldaId, oldfieldbId, fieldbId):
        """
        put field id for field b
        """
        self.isInList([fieldaId, oldfieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = oldfieldbId)
        node['field_id_b'] = fieldbId

        if oldfieldbId != fieldbId :
            self.__freeCouples.remove((fieldaId, fieldbId))
            self.__freeCouples.append((fieldaId, oldfieldbId))
            self.__updatedrag(fieldaId, fieldbId)
            self.__updateWallForce(fieldaId, fieldbId)

        if self.getFieldNature(oldfieldbId) != self.getFieldNature(fieldbId) :
            self.__updatedrag(fieldaId, fieldbId)


    def __updatedrag(self, fieldaId, fieldbId) :
        """
        update drag model after fieldIdChange
        """
        model = ""
        if self.getFieldNature(fieldbId) == "gas" :
            model = self.defaultValues(fieldaId, fieldbId)['gasdisperseddragmodel']
        else :
            model = self.defaultValues(fieldaId, fieldbId)['liquidsoliddisperseddragmodel']
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
        self.isInList([fieldaId, fieldbId], self.getForceList())
        self.isInList(model, self.getAvailableDragModels(fieldaId, fieldbId))

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        ChildNode = node.xmlGetChildNode('drag_model')
        ChildNode['model'] = model
        if model != "none" :
            Variables(self.case).setNewVariableProperty("property", "", self.XMLNodeproperty, fieldbId, "drag_coefficient", "drag_coef"+str(fieldbId))
        else :
            Variables(self.case).removeVariableProperty("property", self.XMLNodeproperty, fieldbId, "drag_coefficient")


    @Variables.noUndo
    def getDragModel(self, fieldaId, fieldbId) :
        """
        get drag model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('drag_model')
        if childNode == None :
            model = ""
            if self.getFieldNature(fieldbId) == "gas" :
                model = self.defaultValues(fieldaId, fieldbId)['gasdisperseddragmodel']
            else :
                model = self.defaultValues(fieldaId, fieldbId)['liquidsoliddisperseddragmodel']
            self.setDragModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('drag_model')['model']
        return model


    @Variables.undoLocal
    def setAddMassModel(self, fieldaId, fieldbId, model) :
        """
        set added mass model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())
        self.isInList(model, self.getAvailableAddedMassModels())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        ChildNode = node.xmlGetChildNode('added_mass')
        ChildNode['model'] = model


    @Variables.noUndo
    def getAddMassModel(self, fieldaId, fieldbId) :
        """
        get added mass model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('added_mass')
        if childNode == None :
            model = self.defaultValues(fieldaId, fieldbId)['addedmassmodel']
            self.setAddMassModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('added_mass')['model']
        return model


    @Variables.undoLocal
    def setLiftModel(self, fieldaId, fieldbId, model) :
        """
        set lift model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())
        self.isInList(model, self.getAvailableLiftModels())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        ChildNode = node.xmlGetChildNode('lift_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getLiftModel(self, fieldaId, fieldbId) :
        """
        get lift model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('lift_model')
        if childNode == None :
            model = self.defaultValues(fieldaId, fieldbId)['liftmodel']
            self.setLiftModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('lift_model')['model']
        return model


    @Variables.undoLocal
    def setWallForceModel(self, fieldaId, fieldbId, model) :
        """
        set wall force model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())
        self.isInList(model, self.getAvailableWallForcesModelList(fieldaId, fieldbId))

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        ChildNode = node.xmlGetChildNode('wall_force_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getWallForceModel(self, fieldaId, fieldbId) :
        """
        get wall force model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('wall_force_model')
        if childNode == None :
            if len(self.getAvailableWallForcesModelList(fieldaId, fieldbId)) > 1:
                model = self.defaultValues(fieldaId, fieldbId)['wallforcemodel']
            else:
                model = self.defaultValues(fieldaId, fieldbId)['nowallforcemodel']
            self.setWallForceModel(fieldaId, fieldbId, model)
        model = node.xmlGetNode('wall_force_model')['model']
        return model


    @Variables.undoLocal
    def setTurbDispModel(self, fieldaId, fieldbId, model) :
        """
        set turbulent dispersion model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())
        self.isInList(model, self.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId))

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        ChildNode = node.xmlGetChildNode('turbulent_dispersion_model')
        ChildNode['model'] = model


    @Variables.noUndo
    def getTurbDispModel(self, fieldaId, fieldbId) :
        """
        get turbulent dispersion model for a couple fieldId
        """
        self.isInList([fieldaId, fieldbId], self.getForceList())

        node = self.XMLInterForce.xmlGetNode('force', field_id_a = fieldaId, field_id_b = fieldbId)
        childNode = node.xmlGetNode('turbulent_dispersion_model')
        if childNode == None :
            model = self.defaultValues(fieldaId, fieldbId)['turbulent_dispersion_model']
            self.setWallForceModel(fieldaId, fieldbId, model)
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

        # update free couples
        self.__freeCouples.append((fieldaId, fieldbId))


    @Variables.undoGlobal
    def setContinuousCouplingModel(self, model) :
        """
        set turbulent dispersion model for a couple fieldId
        """
        self.isInList(model, self.getAvailableContinuousDragModelList())
        ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        if ChildNode == None :
            ChildNode = self.XMLInterForce.xmlInitChildNode('continuous_field_momentum_transfer')
        ChildNode['model'] = model
        if model != "separate_phases" :
            ChildNode2 = ChildNode.xmlGetChildNode('gradP_correction')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()
            ChildNode2 = ChildNode.xmlGetChildNode('gradP_correction_model')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()

        if model != "Large_Interface_Model" :
            ChildNode2 = ChildNode.xmlGetChildNode('BubblesForLIM')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()

        if model == "none" :
            ChildNode2 = ChildNode.xmlGetChildNode('InterfaceSharpening')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()
            ChildNode2 = ChildNode.xmlGetChildNode('UnsharpenedCells')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()
            ChildNode2 = ChildNode.xmlGetChildNode('ITMSurfaceTension')
            if ChildNode2 != None :
                ChildNode2.xmlRemoveNode()


    @Variables.noUndo
    def getContinuousCouplingModel(self) :
        """
        get turbulent dispersion model for a couple fieldId
        """
        ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        if ChildNode == None :
            model = self.defaultValuesContinuous()['continuousdragmodel']
            self.setContinuousCouplingModel(model)
            ChildNode = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        model = ChildNode['model']
        return model


    @Variables.undoLocal
    def setGradPCorrectionStatus(self, status) :
        """
        set gradP correction status for separate phases model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        node.xmlSetData('gradP_correction', status)
        if status == "off" :
            ChildNode = node.xmlGetChildNode('gradP_correction_model')
            if ChildNode != None :
                ChildNode.xmlRemoveNode()
        else :
            # to impose a GradP correction model in XML
            self.getGradPCorrectionModel()


    @Variables.noUndo
    def getGradPCorrectionStatus(self) :
        """
        get gradP correction status for separate phases model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('gradP_correction')
        if ChildNode == None :
           status = self.defaultValuesContinuous()['gradP_correction_status']
           self.setGradPCorrectionStatus(status)
        status = node.xmlGetString('gradP_correction')
        return status


    @Variables.undoLocal
    def setGradPCorrectionModel(self, model) :
        """
        set gradP correction model for separate phases model
        """
        self.isInList(model,('refined_gradient', 'refined_gradient_2_criteria'))
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        node.xmlSetData('gradP_correction_model', model)


    @Variables.noUndo
    def getGradPCorrectionModel(self) :
        """
        get gradP correction model for separate phases model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('gradP_correction_model')
        if ChildNode == None :
           model = self.defaultValuesContinuous()['gradP_correction_model']
           self.setGradPCorrectionModel(model)
        model = node.xmlGetString('gradP_correction_model')
        return model


    @Variables.undoLocal
    def setBubblesForLIMStatus(self, status) :
        """
        set bubbles forces for LIM status for continuous momentum transfert model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        old_status = node.xmlGetString('BubblesForLIM')
        node.xmlSetData('BubblesForLIM', status)

        if status == "off" and old_status:
            node = self.XMLClosure.xmlGetChildNode('interfacial_area_diameter')
            if node:
                node.xmlRemoveNode()


    @Variables.noUndo
    def getBubblesForLIMStatus(self) :
        """
        get bubbles forces for LIM status for continuous momentum transfert model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        if node == None :
            node = self.XMLInterForce.xmlInitChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('BubblesForLIM')
        if ChildNode == None :
           status = self.defaultValuesContinuous()['BubblesForLIM']
           self.setBubblesForLIMStatus(status)
        status = node.xmlGetString('BubblesForLIM')
        return status


    @Variables.undoLocal
    def setInterfaceSharpeningStatus(self, status) :
        """
        set interface sharpening status for continuous momentum transfert model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        node.xmlSetData('InterfaceSharpening', status)


    @Variables.noUndo
    def getInterfaceSharpeningStatus(self) :
        """
        get interface sharpening status for continuous momentum transfert model
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('InterfaceSharpening')
        if ChildNode == None :
           status = self.defaultValuesContinuous()['InterfaceSharpening']
           self.setInterfaceSharpeningStatus(status)
        status = node.xmlGetString('InterfaceSharpening')
        return status


    @Variables.undoLocal
    def setUnsharpenedCellsStatus(self, status) :
        """
        set unsharpened cells status
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        node.xmlSetData('UnsharpenedCells', status)


    @Variables.noUndo
    def getUnsharpenedCellsStatus(self) :
        """
        get unsharpened cells status
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('UnsharpenedCells')
        if ChildNode == None :
           status = self.defaultValuesContinuous()['UnsharpenedCells']
           self.setUnsharpenedCellsStatus(status)
        status = node.xmlGetString('UnsharpenedCells')
        return status


    @Variables.undoLocal
    def setSurfaceTensionStatus(self, status) :
        """
        set
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        node.xmlSetData('ITMSurfaceTension', status)


    @Variables.noUndo
    def getSurfaceTensionStatus(self) :
        """
        get
        """
        node = self.XMLInterForce.xmlGetChildNode('continuous_field_momentum_transfer')
        ChildNode = node.xmlGetChildNode('ITMSurfaceTension')
        if ChildNode == None :
           status = self.defaultValuesContinuous()['SurfaceTension']
           self.setSurfaceTensionStatus(status)
        status = node.xmlGetString('ITMSurfaceTension')
        return status


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


    def checkGetFreeCouples(self):
        """Check whether the InterfacialEnthalpyModel class could get FreeCouples"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getFreeCouples() ==[('1', '3'), ('2', '3')],\
            'Could not get FreeCouples'


    def checkGetFieldIdaList(self):
        """Check whether the InterfacialEnthalpyModel class could get FieldIdaList"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getFieldIdaList('1') == ['1', '2'],\
            'Could not get FieldIdaList'
        assert mdl.getFieldIdaList('2') == ['2', '1'],\
            'Could not get FieldIdaList'
        assert mdl.getFieldIdaList('3') == ['3', '1', '2'],\
            'Could not get FieldIdaList'


    def checkGetFieldIdbList(self):
        """Check whether the InterfacialEnthalpyModel class could get FieldIdbList"""
        MainFieldsModel(self.case).addField()
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getFieldIdbList('1') ==['3'] ,\
            'Could not get FieldIdbList'
        assert mdl.getFieldIdbList('2') == ['3'] ,\
            'Could not get FieldIdbList'
        assert mdl.getFieldIdbList('3') == [],\
            'Could not get FieldIdbList'


    def checkGetAvailableContinuousDragModelList(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableContinuousDragModelList"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableContinuousDragModelList() == ["none", "separate_phases" ,"Large_Interface_Model", "Large_Bubble_Model"],\
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
        assert mdl.getAvailableDragModels('1','2') == ["none", "ishii"],\
            'Could not get AvailableDragModels'
        MainFieldsModel(self.case).setFieldNature('2','solid')
        assert mdl.getAvailableDragModels('1','2') == ["none", "inclusions", "Wen_Yu"],\
            'Could not get AvailableDragModels'


    def checkGetAvailableAddedMassModels(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableAddedMassModels"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableAddedMassModels() == ["none", "standard", "zuber"],\
            'Could not get AvailableAddedMassModels'


    def checkGetAvailableLiftModels(self):
        """Check whether the InterfacialEnthalpyModel class could get the AvailableLiftModels"""
        mdl = InterfacialForcesModel(self.case)
        assert mdl.getAvailableLiftModels() == ["none", "coef_cst", "Tomiyama_SMD"],\
            'Could not get AvailableLiftModels'


    def checkGetForceList(self):
        """Check whether the InterfacialEnthalpyModel class could get force list"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.addForce()
        assert mdl.getForceList() == [['1', '2']],\
            'Could not get force list'


    def checkaddForce(self):
        """Check whether the InterfacialEnthalpyModel class could add force"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.addForce()
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
            'Could not add force'


    def checksetFielda(self):
        """Check whether the InterfacialEnthalpyModel class could setFielda"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'continuous', 'liquid', 'on', 'on', 'off', 3)
        mdl = InterfacialForcesModel(self.case)
        mdl.addForce()
        mdl.setFielda('1','2','3')
        doc = '''<interfacial_forces>
                         <force field_id_a="3" field_id_b="2">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set Fielda'


    def checksetFieldb(self):
        """Check whether the InterfacialEnthalpyModel class could setFieldb"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        MainFieldsModel(self.case).addDefinedField('3', 'field3', 'dispersed', 'gas', 'on', 'on', 'off', 3)
        mdl = InterfacialForcesModel(self.case)
        mdl.addForce()
        mdl.setFieldb('1','2','3')
        doc = '''<interfacial_forces>
                         <force field_id_a="1" field_id_b="3">
                                 <drag_model model="ishii"/>
                                 <added_mass model="none"/>
                                 <lift_model model="none"/>
                                 <wall_force_model model="none"/>
                                 <turbulent_dispersion_model model="none"/>
                         </force>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set Fieldb'


    def checkGetandSetDragModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get DragModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'dispersed', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.addForce()
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
        mdl.addForce()
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
        mdl.addForce()
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
        mdl.addForce()
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
        mdl.addForce()
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
        mdl.addForce()
        mdl.deleteForce('1','2')
        doc = '''<interfacial_forces/>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not delete Force'


    def checkGetandSetContinuousCouplingModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get ContinuousCouplingModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('Large_Interface_Model')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="Large_Interface_Model"/>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set ContinuousCouplingModel'
        assert mdl.getContinuousCouplingModel() == 'Large_Interface_Model',\
            'Could not get ContinuousCouplingModel'


    def checkGetandSetGradPCorrectionStatus(self):
        """Check whether the InterfacialEnthalpyModel class could set and get GradPCorrectionStatus"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('separate_phases')
        mdl.setGradPCorrectionStatus('on')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="separate_phases">
                                 <gradP_correction>
                                         on
                                 </gradP_correction>
                                 <gradP_correction_model>
                                         refined_gradient
                                 </gradP_correction_model>
                         </continuous_field_momentum_transfer>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set GradPCorrectionStatus'
        assert mdl.getGradPCorrectionStatus() == 'on',\
            'Could not set GradPCorrectionStatus'


    def checkGetandSetGradPCorrectionModel(self):
        """Check whether the InterfacialEnthalpyModel class could set and get GradPCorrectionModel"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('separate_phases')
        mdl.setGradPCorrectionModel('refined_gradient_2_criteria')
        doc = '''<interfacial_forces>
                         <continuous_field_momentum_transfer model="separate_phases">
                                 <gradP_correction_model>
                                         refined_gradient_2_criteria
                                 </gradP_correction_model>
                         </continuous_field_momentum_transfer>
                 </interfacial_forces>'''
        assert mdl.XMLInterForce == self.xmlNodeFromString(doc),\
            'Could not set GradPCorrectionModel'
        assert mdl.getGradPCorrectionModel() == 'refined_gradient_2_criteria',\
            'Could not set GradPCorrectionModel'


    def checkGetandSetTurbulenceEffectStatus(self):
        """Check whether the InterfacialEnthalpyModel class could set and get TurbulenceEffectStatus"""
        MainFieldsModel(self.case).addDefinedField('1', 'field1', 'continuous', 'liquid', 'on', 'on', 'off', 1)
        MainFieldsModel(self.case).addDefinedField('2', 'field2', 'continuous', 'gas', 'on', 'on', 'off', 2)
        mdl = InterfacialForcesModel(self.case)
        mdl.setContinuousCouplingModel('Large_Interface_Model')
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
        mdl.setContinuousCouplingModel('Large_Interface_Model')
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
