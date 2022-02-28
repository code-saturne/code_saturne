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
This module defines the 'Interfacial enthalpy transfer' page.

This module contains the following classes:
- FieldDelegate
- StandardItemModelInterfacialEnthalpy
- InterfacialEnthalpyView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator, BasicTableModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.InterfacialEnthalpy import Ui_InterfacialEnthalpy
from code_saturne.model.InterfacialEnthalpyModel import InterfacialEnthalpyModel
from code_saturne.model.NonCondensableModel import NonCondensableModel
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel
from code_saturne.model.NucleateBoilingModel import NucleateBoilingModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialEnthalpyView")
log.setLevel(GuiParam.DEBUG)

# -------------------------------------------------------------------------------
# InterfacialEnthalpyView class
# -------------------------------------------------------------------------------

class InterfacialEnthalpyView(QWidget, Ui_InterfacialEnthalpy):
    """
    InterfacialEnthalpyView layout.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InterfacialEnthalpy.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = InterfacialEnthalpyModel(self.case)

        # Combo models
        liquid_vapor_couples = self.mdl.getLiquidVaporCouples()
        self.modelLiquidVaporFields = ComboModel(self.comboBoxLiquidVaporFields, 1, 1)
        self.modelLiquidVaporFields.addItem("None", "none")
        for id_pair in liquid_vapor_couples:
            field_descriptions = []
            for id in id_pair:
                name = self.mdl.getLabel(id)
                phase = self.mdl.getFieldNature(id)
                criterion = self.mdl.getCriterion(id)
                field_descriptions.append("{0} ({1} {2})".format(name, criterion, phase))
            # Display fields as "Field1 (continuous liquid) / Field2 (dispersed gas)"
            view = self.tr(" / ".join(field_descriptions))
            model = "{0}_{1}".format(*id_pair)
            self.modelLiquidVaporFields.addItem(view, model)

        self.modelPonderationCoefFielda = ComboModel(self.comboBoxPonderationCoefFielda, 3, 1)
        self.modelPonderationCoefFielda.addItem(self.tr("alp1"), "alp1")
        self.modelPonderationCoefFielda.addItem(self.tr("alp2"), "alp2")
        self.modelPonderationCoefFielda.addItem(self.tr("alp1*alp2"), "alp1_alp2")

        self.modelPonderationCoefFieldb = ComboModel(self.comboBoxPonderationCoefFieldb, 3, 1)
        self.modelPonderationCoefFieldb.addItem(self.tr("alp1"), "alp1")
        self.modelPonderationCoefFieldb.addItem(self.tr("alp2"), "alp2")
        self.modelPonderationCoefFieldb.addItem(self.tr("alp1*alp2"), "alp1_alp2")

        self.modelSolidEnergyTransfer = ComboModel(self.comboBoxSolidEnergyTransfer, 2, 1)
        self.modelSolidEnergyTransfer.addItem(self.tr("none"), "none")
        self.modelSolidEnergyTransfer.addItem(self.tr("gas-particle"), "gas_particule")

        self.modelFieldaModel = ComboModel(self.comboBoxFieldaModel, 1, 1)
        self.modelFieldaModel.addItem(self.tr("No source term"), "no_source_term")
        self.modelFieldaModel.addItem(self.tr("Relaxation time : alphk.rok.cpk.(Ts-Tk)/tauk"), "relaxation_time")
        self.modelFieldbModel = ComboModel(self.comboBoxFieldbModel, 1, 1)
        self.modelFieldbModel.addItem(self.tr("No source term"), "no_source_term")
        self.modelFieldbModel.addItem(self.tr("Relaxation time : alphk.rok.cpk.(Ts-Tk)/tauk"), "relaxation_time")

        self.groupBoxLiquidVaporModel.hide()
        self.groupBoxSolidEnergyTransfer.hide()
        if (len(self.mdl.getSolidFieldIdList()) > 0):
            model = self.mdl.getSolidEnergyTransfer()
            self.modelSolidEnergyTransfer.setItem(str_model=model)
            self.groupBoxSolidEnergyTransfer.show()
            self.groupBoxEnergyTransfer.hide()

        self.setValidators()
        self.setConnections()

        # Initial state of Pool boiling model
        if self.mdl.getPoolBoiling() == 'on':
            self.checkBoxActivatePool.setChecked(True)

        # Initialize pair of fields
        predefined_flow = self.mdl.getPredefinedFlow()
        if predefined_flow == "None":
            if self.mdl.getEnthalpyCoupleFieldId() is None:
                self.modelLiquidVaporFields.setItem(str_model="none")
            else:
                field_id_a, field_id_b = self.mdl.getEnthalpyCoupleFieldId()
                model = "{0}_{1}".format(field_id_a, field_id_b)
                self.modelLiquidVaporFields.setItem(str_model=model)
        elif predefined_flow == "free_surface":
            self.lockFreeSurfaceOptions()
        elif predefined_flow == "boiling_flow":
            self.lockBubblyFlowOptions()
        elif predefined_flow == "droplet_flow":
            self.lockDropletFlowOptions()
        elif predefined_flow == "particles_flow":
            # self.lockParticlesFlowOptions()
            # TODO check if options should be locked or not
            pass
        elif predefined_flow == "multiregime":
            self.lockMultiregimeFlowOptions()

        self.case.undoStartGlobal()

    def lockMultiregimeFlowOptions(self):
        self.modelLiquidVaporFields.setItem(str_model="1_2")
        self.modelFieldaModel.setItem(str_model="wall_law_type_model")
        self.modelFieldbModel.setItem(str_model="sublayer_LI3C")
        self.comboBoxLiquidVaporFields.setEnabled(False)
        self.comboBoxFieldaModel.setEnabled(False)
        self.comboBoxFieldbModel.setEnabled(False)

    def lockParticlesFlowOptions(self):
        self.modelSolidEnergyTransfer.setItem(str_model="gas_particule")
        self.comboBoxSolidEnergyTransfer.setEnabled(False)

    def lockDropletFlowOptions(self):
        self.modelLiquidVaporFields.setItem(str_model="1_2")
        self.modelFieldaModel.setItem(str_model="droplet_model_for_vapour")
        self.modelFieldbModel.setItem(str_model="droplet_model_for_liquid")
        self.comboBoxLiquidVaporFields.setEnabled(False)
        self.comboBoxFieldaModel.setEnabled(False)
        self.comboBoxFieldbModel.setEnabled(False)

    def lockBubblyFlowOptions(self):
        self.modelLiquidVaporFields.setItem(str_model="1_2")
        self.modelFieldaModel.setItem(str_model="bubble_model_for_liquid")
        self.modelFieldbModel.setItem(str_model="relaxation_time_subcooled")
        self.comboBoxLiquidVaporFields.setEnabled(False)
        self.comboBoxFieldaModel.setEnabled(False)
        self.comboBoxFieldbModel.setEnabled(False)

    def lockFreeSurfaceOptions(self):
        self.modelLiquidVaporFields.setItem(str_model="1_2")
        self.modelFieldaModel.setItem(str_model="wall_law_type_model")
        self.modelFieldbModel.setItem(str_model="sublayer_LI3C")
        self.comboBoxLiquidVaporFields.setEnabled(False)
        self.comboBoxFieldaModel.setEnabled(False)
        self.comboBoxFieldbModel.setEnabled(False)

    def setValidators(self):
        validatorRelaxa = DoubleValidator(self.lineEditRelaxationTimeFielda, min=0.0)
        validatorRelaxb = DoubleValidator(self.lineEditRelaxationTimeFieldb, min=0.0)
        self.lineEditRelaxationTimeFielda.setValidator(validatorRelaxa)
        self.lineEditRelaxationTimeFieldb.setValidator(validatorRelaxb)

    def setConnections(self):
        self.comboBoxLiquidVaporFields.currentTextChanged[str].connect(self.slotSelectInteraction)
        self.comboBoxSolidEnergyTransfer.currentTextChanged[str].connect(self.slotSolidEnergyTransfer)
        self.comboBoxFieldaModel.currentTextChanged[str].connect(self.slotFieldaModel)
        self.comboBoxPonderationCoefFielda.activated[str].connect(self.slotPonderationCoefFielda)
        self.comboBoxFieldbModel.currentTextChanged[str].connect(self.slotFieldbModel)
        self.comboBoxPonderationCoefFieldb.activated[str].connect(self.slotPonderationCoefFieldb)
        self.lineEditRelaxationTimeFielda.textChanged[str].connect(self.slotRelaxationTimeFielda)
        self.lineEditRelaxationTimeFieldb.textChanged[str].connect(self.slotRelaxationTimeFieldb)
        self.checkBoxActivatePool.stateChanged.connect(self.slotPoolBoilingModel)

    @pyqtSlot(str)
    def slotSelectInteraction(self, value):
        """
        Select a Field in the QTable
        """
        selection = self.modelLiquidVaporFields.dicoV2M[value]
        if selection == "none":
            self.groupBoxLiquidVaporModel.hide()
            self.mdl.deleteLiquidVaporEnthalpyTransfer()
            return
        fieldaId, fieldbId = selection.split("_")
        if self.mdl.getEnthalpyCoupleFieldId() is None:
            self.mdl.addLiquidVaporEnthalpyTransfer(fieldaId, fieldbId)
        self.mdl.setEnthalpyCoupleFieldId(fieldaId, fieldbId)
        self.updateLiquidVaporModel()
        self.groupBoxLiquidVaporModel.show()
        ifm = InterfacialForcesModel(self.case)

        if len(NonCondensableModel(self.case).getNonCondensableLabelList()) > 0 \
                and ifm.getContinuousCouplingModel(fieldaId, fieldbId) in \
                ifm.getAvailableContinuousDragModelList():
            self.checkBoxActivatePool.show()
        else:
            self.checkBoxActivatePool.hide()

    def updateLiquidVaporModel(self):
        """
        update if necessary
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        self.fillLiquidVaporModels(fieldIda, fieldIdb)

        model = self.mdl.getFieldModel(fieldIda)
        if model is None:
            model = "no_source_term"
        self.modelFieldaModel.setItem(str_model=model)

        if model == 'relaxation_time':
            model = self.mdl.getPonderationCoef(fieldIda)
            self.modelPonderationCoefFielda.setItem(str_model=model)
            value = self.mdl.getRelaxationTime(fieldIda)
            self.lineEditRelaxationTimeFielda.setText(str(value))

            self.comboBoxPonderationCoefFielda.show()
            self.labelPonderationCoefFielda.show()
            self.lineEditRelaxationTimeFielda.show()
            self.labelRelaxationTimeFielda.show()
        else :
            self.comboBoxPonderationCoefFielda.hide()
            self.labelPonderationCoefFielda.hide()
            self.lineEditRelaxationTimeFielda.hide()
            self.labelRelaxationTimeFielda.hide()

        model = self.mdl.getFieldModel(fieldIdb)
        self.modelFieldbModel.setItem(str_model = model)

        if model == 'relaxation_time' :
            model = self.mdl.getPonderationCoef(fieldIdb)
            self.modelPonderationCoefFieldb.setItem(str_model = model)
            value = self.mdl.getRelaxationTime(fieldIdb)
            self.lineEditRelaxationTimeFieldb.setText(str(value))

            self.comboBoxPonderationCoefFieldb.show()
            self.labelPonderationCoefFieldb.show()
            self.lineEditRelaxationTimeFieldb.show()
            self.labelRelaxationTimeFieldb.show()
        else:
            self.comboBoxPonderationCoefFieldb.hide()
            self.labelPonderationCoefFieldb.hide()
            self.lineEditRelaxationTimeFieldb.hide()
            self.labelRelaxationTimeFieldb.hide()

    def fillLiquidVaporModels(self, fieldIda, fieldIdb):
        # fieldIda == continuous by construction and nature(fieldIda) != nature(fieldIdb)
        if self.mdl.getFieldNature(fieldIda) == "liquid":
            if self.mdl.getCriterion(fieldIdb) == "continuous":
                self.modelFieldaModel.addItem(self.tr("Coste-Lavieville NURETH13 model, Wall Law Type Model"),
                                              "wall_law_type_model")
                self.modelFieldbModel.addItem(self.tr("Interfacial Sublayer Model for LI3C"), "sublayer_LI3C")
            else:
                self.modelFieldaModel.addItem(self.tr("Bulk model (Ranz-Marshall)"), "bulk")
                self.modelFieldaModel.addItem(self.tr("Flashing (Cathare)"), "flashing")
                self.modelFieldaModel.addItem(self.tr("Bubble model for liquid (Manon-Berne)"),
                                              "bubble_model_for_liquid")
                # suppression temporaire le temps que le modele soit au point
                # self.modelFieldbModel.addItem(self.tr("Bubble model for vapour"),"bubble_model_for_vapour")
                self.modelFieldbModel.addItem(self.tr("Relaxation time + subcooled gas treatment"),
                                              "relaxation_time_subcooled")
        else:
            if self.mdl.getCriterion(fieldIdb) == "continuous":
                self.modelFieldaModel.addItem(self.tr("Interfacial Sublayer Model for LI3C"), "sublayer_LI3C")
                self.modelFieldbModel.addItem(self.tr("Coste-Lavieville NURETH13 model, Wall Law Type Model"),
                                              "wall_law_type_model")
            else:
                self.modelFieldaModel.addItem(self.tr("Bulk model (Ranz-Marshall)"), "bulk")
                self.modelFieldaModel.addItem(self.tr("Droplet model for vapour"), "droplet_model_for_vapour")
                self.modelFieldbModel.addItem(self.tr("Droplet model for liquid"), "droplet_model_for_liquid")

    @pyqtSlot(str)
    def slotSolidEnergyTransfer(self, text):
        """
        set model for solid enthalpy transfer
        """
        choice = self.modelSolidEnergyTransfer.dicoV2M[text]
        self.mdl.setSolidEnergyTransfer(choice)

    @pyqtSlot(str)
    def slotFieldaModel(self, text):
        """
        set model for field a
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        choice = self.modelFieldaModel.dicoV2M[text]
        self.mdl.setFieldModel(fieldIda, choice)
        self.updateLiquidVaporModel()


    @pyqtSlot(str)
    def slotPonderationCoefFielda(self, text):
        """
        set ponderation coefficient for field a
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        choice = self.modelPonderationCoefFielda.dicoV2M[text]
        self.mdl.setPonderationCoef(fieldIda, choice)


    @pyqtSlot(str)
    def slotFieldbModel(self, text):
        """
        set model for field b
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        choice = self.modelFieldbModel.dicoV2M[text]
        self.mdl.setFieldModel(fieldIdb, choice)
        self.updateLiquidVaporModel()


    @pyqtSlot(str)
    def slotPonderationCoefFieldb(self, text):
        """
        set ponderation coefficient for field b
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        choice = self.modelPonderationCoefFieldb.dicoV2M[text]
        self.mdl.setPonderationCoef(fieldIdb, choice)


    @pyqtSlot(str)
    def slotRelaxationTimeFielda(self, text):
        """
        Update the relaxation time for field a
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        if self.lineEditRelaxationTimeFielda.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setRelaxationTime(fieldIda, value)


    @pyqtSlot(str)
    def slotRelaxationTimeFieldb(self, text):
        """
        Update the relaxation time for field b
        """
        selection = self.modelLiquidVaporFields.dicoV2M[self.comboBoxLiquidVaporFields.currentText()]
        fieldIda, fieldIdb = selection.split("_")
        if self.lineEditRelaxationTimeFieldb.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setRelaxationTime(fieldIdb, value)

    @pyqtSlot()
    def slotPoolBoilingModel(self):
        """
        Activate or deactivate the pool boiling model
        """
        pool_boiling_state = self.checkBoxActivatePool.isChecked()
        self.mdl.setPoolBoiling(state = pool_boiling_state)
        wall_model = NucleateBoilingModel(self.case)
        if pool_boiling_state:
            wall_model.setWallFunctionModel("standard")
        else:
            if self.mdl.getPredefinedFlow() != "None":
                wall_model.setWallFunctionModel(wall_model.defaultValues()["wallfunction"])
