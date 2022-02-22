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
This module defines the 'Interfacial forces' page.

This module contains the following classes:
- FieldDelegate
- DragDelegate
- AddedMassDelegate
- LiftDelegate
- DispersionTurbulentDelegate
- WallForceDelegate
- StandardItemModelInterfacialForces
- InterfacialForcesView
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

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, from_qvariant, to_text_string, BasicTableModel, DoubleValidator
from code_saturne.gui.case.InterfacialForces import Ui_InterfacialForces
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel

from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialForcesView")
log.setLevel(GuiParam.DEBUG)

# -------------------------------------------------------------------------------
# AbstractTableModelInteractions class
# -------------------------------------------------------------------------------

class InteractionsTableModel(BasicTableModel):
    """
    Display the list of interactions between pairs of phases
    """

    def __init__(self, parent, data, xml_model):
        super(InteractionsTableModel, self).__init__(parent, xml_model, data,
                                                     ["Field A", "Field B", "Type of interaction"],
                                                     ["A", "B", "T"])


# -------------------------------------------------------------------------------
#  class InterfacialForces
# -------------------------------------------------------------------------------

class InterfacialForcesView(QWidget, Ui_InterfacialForces):
    """
    Main fields layout.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InterfacialForces.__init__(self)
        self.setupUi(self)

        self.case = case
        self.browser = tree
        self.case.undoStopGlobal()
        self.mdl = InterfacialForcesModel(self.case)
        self.field_id_a = None
        self.field_id_b = None

        # Dico
        self.dicoM2V = {"none": 'None',
                        "LLB_model": 'LLB model',
                        "GTD_model": 'GTD model',
                        "antal": 'Antal',
                        "tomiyama": 'Tomiyama',
                        "ishii": 'Ishii',
                        "Gobin": 'Gobin et al.',
                        "Wen_Yu": 'Wen and Yu',
                        "standard": 'Standard',
                        "zuber": 'Zuber',
                        "coef_cst": 'Constant coefficient',
                        "Tomiyama_SMD": 'Tomiyama SMD',
                        "Zeng_Baalbaki": "Zeng and Baalbaki",
                        "Large_Interface_Model": "Large Interface Model",
                        "Large_Bubble_Model": "Large Bubble Model",
                        "G_Large_Interface_Model": "Generalized Large Interface Model"}

        self.dicoV2M = {value: key for key, value in self.dicoM2V.items()}

        allCouples = self.mdl.getAllCouples()

        self.groupBoxDispersedMomentumTransfer.hide()
        self.groupBoxContinuousMomentumTransfer.hide()

        self.tableModelInteractions = InteractionsTableModel(self, allCouples, self.mdl)

        # Dispersed models
        self.modelDispersedDrag = ComboModel(self.comboBoxDispersedDrag, 1, 1)
        self.modelLift = ComboModel(self.comboBoxLift, 1, 1)
        self.modelAddedMass = ComboModel(self.comboBoxAddedMass, 1, 1)
        self.modelTurbulenceDispersion = ComboModel(self.comboBoxTurbulenceDispersion, 1, 1)
        self.modelWallForce = ComboModel(self.comboBoxWallForce, 1, 1)

        # Continuous models
        self.modelContinuousMomentumTransfer = ComboModel(self.comboBoxContinuousMomentumTransfer, 3, 1)
        self.modelInterfaceSharpening = ComboModel(self.comboBoxInterfaceSharpening, 1, 1)
        self.modelSurfaceTension = ComboModel(self.comboBoxSurfaceTension, 1, 1)
        self._fillContinuousComboBoxes()

        self._connectSignalsToSlots()
        self._initializeInteractionsTable()

        self.case.undoStartGlobal()

    def _initializeInteractionsTable(self):
        self.tableViewInteractions.setModel(self.tableModelInteractions)

        self.tableViewInteractions.resizeColumnsToContents()
        self.tableViewInteractions.resizeRowsToContents()
        self.tableViewInteractions.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewInteractions.setSelectionMode(QAbstractItemView.SingleSelection)
        if QT_API == "PYQT4":
            self.tableViewInteractions.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInteractions.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInteractions.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewInteractions.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInteractions.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInteractions.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

    def _connectSignalsToSlots(self):
        self.tableModelInteractions.dataChanged.connect(self.dataChanged)
        self.tableViewInteractions.clicked[QModelIndex].connect(self.slotSelectInteraction)

        self.comboBoxContinuousMomentumTransfer.currentTextChanged[str].connect(self.slotContinuousMomentumTransfer)
        self.comboBoxInterfaceSharpening.activated[str].connect(self.slotInterfaceSharpening)

        self.comboBoxDispersedDrag.activated[str].connect(self.slotDispersedDrag)
        self.comboBoxLift.activated[str].connect(self.slotLift)
        self.comboBoxAddedMass.activated[str].connect(self.slotAddedMass)
        self.comboBoxTurbulenceDispersion.activated[str].connect(self.slotTurbulenceDispersion)
        self.comboBoxWallForce.activated[str].connect(self.slotWallForce)

        self.comboBoxSurfaceTension.currentTextChanged[str].connect(self.slotSurfaceTensionModel)

    def _displayDispersedModels(self):
        self.groupBoxDispersedMomentumTransfer.show()
        self.groupBoxContinuousMomentumTransfer.hide()
        self._fillDispersedComboBoxes()

    def _fillDispersedComboBoxes(self):
        self.modelDispersedDrag.flushItems()
        self.modelLift.flushItems()
        self.modelAddedMass.flushItems()
        self.modelTurbulenceDispersion.flushItems()
        self.modelWallForce.flushItems()
        self.modelDispersedDrag.addItemList([[self.dicoM2V[model], model]
                                             for model in
                                             self.mdl.getAvailableDragModels(self.field_id_a, self.field_id_b)])
        self.modelLift.addItemList([[self.dicoM2V[model], model]
                                    for model in self.mdl.getAvailableLiftModels()])
        self.modelAddedMass.addItemList([[self.dicoM2V[model], model]
                                         for model in self.mdl.getAvailableAddedMassModels()])
        self.modelTurbulenceDispersion.addItemList([[self.dicoM2V[model], model]
                                                    for model in
                                                    self.mdl.getAvailableTurbulenteDispersionModelList(self.field_id_a,
                                                                                                       self.field_id_b)])
        self.modelWallForce.addItemList([[self.dicoM2V[model], model]
                                         for model in
                                         self.mdl.getAvailableWallForcesModelList(self.field_id_a, self.field_id_b)])

        # Read data from XML and update comboBoxes accordingly
        self.modelDispersedDrag.setItem(str_model=self.mdl.getDragModel(self.field_id_a, self.field_id_b))
        self.modelLift.setItem(str_model=self.mdl.getLiftModel(self.field_id_a, self.field_id_b))
        self.modelAddedMass.setItem(str_model=self.mdl.getAddMassModel(self.field_id_a, self.field_id_b))
        self.modelTurbulenceDispersion.setItem(str_model=self.mdl.getTurbDispModel(self.field_id_a, self.field_id_b))
        self.modelWallForce.setItem(str_model=self.mdl.getWallForceModel(self.field_id_a, self.field_id_b))

    def _displayContinuousModels(self):
        self.groupBoxContinuousMomentumTransfer.show()
        self.groupBoxDispersedMomentumTransfer.hide()

        model = self.mdl.getContinuousCouplingModel(self.field_id_a, self.field_id_b)
        self.comboBoxContinuousMomentumTransfer.setEnabled(True)
        self.modelContinuousMomentumTransfer.setItem(str_model=model)
        sharpening = self.mdl.getInterfaceSharpeningModel(self.field_id_a, self.field_id_b)
        self.modelInterfaceSharpening.setItem(str_model=sharpening)
        surface_tension = self.mdl.getSurfaceTensionModel(self.field_id_a, self.field_id_b)
        self.modelSurfaceTension.setItem(str_model=surface_tension)

    def _fillContinuousComboBoxes(self):
        self.modelContinuousMomentumTransfer.addItemList([[self.dicoM2V[model], model]
                                                          for model in
                                                          self.mdl.getAvailableContinuousDragModelList()])
        # TODO : move interface sharpening and surface tension options to InterfacialForcesModel.py
        self.modelInterfaceSharpening.addItemList(
            [[self.tr('None'), 'none'],
             [self.tr('Olsson Interface Sharpening'), 'Olsson_Interface_Sharpening'],
             [self.tr('Olsson Partial Interface Sharpening'), "Olsson_Partial_Interface_Sharpening"],
             [self.tr('Conservative Interface Sharpening (Lavieville)'), 'Conservative_Interface_Sharpening']]
        )
        self.modelSurfaceTension.addItemList(
            [[self.tr('None'), 'none'],
             [self.tr('Brackbill'), 'Brackbill']]
        )

    def lockFreeSurfaceOptions(self):
        default_free_surface = "Large_Interface_Model"
        self.comboBoxContinuousMomentumTransfer.setEnabled(False)

    def lockBubblyFlowOptions(self):
        GTD_condition_1 = (TurbulenceModel(self.case).getTurbulenceModel("1") in
                           ["k-epsilon",
                            "k-epsilon_linear_production",
                            "rij-epsilon_ssg",
                            "rij-epsilon_ebrsm"])

        GTD_condition_2 = (TurbulenceModel(self.case).getTurbulenceModel("2") == "none")

        if GTD_condition_1 and GTD_condition_2:
            self.modelTurbulenceDispersion.setItem(str_model="GTD_model")
            self.comboBoxTurbulenceDispersion.activated[str].emit("GTD model")
        self.comboBoxDispersedDrag.setEnabled(False)
        self.comboBoxLift.setEnabled(False)
        self.comboBoxAddedMass.setEnabled(False)
        self.comboBoxTurbulenceDispersion.setEnabled(False)
        self.comboBoxWallForce.setEnabled(False)

    def lockDropletFlowOptions(self):
        self.comboBoxDispersedDrag.setEnabled(False)

        self.modelLift.disableItem(str_model="Tomiyama_SMD")
        self.modelLift.disableItem(str_model="coef_cst")
        self.comboBoxLift.setEnabled(True)

        self.comboBoxAddedMass.setEnabled(False)
        self.comboBoxTurbulenceDispersion.setEnabled(False)
        self.comboBoxWallForce.setEnabled(False)

    def lockMultiregimeFlowOptions(self):
        default_free_surface = "G_Large_Interface_Model"
        self.comboBoxContinuousMomentumTransfer.setEnabled(False)

    def dataChanged(self, topLeft, bottomRight):
        self.tableViewInteractions.resizeColumnsToContents()
        self.tableViewInteractions.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewInteractions.horizontalHeader().setResizeMode(0, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewInteractions.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)

    @pyqtSlot("QModelIndex")
    def slotSelectInteraction(self, index):
        if not(index.isValid()):
            self.groupBoxContinuousMomentumTransfer.hide()
            self.groupBoxDispersedMomentumTransfer.hide()
            return

        field_a, field_b, interaction_type = self.tableModelInteractions.data_table[index.row()]
        self.field_id_a = self.mdl.getFieldId(field_a)
        self.field_id_b = self.mdl.getFieldId(field_b)

        if interaction_type == "continuous":
            self._displayContinuousModels()
        elif interaction_type == "dispersed":
            self._displayDispersedModels()
        else:
            raise TypeError("Unknown interaction type : {0}".format(interaction_type))
        predefined_flow = self.mdl.getPredefinedFlow()
        if predefined_flow == "free_surface":
            self.lockFreeSurfaceOptions()
        elif predefined_flow == "boiling_flow":
            self.lockBubblyFlowOptions()
        elif predefined_flow == "droplet_flow":
            self.lockDropletFlowOptions()
        elif predefined_flow == "multiregime":
            self.lockMultiregimeFlowOptions()

    @pyqtSlot(str)
    def slotContinuousMomentumTransfer(self, text):
        """
        configure momentum transfer for continuous phases
        """
        value = self.modelContinuousMomentumTransfer.dicoV2M[text]
        log.debug("slotContinuousMomentumTransfer -> %s" % value)

        self.mdl.setContinuousCouplingModel(self.field_id_a, self.field_id_b, value)

    @pyqtSlot(str)
    def slotInterfaceSharpening(self, text):
        model = self.modelInterfaceSharpening.dicoV2M[text]
        log.debug("slotInterfaceSharpening -> %s" % model)
        self.mdl.setInterfaceSharpeningModel(self.field_id_a, self.field_id_b, model)

    @pyqtSlot(str)
    def slotSurfaceTensionModel(self, text):
        model = self.modelSurfaceTension.dicoV2M[text]
        log.debug("slotSurfaceTension -> %s" % model)
        self.mdl.setSurfaceTensionModel(self.field_id_a, self.field_id_b, model)


    @pyqtSlot(str)
    def slotDispersedDrag(self, text):
        model = self.modelDispersedDrag.dicoV2M[text]
        log.debug("slotDispersedDrag -> %s" % model)
        self.mdl.setDragModel(self.field_id_a, self.field_id_b, model)

    @pyqtSlot(str)
    def slotLift(self, text):
        model = self.modelLift.dicoV2M[text]
        log.debug("slotLift -> %s" % model)
        self.mdl.setLiftModel(self.field_id_a, self.field_id_b, model)

    @pyqtSlot(str)
    def slotAddedMass(self, text):
        model = self.modelAddedMass.dicoV2M[text]
        log.debug("slotAddedMass -> %s" % model)
        self.mdl.setAddMassModel(self.field_id_a, self.field_id_b, model)

    @pyqtSlot(str)
    def slotTurbulenceDispersion(self, text):
        model = self.modelTurbulenceDispersion.dicoV2M[text]
        log.debug("slotTurbulenceDispersion -> %s" % model)
        self.mdl.setTurbDispModel(self.field_id_a, self.field_id_b, model)

    @pyqtSlot(str)
    def slotWallForce(self, text):
        model = self.modelWallForce.dicoV2M[text]
        log.debug("slotWallForce -> %s" % model)
        self.mdl.setWallForceModel(self.field_id_a, self.field_id_b, model)
