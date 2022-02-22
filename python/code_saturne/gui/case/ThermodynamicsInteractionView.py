# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
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

# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

import logging

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *
from code_saturne.gui.base.QtPage import ComboModel, from_qvariant, to_text_string, BasicTableModel, DoubleValidator

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------
from code_saturne.model.Common import GuiParam
from code_saturne.model.ThermodynamicsModel import ThermodynamicsInteractionModel
from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel

from code_saturne.gui.case.InterfacialForcesView import InteractionsTableModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

from code_saturne.gui.case.ThermodynamicsInteraction import Ui_ThermodynamicsInteraction

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialForcesView")
log.setLevel(GuiParam.DEBUG)


class ThermodynamicsInteractionView(QWidget, Ui_ThermodynamicsInteraction):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        Ui_ThermodynamicsInteraction.__init__(self)
        self.setupUi(self)

        self.case      = None
        self.zone      = None
        self.zone_name = None
        self.zone_id   = None
        self.model     = None

    def setup(self, case, zone_name):
        self.case = case

        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_name = zone.getLabel()
                self.zone_id   = zone.getCodeNumber()

        self.model = ThermodynamicsInteractionModel(self.case)

        self.groupBoxGeneral.hide()

        allCouples = InterfacialForcesModel(self.case).getAllCouples()
        self.tableModelInteractions = InteractionsTableModel(self, allCouples, self.model)
        self.tableViewInteractions.setModel(self.tableModelInteractions)
        self._initializeInteractionsTable()

        self.modelSurfaceTensionValue = ComboModel(self.comboBoxSurfaceTensionValue, 1, 1)
        self.modelSurfaceTensionValue.addItem(self.tr('constant'), 'constant')
        self.modelSurfaceTensionValue.addItem(self.tr('user law'), 'user_law')
        self.modelSurfaceTensionValue.addItem(self.tr('eos'), 'eos')
        self.modelSurfaceTensionValue.disableItem(str_model="user_law")
        validatorSurfaceTension = DoubleValidator(self.lineEditSurfaceTensionValue, min=0.0)
        self.lineEditSurfaceTensionValue.setValidator(validatorSurfaceTension)

        self.tableModelInteractions.dataChanged.connect(self.dataChanged)
        self.tableViewInteractions.clicked[QModelIndex].connect(self.slotSelectInteraction)
        self.comboBoxSurfaceTensionValue.currentTextChanged[str].connect(self.slotSurfaceTensionType)
        self.lineEditSurfaceTensionValue.textChanged[str].connect(self.slotSurfaceTensionValue)
        self.pushButtonSurfaceTension.clicked.connect(self.slotFormulaSt)

    def refresh(self):
        self._initializeInteractionsTable()
        self.tableViewInteractions.clearSelection()
        self.groupBoxGeneral.hide()

    def _initializeInteractionsTable(self):
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

    def dataChanged(self, topLeft, bottomRight):
        self.tableViewInteractions.resizeColumnsToContents()
        self.tableViewInteractions.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewInteractions.horizontalHeader().setResizeMode(0, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewInteractions.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)

    @pyqtSlot("QModelIndex")
    def slotSelectInteraction(self, index):
        self.groupBoxGeneral.show()

        field_a, field_b, interaction_type = self.tableModelInteractions.data_table[index.row()]
        self.field_id_a = self.model.getFieldId(field_a)
        self.field_id_b = self.model.getFieldId(field_b)

        tag = "surface_tension"
        choice = self.model.getPropertyMode(self.field_id_a, self.field_id_b, tag)
        thermo_a = self.model.getMethod(self.field_id_a)
        thermo_b = self.model.getMethod(self.field_id_b)
        if thermo_a != "user_properties" and thermo_b != "user_properties":
            choice = "eos"
            self.comboBoxSurfaceTensionValue.setEnabled(False)
        else:
            self.modelSurfaceTensionValue.disableItem(str_model="eos")
            if choice == "eos":
                choice = "constant"

        text = self.modelSurfaceTensionValue.dicoM2V[choice]
        self.modelSurfaceTensionValue.setItem(str_model=choice)
        self.comboBoxSurfaceTensionValue.currentTextChanged[str].emit(text)

    @pyqtSlot(str)
    def slotSurfaceTensionType(self, text):
        choice = self.modelSurfaceTensionValue.dicoV2M[text]
        tag = "surface_tension"
        self.model.setPropertyMode(self.field_id_a, self.field_id_b, tag, choice)
        if text == "eos":
            self.lineEditSurfaceTensionValue.setText("")
            self.lineEditSurfaceTensionValue.setEnabled(False)
            self.pushButtonSurfaceTension.setStyleSheet("background-color: None")
            self.pushButtonSurfaceTension.setEnabled(False)
        elif text == "constant":
            surface_tension_value = self.model.getInitialValue(self.field_id_a, self.field_id_b, tag)
            self.lineEditSurfaceTensionValue.setText(str(surface_tension_value))
            self.lineEditSurfaceTensionValue.setEnabled(True)
            self.pushButtonSurfaceTension.setStyleSheet("background-color: None")
            self.pushButtonSurfaceTension.setEnabled(False)
        elif text == "user law":
            self.lineEditSurfaceTensionValue.setText("")
            self.lineEditSurfaceTensionValue.setEnabled(False)
            exp = self.model.getFormula(self.field_id_a, self.field_id_b, tag)
            if exp:
                self.pushButtonSurfaceTension.setStyleSheet("background-color: green")
                self.pushButtonSurfaceTension.setToolTip(exp)
            else:
                self.pushButtonSurfaceTension.setStyleSheet("background-color: red")
            self.pushButtonSurfaceTension.setEnabled(True)
        return

    @pyqtSlot(str)
    def slotSurfaceTensionValue(self, value):
        """
        Update the surface tension
        """
        if self.lineEditSurfaceTensionValue.validator().state == QValidator.Acceptable:
            log.debug("slotSurfaceTensionValue -> %s" % float(value))
            self.model.setInitialValueTens(self.field_id_a, self.field_id_b, float(value))

    @pyqtSlot()
    def slotFormulaSt(self):
        """
        User formula for surface tension
        """
        # TODO getFormulaStComponents should probably be transferred to InterfacialForcesModel
        exp, req, sca, symbols_st = self.model.getFormulaStComponents(self.field_id_a, self.field_id_b, self.zone_id)

        exa = """# water-air at 20°C
                    sigma = 0.075;

            """

        vname = "SurfaceTension"
        dialog = QMegEditorView(parent=self,
                                function_type='vol',
                                zone_name=self.zone_name,
                                variable_name=vname,
                                expression=exp,
                                required=req,
                                symbols=symbols_st,
                                known_fields=sca,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaSt -> %s" % str(result))
            self.model.setFormula(self.field_id_a, self.field_id_b, "surface_tension", result, self.zone_id)
            self.pushButtonSurfaceTension.setStyleSheet("background-color: green")
            self.pushButtonSurfaceTension.setToolTip(exp)
