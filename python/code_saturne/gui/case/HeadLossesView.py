# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

"""
This module defines the HeadLosses model data management.

This module contains the following classes:
- HeadLosses
- HeadLossesView
"""

# -------------------------------------------------------------------------------
# Library modules import
# -------------------------------------------------------------------------------

import sys, logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.HeadLossesForm import Ui_HeadLossesForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.HeadLossesModel import HeadLossesModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("HeadLossesView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main view class
# -------------------------------------------------------------------------------

class HeadLossesView(QWidget, Ui_HeadLossesForm):
    # TODO : unify redundant methods slotKxx, ... slotKzz, slotA11, ..., slotA33

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_HeadLossesForm.__init__(self)
        self.setupUi(self)

        self.case = None
        self.zone = None
        self.model = None

    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()
        localization_model = LocalizationModel("VolumicZone", self.case)
        for zone in localization_model.getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
        if self.zone.isNatureActivated("head_losses"):
            self.model = HeadLossesModel(self.case)
            self.selectHeadLossesZones()
            self.setConnections()
            self.setValidators()
        else:  # TODO check if content of tab should remain visible
            self.displayDefaultView()
            self.setEnabled(False)
        self.case.undoStartGlobal()

    def setValidators(self):
        validator_kxx = DoubleValidator(self.lineEdit_kxx, min=0.0)
        validator_kyy = DoubleValidator(self.lineEdit_kyy, min=0.0)
        validator_kzz = DoubleValidator(self.lineEdit_kzz, min=0.0)
        validator_a11 = DoubleValidator(self.lineEdit_a11)
        validator_a12 = DoubleValidator(self.lineEdit_a12)
        validator_a13 = DoubleValidator(self.lineEdit_a13)
        validator_a21 = DoubleValidator(self.lineEdit_a21)
        validator_a22 = DoubleValidator(self.lineEdit_a22)
        validator_a23 = DoubleValidator(self.lineEdit_a23)
        validator_a31 = DoubleValidator(self.lineEdit_a31)
        validator_a32 = DoubleValidator(self.lineEdit_a32)
        validator_a33 = DoubleValidator(self.lineEdit_a33)
        # Apply validators
        self.lineEdit_kxx.setValidator(validator_kxx)
        self.lineEdit_kyy.setValidator(validator_kyy)
        self.lineEdit_kzz.setValidator(validator_kzz)
        self.lineEdit_a11.setValidator(validator_a11)
        self.lineEdit_a12.setValidator(validator_a12)
        self.lineEdit_a13.setValidator(validator_a13)
        self.lineEdit_a21.setValidator(validator_a21)
        self.lineEdit_a22.setValidator(validator_a22)
        self.lineEdit_a23.setValidator(validator_a23)
        self.lineEdit_a31.setValidator(validator_a31)
        self.lineEdit_a32.setValidator(validator_a32)
        self.lineEdit_a33.setValidator(validator_a33)

    def setConnections(self):
        self.groupBox_3.clicked[bool].connect(self.slotTransfoMatrix)
        self.lineEdit_kxx.textChanged[str].connect(self.slotKxx)
        self.lineEdit_kyy.textChanged[str].connect(self.slotKyy)
        self.lineEdit_kzz.textChanged[str].connect(self.slotKzz)
        self.lineEdit_a11.textChanged[str].connect(self.slotA11)
        self.lineEdit_a12.textChanged[str].connect(self.slotA12)
        self.lineEdit_a13.textChanged[str].connect(self.slotA13)
        self.lineEdit_a21.textChanged[str].connect(self.slotA21)
        self.lineEdit_a22.textChanged[str].connect(self.slotA22)
        self.lineEdit_a23.textChanged[str].connect(self.slotA23)
        self.lineEdit_a31.textChanged[str].connect(self.slotA31)
        self.lineEdit_a32.textChanged[str].connect(self.slotA32)
        self.lineEdit_a33.textChanged[str].connect(self.slotA33)

    def selectHeadLossesZones(self):
        label = self.zone.getLabel()
        name = self.zone.getCodeNumber()

        if hasattr(self, "modelScalars"): del self.modelScalars
        log.debug("slotSelectHeadLossesZones label %s " % label)
        self.groupBoxDef.show()
        self.groupBox_3.show()
        kxx, kyy, kzz = self.model.getKCoefficients(name)
        self.lineEdit_kxx.setText(str(kxx))
        self.lineEdit_kyy.setText(str(kyy))
        self.lineEdit_kzz.setText(str(kzz))

        if self.model.getMatrixChoice(name, 'choice') == 'on':
            self.groupBox_3.setChecked(True)
            checked = True
        else:
            self.groupBox_3.setChecked(False)
            checked = False

        self.slotTransfoMatrix(checked)

    def displayDefaultView(self):
        self.groupBoxDef.show()
        self.groupBox_3.setChecked(False)
        self.groupBox_3.hide()

    @pyqtSlot(str)
    def slotKxx(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_kxx.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setCoefficient(zone_id, 'kxx', value)

    @pyqtSlot(str)
    def slotKyy(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_kyy.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setCoefficient(zone_id, 'kyy', value)

    @pyqtSlot(str)
    def slotKzz(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_kzz.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setCoefficient(zone_id, 'kzz', value)

    @pyqtSlot(bool)
    def slotTransfoMatrix(self, checked):
        self.groupBox_3.setFlat(not checked)
        zone_id = self.zone.getCodeNumber()

        if checked:
            self.model.setMatrixChoice(zone_id, 'choice', 'on')
            self.groupBox_3.setChecked(True)
            self.frameTransfo.show()
            a11, a12, a13, a21, a22, a23, a31, a32, a33 = self.model.getMatrix(zone_id)
        else:
            self.model.setMatrixChoice(zone_id, 'choice', 'off')
            self.frameTransfo.hide()
            a11, a12, a13, a21, a22, a23, a31, a32, a33 = 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0
            self.groupBox_3.setChecked(False)

        self.lineEdit_a11.setText(str(a11))
        self.lineEdit_a12.setText(str(a12))
        self.lineEdit_a13.setText(str(a13))
        self.lineEdit_a21.setText(str(a21))
        self.lineEdit_a22.setText(str(a22))
        self.lineEdit_a23.setText(str(a23))
        self.lineEdit_a31.setText(str(a31))
        self.lineEdit_a32.setText(str(a32))
        self.lineEdit_a33.setText(str(a33))

    @pyqtSlot(str)
    def slotA11(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a11.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a11', value)

    @pyqtSlot(str)
    def slotA12(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a12.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a12', value)

    @pyqtSlot(str)
    def slotA13(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a13.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a13', value)

    @pyqtSlot(str)
    def slotA21(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a21.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a21', value)

    @pyqtSlot(str)
    def slotA22(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a22.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a22', value)

    @pyqtSlot(str)
    def slotA23(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a23.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a23', value)

    @pyqtSlot(str)
    def slotA31(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a31.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a31', value)

    @pyqtSlot(str)
    def slotA32(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a32.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a32', value)

    @pyqtSlot(str)
    def slotA33(self, text):
        zone_id = self.zone.getCodeNumber()
        model = HeadLossesModel(self.case)
        if self.lineEdit_a33.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            model.setMatrixComponent(zone_id, 'a33', value)


# -------------------------------------------------------------------------------
# Testing part
# -------------------------------------------------------------------------------


if __name__ == "__main__":
    pass

# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
