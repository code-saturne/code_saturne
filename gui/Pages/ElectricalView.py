# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
This module contains the following classes:
- ElectricalView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, RegExpValidator
from code_saturne.Base.QtPage import from_qvariant

from code_saturne.Pages.ElectricalForm import Ui_ElectricalForm
from code_saturne.Pages.ElectricalModel import ElectricalModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ElectricalView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ElectricalView(QWidget, Ui_ElectricalForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ElectricalForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.stbar = stbar
        self.case.undoStopGlobal()

        self.model = ElectricalModel(self.case)

        # Combo model
        self.modelJoule = ComboModel(self.comboBoxJouleModel, 4, 1)
        self.modelJoule.addItem(self.tr("AC/DC"), "AC/DC")
        self.modelJoule.addItem(self.tr("three-phase"), "three-phase")
        self.modelJoule.addItem(self.tr("AC/DC with Transformer coupling"), "AC/DC+Transformer")
        self.modelJoule.addItem(self.tr("three-phase with Transformer coupling"), "three-phase+Transformer")
        self.modelJoule.disableItem(str_model="AC/DC+Transformer")
        self.modelJoule.disableItem(str_model="three-phase+Transformer")

        self.modelScaling = ComboModel(self.comboBoxScalingModel, 3, 1)
        self.modelScaling.addItem(self.tr("general case"), "general_case")
        self.modelScaling.addItem(self.tr("plane define"), "plane_define")
        self.modelScaling.addItem(self.tr("user define"), "user")

        self.modelDirection = ComboModel(self.comboBoxDirection, 3, 1)
        self.modelDirection.addItem(self.tr("X"), "X")
        self.modelDirection.addItem(self.tr("Y"), "Y")
        self.modelDirection.addItem(self.tr("Z"), "Z")

        # Connections
        self.lineEditPower.textChanged[str].connect(self.slotPower)
        self.lineEditCurrent.textChanged[str].connect(self.slotCurrent)
        self.checkBoxScaling.clicked.connect(self.slotScaling)
        self.comboBoxJouleModel.activated[str].connect(self.slotJouleModel)
        self.comboBoxScalingModel.activated[str].connect(self.slotScalingModel)
        self.comboBoxDirection.activated[str].connect(self.slotDirection)
        self.lineEditPlaneDefinitionA.textChanged[str].connect(self.slotPlaneDefA)
        self.lineEditPlaneDefinitionB.textChanged[str].connect(self.slotPlaneDefB)
        self.lineEditPlaneDefinitionC.textChanged[str].connect(self.slotPlaneDefC)
        self.lineEditPlaneDefinitionD.textChanged[str].connect(self.slotPlaneDefD)
        self.lineEditEpsilon.textChanged[str].connect(self.slotPlaneDefEpsilon)

        # Validators
        validatorPower = DoubleValidator(self.lineEditPower, min=0.0)
        validatorPower.setExclusiveMin(False)
        validatorCurrent = DoubleValidator(self.lineEditCurrent, min=0.0)
        validatorCurrent.setExclusiveMin(False)
        validatorDefinitionA = DoubleValidator(self.lineEditPlaneDefinitionA)
        validatorDefinitionB = DoubleValidator(self.lineEditPlaneDefinitionB)
        validatorDefinitionC = DoubleValidator(self.lineEditPlaneDefinitionC)
        validatorDefinitionD = DoubleValidator(self.lineEditPlaneDefinitionD)
        validatorEpsilon     = DoubleValidator(self.lineEditEpsilon)
        self.lineEditPower.setValidator(validatorPower)
        self.lineEditCurrent.setValidator(validatorCurrent)
        self.lineEditPlaneDefinitionA.setValidator(validatorDefinitionA)
        self.lineEditPlaneDefinitionB.setValidator(validatorDefinitionB)
        self.lineEditPlaneDefinitionC.setValidator(validatorDefinitionC)
        self.lineEditPlaneDefinitionD.setValidator(validatorDefinitionD)
        self.lineEditEpsilon.setValidator(validatorEpsilon)

        # Initialize widget
        self.__initializeWidget()

        self.case.undoStartGlobal()


    @pyqtSlot()
    def __initializeWidget(self):
        """
        Initialize widget
        """
        self.groupBoxRecalage.hide()

        if self.model.getScaling() == 'on':
            self.checkBoxScaling.setChecked(True)
            self.labelScalingModel.show()
            self.comboBoxScalingModel.show()
        else:
            self.checkBoxScaling.setChecked(False)
            self.labelScalingModel.hide()
            self.comboBoxScalingModel.hide()

        if self.model.getElectricalModel() == "joule":
            self.groupBoxJoule.show()
            self.groupBoxElectricArc.hide()

            model = self.model.getJouleModel()
            self.modelJoule.setItem(str_model=str(model))
            power = self.model.getPower()
            self.lineEditPower.setText(str(power))

        elif self.model.getElectricalModel() == "arc":
            self.groupBoxJoule.hide()
            self.groupBoxElectricArc.show()

            current = self.model.getCurrent()
            self.lineEditCurrent.setText(str(current))

            if self.model.getScaling() == 'on':
                model = self.model.getScalingModel()
                self.modelScaling.setItem(str_model=str(model))
                if model == 'plane_define':
                    self.groupBoxRecalage.show()
                    direction = self.model.getDirection()
                    self.modelDirection.setItem(str_model=str(direction))
                    definition = self.model.getPlaneDefinition("A")
                    self.lineEditPlaneDefinitionA.setText(str(definition))
                    definition = self.model.getPlaneDefinition("B")
                    self.lineEditPlaneDefinitionB.setText(str(definition))
                    definition = self.model.getPlaneDefinition("C")
                    self.lineEditPlaneDefinitionC.setText(str(definition))
                    definition = self.model.getPlaneDefinition("D")
                    self.lineEditPlaneDefinitionD.setText(str(definition))
                    definition = self.model.getPlaneDefinition("epsilon")
                    self.lineEditEpsilon.setText(str(definition))


    @pyqtSlot(str)
    def slotPower(self, text):
        """
        Input Imposed Power
        """
        if self.lineEditPower.validator().state == QValidator.Acceptable:
            power = from_qvariant(text, float)
            self.model.setPower(power)


    @pyqtSlot(str)
    def slotCurrent(self, text):
        """
        Input Imposed current intensity
        """
        if self.lineEditCurrent.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setCurrent(current)


    @pyqtSlot(str)
    def slotJouleModel(self, text):
        """
        Input Joule model.
        """
        model = self.modelJoule.dicoV2M[str(text)]
        self.model.setJouleModel(model)


    @pyqtSlot()
    def slotScaling(self):
        """
        Input "Electric variables" scaling.
        """
        if self.checkBoxScaling.isChecked():
            self.model.setScaling("on")
        else:
            self.model.setScaling("off")

        self.__initializeWidget()


    @pyqtSlot(str)
    def slotScalingModel(self, text):
        """
        Input scaling model.
        """
        model = self.modelScaling.dicoV2M[str(text)]
        self.model.setScalingModel(model)
        self.__initializeWidget()


    @pyqtSlot(str)
    def slotDirection(self, text):
        """
        Input current density direction for scaling.
        """
        direction = self.modelDirection.dicoV2M[str(text)]
        self.model.setDirection(direction)


    @pyqtSlot(str)
    def slotPlaneDefA(self, text):
        """
        Input define plane
        """
        if self.lineEditPlaneDefinitionA.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setPlaneDefinition("A", current)


    @pyqtSlot(str)
    def slotPlaneDefB(self, text):
        """
        Input define plane
        """
        if self.lineEditPlaneDefinitionB.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setPlaneDefinition("B", current)


    @pyqtSlot(str)
    def slotPlaneDefC(self, text):
        """
        Input define plane
        """
        if self.lineEditPlaneDefinitionC.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setPlaneDefinition("C", current)


    @pyqtSlot(str)
    def slotPlaneDefD(self, text):
        """
        Input define plane
        """
        if self.lineEditPlaneDefinitionD.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setPlaneDefinition("D", current)


    @pyqtSlot(str)
    def slotPlaneDefEpsilon(self, text):
        """
        Input define plane
        """
        if self.lineEditEpsilon.validator().state == QValidator.Acceptable:
            current = from_qvariant(text, float)
            self.model.setPlaneDefinition("epsilon", current)


    def tr(self, text):
        """
        Translation.
        """
        return text


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
