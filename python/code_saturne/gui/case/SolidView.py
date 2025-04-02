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

"""
This module defines the 'Particles interactions' page.

This module contains the following classes:
- SolidView
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator, from_qvariant
from Solid import Ui_Solid
from code_saturne.model.SolidModel import SolidModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SolidView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# SolidView class
#-------------------------------------------------------------------------------

class SolidView(QWidget, Ui_Solid):
    """
    Particles interaction layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_Solid.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = SolidModel(self.case)

        if self.mdl.mainFieldsModel.getSolidPhaseList() == []:
            self.groupBoxInteractions.hide()
            self.groupBoxGeneral.hide()
            self.labelNoParticles.show()
            return

        # Combo box models

        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for field in self.mdl.mainFieldsModel.getSolidPhaseList():
            label = field.label
            name = field.f_id
            self.modelField.addItem(self.tr(label), name)

        self.currentid = -1

        if len(self.mdl.mainFieldsModel.getSolidPhaseList()) > 0 :
            self.currentid = self.mdl.mainFieldsModel.getSolidPhaseList()[0].f_id
            self.modelField.setItem(str_model = self.currentid)

        self.modelFriction = ComboModel(self.comboBoxFriction, 3, 1)
        self.modelFriction.addItem(self.tr("none"), "none")
        self.modelFriction.addItem(self.tr("pressure"), "pressure")
        self.modelFriction.addItem(self.tr("fluxes"), "fluxes")

        self.modelGranular = ComboModel(self.comboBoxGranular, 3, 1)
        self.modelGranular.addItem(self.tr("none"), "none")
        self.modelGranular.addItem(self.tr("pressure"), "pressure")
        self.modelGranular.addItem(self.tr("fluxes"), "fluxes")

        self.modelKinetic = ComboModel(self.comboBoxKinetic, 3, 1)
        self.modelKinetic.addItem(self.tr("none"), "none")
        self.modelKinetic.addItem(self.tr("uncorrelated collision"), "uncorrelate_collision")
        self.modelKinetic.addItem(self.tr("correlated collision"), "correlate_collision")

        # Validators

        validatorComp = DoubleValidator(self.lineEditCompaction, min = 0.0)
        validatorComp.setExclusiveMin(True)
        self.lineEditCompaction.setValidator(validatorComp)

        validatorFricThres = DoubleValidator(self.lineEditFrictonalThres, min = 0.0)
        validatorFricThres.setExclusiveMin(True)
        self.lineEditFrictonalThres.setValidator(validatorFricThres)

        validatorElast = DoubleValidator(self.lineEditElastCoef, min = 0.0)
        validatorElast.setExclusiveMin(False)
        self.lineEditElastCoef.setValidator(validatorElast)

        # Connect signals to slots
        self.comboBoxField.activated[str].connect(self.slotField)
        self.comboBoxFriction.activated[str].connect(self.slotFriction)
        self.comboBoxGranular.activated[str].connect(self.slotGranular)
        self.comboBoxKinetic.activated[str].connect(self.slotKinetic)
        self.lineEditCompaction.textChanged[str].connect(self.slotCompaction)
        self.lineEditFrictonalThres.textChanged[str].connect(self.slotFrictionalThreshold)
        self.lineEditElastCoef.textChanged[str].connect(self.slotSetElasticity)
        self.checkBoxCoupling.clicked[bool].connect(self.slotCoupling)

        # Show / hide polydispersed parameter
        has_qp_qfp = False
        for field in self.mdl.mainFieldsModel.getSolidPhaseList():
            if self.mdl.getTurbulenceModel(field.f_id) == "q2-q12":
                has_qp_qfp = True

        self.labelCoupling.setVisible(has_qp_qfp)
        self.checkBoxCoupling.setVisible(has_qp_qfp)
        if has_qp_qfp:
            isCoupling = self.mdl.getCouplingStatus() == "on"
            self.checkBoxCoupling.setChecked(isCoupling)


        # Initialize widget
        self.initializeVariables(self.currentid)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotField(self, text):
        """
        INPUT label for choice of field
        """
        self.currentid = self.modelField.dicoV2M[text]
        self.initializeVariables(self.currentid)


    @pyqtSlot(str)
    def slotFriction(self, text):
        """
        INPUT type for choice of friction model
        """
        model = self.modelFriction.dicoV2M[text]
        self.mdl.setFrictionModel(self.currentid, model)


    @pyqtSlot(str)
    def slotGranular(self, text):
        """
        INPUT type for choice of granular model
        """
        model = self.modelGranular.dicoV2M[text]
        self.mdl.setGranularModel(self.currentid, model)


    @pyqtSlot(str)
    def slotKinetic(self, text):
        """
        INPUT type for choice of kinetic model
        """
        model = self.modelKinetic.dicoV2M[text]
        self.mdl.setKineticModel(self.currentid, model)


    @pyqtSlot(str)
    def slotCompaction(self, var):
        """
        """
        if self.lineEditCompaction.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setCompaction(value)


    @pyqtSlot(str)
    def slotFrictionalThreshold(self, var):
        """
        Setter slot.
        """
        if self.lineEditFrictonalThres.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMinFrictionalThreshold(value)


    @pyqtSlot(str)
    def slotSetElasticity(self, var):
        """
        Set elasiticity coefficient
        """
        if self.lineEditElastCoef.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setElastCoeff(value, self.currentid)


    def initializeVariables(self, fieldId):
        """
        Initialize variables when a new fieldId is choosen
        """
        self.labelNoParticles.hide()

        value = self.mdl.getCompaction()
        self.lineEditCompaction.setText(str(value))

        value_min_fric = self.mdl.getMinFrictionalThreshold()
        self.lineEditFrictonalThres.setText(str(value_min_fric))

        self.lineEditElastCoef.setText(str(self.mdl.getElastCoeff(fieldId)))

        model = self.mdl.getFrictionModel(fieldId)
        self.modelFriction.setItem(str_model=model)

        model = self.mdl.getGranularModel(fieldId)
        self.modelGranular.setItem(str_model=model)

        model = self.mdl.getKineticModel(fieldId)
        self.modelKinetic.setItem(str_model = model)


    @pyqtSlot(bool)
    def slotCoupling(self, checked):
        """
        check box for polydispersed coupling
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setCouplingStatus(status)

