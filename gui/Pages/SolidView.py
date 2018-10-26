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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from Solid import Ui_Solid
from SolidModel import SolidModel

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

        # Combo box models

        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for fieldId in self.mdl.getSolidFieldIdList() :
            label = self.mdl.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)

        self.currentid = -1

        if len(self.mdl.getSolidFieldIdList()) > 0 :
            self.currentid = self.mdl.getSolidFieldIdList()[0]
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
        self.modelKinetic.addItem(self.tr("uncorrelate collision"), "uncorrelate_collision")
        self.modelKinetic.addItem(self.tr("correlate collision"), "correlate_collision")

        # Validators

        validatorComp = DoubleValidator(self.lineEditCompaction, min = 0.0)
        validatorComp.setExclusiveMin(True)
        self.lineEditCompaction.setValidator(validatorComp)

        # Connect signals to slots
        self.comboBoxField.activated[str].connect(self.slotField)
        self.comboBoxFriction.activated[str].connect(self.slotFriction)
        self.comboBoxGranular.activated[str].connect(self.slotGranular)
        self.comboBoxKinetic.activated[str].connect(self.slotKinetic)
        self.lineEditCompaction.textChanged[str].connect(self.slotCompaction)
        self.checkBoxCoupling.clicked[bool].connect(self.slotCoupling)

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


    def initializeVariables(self, fieldId):
        """
        Initialize variables when a new fieldId is choosen
        """
        value = self.mdl.getCompaction()
        self.lineEditCompaction.setText(str(value))

        model = self.mdl.getFrictionModel(fieldId)
        self.modelFriction.setItem(str_model = model)

        model = self.mdl.getGranularModel(fieldId)
        self.modelGranular.setItem(str_model = model)

        model = self.mdl.getKineticModel(fieldId)
        self.modelKinetic.setItem(str_model = model)

        if self.mdl.getTurbulenceModel(fieldId) == "q2-q12" :
            self.labelCoupling.show()
            self.checkBoxCoupling.show()

            isCoupling = self.mdl.getCouplingStatus(fieldId) == "on"
            self.checkBoxCoupling.setChecked(isCoupling)

        else :
            self.labelCoupling.hide()
            self.checkBoxCoupling.hide()


    @pyqtSlot(bool)
    def slotCoupling(self, checked):
        """
        check box for polydispersed coupling
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setCouplingStatus(self.currentid, status)

