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

"""
This module defines the 'Interfacial area' page.

This module contains the following classes:
- InterfacialAreaView
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
from InterfacialArea import Ui_InterfacialArea
from InterfacialAreaModel import InterfacialAreaModel
from code_saturne.Pages.MainFieldsModel import MainFieldsModel
from code_saturne.Pages.InterfacialForcesModel import InterfacialForcesModel
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialAreaView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# SolidView class
#-------------------------------------------------------------------------------

class InterfacialAreaView(QWidget, Ui_InterfacialArea):
    """
    InterfacialAreaView layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InterfacialArea.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = InterfacialAreaModel(self.case)

        # Combo box models
        id_to_set = -1
        self.modelField = ComboModel(self.comboBoxField, 1, 1)

        # For consistency with the previous pages, the second phase of the
        # Large Interface Model is set before the dispersed fields
        if InterfacialForcesModel(self.case).getBubblesForLIMStatus() == 'on':
            id_to_set = 2
            label = self.mdl.getLabel(id_to_set)
            name = str(id_to_set)
            self.modelField.addItem(self.tr(label), name)
            self.modelField.setItem(str_model = name)

        for fieldId in self.mdl.getDispersedFieldList() :
            label = self.mdl.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)

        if len(self.mdl.getDispersedFieldList()) > 0 and id_to_set == -1:
            id_to_set = self.mdl.getDispersedFieldList()[0]
            self.modelField.setItem(str_model = id_to_set)

        # case no field
        self.currentid = id_to_set

        self.modelModel = ComboModel(self.comboBoxModel, 2, 1)
        self.modelModel.addItem(self.tr("constant"),"constant")
        self.modelModel.addItem(self.tr("interfacial area transport"),"interfacial_area_transport")

        self.modelSourceTerm = ComboModel(self.comboBoxSourceTerm, 4, 1)

        self.modelSourceTerm.addItem(self.tr("No coalescence, no fragmentation"),"no_coalescence_no_fragmentation")
        self.modelSourceTerm.addItem(self.tr("Yao & Morel"),"wei_yao")
        self.modelSourceTerm.addItem(self.tr("Kamp & Colin"),"kamp_colin")
        self.modelSourceTerm.addItem(self.tr("Ruyer & Seiler"),"ruyer_seiler")
        self.modelSourceTerm.disableItem(2)

        self.modelSolutionMethod = ComboModel(self.comboBoxSolutionMethod, 2, 1)
        self.modelSolutionMethod.addItem(self.tr("Uncoupled (after navsto)"),"uncoupled")
        self.modelSolutionMethod.addItem(self.tr("Coupled (alpha-P-H cycle)"),"coupled")

        # Validators
        validatorDefDiam = DoubleValidator(self.lineEditDefaultDiameter, min = 0.0)
        validatorMinDiam = DoubleValidator(self.lineEditMinDiameter, min = 0.0)
        validatorMaxDiam = DoubleValidator(self.lineEditMaxDiameter, min = 0.0)

        validatorDefDiam.setExclusiveMin(True)
        validatorMinDiam.setExclusiveMin(True)
        validatorMaxDiam.setExclusiveMin(True)

        self.lineEditDefaultDiameter.setValidator(validatorDefDiam)
        self.lineEditMinDiameter.setValidator(validatorMinDiam)
        self.lineEditMaxDiameter.setValidator(validatorMaxDiam)

        # Connect signals to slots
        self.comboBoxField.activated[str].connect(self.slotField)
        self.comboBoxModel.activated[str].connect(self.slotModel)
        self.comboBoxSourceTerm.activated[str].connect(self.slotSourceTerm)
        self.comboBoxSolutionMethod.activated[str].connect(self.slotSolutionMethod)
        self.lineEditDefaultDiameter.textChanged[str].connect(self.slotDefaultDiameter)
        self.lineEditMinDiameter.textChanged[str].connect(self.slotMinDiameter)
        self.lineEditMaxDiameter.textChanged[str].connect(self.slotMaxDiameter)

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

        if self.mdl.getFieldNature(self.currentid) == "gas" :
            self.modelSourceTerm.enableItem(0)
        else :
            self.modelSourceTerm.disableItem(0)


    @pyqtSlot(str)
    def slotModel(self, text):
        """
        INPUT type for choice of model
        """
        model = self.modelModel.dicoV2M[text]
        self.mdl.setAreaModel(self.currentid, model)
        self.initializeVariables(self.currentid)


    @pyqtSlot(str)
    def slotSourceTerm(self, text):
        """
        INPUT type for choice of model source term
        """
        model = self.modelSourceTerm.dicoV2M[text]
        self.mdl.setSourceTerm(self.currentid, model)


    @pyqtSlot(str)
    def slotSolutionMethod(self, text):
        """
        INPUT type for choice of solution method
        """
        model = self.modelSolutionMethod.dicoV2M[text]
        self.mdl.setSolutionMethod(self.currentid, model)


    @pyqtSlot(str)
    def slotDefaultDiameter(self, var):
        """
        """
        if self.lineEditDefaultDiameter.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setInitialDiameter(self.currentid, value)


    @pyqtSlot(str)
    def slotMinDiameter(self, var):
        """
        """
        if self.lineEditMinDiameter.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMinDiameter(value)


    @pyqtSlot(str)
    def slotMaxDiameter(self, var):
        """
        """
        if self.lineEditMaxDiameter.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.mdl.setMaxDiameter(value)


    def initializeVariables(self, fieldId):
        """
        Initialize variables when a new fieldId is choosen
        """
        model = self.mdl.getAreaModel(fieldId)
        self.modelModel.setItem(str_model = model)

        value = self.mdl.getInitialDiameter(self.currentid)
        self.lineEditDefaultDiameter.setText(str(value))

        if self.mdl.getAreaModel(fieldId) == "constant" :
            self.groupBoxAreaTransport.hide()
            self.groupBoxMinMaxDiameter.hide()
        else :
            self.groupBoxAreaTransport.show()
            model = self.mdl.getSourceTerm(fieldId)
            self.modelSourceTerm.setItem(str_model = model)

            model = self.mdl.getSolutionMethod(fieldId)
            self.modelSolutionMethod.setItem(str_model = model)

            self.groupBoxMinMaxDiameter.show()

            value = self.mdl.getMinDiameter()
            self.lineEditMinDiameter.setText(str(value))

            value = self.mdl.getMaxDiameter()
            self.lineEditMaxDiameter.setText(str(value))

            if MainFieldsModel(self.case).getFieldNature(fieldId) != 'gas' :
                self.modelSourceTerm.disableItem(1)
                self.modelSourceTerm.disableItem(2)
                self.modelSourceTerm.disableItem(3)
            else :
                self.modelSourceTerm.enableItem(1)
                self.modelSourceTerm.enableItem(2)
                self.modelSourceTerm.enableItem(3)

