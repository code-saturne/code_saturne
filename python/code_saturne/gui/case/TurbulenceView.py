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
This module defines the turbulence model data management.

This module contains the following classes:
- TurbulenceAdvancedOptionsDialogView
- TurbulenceView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, logging

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
from code_saturne.gui.case.TurbulenceForm import Ui_TurbulenceForm
from code_saturne.model.TurbulenceModel import TurbulenceModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TurbulenceView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Main view class
#-------------------------------------------------------------------------------

class TurbulenceView(QWidget, Ui_TurbulenceForm):
    """
    Class to open Turbulence Page.
    """
    def __init__(self, parent=None, case=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_TurbulenceForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = TurbulenceModel(self.case)

        # Combo model

        self.modelTurbModel = ComboModel(self.comboBoxTurbModel,10,1)

        self.modelTurbModel.addItem(self.tr("No model (i.e. laminar flow)"), "off")

        # RANS - Algebraic
        self.modelTurbModel.addItemGroup(self.tr("RANS - Algebraic"))
        self.modelTurbModel.addItem(self.tr("Mixing length"), "mixing_length", groupName="RANS - Algebraic")

        # RANS - 1st order
        self.modelTurbModel.addItemGroup(self.tr("RANS - 1st order"))

        e = {"k-epsilon-PL": "k-\u03B5 Linear Production",
             "v2f-BL-v2/k": "v\u00B2-f BL-v\u00B2/k",
             "k-omega-SST": "k-\u03C9 SST",
             "Spalart-Allmaras": "Spalart-Allmaras"}

        if QT_API == "PYQT4":
            e["k-epsilon-PL"] = "k-epsilon Linear Production"
            e["v2f-BL-v2/k"] = "v2f BL-v2/k"
            e["k-omega-SST"] = "k-omega SST"

        for k in ("k-epsilon-PL", "v2f-BL-v2/k",
                  "k-omega-SST", "Spalart-Allmaras"):
            self.modelTurbModel.addItem(self.tr(e[k]), k,
                                        groupName="RANS - 1st order")

        # RANS - 2nd order
        self.modelTurbModel.addItemGroup(self.tr("RANS - 2nd order"))

        e = {"Rij-SSG": "R\u1D62\u2C7C-\u03B5 SSG",
             "Rij-EBRSM": "R\u1D62\u2C7C-\u03B5 EBRSM"}

        if QT_API == "PYQT4":
            e["Rij-SSG"] = "Rij-SSG"
            e["Rij-EBRSM"] = "Rij-EBRSM"

        for k in ("Rij-SSG", "Rij-EBRSM"):
            self.modelTurbModel.addItem(self.tr(e[k]), k,
                                        groupName="RANS - 2nd order")

        # LES
        self.modelTurbModel.addItemGroup(self.tr("LES"))

        self.modelTurbModel.addItem(self.tr("Smagorinsky"),
                                    "LES_Smagorinsky",
                                    groupName="LES")
        self.modelTurbModel.addItem(self.tr("Standard dynamic model"),
                                    "LES_dynamique",
                                    groupName="LES")
        self.modelTurbModel.addItem(self.tr("WALE"),
                                    "LES_WALE",
                                    groupName="LES")

        # Others
        self.modelTurbModel.addItemGroup(self.tr(""))
        self.modelTurbModel.addItemGroup(self.tr("Others"))

        e = {"Rij-epsilon": "R\u1D62\u2C7C-\u03B5 LRR",
             "k-epsilon": "k-\u03B5"}

        if QT_API == "PYQT4":
            e["Rij-epsilon"] = "Rij-LRR"
            e["k-epsilon"] = "k-epsilon"

        for k in ("Rij-epsilon", "k-epsilon"):
            self.modelTurbModel.addItem(self.tr(e[k]), k,
                                        groupName="Others")

        self.modelLength = ComboModel(self.comboBoxLength,2,1)
        self.modelLength.addItem(self.tr("Automatic"), 'automatic')
        self.modelLength.addItem(self.tr("Prescribed"), 'prescribed')
        self.comboBoxLength.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        self.labelTurbDiff.hide()
        self.comboBoxTurbDiff.hide()
        self.turbDiff = None

        # Connections

        self.comboBoxTurbModel.activated[str].connect(self.slotTurbulenceModel)
        self.comboBoxWallFunctions.activated[str].connect(self.slotWallFunction)
        self.comboBoxTurbDiff.activated[str].connect(self.slotTurbDiff)
        self.checkBoxGravity.clicked.connect(self.slotGravity)
        self.checkBoxRijCoupled.clicked.connect(self.slotRijCoupled)
        self.lineEditLength.textChanged[str].connect(self.slotLengthScale)

        self.lineEditV0.textChanged[str].connect(self.slotVelocity)
        self.comboBoxLength.activated[str].connect(self.slotLengthChoice)
        self.lineEditL0.textChanged[str].connect(self.slotLength)

        # Frames display

        self.groupBoxAdvanced.hide()
        self.frameLength.hide()

        # Validators

        validator = DoubleValidator(self.lineEditLength, min=0.0)
        validator.setExclusiveMin(True)
        self.lineEditLength.setValidator(validator)

        validatorV0 = DoubleValidator(self.lineEditV0, min=0.0)
        self.lineEditV0.setValidator(validatorV0)

        validatorL0 = DoubleValidator(self.lineEditL0, min=0.0)
        self.lineEditL0.setValidator(validatorL0)

        # Update the turbulence models list with the calculation features

        for turb in self.model.turbulenceModels():
            if turb not in self.model.turbulenceModelsList():
                self.modelTurbModel.disableItem(str_model=turb)

        # Select the turbulence model

        model = self.model.getTurbulenceModel()
        self.modelTurbModel.setItem(str_model=model)
        self.__initializeView()

        # Length scale

        l_scale = self.model.getLengthScale()
        self.lineEditLength.setText(str(l_scale))

        # Initialization

        v = self.model.getVelocity()
        self.lineEditV0.setText(str(v))

        init_length_choice = self.model.getLengthChoice()
        self.modelLength.setItem(str_model=init_length_choice)
        if init_length_choice == 'automatic':
            self.lineEditL0.setText(str())
            self.lineEditL0.hide()
            self.labelUnitL0.hide()
        else:
            self.lineEditL0.show()
            self.labelUnitL0.show()
            l = self.model.getLength()
            self.lineEditL0.setText(str(l))

        self.case.undoStartGlobal()


    def __initializeView(self):
        """
        Private Method.
        initalize view for a turbulence model
        """
        turb_model = self.model.getTurbulenceModel()

        self.groupBoxAdvanced.hide()
        self.frameLength.hide()
        self.groupBoxReferenceValues.hide()

        if turb_model not in ('off', 'mixing_length',
                              'LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            self.groupBoxReferenceValues.show()

        if turb_model == 'mixing_length':
            self.frameLength.show()
            self.model.getLengthScale()

        elif turb_model not in ('off',
                                'LES_Smagorinsky', 'LES_dynamique', 'LES_WALE',
                                'Spalart-Allmaras'):
            self.frameLength.hide()
            self.__init_advanced__(turb_model)
            self.groupBoxAdvanced.show()


    def __init_advanced__(self, turb_model):
        """
        Update advanced options, rebuilding combobox widgets.
        """

        self.labelTurbDiff.hide()
        self.labelWallFunctionsDesc.hide()
        self.comboBoxTurbDiff.hide()
        self.checkBoxRijCoupled.hide()

        if turb_model in ('Rij-epsilon', 'Rij-SSG', 'Rij-EBRSM'):
            turb_diff = self.model.getTurbDiffModel()
            self.labelTurbDiff.show()
            self.comboBoxTurbDiff.show()
            self.turbDiff = ComboModel(self.comboBoxTurbDiff, 2, 1)
            self.turbDiff.addItem(self.tr("Scalar diffusivity (Shir model)"), 'shir')
            self.turbDiff.addItem(self.tr("Tensorial diffusivity (Daly and Harlow model)"), 'daly_harlow')
            self.turbDiff.setItem(str_model=str(turb_diff))

            self.checkBoxRijCoupled.show()

        self.checkBoxGravity.setEnabled(True)
        self.comboBoxWallFunctions.setEnabled(True)

        wall_f_types, wall_f_default = self.model.wall_function_types()
        self.default_desc = wall_f_types[wall_f_default]

        # Initialization of wall function model

        n_wfo = len(wall_f_types)

        self.wallFunctions = ComboModel(self.comboBoxWallFunctions, n_wfo, 1)

        for o in wall_f_types:
            descr = wall_f_types[o]
            self.wallFunctions.addItem(self.tr(descr), o)

        wall_function = self.model.getWallFunction()
        self.wallFunctions.setItem(str_model=str(wall_function))

        if wall_function == -1:
            self.labelWallFunctionsDesc.setText(self.default_desc)
            self.labelWallFunctionsDesc.show()

        # Initialization of gravity terms
        if self.model.getGravity() == 'on':
            self.checkBoxGravity.setChecked(True)
        else:
            self.checkBoxGravity.setChecked(False)

        # Initialization of coupled component option
        if self.model.getRijCoupled() == 'on':
            self.checkBoxRijCoupled.setChecked(True)
        else:
            self.checkBoxRijCoupled.setChecked(False)


    @pyqtSlot(str)
    def slotLengthScale(self, text):
        """
        Private slot.
        Input XLOMLG.
        """
        if self.lineEditLength.validator().state == QValidator.Acceptable:
            l_scale = from_qvariant(text, float)
            self.model.setLengthScale(l_scale)


    @pyqtSlot(str)
    def slotTurbulenceModel(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.modelTurbModel.dicoV2M[str(text)]
        self.model.setTurbulenceModel(model)
        self.__initializeView()


    @pyqtSlot(str)
    def slotWallFunction(self, text):
        """
        Private slot.
        Input iwallf.
        """
        m = self.wallFunctions.dicoV2M[str(self.comboBoxWallFunctions.currentText())]
        self.model.setWallFunction(m)
        if m == '-1':
            self.labelWallFunctionsDesc.show()
        else:
            self.labelWallFunctionsDesc.hide()


    @pyqtSlot(str)
    def slotTurbDiff(self, text):
        """
        Private slot.
        Input iwallf.
        """
        m = self.turbDiff.dicoV2M[str(self.comboBoxTurbDiff.currentText())]
        self.model.setTurbDiffModel(m)


    @pyqtSlot(str)
    def slotVelocity(self,  text):
        """
        Private slot.
        Input reference velocity.
        """
        if self.lineEditV0.validator().state == QValidator.Acceptable:
            v = from_qvariant(text, float)
            self.model.setVelocity(v)


    @pyqtSlot(str)
    def slotLengthChoice(self,text):
        """
        Private slot.
        Input mode for reference length.
        """
        choice = self.modelLength.dicoV2M[str(text)]
        self.model.setLengthChoice(choice)
        if choice == 'automatic':
            self.lineEditL0.setText(str())
            self.lineEditL0.hide()
            self.labelUnitL0.hide()
        else:
            self.lineEditL0.show()
            self.labelUnitL0.show()
            value = self.model.getLength()
            self.lineEditL0.setText(str(value))
        log.debug("slotlengthchoice-> %s" % choice)


    @pyqtSlot(str)
    def slotLength(self,  text):
        """
        Private slot.
        Input reference length.
        """
        if self.lineEditL0.validator().state == QValidator.Acceptable:
            l = from_qvariant(text, float)
            self.model.setLength(l)


    @pyqtSlot()
    def slotGravity(self):
        """
        Activate or deactivate gravity source terms
        """
        if self.checkBoxGravity.isChecked():
            self.model.setGravity('on')
        else:
            self.model.setGravity('off')


    @pyqtSlot()
    def slotRijCoupled(self):
        """
        Activate or deactivate Rij component coupling
        """
        if self.checkBoxRijCoupled.isChecked():
            self.model.setRijCoupled('on')
        else:
            self.model.setRijCoupled('off')


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
