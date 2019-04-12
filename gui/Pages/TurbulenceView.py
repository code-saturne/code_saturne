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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from code_saturne.Pages.TurbulenceForm import Ui_TurbulenceForm
from code_saturne.Pages.TurbulenceAdvancedOptionsDialogForm import Ui_TurbulenceAdvancedOptionsDialogForm
from code_saturne.model.TurbulenceModel import TurbulenceModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TurbulenceView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Advanced dialog
#-------------------------------------------------------------------------------

class TurbulenceAdvancedOptionsDialogView(QDialog, Ui_TurbulenceAdvancedOptionsDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_TurbulenceAdvancedOptionsDialogForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        if default['model'] in ('k-epsilon', 'k-epsilon-PL'):
            title = self.tr("Options for k-epsilon model")
        elif default['model'] in ('Launder-Sharma'):
            title = self.tr("Options for Launder-Sharma k-epsilon model")
        elif default['model'] in ('Rij-epsilon', 'Rij-SSG', 'Rij-EBRSM'):
            title = self.tr("Options for Rij-epsilon model")
        elif default['model'] == 'k-omega-SST':
            title = self.tr("Options for k-omega-SST model")
        elif default['model'] == 'v2f-BL-v2/k':
            title = self.tr("Options for v2f-BL-v2/k model")
        elif default['model'] == 'Spalart-Allmaras':
            title = self.tr("Options for Spalart-Allmaras model")
        self.setWindowTitle(title)
        self.default = default
        self.result  = self.default.copy()

        self.checkBoxGravity.setEnabled(True)
        self.comboBoxWallFunctions.setEnabled(True)

        if default['model'] == 'Rij-EBRSM':
            # Combo - piecewise laws (iwallf=2,3) unavailable through the GUI
            self.wallFunctions = ComboModel(self.comboBoxWallFunctions, 2, 1)
            self.wallFunctions.addItem(self.tr("No wall function"), '0')
            self.wallFunctions.addItem(self.tr("2-scale model (all y+)"), '7')

            # Initialization of wall function model
            self.wallFunctions.setItem(str_model=str(self.result['wall_function']))

        elif default['model'] == 'v2f-BL-v2/k':
            self.wallFunctions = ComboModel(self.comboBoxWallFunctions, 1, 1)
            self.wallFunctions.addItem(self.tr("No wall function"), '0')
            self.comboBoxWallFunctions.setEnabled(False)
        elif default['model'] == 'Spalart-Allmaras':
            self.wallFunctions = ComboModel(self.comboBoxWallFunctions, 1, 1)
            self.wallFunctions.addItem(self.tr("One scale model (log law)"), '2')
            self.comboBoxWallFunctions.setEnabled(False)
        elif default['model'] == 'k-omega-SST':
            # Combo - power law (iwallf=1) unavailable through the GUI
            self.wallFunctions = ComboModel(self.comboBoxWallFunctions, 5, 1)
            self.wallFunctions.addItem(self.tr("No wall function"), '0')
            self.wallFunctions.addItem(self.tr("1-scale model (log law)"), '2')
            self.wallFunctions.addItem(self.tr("2-scale model (log law)"), '3')
            self.wallFunctions.addItem(self.tr("2-scale model (all y+)"), '7')
            self.wallFunctions.addItem(self.tr("Scalable 2-scale model (log law)"), '4')

            # Initialization of wall function model
            self.wallFunctions.setItem(str_model=str(self.result['wall_function']))
        else:
            # Combo - power law (iwallf=1) unavailable through the GUI
            self.wallFunctions = ComboModel(self.comboBoxWallFunctions, 4, 1)
            self.wallFunctions.addItem(self.tr("No wall function"), '0')
            self.wallFunctions.addItem(self.tr("1-scale model (log law)"), '2')
            self.wallFunctions.addItem(self.tr("2-scale model (log law)"), '3')
            self.wallFunctions.addItem(self.tr("Scalable 2-scale model (log law)"), '4')

            # Initialization of wall function model
            self.wallFunctions.setItem(str_model=str(self.result['wall_function']))

        # Initialization gravity terms
        if self.result['gravity_terms'] == 'on':
            self.checkBoxGravity.setChecked(True)
        else:
            self.checkBoxGravity.setChecked(False)

        self.case.undoStartGlobal()


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        if self.checkBoxGravity.isChecked():
            self.result['gravity_terms'] = "on"
        else:
            self.result['gravity_terms'] = "off"
        self.result['wall_function'] = \
        int(self.wallFunctions.dicoV2M[str(self.comboBoxWallFunctions.currentText())])

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text

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
        self.modelTurbModel.addItem(self.tr("Mixing length"), "mixing_length")
        self.modelTurbModel.addItem(self.tr("k-epsilon"), "k-epsilon")
        self.modelTurbModel.addItem(self.tr("k-epsilon Linear Production"), "k-epsilon-PL")
        self.modelTurbModel.addItem(self.tr("Launder-Sharma"), "Launder-Sharma")
        self.modelTurbModel.addItem(self.tr("Rij-epsilon LRR"), "Rij-epsilon")
        self.modelTurbModel.addItem(self.tr("Rij-epsilon SSG"), "Rij-SSG")
        self.modelTurbModel.addItem(self.tr("Rij-epsilon EBRSM"), "Rij-EBRSM")
        self.modelTurbModel.addItem(self.tr("v2f BL-v2/k"), "v2f-BL-v2/k")
        self.modelTurbModel.addItem(self.tr("k-omega SST"), "k-omega-SST")
        self.modelTurbModel.addItem(self.tr("Spalart-Allmaras"), "Spalart-Allmaras")
        self.modelTurbModel.addItem(self.tr("LES (Smagorinsky)"), "LES_Smagorinsky")
        self.modelTurbModel.addItem(self.tr("LES (classical dynamic model)"), "LES_dynamique")
        self.modelTurbModel.addItem(self.tr("LES (WALE)"), "LES_WALE")

        self.modelLength = ComboModel(self.comboBoxLength,2,1)
        self.modelLength.addItem(self.tr("Automatic"), 'automatic')
        self.modelLength.addItem(self.tr("Prescribed"), 'prescribed')
        self.comboBoxLength.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        # Connections

        self.comboBoxTurbModel.activated[str].connect(self.slotTurbulenceModel)
        self.pushButtonAdvanced.clicked.connect(self.slotAdvancedOptions)
        self.lineEditLength.textChanged[str].connect(self.slotLengthScale)

        self.lineEditV0.textChanged[str].connect(self.slotVelocity)
        self.comboBoxLength.activated[str].connect(self.slotLengthChoice)
        self.lineEditL0.textChanged[str].connect(self.slotLength)

        # Frames display

        self.frameAdvanced.hide()
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
        model = self.model.getTurbulenceModel()

        self.frameAdvanced.hide()
        self.frameLength.hide()
        self.groupBoxReferenceValues.hide()

        if model not in ('off', 'mixing_length', 'LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            self.groupBoxReferenceValues.show()

        if model == 'mixing_length':
            self.frameLength.show()
            self.frameAdvanced.hide()
            self.model.getLengthScale()
        elif model not in ('off', 'LES_Smagorinsky', 'LES_dynamique', 'LES_WALE', 'Spalart-Allmaras'):
            self.frameLength.hide()
            self.frameAdvanced.show()

        if model in ('off', 'LES_Smagorinsky', 'LES_dynamique', 'LES_WALE', 'Spalart-Allmaras'):
            self.line.hide()
        else:
            self.line.show()


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
    def slotAdvancedOptions(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        default = {}
        default['model']         = self.model.getTurbulenceModel()
        default['wall_function'] = self.model.getWallFunction()
        default['gravity_terms'] = self.model.getGravity()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = TurbulenceAdvancedOptionsDialogView(self, self.case, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.model.setTurbulenceModel(result['model'])
            self.model.setWallFunction(result['wall_function'])
            self.model.setGravity(result['gravity_terms'])


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
