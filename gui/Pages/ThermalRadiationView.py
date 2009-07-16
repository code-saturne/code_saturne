# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module defines the Thermal Radiation model

This module contains the following classes and function:
- ThermalRadiationView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from ThermalRadiationForm import Ui_ThermalRadiationForm
from ThermalRadiationAdvancedDialogForm import Ui_ThermalRadiationAdvancedDialogForm
import Base.QtPage as QtPage
from Pages.ThermalRadiationModel import ThermalRadiationModel
from Pages.OutputControlModel import OutputControlModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermalRadiationView")
log.setLevel(GuiParam.DEBUG)

#--------------------------------------------------------------------------------
# Popup Class
#--------------------------------------------------------------------------------

class ThermalRadiationAdvancedDialogView(QDialog, Ui_ThermalRadiationAdvancedDialogForm):
    """
    Building of popup window for advanced options.
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_ThermalRadiationAdvancedDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result  = self.default.copy()

        # Combo models

        self.modelTSRay  = QtPage.ComboModel(self.comboBoxTSRay, 3, 1)
        self.modelPrintT = QtPage.ComboModel(self.comboBoxPrintT, 3, 1)
        self.modelPrintL = QtPage.ComboModel(self.comboBoxPrintL, 3, 1)

        self.modelTSRay.addItem('0', '0')
        self.modelTSRay.addItem('1', '1')
        self.modelTSRay.addItem('2', '2')

        self.modelPrintT.addItem('0', '0')
        self.modelPrintT.addItem('1', '1')
        self.modelPrintT.addItem('2', '2')

        self.modelPrintL.addItem('0', '0')
        self.modelPrintL.addItem('1', '1')
        self.modelPrintL.addItem('2', '2')

        self.frequ     = self.default['frequency']
        self.tsr       = self.default['idiver']
        self.printTemp = self.default['tempP']
        self.printLum  = self.default['intensity']
        model          = self.default['model']

        # Initialization

        self.lineEditFreq.setText(QString(str(self.frequ)))
        self.modelTSRay.setItem(str_model=str(self.tsr))
        self.modelPrintT.setItem(str_model=str(self.printTemp))
        self.modelPrintL.setItem(str_model=str(self.printLum))

        if model == 'dom':
            self.labelPrintL.show()
            self.comboBoxPrintL.show()
        else:
            self.labelPrintL.hide()
            self.comboBoxPrintL.hide()

        # Validator

        validatorFreq = QtPage.IntValidator(self.lineEditFreq, min=1)
        self.lineEditFreq.setValidator(validatorFreq)


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """
        if self.lineEditFreq.validator().state == QValidator.Acceptable:
            self.result['frequency'], ok = self.lineEditFreq.text().toInt()
        self.result['idiver'], ok    = self.comboBoxTSRay.currentText().toInt()
        self.result['tempP'], ok     = self.comboBoxPrintT.currentText().toInt()
        self.result['intensity'], ok = self.comboBoxPrintL.currentText().toInt()

        QDialog.accept(self)


    def reject(self):
        """
        Method called when 'Cancel' button is clicked.
        """
        QDialog.reject(self)


    def get_result(self):
        """
        Method to get the result.
        """
        return self.result


    def tr(self, text):
        """
        Translation.
        """
        return text


#--------------------------------------------------------------------------------
# Main class
#--------------------------------------------------------------------------------

class ThermalRadiationView(QWidget, Ui_ThermalRadiationForm):
    """
    Class to open Thermal Scalar Transport Page.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ThermalRadiationForm.__init__(self)
        self.setupUi(self)

        self.browser = tree
        self.case = case
        self.mdl = ThermalRadiationModel(self.case)

        # Combo models

        self.modelRadModel   = QtPage.ComboModel(self.comboBoxRadModel, 3, 1)
        self.modelDirection  = QtPage.ComboModel(self.comboBoxDirection, 2, 1)
        self.modelAbsorption = QtPage.ComboModel(self.comboBoxAbsorption, 3, 1)

        self.modelRadModel.addItem("No radiative transfers", 'off')
        self.modelRadModel.addItem("Discrete ordinates method", 'dom')
        self.modelRadModel.addItem("P-1 Model", 'p-1')

        self.modelDirection.addItem("32")
        self.modelDirection.addItem("128")

        # Connections

        self.connect(self.comboBoxRadModel,
                     SIGNAL("activated(const QString&)"),
                     self.slotRadiativeTransfer)
        self.connect(self.radioButtonOn,
                     SIGNAL("clicked()"),
                     self.slotStartRestart)
        self.connect(self.radioButtonOff,
                     SIGNAL("clicked()"),
                     self.slotStartRestart)
        self.connect(self.comboBoxDirection,
                     SIGNAL("activated(const QString&)"),
                     self.slotDirection)
        self.connect(self.comboBoxAbsorption,
                     SIGNAL("activated(const QString&)"),
                     self.slotTypeCoefficient)
        self.connect(self.lineEditCoeff,
                     SIGNAL("textChanged(const QString &)"),
                     self.slotAbsorptionCoefficient)
        self.connect(self.toolButtonAdvanced,
                     SIGNAL("clicked()"), 
                     self.slotAdvancedOptions)

        # Validator

        validatorCoeff = QtPage.DoubleValidator(self.lineEditCoeff, min=0.0)
        self.lineEditCoeff.setValidator(validatorCoeff)

        self.modelAbsorption.addItem('constant',                   'constant')
        self.modelAbsorption.addItem('user subroutine (usray3)',   'variable')
        self.modelAbsorption.addItem('user law',                   'formula')
        self.modelAbsorption.addItem('H2O and CO2 mixing (Modak)', 'modak')

        if self.mdl.isCoalCombustion():
            self.modelAbsorption.disableItem(str_model='variable')
            self.modelAbsorption.enableItem(str_model='modak')
        else:
            self.modelAbsorption.disableItem(str_model='modak')
            self.modelAbsorption.enableItem(str_model='variable')

        self.modelAbsorption.disableItem(str_model='formula')

        # Initialization

        self.modelRadModel.setItem(str_model=self.mdl.getRadiativeModel())
        self.slotRadiativeTransfer()

        if self.mdl.getRestart() == 'on':
            self.radioButtonOn.setChecked(True)
            self.radioButtonOff.setChecked(False)
        else:
            self.radioButtonOn.setChecked(False)
            self.radioButtonOff.setChecked(True)

        value = self.mdl.getTypeCoeff()
        self.modelAbsorption.setItem(str_model=value)
        self.slotTypeCoefficient(self.modelAbsorption.dicoM2V[value])

        self.pushButtonCoeffFormula.setEnabled(False)

        self.lineEditCoeff.setText(QString(str(self.mdl.getAbsorCoeff())))


    @pyqtSignature("const QString &")
    def slotRadiativeTransfer(self):
        """
        """
        model = self.modelRadModel.dicoV2M[str(self.comboBoxRadModel.currentText())]
        self.mdl.setRadiativeModel(model)
        if model == 'off':
            self.frameOptions.hide()
            self.line.hide()
            OutputControlModel(self.case).setDomainBoundaryPostProStatus('off')
        else:
            self.frameOptions.show()
            self.line.show()
            OutputControlModel(self.case).setDomainBoundaryPostProStatus('on')

            if model == 'p-1':
                self.frameDirection.hide()
            elif model == 'dom':
                self.frameDirection.show()
                self.modelDirection.setItem(str_model=str(self.mdl.getNbDir()))

        self.browser.configureTree(self.case)


    @pyqtSignature("")
    def slotStartRestart(self):
        """
        """
        if self.radioButtonOn.isChecked():
            self.mdl.setRestart("on")
        else:
            self.mdl.setRestart("off")


    @pyqtSignature("const QString &")
    def slotDirection(self, text):
        """
        """
        n, ok = text.toInt()
        self.mdl.setNbDir(n)


    @pyqtSignature("const QString &")
    def slotTypeCoefficient(self, text):
        """
        """
        typeCoeff = self.modelAbsorption.dicoV2M[str(text)]
        self.mdl.setTypeCoeff(typeCoeff)

        if typeCoeff == 'constant':
            self.lineEditCoeff.setEnabled(True)
        elif typeCoeff == 'modak':
            self.lineEditCoeff.setDisabled(True)
        else:
            self.lineEditCoeff.setDisabled(True)


    @pyqtSignature("const QString &")
    def slotAbsorptionCoefficient(self, text):
        """
        """
        c, ok  = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setAbsorCoeff(c)


    @pyqtSignature("")
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        default = {}
        default['frequency'] = self.mdl.getFrequency()
        default['idiver']    = self.mdl.getTrs()
        default['tempP']     = self.mdl.getTemperatureListing()
        default['intensity'] = self.mdl.getIntensityResolution()
        default['model']     = self.mdl.getRadiativeModel()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = ThermalRadiationAdvancedDialogView(self, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.mdl.setFrequency(result['frequency'])
            self.mdl.setTrs(result['idiver'])
            self.mdl.setTemperatureListing(result['tempP'])
            self.mdl.setIntensityResolution(result['intensity'])


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
