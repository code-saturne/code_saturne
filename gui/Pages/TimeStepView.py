# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
This module contains the following classes and function:
- TimeStepView
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
from Pages.TimeStepForm import Ui_TimeStepForm
import Base.QtPage as QtPage
from Pages.TimeStepModel import TimeStepModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TimeStepView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class TimeStepView(QWidget, Ui_TimeStepForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_TimeStepForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = TimeStepModel(self.case)

       # Combo model

        self.modelTimeOptions = QtPage.ComboModel(self.comboBoxOptions,3,1)

        self.modelTimeOptions.addItem(self.tr("Uniform and constant"), '0')
        self.modelTimeOptions.addItem(self.tr("Variable in time and uniform in space"), '1')
        self.modelTimeOptions.addItem(self.tr("Variable in time and in space"), '2')

        # Connections
        self.connect(self.comboBoxOptions, SIGNAL("activated(const QString&)"), self.slotTimePassing)
        self.connect(self.lineEditDTREF, SIGNAL("textChanged(const QString &)"), self.slotTimeStep)
        self.connect(self.lineEditNTMABS, SIGNAL("textChanged(const QString &)"), self.slotIter)
        self.connect(self.lineEditCOUMAX, SIGNAL("textChanged(const QString &)"), self.slotTimeOptionCOUMAX)
        self.connect(self.lineEditFOUMAX, SIGNAL("textChanged(const QString &)"), self.slotTimeOptionFOUMAX)
        self.connect(self.lineEditCDTMIN, SIGNAL("textChanged(const QString &)"), self.slotTimeOptionCDTMIN)
        self.connect(self.lineEditCDTMAX, SIGNAL("textChanged(const QString &)"), self.slotTimeOptionCDTMAX)
        self.connect(self.lineEditVARRDT, SIGNAL("textChanged(const QString &)"), self.slotTimeOptionVARRDT)
        self.connect(self.checkBoxIPTLRO, SIGNAL("clicked()"), self.slotThermalTimeStep)
        self.connect(self.checkBoxINPDT0, SIGNAL("clicked()"), self.slotZeroTimeStep)

        # Validators

        validatorDTREF = QtPage.DoubleValidator(self.lineEditDTREF, min=0.0)
        validatorDTREF.setExclusiveMin(True)
        validatorNTMABS = QtPage.IntValidator(self.lineEditNTMABS, min=1)
        validatorCOUMAX = QtPage.DoubleValidator(self.lineEditCOUMAX, min=0.0)
        validatorCOUMAX.setExclusiveMin(True)
        validatorFOUMAX = QtPage.DoubleValidator(self.lineEditFOUMAX, min=0.0)
        validatorFOUMAX.setExclusiveMin(True)
        validatorCDTMIN = QtPage.DoubleValidator(self.lineEditCDTMIN, min=0.0, max=1.0)
        validatorCDTMIN.setExclusiveMin(True)
        validatorCDTMAX = QtPage.DoubleValidator(self.lineEditCDTMAX, min=1.0)
        validatorVARRDT = QtPage.DoubleValidator(self.lineEditVARRDT, min=0.0, max=1.0)
        validatorVARRDT.setExclusiveMin(True)

        self.lineEditDTREF.setValidator(validatorDTREF)
        self.lineEditNTMABS.setValidator(validatorNTMABS)
        self.lineEditCOUMAX.setValidator(validatorCOUMAX)
        self.lineEditFOUMAX.setValidator(validatorFOUMAX)
        self.lineEditCDTMIN.setValidator(validatorCDTMIN)
        self.lineEditCDTMAX.setValidator(validatorCDTMAX)
        self.lineEditVARRDT.setValidator(validatorVARRDT)

        # Initialization

        idtvar = self.mdl.getTimePassing()
        self.modelTimeOptions.setItem(str_model=str(idtvar))

        from Pages.TurbulenceModel import TurbulenceModel
        model = TurbulenceModel(self.case).getTurbulenceModel()
        del TurbulenceModel

        if model in ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            idtvar = 0
            self.modelTimeOptions.setItem(str_model=str(idtvar))
            self.modelTimeOptions.disableItem(str_model='0')
            self.modelTimeOptions.disableItem(str_model='1')
            self.modelTimeOptions.disableItem(str_model='2')

        text = self.comboBoxOptions.currentText()
        self.slotTimePassing(text)

        dtref = self.mdl.getTimeStep()
        self.lineEditDTREF.setText(QString(str(dtref)))

        ntmabs = self.mdl.getIterationsNumber()
        self.lineEditNTMABS.setText(QString(str(ntmabs)))

        if self.mdl.thermalCase():
            if self.mdl.getThermalTimeStep() == 'on':
                self.checkBoxIPTLRO.setChecked(True)
            else:
                self.checkBoxIPTLRO.setChecked(False)
        else:
            self.lineIPTLRO.hide()
            self.labelIPTLRO.hide()
            self.checkBoxIPTLRO.hide()
            self.mdl.RemoveThermalTimeStepNode()

        if self.mdl.getZeroTimeStep() == 'on':
            self.checkBoxINPDT0.setChecked(True)
        else:
            self.checkBoxINPDT0.setChecked(False)


    @pyqtSignature("")
    def slotTimePassing(self, text):
        """
        Input IDTVAR.
        """
        idtvar = int(self.modelTimeOptions.dicoV2M[str(text)])

        self.mdl.setTimePassing(idtvar)

        if idtvar in (1, 2):
            courant_max   = self.mdl.getOptions('max_courant_num')
            fourier_max   = self.mdl.getOptions('max_fourier_num')
            time_step_min_factor = self.mdl.getOptions('time_step_min_factor')
            time_step_max_factor = self.mdl.getOptions('time_step_max_factor')
            time_step_var = self.mdl.getOptions('time_step_var')

            self.lineEditCOUMAX.setText(QString(str(courant_max)))
            self.lineEditFOUMAX.setText(QString(str(fourier_max)))
            self.lineEditCDTMIN.setText(QString(str(time_step_min_factor)))
            self.lineEditCDTMAX.setText(QString(str(time_step_max_factor)))
            self.lineEditVARRDT.setText(QString(str(time_step_var)))

            self.groupBoxLabels.show()
        else:
            self.groupBoxLabels.hide()


    @pyqtSignature("const QString &")
    def slotTimeStep(self, text):
        """
        Input DTREF.
        """
        time_step, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setTimeStep(time_step)


    @pyqtSignature("const QString &")
    def slotIter(self, text):
        """
        Input NTMABS.
        """
        iteration, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setIterationsNumber(iteration)


    @pyqtSignature("const QString &")
    def slotTimeOptionCOUMAX(self, text):
        """
        Input COUMAX.
        """
        courant_max, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOptions('max_courant_num', courant_max)


    @pyqtSignature("const QString &")
    def slotTimeOptionFOUMAX(self, text):
        """
        Input FOUMAX.
        """
        fourier_max, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOptions('max_fourier_num', fourier_max)


    @pyqtSignature("const QString &")
    def slotTimeOptionCDTMIN(self, text):
        """
        Input CDTMIN.
        """
        time_step_min_factor, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOptions('time_step_min_factor', time_step_min_factor)


    @pyqtSignature("const QString &")
    def slotTimeOptionCDTMAX(self, text):
        """
        Input CDTMAX.
        """
        time_step_max_factor, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOptions('time_step_max_factor', time_step_max_factor)


    @pyqtSignature("const QString &")
    def slotTimeOptionVARRDT(self, text):
        """
        Input VARRDT.
        """
        time_step_var, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOptions('time_step_var', time_step_var)


    @pyqtSignature("")
    def slotThermalTimeStep(self):
        """
        Input IPTLRO.
        """
        if self.checkBoxIPTLRO.isChecked():
            self.mdl.setThermalTimeStep("on")
        else:
            self.mdl.setThermalTimeStep("off")


    @pyqtSignature("")
    def slotZeroTimeStep(self):
        """
        Input INPDT0.
        """
        if self.checkBoxINPDT0.isChecked():
            self.mdl.setZeroTimeStep("on")
        else:
            self.mdl.setZeroTimeStep("off")


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
