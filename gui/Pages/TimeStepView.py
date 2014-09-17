# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator, from_qvariant
from code_saturne.Pages.TimeStepForm import Ui_TimeStepForm
from code_saturne.Pages.TimeStepModel import TimeStepModel
from code_saturne.Pages.SteadyManagementModel import SteadyManagementModel

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
        self.case.undoStopGlobal()
        self.mdl = TimeStepModel(self.case)

       # Combo model

        self.modelTimeOptions = ComboModel(self.comboBoxOptions,2,1)
        self.modelTimeOptions.addItem(self.tr("Constant"), '0')
        self.modelTimeOptions.addItem(self.tr("Variable"), '1')

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

        validatorDTREF = DoubleValidator(self.lineEditDTREF, min=0.0)
        validatorDTREF.setExclusiveMin(True)
        validatorNTMABS = IntValidator(self.lineEditNTMABS, min=1)
        validatorCOUMAX = DoubleValidator(self.lineEditCOUMAX, min=0.0)
        validatorCOUMAX.setExclusiveMin(True)
        validatorFOUMAX = DoubleValidator(self.lineEditFOUMAX, min=0.0)
        validatorFOUMAX.setExclusiveMin(True)
        validatorCDTMIN = DoubleValidator(self.lineEditCDTMIN, min=0.0, max=1.0)
        validatorCDTMIN.setExclusiveMin(True)
        validatorCDTMAX = DoubleValidator(self.lineEditCDTMAX, min=1.0)
        validatorVARRDT = DoubleValidator(self.lineEditVARRDT, min=0.0, max=1.0)
        validatorVARRDT.setExclusiveMin(True)

        self.lineEditDTREF.setValidator(validatorDTREF)
        self.lineEditNTMABS.setValidator(validatorNTMABS)
        self.lineEditCOUMAX.setValidator(validatorCOUMAX)
        self.lineEditFOUMAX.setValidator(validatorFOUMAX)
        self.lineEditCDTMIN.setValidator(validatorCDTMIN)
        self.lineEditCDTMAX.setValidator(validatorCDTMAX)
        self.lineEditVARRDT.setValidator(validatorVARRDT)

        # Initialization

        status = SteadyManagementModel(self.case).getSteadyFlowManagement()
        if status == 'on':
            self.comboBoxOptions.hide()

            self.mdl.setTimePassing(2)

            courant_max   = self.mdl.getOptions('max_courant_num')
            fourier_max   = self.mdl.getOptions('max_fourier_num')
            time_step_min_factor = self.mdl.getOptions('time_step_min_factor')
            time_step_max_factor = self.mdl.getOptions('time_step_max_factor')
            time_step_var = self.mdl.getOptions('time_step_var')

            self.lineEditCOUMAX.setText(str(courant_max))
            self.lineEditFOUMAX.setText(str(fourier_max))
            self.lineEditCDTMIN.setText(str(time_step_min_factor))
            self.lineEditCDTMAX.setText(str(time_step_max_factor))
            self.lineEditVARRDT.setText(str(time_step_var))

            self.groupBoxLabels.show()

        else:
            self.comboBoxOptions.show()

            idtvar = self.mdl.getTimePassing()
            self.modelTimeOptions.setItem(str_model=str(idtvar))

            from code_saturne.Pages.TurbulenceModel import TurbulenceModel
            model = TurbulenceModel(self.case).getTurbulenceModel()
            del TurbulenceModel

            if model in ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
                idtvar = 0
                self.modelTimeOptions.setItem(str_model=str(idtvar))
                self.modelTimeOptions.disableItem(str_model='0')
                self.modelTimeOptions.disableItem(str_model='1')

            text = self.comboBoxOptions.currentText()
            self.slotTimePassing(text)

        dtref = self.mdl.getTimeStep()
        self.lineEditDTREF.setText(str(dtref))

        ntmabs = self.mdl.getIterationsNumber()
        self.lineEditNTMABS.setText(str(ntmabs))

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

        self.case.undoStartGlobal()


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

            self.lineEditCOUMAX.setText(str(courant_max))
            self.lineEditFOUMAX.setText(str(fourier_max))
            self.lineEditCDTMIN.setText(str(time_step_min_factor))
            self.lineEditCDTMAX.setText(str(time_step_max_factor))
            self.lineEditVARRDT.setText(str(time_step_var))

            self.groupBoxLabels.show()
        else:
            self.groupBoxLabels.hide()


    @pyqtSignature("const QString &")
    def slotTimeStep(self, text):
        """
        Input DTREF.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            time_step = from_qvariant(text, float)
            self.mdl.setTimeStep(time_step)


    @pyqtSignature("const QString &")
    def slotIter(self, text):
        """
        Input NTMABS.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            iteration = from_qvariant(text, int)
            self.mdl.setIterationsNumber(iteration)


    @pyqtSignature("const QString &")
    def slotTimeOptionCOUMAX(self, text):
        """
        Input COUMAX.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            courant_max = from_qvariant(text, float)
            self.mdl.setOptions('max_courant_num', courant_max)


    @pyqtSignature("const QString &")
    def slotTimeOptionFOUMAX(self, text):
        """
        Input FOUMAX.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            fourier_max = from_qvariant(text, float)
            self.mdl.setOptions('max_fourier_num', fourier_max)


    @pyqtSignature("const QString &")
    def slotTimeOptionCDTMIN(self, text):
        """
        Input CDTMIN.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            time_step_min_factor = from_qvariant(text, float)
            self.mdl.setOptions('time_step_min_factor', time_step_min_factor)


    @pyqtSignature("const QString &")
    def slotTimeOptionCDTMAX(self, text):
        """
        Input CDTMAX.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            time_step_max_factor = from_qvariant(text, float)
            self.mdl.setOptions('time_step_max_factor', time_step_max_factor)


    @pyqtSignature("const QString &")
    def slotTimeOptionVARRDT(self, text):
        """
        Input VARRDT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            time_step_var = from_qvariant(text, float)
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
