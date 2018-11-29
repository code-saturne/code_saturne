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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator, from_qvariant
from code_saturne.Pages.TimeStepForm import Ui_TimeStepForm
from code_saturne.Pages.TimeStepModel import TimeStepModel

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
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_TimeStepForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = TimeStepModel(self.case)
        self.browser = tree

       # Combo model

        self.modelTimeOptions = ComboModel(self.comboBoxOptions,4,1)
        self.modelTimeOptions.addItem(self.tr("Constant"), '0')
        self.modelTimeOptions.addItem(self.tr("Time varying (adaptive)"), '1')
        self.modelTimeOptions.addItem(self.tr("Space & time varying (pseudo-steady)"), '2')
        self.modelTimeOptions.addItem(self.tr("Steady (constant relaxation coefficient)"), '-1')

        self.modelNTERUP = ComboModel(self.comboBoxNTERUP,3,1)
        self.modelNTERUP.addItem(self.tr("SIMPLE"), 'simple')
        self.modelNTERUP.addItem(self.tr("SIMPLEC"), 'simplec')
        self.modelNTERUP.addItem(self.tr("PISO"), 'piso')
        self.comboBoxNTERUP.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        # Connections
        self.comboBoxOptions.activated[str].connect(self.slotTimePassing)
        self.lineEditDTREF.textChanged[str].connect(self.slotTimeStep)
        self.lineEditRELXST.textChanged[str].connect(self.slotRelaxCoef)
        self.lineEditNTMABS.textChanged[str].connect(self.slotIter)
        self.lineEditCOUMAX.textChanged[str].connect(self.slotTimeOptionCOUMAX)
        self.lineEditFOUMAX.textChanged[str].connect(self.slotTimeOptionFOUMAX)
        self.lineEditCDTMIN.textChanged[str].connect(self.slotTimeOptionCDTMIN)
        self.lineEditCDTMAX.textChanged[str].connect(self.slotTimeOptionCDTMAX)
        self.lineEditVARRDT.textChanged[str].connect(self.slotTimeOptionVARRDT)
        self.checkBoxIPTLRO.clicked.connect(self.slotThermalTimeStep)
        self.checkBoxINPDT0.clicked.connect(self.slotZeroTimeStep)
        self.comboBoxNTERUP.activated[str].connect(self.slotNTERUP)
        self.spinBoxNTERUP.valueChanged[int].connect(self.slotNTERUP2)

        # Validators

        validatorDTREF = DoubleValidator(self.lineEditDTREF, min=0.0)
        validatorDTREF.setExclusiveMin(True)
        validatorRELXST = DoubleValidator(self.lineEditRELXST, min=0.0, max=1.0)
        validatorRELXST.setExclusiveMin(True)
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
        self.lineEditRELXST.setValidator(validatorRELXST)
        self.lineEditNTMABS.setValidator(validatorNTMABS)
        self.lineEditCOUMAX.setValidator(validatorCOUMAX)
        self.lineEditFOUMAX.setValidator(validatorFOUMAX)
        self.lineEditCDTMIN.setValidator(validatorCDTMIN)
        self.lineEditCDTMAX.setValidator(validatorCDTMAX)
        self.lineEditVARRDT.setValidator(validatorVARRDT)

        # Initialization

        idtvar = self.mdl.getTimePassing()
        idtvar_p = idtvar

        # Constraints on time step from Turbulence model

        from code_saturne.Pages.TurbulenceModel import TurbulenceModel
        model = TurbulenceModel(self.case).getTurbulenceModel()
        del TurbulenceModel

        if model in ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            idtvar = 0
            self.modelTimeOptions.disableItem(str_model='1')
            self.modelTimeOptions.disableItem(str_model='2')
            self.modelTimeOptions.disableItem(str_model='-1')

        # Constraints on time step from Lagrangian model

        from code_saturne.Pages.LagrangianModel import LagrangianModel
        model = LagrangianModel(self.case).getLagrangianModel()
        if model in ['one_way', 'two_way']:
            if idtvar not in [0, 1]:
                idtvar = 0
            self.modelTimeOptions.disableItem(str_model='2')
            self.modelTimeOptions.disableItem(str_model='-1')
            if model == 'two_way':
                idtvar = 0
                self.modelTimeOptions.disableItem(str_model='1')

        # Constraints on time step from compressible model

        from code_saturne.Pages.CompressibleModel import CompressibleModel
        model = CompressibleModel(self.case).getCompressibleModel()
        if model != 'off':
            if idtvar not in [0, 1]:
                idtvar = 0
            self.modelTimeOptions.disableItem(str_model='2')
            self.modelTimeOptions.disableItem(str_model='-1')
            self.labelNTERUP.setText("Velocity-Pressure algorithm\nsub-iterations on Navier-Stokes")
            self.comboBoxNTERUP.hide()
            self.spinBoxNTERUP.show()

        # Constraints on time step from groundwater model

        from code_saturne.Pages.GroundwaterModel import GroundwaterModel
        model = GroundwaterModel(self.case).getGroundwaterModel()
        if model != 'off':
            if idtvar not in [0, 1]:
                idtvar = 0
            self.modelTimeOptions.disableItem(str_model='1')
            self.modelTimeOptions.disableItem(str_model='2')
            self.modelTimeOptions.disableItem(str_model='-1')

        # Change time step option if required by model constraints

        if idtvar_p != idtvar:
            self.mdl.setTimePassing(idtvar)

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

        self.__setTimePassingDisplay(idtvar)

        self.case.undoStartGlobal()


    def __setTimePassingDisplay(self, idtvar):
        """
        Choices based on IDTVAR.
        """

        self.modelTimeOptions.setItem(str_model=str(idtvar))

        if idtvar == -1:
            self.modelNTERUP.enableItem(str_model = 'simple')
            self.modelNTERUP.disableItem(str_model = 'simplec')
            self.modelNTERUP.disableItem(str_model = 'piso')
        elif idtvar == 2:
            self.modelNTERUP.disableItem(str_model = 'simple')
            self.modelNTERUP.enableItem(str_model = 'simplec')
            self.modelNTERUP.disableItem(str_model = 'piso')
        else:
            self.modelNTERUP.disableItem(str_model = 'simple')
            self.modelNTERUP.enableItem(str_model = 'simplec')
            self.modelNTERUP.enableItem(str_model = 'piso')

        algo = self.mdl.getVelocityPressureAlgorithm()
        self.modelNTERUP.setItem(str_model=algo)
        if algo == 'piso':
            self.spinBoxNTERUP.show()
        else:
            self.spinBoxNTERUP.hide()

        value = self.mdl.getPisoSweepNumber()
        self.spinBoxNTERUP.setValue(value)

        if idtvar == -1:
            self.labelRELXST.show()
            self.lineEditRELXST.show()
            self.labelDTREF.hide()
            self.labelDTREFunit.hide()
            self.lineEditDTREF.hide()
            relax_coef = self.mdl.getRelaxCoefficient()
            self.lineEditRELXST.setText(str(relax_coef))

        else:
            self.labelRELXST.hide()
            self.lineEditRELXST.hide()
            self.labelDTREF.show()
            self.labelDTREFunit.show()
            self.lineEditDTREF.show()
            dtref = self.mdl.getTimeStep()
            self.lineEditDTREF.setText(str(dtref))

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


    @pyqtSlot(str)
    def slotTimePassing(self, text):
        """
        Input IDTVAR.
        """
        idtvar = int(self.modelTimeOptions.dicoV2M[str(text)])

        self.mdl.setTimePassing(idtvar)

        self.__setTimePassingDisplay(idtvar)


    @pyqtSlot(str)
    def slotTimeStep(self, text):
        """
        Input DTREF.
        """
        if self.lineEditDTREF.validator().state == QValidator.Acceptable:
            time_step = from_qvariant(text, float)
            self.mdl.setTimeStep(time_step)


    @pyqtSlot(str)
    def slotNTERUP(self,text):
        """
        Set value for parameterNTERUP
        """
        NTERUP = self.modelNTERUP.dicoV2M[str(text)]
        self.mdl.setVelocityPressureAlgorithm(NTERUP)
        if NTERUP == 'piso':
            self.spinBoxNTERUP.show()
            value = self.mdl.getPisoSweepNumber()
            self.spinBoxNTERUP.setValue(value)
        else:
            self.spinBoxNTERUP.hide()
        self.browser.configureTree(self.case)
        log.debug("slotNTERUP-> %s" % NTERUP)


    @pyqtSlot(int)
    def slotNTERUP2(self, var):
        """
        Set value for parameter piso sweep number
        """
        self.mdl.setPisoSweepNumber(var)
        log.debug("slotNTERUP2-> %s" % var)


    @pyqtSlot(str)
    def slotRelaxCoef(self, text):
        """
        Input relaxation coefficient.
        """
        if self. lineEditRELXST.validator().state == QValidator.Acceptable:
            relax_coef = from_qvariant(text, float)
            self.mdl.setRelaxCoefficient(relax_coef)


    @pyqtSlot(str)
    def slotIter(self, text):
        """
        Input NTMABS.
        """
        if self.lineEditNTMABS.validator().state == QValidator.Acceptable:
            iteration = from_qvariant(text, int)
            self.mdl.setIterationsNumber(iteration)


    @pyqtSlot(str)
    def slotTimeOptionCOUMAX(self, text):
        """
        Input COUMAX.
        """
        if self.lineEditCOUMAX.validator().state == QValidator.Acceptable:
            courant_max = from_qvariant(text, float)
            self.mdl.setOptions('max_courant_num', courant_max)


    @pyqtSlot(str)
    def slotTimeOptionFOUMAX(self, text):
        """
        Input FOUMAX.
        """
        if self.lineEditFOUMAX.validator().state == QValidator.Acceptable:
            fourier_max = from_qvariant(text, float)
            self.mdl.setOptions('max_fourier_num', fourier_max)


    @pyqtSlot(str)
    def slotTimeOptionCDTMIN(self, text):
        """
        Input CDTMIN.
        """
        if self.lineEditCDTMIN.validator().state == QValidator.Acceptable:
            time_step_min_factor = from_qvariant(text, float)
            self.mdl.setOptions('time_step_min_factor', time_step_min_factor)


    @pyqtSlot(str)
    def slotTimeOptionCDTMAX(self, text):
        """
        Input CDTMAX.
        """
        if self.lineEditCDTMAX.validator().state == QValidator.Acceptable:
            time_step_max_factor = from_qvariant(text, float)
            self.mdl.setOptions('time_step_max_factor', time_step_max_factor)


    @pyqtSlot(str)
    def slotTimeOptionVARRDT(self, text):
        """
        Input VARRDT.
        """
        if self.lineEditVARRDT.validator().state == QValidator.Acceptable:
            time_step_var = from_qvariant(text, float)
            self.mdl.setOptions('time_step_var', time_step_var)


    @pyqtSlot()
    def slotThermalTimeStep(self):
        """
        Input IPTLRO.
        """
        if self.checkBoxIPTLRO.isChecked():
            self.mdl.setThermalTimeStep("on")
        else:
            self.mdl.setThermalTimeStep("off")


    @pyqtSlot()
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
