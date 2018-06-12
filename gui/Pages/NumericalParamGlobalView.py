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
This module contains the following class:
- NumericalParamGlobalView
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
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from code_saturne.Pages.NumericalParamGlobalForm import Ui_NumericalParamGlobalForm
from code_saturne.Pages.NumericalParamGlobalModel import NumericalParamGlobalModel
from code_saturne.Pages.SteadyManagementModel import SteadyManagementModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NumericalParamGlobalView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class NumericalParamGlobalView(QWidget, Ui_NumericalParamGlobalForm):
    """
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NumericalParamGlobalForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = NumericalParamGlobalModel(self.case)
        self.browser = tree

        self.labelSRROM.hide()
        self.lineEditSRROM.hide()
        self.line_5.hide()

        # Combo models
        self.modelEXTRAG = ComboModel(self.comboBoxEXTRAG,2,1)
        self.modelIMRGRA = ComboModel(self.comboBoxIMRGRA,5,1)
        self.modelNTERUP = ComboModel(self.comboBoxNTERUP,3,1)

        self.modelEXTRAG.addItem(self.tr("Neumann 1st order"), 'neumann')
        self.modelEXTRAG.addItem(self.tr("Extrapolation"), 'extrapolation')

        self.modelIMRGRA.addItem(self.tr("Iterative handling of non-orthogonalities"),'0')
        self.modelIMRGRA.addItem(self.tr("Least squares method over neighboring cells"),'1')
        self.modelIMRGRA.addItem(self.tr("Least squares method over extended cell neighborhood"),'2')
        self.modelIMRGRA.addItem(self.tr("Least squares method over partial extended cell neighborhood"),'3')
        self.modelIMRGRA.addItem(self.tr("Iterative method with least squares initialization"),'4')

        self.modelNTERUP.addItem(self.tr("SIMPLE"), 'simple')
        self.modelNTERUP.addItem(self.tr("SIMPLEC"),'simplec')
        self.modelNTERUP.addItem(self.tr("PISO"),'piso')

        self.comboBoxEXTRAG.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.comboBoxNTERUP.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        # Connections
        self.checkBoxIVISSE.clicked.connect(self.slotIVISSE)
        self.checkBoxIPUCOU.clicked.connect(self.slotIPUCOU)
        self.checkBoxICFGRP.clicked.connect(self.slotICFGRP)
        self.checkBoxImprovedPressure.clicked.connect(self.slotImprovedPressure)
        self.comboBoxEXTRAG.activated[str].connect(self.slotEXTRAG)
        self.lineEditRELAXP.textChanged[str].connect(self.slotRELAXP)
        self.comboBoxIMRGRA.activated[str].connect(self.slotIMRGRA)
        self.lineEditSRROM.textChanged[str].connect(self.slotSRROM)
        self.comboBoxNTERUP.activated[str].connect(self.slotNTERUP)
        self.spinBoxNTERUP.valueChanged[int].connect(self.slotNTERUP2)

        # Validators
        validatorRELAXP = DoubleValidator(self.lineEditRELAXP, min=0., max=1.)
        validatorRELAXP.setExclusiveMin(True)
        validatorSRROM = DoubleValidator(self.lineEditSRROM, min=0., max=1.)
        validatorSRROM.setExclusiveMax(True)
        self.lineEditRELAXP.setValidator(validatorRELAXP)
        self.lineEditSRROM.setValidator(validatorSRROM)

        if self.model.getTransposedGradient() == 'on':
            self.checkBoxIVISSE.setChecked(True)
        else:
            self.checkBoxIVISSE.setChecked(False)

        if self.model.getVelocityPressureCoupling() == 'on':
            self.checkBoxIPUCOU.setChecked(True)
        else:
            self.checkBoxIPUCOU.setChecked(False)

        import code_saturne.Pages.FluidCharacteristicsModel as FluidCharacteristics
        fluid = FluidCharacteristics.FluidCharacteristicsModel(self.case)
        modl_atmo, modl_joul, modl_thermo, modl_gas, modl_coal, modl_comp = fluid.getThermoPhysicalModel()

        if self.model.getHydrostaticPressure() == 'on':
            self.checkBoxImprovedPressure.setChecked(True)
        else:
            self.checkBoxImprovedPressure.setChecked(False)

        self.lineEditRELAXP.setText(str(self.model.getPressureRelaxation()))
        self.modelEXTRAG.setItem(str_model=self.model.getWallPressureExtrapolation())
        self.modelIMRGRA.setItem(str_model=str(self.model.getGradientReconstruction()))

        if modl_joul != 'off' or modl_gas != 'off' or modl_coal != 'off':
            self.labelSRROM.show()
            self.lineEditSRROM.show()
            self.lineEditSRROM.setText(str(self.model.getDensityRelaxation()))
            self.line_5.show()

        algo = self.model.getVelocityPressureAlgorithm()

        from code_saturne.Pages.TimeStepModel import TimeStepModel
        idtvar = TimeStepModel(self.case).getTimePassing()
        if idtvar in [-1, 2]:
            self.modelNTERUP.enableItem(str_model = 'simple')
            self.modelNTERUP.disableItem(str_model = 'piso')
        else:
            self.modelNTERUP.disableItem(str_model = 'simple')
            self.modelNTERUP.enableItem(str_model = 'piso')

        self.modelNTERUP.setItem(str_model=algo)

        if algo == 'piso':
            self.spinBoxNTERUP.show()
        else:
            self.spinBoxNTERUP.hide()

        if modl_comp != 'off':
            self.labelICFGRP.show()
            self.checkBoxICFGRP.show()
            self.line_4.show()
            if self.model.getHydrostaticEquilibrium() == 'on':
                self.checkBoxICFGRP.setChecked(True)
            else:
                self.checkBoxICFGRP.setChecked(False)
            self.checkBoxIPUCOU.hide()
            self.labelIPUCOU.hide()
            self.lineEditRELAXP.hide()
            self.labelRELAXP.hide()
            self.checkBoxImprovedPressure.hide()
            self.labelImprovedPressure.hide()
            self.line_2.hide()
            self.line_5.hide()
            self.line_7.hide()
            self.line_8.hide()
            self.labelNTERUP.setText("Velocity-Pressure algorithm\nsub-iterations on Navier-Stokes")
            self.comboBoxNTERUP.hide()
            self.spinBoxNTERUP.show()
        else:
            self.labelICFGRP.hide()
            self.checkBoxICFGRP.hide()
            self.line_4.hide()
            self.checkBoxIPUCOU.show()
            self.labelIPUCOU.show()
            self.lineEditRELAXP.show()
            self.labelRELAXP.show()
            self.checkBoxImprovedPressure.show()
            self.labelImprovedPressure.show()
            self.line_2.show()
            self.line_5.show()
            self.line_7.show()
            self.line_8.show()

        value = self.model.getPisoSweepNumber()
        self.spinBoxNTERUP.setValue(value)

        # Update the Tree files and folders
        self.browser.configureTree(self.case)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotIVISSE(self):
        """
        Set value for parameter IVISSE
        """
        if self.checkBoxIVISSE.isChecked():
            self.model.setTransposedGradient("on")
        else:
            self.model.setTransposedGradient("off")


    @pyqtSlot()
    def slotIPUCOU(self):
        """
        Set value for parameter IPUCOU
        """
        if self.checkBoxIPUCOU.isChecked():
            self.model.setVelocityPressureCoupling("on")
        else:
            self.model.setVelocityPressureCoupling("off")


    @pyqtSlot()
    def slotICFGRP(self):
        """
        Set value for parameter IPUCOU
        """
        if self.checkBoxICFGRP.isChecked():
            self.model.setHydrostaticEquilibrium("on")
        else:
            self.model.setHydrostaticEquilibrium("off")


    @pyqtSlot()
    def slotImprovedPressure(self):
        """
        Input IHYDPR.
        """
        if self.checkBoxImprovedPressure.isChecked():
            self.model.setHydrostaticPressure("on")
        else:
            self.model.setHydrostaticPressure("off")


    @pyqtSlot(str)
    def slotEXTRAG(self, text):
        """
        Set value for parameter EXTRAG
        """
        extrag = self.modelEXTRAG.dicoV2M[str(text)]
        self.model.setWallPressureExtrapolation(extrag)
        log.debug("slotEXTRAG-> %s" % extrag)


    @pyqtSlot(str)
    def slotRELAXP(self, text):
        """
        Set value for parameter RELAXP
        """
        if self.lineEditRELAXP.validator().state == QValidator.Acceptable:
            relaxp = from_qvariant(text, float)
            self.model.setPressureRelaxation(relaxp)
            log.debug("slotRELAXP-> %s" % relaxp)


    @pyqtSlot(str)
    def slotSRROM(self, text):
        """
        Set value for parameter SRROM
        """
        if self.lineEditSRROM.validator().state == QValidator.Acceptable:
            srrom = from_qvariant(text, float)
            self.model.setDensityRelaxation(srrom)
            log.debug("slotSRROM-> %s" % srrom)


    @pyqtSlot(str)
    def slotIMRGRA(self, text):
        """
        Set value for parameter IMRGRA
        """
        imrgra = self.modelIMRGRA.dicoV2M[str(text)]
        self.model.setGradientReconstruction(int(imrgra))
        log.debug("slotIMRGRA-> %s" % imrgra)


    @pyqtSlot(str)
    def slotNTERUP(self,text):
        """
        Set value for parameterNTERUP
        """
        NTERUP = self.modelNTERUP.dicoV2M[str(text)]
        self.model.setVelocityPressureAlgorithm(NTERUP)
        if NTERUP == 'piso':
            self.spinBoxNTERUP.show()
            value = self.model.getPisoSweepNumber()
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
        self.model.setPisoSweepNumber(var)
        log.debug("slotNTERUP2-> %s" % var)


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
