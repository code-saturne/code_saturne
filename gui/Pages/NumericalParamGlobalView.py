# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Pages.NumericalParamGlobalForm import Ui_NumericalParamGlobalForm
import Base.QtPage as QtPage
from Pages.NumericalParamGlobalModel import NumericalParamGlobalModel

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
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NumericalParamGlobalForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.model = NumericalParamGlobalModel(self.case)

        self.labelSRROM.hide()
        self.lineEditSRROM.hide()
        self.line_5.hide()

        # Combo models
        self.modelEXTRAG = QtPage.ComboModel(self.comboBoxEXTRAG,2,1)
        self.modelIMRGRA = QtPage.ComboModel(self.comboBoxIMRGRA,5,1)

        self.modelEXTRAG.addItem(self.tr("Neumann 1st order"), 'neumann')
        self.modelEXTRAG.addItem(self.tr("Extrapolation"), 'extrapolation')

        self.modelIMRGRA.addItem(self.tr("Iterative handling of non-orthogonalities"),'0')
        self.modelIMRGRA.addItem(self.tr("Least squares method over neighboring cells"),'1')
        self.modelIMRGRA.addItem(self.tr("Least squares method over extended cell neighborhood"),'2')
        self.modelIMRGRA.addItem(self.tr("Least squares method over partial extended cell neighborhood"),'3')
        self.modelIMRGRA.addItem(self.tr("Iterative method with least squares initialization"),'4')

        self.comboBoxEXTRAG.setSizeAdjustPolicy(QComboBox.AdjustToContents)

        # Connections
        self.connect(self.checkBoxIVISSE, SIGNAL("clicked()"), self.slotIVISSE)
        self.connect(self.checkBoxIMGR,   SIGNAL("clicked()"), self.slotIMGR)
        self.connect(self.checkBoxIPUCOU, SIGNAL("clicked()"), self.slotIPUCOU)
        self.connect(self.comboBoxEXTRAG, SIGNAL("activated(const QString&)"), self.slotEXTRAG)
        self.connect(self.lineEditRELAXP, SIGNAL("textChanged(const QString &)"), self.slotRELAXP)
        self.connect(self.comboBoxIMRGRA, SIGNAL("activated(const QString&)"), self.slotIMRGRA)
        self.connect(self.lineEditSRROM,  SIGNAL("textChanged(const QString &)"), self.slotSRROM)

        # Validators
        validatorRELAXP = QtPage.DoubleValidator(self.lineEditRELAXP, min=0., max=1.)
        validatorRELAXP.setExclusiveMin(True)
        validatorSRROM = QtPage.DoubleValidator(self.lineEditSRROM, min=0., max=1.)
        validatorSRROM.setExclusiveMin(True)
        self.lineEditRELAXP.setValidator(validatorRELAXP)
        self.lineEditSRROM.setValidator(validatorSRROM)

        if self.model.getMultigrid() == 'on':
            self.checkBoxIMGR.setChecked(True)
        else:
            self.checkBoxIMGR.setChecked(False)

        if self.model.getTransposedGradient() == 'on':
            self.checkBoxIVISSE.setChecked(True)
        else:
            self.checkBoxIVISSE.setChecked(False)

        if self.model.getVelocityPressureCoupling() == 'on':
            self.checkBoxIPUCOU.setChecked(True)
        else:
            self.checkBoxIPUCOU.setChecked(False)

        self.lineEditRELAXP.setText(str(self.model.getPressureRelaxation()))
        self.modelEXTRAG.setItem(str_model=self.model.getWallPressureExtrapolation())
        self.modelIMRGRA.setItem(str_model=self.model.getGradientReconstruction())

        import Pages.FluidCharacteristicsModel as FluidCharacteristics
        fluid = FluidCharacteristics.FluidCharacteristicsModel(self.case)
        modl_atmo, modl_joul, modl_thermo, modl_gas, modl_coal = fluid.getThermoPhysicalModel()
        if modl_joul != 'off' or modl_gas != 'off' or modl_coal != 'off':
            self.labelSRROM.show()
            self.lineEditSRROM.show()
            self.lineEditSRROM.setText(str(self.model.getDensityRelaxation()))
            self.line_5.show()


    @pyqtSignature("")
    def slotIMGR(self):
        """
        Set value for parameter IMGR
        """
        if self.checkBoxIMGR.isChecked():
            self.model.setMultigrid("on")
        else:
            self.model.setMultigrid("off")


    @pyqtSignature("")
    def slotIVISSE(self):
        """
        Set value for parameter IVISSE
        """
        if self.checkBoxIVISSE.isChecked():
            self.model.setTransposedGradient("on")
        else:
            self.model.setTransposedGradient("off")


    @pyqtSignature("")
    def slotIPUCOU(self):
        """
        Set value for parameter IPUCOU
        """
        if self.checkBoxIPUCOU.isChecked():
            self.model.setVelocityPressureCoupling("on")
        else:
            self.model.setVelocityPressureCoupling("off")


    @pyqtSignature("const QString &")
    def slotEXTRAG(self, text):
        """
        Set value for parameter EXTRAG
        """
        extrag = self.modelEXTRAG.dicoV2M[str(text)]
        self.model.setWallPressureExtrapolation(extrag)
        log.debug("slotEXTRAG-> %s" % extrag)


    @pyqtSignature("const QString &")
    def slotRELAXP(self, text):
        """
        Set value for parameter RELAXP
        """
        relaxp, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.model.setPressureRelaxation(relaxp)
            log.debug("slotRELAXP-> %s" % relaxp)


    @pyqtSignature("const QString &")
    def slotSRROM(self, text):
        """
        Set value for parameter SRROM
        """
        srrom, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.model.setDensityRelaxation(srrom)
            log.debug("slotSRROM-> %s" % srrom)


    @pyqtSignature("const QString &")
    def slotIMRGRA(self, text):
        """
        Set value for parameter IMRGRA
        """
        imrgra = self.modelIMRGRA.getIndex(str_view=str(text))
        self.model.setGradientReconstruction(imrgra)
        log.debug("slotIMRGRA-> %s" % imrgra)


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
