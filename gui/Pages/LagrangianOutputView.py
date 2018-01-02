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
- LagrangianOutputView
"""

#-------------------------------------------------------------------------------
# Standard modules
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
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator
from code_saturne.Pages.LagrangianOutputForm import Ui_LagrangianOutputForm
from code_saturne.Pages.LagrangianOutputModel import LagrangianOutputModel
from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianOutputView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class LagrangianOutputView(QWidget, Ui_LagrangianOutputForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LagrangianOutputForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = LagrangianOutputModel(self.case)

        # Connections
        self.checkBoxIVISV1.clicked.connect(self.slotIVISV1)
        self.checkBoxIVISV2.clicked.connect(self.slotIVISV2)
        self.checkBoxIVISTP.clicked.connect(self.slotIVISTP)
        self.checkBoxIVISDM.clicked.connect(self.slotIVISDM)
        self.checkBoxIVISTE.clicked.connect(self.slotIVISTE)
        self.checkBoxIVISMP.clicked.connect(self.slotIVISMP)
        self.checkBoxIVISDK.clicked.connect(self.slotIVISDK)
        self.checkBoxIVISCH.clicked.connect(self.slotIVISCH)
        self.checkBoxIVISCK.clicked.connect(self.slotIVISCK)
        self.checkBoxMoisture.clicked.connect(self.slotMoisture)

        # initialize Widgets
        status = self.model.getFluidVelocityStatus()
        if status == "on":
            self.checkBoxIVISV1.setChecked(True)
        else:
            self.checkBoxIVISV1.setChecked(False)

        status = self.model.getParticlesVelocityStatus()
        if status == "on":
            self.checkBoxIVISV2.setChecked(True)
        else:
            self.checkBoxIVISV2.setChecked(False)

        status = self.model.getResidentTimeStatus()
        if status == "on":
            self.checkBoxIVISTP.setChecked(True)
        else:
            self.checkBoxIVISTP.setChecked(False)

        status = self.model.getParticleDiameterStatus()
        if status == "on":
            self.checkBoxIVISDM.setChecked(True)
        else:
            self.checkBoxIVISDM.setChecked(False)

        status = self.model.getParticleTemperatureStatus()
        if status == "on":
            self.checkBoxIVISTE.setChecked(True)
        else:
            self.checkBoxIVISTE.setChecked(False)

        status = self.model.getParticleMassStatus()
        if status == "on":
            self.checkBoxIVISMP.setChecked(True)
        else:
            self.checkBoxIVISMP.setChecked(False)

        if CoalCombustionModel(self.case).getCoalCombustionModel("only") == \
           "homogeneous_fuel_moisture":
            status = self.model.getCoalParticleDiameterStatus()
            if status == "on":
                self.checkBoxIVISDK.setChecked(True)
            else:
                self.checkBoxIVISDK.setChecked(False)

            status = self.model.getCoalParticleMassStatus()
            if status == "on":
                self.checkBoxIVISCH.setChecked(True)
            else:
                self.checkBoxIVISCH.setChecked(False)

            status = self.model.getCokeParticleMassStatus()
            if status == "on":
                self.checkBoxIVISCK.setChecked(True)
            else:
                self.checkBoxIVISCK.setChecked(False)

            status = self.model.getMoistureMassStatus()
            if status == "on":
                self.checkBoxMoisture.setChecked(True)
            else:
                self.checkBoxMoisture.setChecked(False)
        else:
            self.checkBoxIVISDK.hide()
            self.labelIVISDK.hide()
            self.checkBoxIVISCH.hide()
            self.labelIVISCH.hide()
            self.checkBoxIVISCK.hide()
            self.labelIVISCK.hide()
            self.checkBoxMoisture.hide()
            self.labelIMoisture.hide()

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotIVISV1(self):
        """
        Input IVISV1.
        """
        if self.checkBoxIVISV1.isChecked():
            self.model.setFluidVelocityStatus("on")
        else:
            self.model.setFluidVelocityStatus("off")


    @pyqtSlot()
    def slotIVISV2(self):
        """
        Input IVISV2.
        """
        if self.checkBoxIVISV2.isChecked():
            self.model.setParticlesVelocityStatus("on")
        else:
            self.model.setParticlesVelocityStatus("off")


    @pyqtSlot()
    def slotIVISTP(self):
        """
        Input IVISTP.
        """
        if self.checkBoxIVISTP.isChecked():
            self.model.setResidentTimeStatus("on")
        else:
            self.model.setResidentTimeStatus("off")


    @pyqtSlot()
    def slotIVISDM(self):
        """
        Input IVISDM.
        """
        if self.checkBoxIVISDM.isChecked():
            self.model.setParticleDiameterStatus("on")
        else:
            self.model.setParticleDiameterStatus("off")


    @pyqtSlot()
    def slotIVISTE(self):
        """
        Input IVISTE.
        """
        if self.checkBoxIVISTE.isChecked():
            self.model.setParticleTemperatureStatus("on")
        else:
            self.model.setParticleTemperatureStatus("off")


    @pyqtSlot()
    def slotIVISMP(self):
        """
        Input IVISMP.
        """
        if self.checkBoxIVISMP.isChecked():
            self.model.setParticleMassStatus("on")
        else:
            self.model.setParticleMassStatus("off")


    @pyqtSlot()
    def slotIVISDK(self):
        """
        Input IVISDK.
        """
        if self.checkBoxIVISDK.isChecked():
            self.model.setCoalParticleDiameterStatus("on")
        else:
            self.model.setCoalParticleDiameterStatus("off")


    @pyqtSlot()
    def slotIVISCH(self):
        """
        Input IVISCH.
        """
        if self.checkBoxIVISCH.isChecked():
            self.model.setCoalParticleMassStatus("on")
        else:
            self.model.setCoalParticleMassStatus("off")


    @pyqtSlot()
    def slotIVISCK(self):
        """
        Input IVISCK.
        """
        if self.checkBoxIVISCK.isChecked():
            self.model.setCokeParticleMassStatus("on")
        else:
            self.model.setCokeParticleMassStatus("off")


    @pyqtSlot()
    def slotMoisture(self):
        """
        Input IVISCK.
        """
        if self.checkBoxMoisture.isChecked():
            self.model.setMoistureMassStatus("on")
        else:
            self.model.setMoistureMassStatus("off")


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
