# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.LagrangianOutputForm import Ui_LagrangianOutputForm
from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, IntValidator, DoubleValidator
from Pages.LagrangianOutputModel import LagrangianOutputModel

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

        # Combo model
        self.modelNTLAL = ComboModel(self.comboBoxNTLAL,3,1)
        self.modelNTLAL.addItem(self.tr("No output"), 'None')
        self.modelNTLAL.addItem(self.tr("Output listing at each time step"), 'At each step')
        self.modelNTLAL.addItem(self.tr("Output every 'n' time steps"), 'Frequency_l')

        # Connections
        self.connect(self.checkBoxIENSI1, SIGNAL("clicked()"),    self.slotIENSI1)
        self.connect(self.checkBoxIENSI2, SIGNAL("clicked()"),    self.slotIENSI2)
        self.connect(self.lineEditNBVIS,  SIGNAL("textChanged(const QString &)"), self.slotNBVIS)
        self.connect(self.lineEditNVISLA, SIGNAL("textChanged(const QString &)"), self.slotNVISLA)
        self.connect(self.comboBoxNTLAL,  SIGNAL("activated(const QString&)"),    self.slotChoiceNTLAL)
        self.connect(self.lineEditNTLAL,  SIGNAL("textChanged(const QString &)"), self.slotNTLAL)
        self.connect(self.checkBoxIVISV1, SIGNAL("clicked()"),    self.slotIVISV1)
        self.connect(self.checkBoxIVISV2, SIGNAL("clicked()"),    self.slotIVISV2)
        self.connect(self.checkBoxIVISTP, SIGNAL("clicked()"),    self.slotIVISTP)
        self.connect(self.checkBoxIVISDM, SIGNAL("clicked()"),    self.slotIVISDM)
        self.connect(self.checkBoxIVISTE, SIGNAL("clicked()"),    self.slotIVISTE)
        self.connect(self.checkBoxIVISMP, SIGNAL("clicked()"),    self.slotIVISMP)
        self.connect(self.checkBoxIVISHP, SIGNAL("clicked()"),    self.slotIVISHP)
        self.connect(self.checkBoxIVISDK, SIGNAL("clicked()"),    self.slotIVISDK)
        self.connect(self.checkBoxIVISCH, SIGNAL("clicked()"),    self.slotIVISCH)
        self.connect(self.checkBoxIVISCK, SIGNAL("clicked()"),    self.slotIVISCK)

        validatorNBVIS  = IntValidator(self.lineEditNBVIS, min=0)
        self.lineEditNBVIS.setValidator(validatorNBVIS)

        validatorNVISLA = IntValidator(self.lineEditNVISLA, min=0)
        #setExclusive
        self.lineEditNVISLA.setValidator(validatorNVISLA)

        validatorNTLAL = IntValidator(self.lineEditNTLAL)
        self.lineEditNTLAL.setValidator(validatorNTLAL)

        # initialize Widgets

        # post processing info to display
        status = self.model.getTrajectoryStatus()
        if status == "on":
            self.checkBoxIENSI1.setChecked(True)
        else:
            self.checkBoxIENSI1.setChecked(False)

        status = self.model.getParticlesStatus()
        if status == "on":
            self.checkBoxIENSI2.setChecked(True)
        else:
            self.checkBoxIENSI2.setChecked(False)

        self.modelFormat = ComboModel(self.comboBoxFormat,1,1)
        format = self.model.getPostProcessingFormat()
        self.modelFormat.addItem(format)
        self.modelFormat.setItem(str_model=format)
        self.comboBoxFormat.setDisabled(True)

        self.modelOption = ComboModel(self.comboBoxOptions,1,1)
        option = self.model.getPostProcessingOption()
        self.modelOption.addItem(option)
        self.modelOption.setItem(str_model=option)
        self.comboBoxOptions.setDisabled(True)

        npart = self.model.getDisplayParticlesValue()
        self.lineEditNBVIS.setText(QString(str(npart)))

        period = self.model.getPostProcessingFrequency()
        self.lineEditNVISLA.setText(QString(str(period)))

        period = self.model.getListingFrequency()
        if period == -1:
            m = "None"
        elif period == 1:
            m = "At each step"
        else:
            m = "Frequency_l"
        self.lineEditNTLAL.setText(QString(str(period)))
        t = self.modelNTLAL.dicoM2V[m]
        self.modelNTLAL.setItem(str_model = m)
        self.slotChoiceNTLAL(t)

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

        # FIXME
        # check if coal model is activated
##         coalThermoChModel = CoalThermoChemistry.CoalThermoChemistryModel("dp_FCP", self.case)
##         coals = coalThermoChModel.getCoals()
##         CoalsNumber = coals.getNumber()
##         if CoalsNumber == 0:
##             self.lineEditIVISHP.setDisabled(True)
##             self.checkBoxIVISHP.setDisabled(True)
##             self.lineEditIVISDK.setDisabled(True)
##             self.checkBoxIVISDK.setDisabled(True)
##             self.lineEditIVISCH.setDisabled(True)
##             self.checkBoxIVISCH.setDisabled(True)
##             self.lineEditIVISCK.setDisabled(True)
##             self.checkBoxIVISCK.setDisabled(True)
        status = self.model.getCoalParticleTemperatureStatus()
        if status == "on":
            self.checkBoxIVISHP.setChecked(True)
        else:
            self.checkBoxIVISHP.setChecked(False)

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

        self.case.undoStartGlobal()


    @pyqtSignature("")
    def slotIENSI1(self):
        """
        Input IENSI1.
        """
        if self.checkBoxIENSI1.isChecked():
            self.model.setTrajectoryStatus("on")
        else:
            self.model.setTrajectoryStatus("off")


    @pyqtSignature("")
    def slotIENSI2(self):
        """
        Input IENSI2.
        """
        if self.checkBoxIENSI2.isChecked():
            self.model.setParticlesStatus("on")
        else:
            self.model.setParticlesStatus("off")


    @pyqtSignature("const QString&")
    def slotNBVIS(self, text):
        """
        Input NBVIS.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value, ok = text.toInt()
            log.debug("slotNBVIS value = %i "%value)
            self.model.setDisplayParticlesValue(value)


    @pyqtSignature("const QString&")
    def slotNBVIS(self, text):
        """
        Input NBVIS.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value, ok = text.toInt()
            self.model.setDisplayParticlesValue(value)


    @pyqtSignature("const QString&")
    def slotChoiceNVISLA(self, text):
        """
        Input NVISLA.
        """
        log.debug("slotChoiceNVISLA text = %s " %str(text))


    @pyqtSignature("const QString&")
    def slotNVISLA(self, text):
        """
        Input NVISLA.
        """
        log.debug("slotNVISLA text = %s " %str(text))


    @pyqtSignature("const QString&")
    def slotChoiceNTLAL(self, text):
        """
        Input NTLAL.
        """
        listing = self.modelNTLAL.dicoV2M[str(text)]
        log.debug("slotChoiceNTLAL-> listing = %s" % listing)

        if listing == "None":
            ntlist = -1
            self.model.setListingFrequency(ntlist)
            self.lineEditNTLAL.setText(QString(str(ntlist)))
            self.lineEditNTLAL.setDisabled(True)

        elif listing == "At each step":
            ntlist = 1
            self.model.setListingFrequency(ntlist)
            self.lineEditNTLAL.setText(QString(str(ntlist)))
            self.lineEditNTLAL.setDisabled(True)

        elif listing == "Frequency_l":
            self.lineEditNTLAL.setEnabled(True)
            ntlist, ok = self.lineEditNTLAL.text().toInt()
            if ntlist < 1:
                ntlist = 1
                self.model.setListingFrequency(ntlist)
                self.lineEditNTLAL.setText(QString(str(ntlist)))


    @pyqtSignature("const QString&")
    def slotNTLAL(self, text):
        """
        Input NTLAL.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            period, ok = text.toInt()
            self.model.setListingFrequency(period)


    @pyqtSignature("")
    def slotIVISV1(self):
        """
        Input IVISV1.
        """
        if self.checkBoxIVISV1.isChecked():
            self.model.setFluidVelocityStatus("on")
        else:
            self.model.setFluidVelocityStatus("off")


    @pyqtSignature("")
    def slotIVISV2(self):
        """
        Input IVISV2.
        """
        if self.checkBoxIVISV2.isChecked():
            self.model.setParticlesVelocityStatus("on")
        else:
            self.model.setParticlesVelocityStatus("off")


    @pyqtSignature("")
    def slotIVISTP(self):
        """
        Input IVISTP.
        """
        if self.checkBoxIVISTP.isChecked():
            self.model.setResidentTimeStatus("on")
        else:
            self.model.setResidentTimeStatus("off")


    @pyqtSignature("")
    def slotIVISDM(self):
        """
        Input IVISDM.
        """
        if self.checkBoxIVISDM.isChecked():
            self.model.setParticleDiameterStatus("on")
        else:
            self.model.setParticleDiameterStatus("off")


    @pyqtSignature("")
    def slotIVISTE(self):
        """
        Input IVISTE.
        """
        if self.checkBoxIVISTE.isChecked():
            self.model.setParticleTemperatureStatus("on")
        else:
            self.model.setParticleTemperatureStatus("off")


    @pyqtSignature("")
    def slotIVISMP(self):
        """
        Input IVISMP.
        """
        if self.checkBoxIVISMP.isChecked():
            self.model.setParticleMassStatus("on")
        else:
            self.model.setParticleMassStatus("off")


    @pyqtSignature("")
    def slotIVISHP(self):
        """
        Input IVISHP.
        """
        if self.checkBoxIVISHP.isChecked():
            self.model.setCoalParticleTemperatureStatus("on")
        else:
            self.model.setCoalParticleTemperatureStatus("off")


    @pyqtSignature("")
    def slotIVISDK(self):
        """
        Input IVISDK.
        """
        if self.checkBoxIVISDK.isChecked():
            self.model.setCoalParticleDiameterStatus("on")
        else:
            self.model.setCoalParticleDiameterStatus("off")


    @pyqtSignature("")
    def slotIVISCH(self):
        """
        Input IVISCH.
        """
        if self.checkBoxIVISCH.isChecked():
            self.model.setCoalParticleMassStatus("on")
        else:
            self.model.setCoalParticleMassStatus("off")


    @pyqtSignature("")
    def slotIVISCK(self):
        """
        Input IVISCK.
        """
        if self.checkBoxIVISCK.isChecked():
            self.model.setCokeParticleMassStatus("on")
        else:
            self.model.setCokeParticleMassStatus("off")


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
