# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module contains the following classes:
- AtmosphericFlowsView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging

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
from code_saturne.gui.base.QtPage import DoubleValidator, from_qvariant
from code_saturne.model.AtmosphericFlowsModel  import AtmosphericFlowsModel
from code_saturne.gui.case.AtmosphericFlowsForm   import Ui_AtmosphericFlowsForm

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("AtmosphericFlowsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class AtmosphericFlowsView(QWidget, Ui_AtmosphericFlowsForm):
    """
    Main class
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        # Setup base classes
        QWidget.__init__(self, parent)
        Ui_AtmosphericFlowsForm.__init__(self)
        self.setupUi(self)

        # create model
        model = AtmosphericFlowsModel(case)
        self.__model = model
        self.case = case

        self.case.undoStopGlobal()

        # Define connection
        self.groupBoxMeteoData.clicked[bool].connect(self.__slotGroupBoxMeteoData)
        self.pushButtonMeteoData.pressed.connect(self.__slotSearchMeteoData)

        self.groupBoxLargeScaleMeteo.clicked[bool].connect(self.__slotGroupBoxLargeScaleMeteo)
        #TODO not yet connected
        #self.groupBoxActChemistry.clicked[bool].connect(self.__slotGroupBoxActChemistry)
        self.__slotGroupBoxActChemistry(False)
        self.groupBoxActChemistry.setEnabled(False)

        self.comboBoxUstarOrdLMO.currentIndexChanged[int].connect(self.__slotComboBoxUstarOrDlmo)
        self.comboBoxUrefOrdLMO.currentIndexChanged[int].connect(self.__slotComboBoxUrefOrDlmo)

        self.lineEditLongCenter.textChanged[str].connect(self.slotLongitude)
        self.lineEditLatCenter.textChanged[str].connect(self.slotLatitude)
        self.spinBoxDomainOrientation.valueChanged[int].connect(self.slotDomainOrientation)
        self.dateTimeEdit.dateTimeChanged[QDateTime].connect(self.__slotDateTime)
        self.lineEditLargeScaleRoughness.textChanged[str].connect(self.slotMeteoZ0)
        self.lineEditPressureSeaLevel.textChanged[str].connect(self.slotMeteoPsea)
        self.lineEditTemperatureSeaLevel.textChanged[str].connect(self.slotMeteoT0)
        self.spinBoxWindDir.valueChanged[int].connect(self.slotWindDir)
        self.lineEditUstarOrDlmo.textChanged[str].connect(self.slotMeteoUstarOrDlmo)
        self.lineEditUrefOrDlmo.textChanged[str].connect(self.slotMeteoUrefOrDlmo)

        self.lineEditZref.textChanged[str].connect(self.slotMeteoZref)

        # Initialize the widgets in groupBoxMeteoData
        isMeteoDataChecked = model.getMeteoDataStatus() == 'on'
        self.groupBoxMeteoData.setChecked(isMeteoDataChecked)
        self.labelMeteoFile.setText(str(self.__model.getMeteoDataFileName()))

        # Validate lineEdits in groupBoxLargeScaleMeteo
        validatorLongitude = DoubleValidator(self.lineEditLongCenter)
        validatorLatitude = DoubleValidator(self.lineEditLatCenter)

        validatorZ0 = DoubleValidator(self.lineEditLargeScaleRoughness)
        validatorPsea = DoubleValidator(self.lineEditPressureSeaLevel)
        validatorT0 = DoubleValidator(self.lineEditTemperatureSeaLevel)

        validatorZref = DoubleValidator(self.lineEditZref)
        validatorUstarOrDlmo = DoubleValidator(self.lineEditUstarOrDlmo)
        validatorUrefOrDlmo = DoubleValidator(self.lineEditUrefOrDlmo)

        self.lineEditLongCenter.setValidator(validatorLongitude)
        self.lineEditLatCenter.setValidator(validatorLatitude)

        self.lineEditLargeScaleRoughness.setValidator(validatorZ0)
        self.lineEditPressureSeaLevel.setValidator(validatorPsea)
        self.lineEditTemperatureSeaLevel.setValidator(validatorT0)

        self.lineEditZref.setValidator(validatorZref)
        self.lineEditUstarOrDlmo.setValidator(validatorUstarOrDlmo)
        self.lineEditUrefOrDlmo.setValidator(validatorUrefOrDlmo)

        # Initialize the widgets in groupBoxLargeScaleMeteo
        isLargeScaleMeteoChecked = model.getLargeScaleMeteoStatus() == 'on'
        self.groupBoxLargeScaleMeteo.setChecked(isLargeScaleMeteoChecked)

        tmpVar = model.getLongitude();
        self.lineEditLongCenter.setText(str(tmpVar))
        tmpVar = model.getLatitude();
        self.lineEditLatCenter.setText(str(tmpVar))

        tmpVar = model.getMeteoZ0();
        self.lineEditLargeScaleRoughness.setText(str(tmpVar))

        tmpVar = model.getMeteoPsea();
        self.lineEditPressureSeaLevel.setText(str(tmpVar))

        tmpVar = model.getMeteoT0();
        self.lineEditTemperatureSeaLevel.setText(str(tmpVar))

        tmpVar = model.getMeteoZref();
        self.lineEditZref.setText(str(tmpVar))

        tmpVar = model.getMeteoUstar()
        if tmpVar >= 0:
            self.comboBoxUstarOrdLMO.setCurrentIndex(0)
        else:
            # Get value before overwritting it
            tmpVar = model.getMeteoDlmo()
            self.comboBoxUstarOrdLMO.setCurrentIndex(1)

        self.lineEditUstarOrDlmo.setText(str(tmpVar))

        tmpVar = model.getMeteoUref()
        if tmpVar >= 0:
            self.comboBoxUrefOrdLMO.setCurrentIndex(0)
        else:
            # Get value before overwritting it
            tmpVar = model.getMeteoDlmo()
            self.comboBoxUrefOrdLMO.setCurrentIndex(1)

        self.lineEditUrefOrDlmo.setText(str(tmpVar))

        tmpVar=model.getDomainOrientation()
        self.spinBoxDomainOrientation.setValue(int(tmpVar))

        tmpVar=model.getWindDir()
        self.spinBoxWindDir.setValue(int(tmpVar))

        startTimeStr = model.getStartTime()
        startTime = QDateTime.fromString(startTimeStr, "yyyy-MM-dd HH:mm:ss")
        self.dateTimeEdit.setDateTime(startTime)

        # Initialize the widgets in groupBoxActChemistry
        isChemistryChecked = model.getChemistryStatus() == 'on'
        self.groupBoxActChemistry.setChecked(isChemistryChecked)

        self.case.undoStartGlobal()


    #--------------- Functions for the groupBox LargeScalaMeteData--------------
    @pyqtSlot(QDateTime)
    def __slotDateTime(self, startTime):
            self.__model.setStartTime(startTime.toPyDateTime())

    @pyqtSlot(int)
    def __slotComboBoxUrefOrDlmo(self, indCurrent):
        text = self.lineEditUrefOrDlmo.text()
        val = from_qvariant(text, float)
        if indCurrent==0:
            self.labelDimRefVel.setText("m/s")
            self.comboBoxUrefOrdLMO.model().item(1).setEnabled(True)
            self.__model.setMeteoUref(val)
            self.labelReferenceHeight.setEnabled(True)
            self.lineEditZref.setEnabled(True)
            self.labelDimZref.setEnabled(True)
        elif indCurrent==1:
            self.labelDimRefVel.setText("m<sup>-1</sup>")
            self.comboBoxUrefOrdLMO.SelectedIndex = 0;
            self.comboBoxUrefOrdLMO.model().item(1).setEnabled(False)
            self.__model.setMeteoUref(-1.)
            self.__model.setMeteoDlmo(val)
            self.labelReferenceHeight.setEnabled(False)
            self.lineEditZref.setEnabled(False)
            self.labelDimZref.setEnabled(False)

    @pyqtSlot(int)
    def __slotComboBoxUstarOrDlmo(self, indCurrent):
        text = self.lineEditUstarOrDlmo.text()
        val = from_qvariant(text, float)
        if indCurrent==0:
            self.labelDimZRef.setText("m/s")
            self.comboBoxUstarOrdLMO.model().item(1).setEnabled(True)
            self.__model.setMeteoUstar(val)
        elif indCurrent==1:
            self.labelDimZRef.setText("m<sup>-1</sup>")
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;
            self.comboBoxUstarOrdLMO.model().item(1).setEnabled(False)
            self.__model.setMeteoDlmo(val)
            self.__model.setMeteoUstar(-1.)

    @pyqtSlot(bool)
    def __slotGroupBoxLargeScaleMeteo(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'
            self.__slotGroupBoxMeteoData(False)

        self.groupBoxLargeScaleMeteo.setChecked(checked)
        self.__model.setLargeScaleMeteoStatus(status)

    @pyqtSlot(bool)
    def __slotApplyLargeScaleMeteo(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'
            self.__slotGroupBoxMeteoData(False)

        self.groupBoxLargeScaleMeteo.setChecked(checked)
        self.__model.setLargeScaleMeteoStatus(status)

    @pyqtSlot(str)
    def slotLongitude(self, text):
        if self.lineEditLongCenter.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setLongitude(val)

    @pyqtSlot(str)
    def slotLatitude(self, text):
        if self.lineEditLatCenter.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setLatitude(val)

    @pyqtSlot(int)
    def slotDomainOrientation(self, text):
        val = self.spinBoxDomainOrientation.value()
        self.__model.setDomainOrientation(val)

    @pyqtSlot(str)
    def slotMeteoZ0(self, text):
        if self.lineEditLargeScaleRoughness.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoZ0(val)

    @pyqtSlot(str)
    def slotMeteoPsea(self, text):
        if self.lineEditPressureSeaLevel.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoPsea(val)

    @pyqtSlot(str)
    def slotMeteoT0(self, text):
        if self.lineEditTemperatureSeaLevel.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoT0(val)

    @pyqtSlot(int)
    def slotWindDir(self, text):
        val = self.spinBoxWindDir.value()
        self.__model.setWindDir(val)

    @pyqtSlot(str)
    def slotMeteoZref(self, text):
        if self.lineEditZref.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoZref(val)

    @pyqtSlot(str)
    def slotMeteoUstarOrDlmo(self, text):
        if self.lineEditUstarOrDlmo.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            if self.comboBoxUstarOrdLMO.currentIndex() == 0:
                self.__model.setMeteoUstar(val)
            elif self.comboBoxUstarOrdLMO.currentIndex() == 1:
                self.__model.setMeteoDlmo(val)
                self.__model.setMeteoUstar(-1.)

    @pyqtSlot(str)
    def slotMeteoUrefOrDlmo(self, text):
        if self.lineEditUrefOrDlmo.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            if self.comboBoxUrefOrdLMO.currentIndex() == 0:
                self.__model.setMeteoUref(val)
            elif self.comboBoxUrefOrdLMO.currentIndex() == 1:
                self.__model.setMeteoDlmo(val)
                self.__model.setMeteoUref(-1.)

    #--------------- Functions for the groupBox Activate Chemistry--------------
    @pyqtSlot(bool)
    def __slotGroupBoxActChemistry(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        #TODO not yet activated
        self.groupBoxActChemistry.setEnabled(False)
        self.groupBoxActChemistry.setChecked(checked)
        self.__model.setChemistryStatus(status)

    #--------------- Functions for the groupBox  MeteoDataFile-----------------

    @pyqtSlot(bool)
    def __slotGroupBoxMeteoData(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        self.groupBoxMeteoData.setChecked(checked)
        self.__model.setMeteoDataStatus(status)
        if checked:
            self.__slotGroupBoxLargeScaleMeteo(False)


    @pyqtSlot()
    def __slotSearchMeteoData(self):
        """
        Select a meteorological file of data
        """
        data = self.case['data_path']
        title = self.tr("Meteorological file of data.")
        filetypes = self.tr("Meteo data (*meteo*);;All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, data, filetypes)[0]
        file = str(file)
        if not file:
            return
        file = os.path.basename(file)
        if not data:
            data = os.getcwd()
        if file not in os.listdir(data):
            title = self.tr("Warning")
            msg   = self.tr("This selected file is not in the DATA directory")
            QMessageBox.information(self, title, msg)
        else:
            self.labelMeteoFile.setText(str(file))
            self.__model.setMeteoDataFileName(file)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
