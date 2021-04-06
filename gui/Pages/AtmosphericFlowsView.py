# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2021 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, from_qvariant
from code_saturne.Pages.AtmosphericFlowsForm   import Ui_AtmosphericFlowsForm
from code_saturne.model.AtmosphericFlowsModel  import AtmosphericFlowsModel

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
        self.groupBoxActChemistry.clicked[bool].connect(self.__slotGroupBoxActChemistry)
        self.comboBoxUstarOrdLMO.currentIndexChanged[int].connect(self.__slotComboBoxUstarOrDlmo)
        self.comboBoxUrefOrdLMO.currentIndexChanged[int].connect(self.__slotComboBoxUrefOrDlmo)

        # Initialize the widgets in groupBoxMeteoData
        isMeteoDataChecked = model.getMeteoDataStatus() == 'on'
        self.groupBoxMeteoData.setChecked(isMeteoDataChecked)
        self.labelMeteoFile.setText(str(self.__model.getMeteoDataFileName()))

        # Validate lineEdits in groupBoxLargeScaleMeteo
        validatatorLongitude = DoubleValidator(self.lineEditLongCenter)
        validatatorLatitude = DoubleValidator(self.lineEditLatCenter)
        validatatorZ0 = DoubleValidator(self.lineEditLargeScaleRoughness)
        validatatorPsea = DoubleValidator(self.lineEditPressureSeaLevel)
        validatatorZref = DoubleValidator(self.lineEditZref)
        validatatorUstarOrDlmo = DoubleValidator(self.lineEditUstarOrDlmo)
        validatatorUrefOrDlmo = DoubleValidator(self.lineEditUrefOrDlmo)

        self.lineEditLongCenter.setValidator(validatatorLongitude)
        self.lineEditLatCenter.setValidator(validatatorLatitude)
        self.lineEditLargeScaleRoughness.setValidator(validatatorZ0)
        self.lineEditPressureSeaLevel.setValidator(validatatorPsea)
        self.lineEditZref.setValidator(validatatorZref)
        self.lineEditUstarOrDlmo.setValidator(validatatorUstarOrDlmo)
        self.lineEditUrefOrDlmo.setValidator(validatatorUrefOrDlmo)

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

        tmpVar = model.getMeteoZref();
        self.lineEditZref.setText(str(tmpVar))

        meteoUstar = float(model.getMeteoUstar());
        meteoDlmo = float(model.getMeteoDlmo());
        meteoUref = float(model.getMeteoUref());
        if (meteoUstar>0.0) and (meteoUref>0.0):
            self.lineEditUstarOrDlmo.setText(str(meteoUstar))
            self.lineEditUrefOrDlmo.setText(str(meteoUref))
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;
        elif (meteoUstar>0.0) and (meteoUref<=0.0):
            self.lineEditUstarOrDlmo.setText(str(meteoUstar))
            self.lineEditUrefOrDlmo.setText(str(meteoDlmo))
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;
            self.comboBoxUstarOrdLMO.SelectedIndex = 1;
        else:
            self.lineEditUstarOrDlmo.setText(str(meteoDlmo))
            self.lineEditUrefOrDlmo.setText(str(meteoUref))
            self.comboBoxUstarOrdLMO.SelectedIndex = 1;
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;

        tmpVar=model.getDomainOrientation()
        self.spinBoxDomainOrientation.setValue(int(tmpVar))

        tmpVar=model.getWindDir()
        self.spinBoxWindDir.setValue(int(tmpVar))


        startTime = model.getStartTime()
        self.dateTimeEdit.setDateTime(startTime)

        # Initialize the widgets in groupBoxActChemistry
        isChemistryChecked = model.getChemistryStatus() == 'on'
        self.groupBoxActChemistry.setChecked(isChemistryChecked)


        self.case.undoStartGlobal()


    #--------------- Fuctions for the groupBox LargeScalaMeteData--------------

    @pyqtSlot(int)
    def __slotComboBoxUrefOrDlmo(self, indCurrent):
        if indCurrent==0:
            self.labelDimRefVel.setText("m/s")
            self.comboBoxUstarOrdLMO.model().item(1).setEnabled(True)
            self.labelReferenceHeight.setEnabled(True)
            self.lineEditZref.setEnabled(True)
            self.labelDimZref.setEnabled(True)
        elif indCurrent==1:
            self.labelDimRefVel.setText("m<sup>-1</sup>")
            self.comboBoxUstarOrdLMO.SelectedIndex = 0;
            self.comboBoxUstarOrdLMO.model().item(1).setEnabled(False)
            self.labelReferenceHeight.setEnabled(False)
            self.lineEditZref.setEnabled(False)
            self.labelDimZref.setEnabled(False)

    @pyqtSlot(int)
    def __slotComboBoxUstarOrDlmo(self, indCurrent):
        if indCurrent==0:
            self.labelDimZRef.setText("m/s")
            self.comboBoxUrefOrdLMO.model().item(1).setEnabled(True)
        elif indCurrent==1:
            self.labelDimZRef.setText("m<sup>-1</sup>")
            self.comboBoxUrefOrdLMO.SelectedIndex = 0;
            self.comboBoxUrefOrdLMO.model().item(1).setEnabled(False)

    @pyqtSlot(bool)
    def __slotGroupBoxLargeScaleMeteo(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        self.groupBoxLargeScaleMeteo.setChecked(checked)
        self.__model.setLargeScaleMeteoStatus(status)
        if checked:
            self.__slotGroupBoxMeteoData(False)

    @pyqtSlot(bool)
    def __slotApplyLargeScaleMeteo(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        self.groupBoxLargeScaleMeteo.setChecked(checked)
        self.__model.setLargeScaleMeteoStatus(status)
        if checked:
            self.__slotGroupBoxMeteoData(False)

    @pyqtSlot(str)
    def slotLatitude(self, text):
        if self.lineEditLatCenter.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setLatitude(val)

    @pyqtSlot(str)
    def slotLongitude(self, text):
        if self.lineEditLongCenter.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setLongitude(val)

    @pyqtSlot(str)
    def slotMeteZ0(self, text):
        if self.lineEditLargeScaleRoughness.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoZ0(val)

    @pyqtSlot(str)
    def slotMeteoZref(self, text):
        if self.lineEditZref.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.__model.setMeteoZref(val)


    @pyqtSlot(str)
    def slotMeteoUstarOrDlmo(self, text):
        if self.lineEditLongitude.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            if self.comboBoxUstarOrdLMO.currentIndex() == 0:
                self.__model.setMeteoUstar(val)
                self.__model.setMeteoDlmo("0.0")
            elif self.comboBoxUstarOrdLMO.currentIndex() == 1:
                self.__model.setMeteoDlmo(val)
                self.__model.setMeteoUstar("0.0")

    @pyqtSlot(str)
    def slotMeteoUrefOrDlmo(self, text):
        if self.lineEditUrefOrDlmo.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            if self.comboBoxUrefOrdLMO.currentIndex() == 0:
                self.__model.setMeteoUref(val)
                self.__model.setMeteoDlmo("0.0")
            elif self.comboBoxUrefOrdLMO.currentIndex() == 1:
                self.__model.setMeteoDlmo(val)
                self.__model.setMeteoUref("0.0")

    #--------------- Fuctions for the groupBox Activate Chemistry--------------
    @pyqtSlot(bool)
    def __slotGroupBoxActChemistry(self, checked):
        """
        Called when groupBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        self.groupBoxActChemistry.setChecked(checked)

    #--------------- Fuctions for the groupBox  MeteoDataFile-----------------

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
