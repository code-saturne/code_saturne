# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
- StandardItemModelVolumicNames
- StandardItemModelBoundariesNames
- LagrangianStatisticsView
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

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.LagrangianStatisticsForm import Ui_LagrangianStatisticsForm
from code_saturne.model.LagrangianStatisticsModel import LagrangianStatisticsModel
from code_saturne.model.LagrangianModel import LagrangianModel


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------


logging.basicConfig()
log = logging.getLogger("LagrangianStatisticsView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# StandarItemModel for volumic variables names
#-------------------------------------------------------------------------------


class StandardItemModelVolumicNames(QStandardItemModel):
    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Variable Name (Mean value)"),
                        self.tr("Variable Name (Variance)"),
                        self.tr("Active")]
        self.setColumnCount(len(self.headers))
        self.model = model
        self.initData()


    def initData(self):

        self.dataVolumicNames = []
        vnames = self.model.getVariablesNamesVolume()
        for vname in vnames:
            if vname[0:5] == "mean_":
                labelv = "var_" + vname[5:]
                compute = self.model.getComputeVolStatusFromName(vname)
                line = [vname, labelv, compute]
            elif vname[0:4] == "var_":
                compute = self.model.getComputeVolStatusFromName(vname)
                continue # group with mean
            else:
                labelv = ""
                compute = self.model.getComputeVolStatusFromName(vname)
                line = [vname, labelv, compute]

            row = self.rowCount()
            self.setRowCount(row+1)
            self.dataVolumicNames.append(line)


    def data(self, index, role):

        if not index.isValid():
            return

        # ToolTips
        if role == Qt.ToolTipRole:
            if index.column() in [0, 1]:
                return self.tr("field label base")

        # Display
        if role == Qt.DisplayRole:
            if index.column() in [0,1]:
                return self.dataVolumicNames[index.row()][index.column()]

        # CheckState
        elif role == Qt.CheckStateRole:
            if index.column() == 2:
                if self.dataVolumicNames[index.row()][index.column()] == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        return


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == [0,1]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsUserCheckable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return


    def setData(self, index, value, role):
        #
        if index.column() == 0:
            self.dataVolumicNames[index.row()][index.column()] = \
                        str(from_qvariant(value, to_text_string))
            vname = self.dataVolumicNames[index.row()][0]

        elif index.column() == 1:
            labelv = str(from_qvariant(value, to_text_string))
            self.dataVolumicNames[index.row()][index.column()] = labelv
            name = self.dataVolumicNames[index.row()][0]

        elif index.column() == 2:
            v = from_qvariant(value, int)
            if v == Qt.Unchecked:
                status = "off"
                self.dataVolumicNames[index.row()][index.column()] = "off"
            else:
                status = "on"
                self.dataVolumicNames[index.row()][index.column()] = "on"

            vname = self.dataVolumicNames[index.row()][0]
            self.model.setComputeVolStatusFromName(vname, status)
            vname = self.dataVolumicNames[index.row()][1]
            if vname:
                self.model.setComputeVolStatusFromName(vname, status)

        self.dataChanged.emit(index, index)
        return True


#-------------------------------------------------------------------------------
# StandarItemModel for boundaries variables names
#-------------------------------------------------------------------------------


class StandardItemModelBoundariesNames(QStandardItemModel):
    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Variable Name (Mean value)"),
                        self.tr("Variable Name (Variance)"),
                        self.tr("Active")]
        self.setColumnCount(len(self.headers))
        self.model = model
        self.initData()


    def initData(self):

        self.dataBoundariesNames = []
        vnames = self.model.getVariablesNamesBoundary()
        for vname in vnames:
            if vname[0:5] == "mean_":
                labelv = "var_" + vname[5:]
                compute = self.model.getComputeBoundaryStatusFromName(vname)
                line = [vname, labelv, compute]
            elif vname[0:4] == "var_":
                compute = self.model.getComputeBoundaryStatusFromName(vname)
                continue # group with mean
            else:
                labelv = ""
                compute = self.model.getComputeBoundaryStatusFromName(vname)
                line = [vname, labelv, compute]

            row = self.rowCount()
            self.setRowCount(row+1)
            self.dataBoundariesNames.append(line)


    def data(self, index, role):

        if not index.isValid():
            return

        # ToolTips
        if role == Qt.ToolTipRole:
            if index.column() in [0, 1]:
                return self.tr("field label base")

        # Display
        if role == Qt.DisplayRole:
            if index.column() in [0,1]:
                return self.dataBoundariesNames[index.row()][index.column()]

        # CheckState
        elif role == Qt.CheckStateRole:
            if index.column() ==2:
                if self.dataBoundariesNames[index.row()][index.column()] == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        return


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return


    def setData(self, index, value, role):

        if index.column() == 1:
            v = from_qvariant(value, int)
            if v == Qt.Unchecked:
                status = "off"
                self.dataBoundariesNames[index.row()][index.column()] = "off"
            else:
                status = "on"
                self.dataBoundariesNames[index.row()][index.column()] = "on"

            vname = self.dataBoundariesNames[index.row()][0]
            self.model.setComputeBoundaryStatusFromName(vname, status)
            vname = self.dataBoundariesNames[index.row()][1]
            if vname:
                self.model.setComputeVolStatusFromName(vname, status)


        self.dataChanged.emit(index, index)
        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class LagrangianStatisticsView(QWidget, Ui_LagrangianStatisticsForm):
    """
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LagrangianStatisticsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = LagrangianStatisticsModel(self.case)

        self.checkBoxISUIST.clicked.connect(self.slotISUIST)
        self.lineEditNBCLST.textChanged[str].connect(self.slotNBCLST)

        self.groupBoxISTALA.clicked.connect(self.slotISTALA)
        self.lineEditIDSTNT.editingFinished.connect(self.slotIDSTNT)
        self.lineEditNSTIST.editingFinished.connect(self.slotNSTIST)

        self.lineEditSEUIL.textChanged[str].connect(self.slotSEUIL)

        self.groupBoxIENSI3.clicked.connect(self.slotIENSI3)

        validatorNBCLST = IntValidator(self.lineEditNBCLST, min=0) # max=100

        validatorIDSTNT = IntValidator(self.lineEditIDSTNT, min=0)
        validatorIDSTNT.setExclusiveMin(True)
        validatorNSTIST = IntValidator(self.lineEditNSTIST, min=0)
        validatorSEUIL = DoubleValidator(self.lineEditSEUIL, min=0.)

        self.lineEditNBCLST.setValidator(validatorNBCLST)
        self.lineEditIDSTNT.setValidator(validatorIDSTNT)
        self.lineEditNSTIST.setValidator(validatorNSTIST)
        self.lineEditSEUIL.setValidator(validatorSEUIL)

        # initialize Widgets
        # FIXME
        # test if restart lagrangian is on
##         mdl.lagr = LagrangianModel()
##         is_restart = mdl_lagr.getRestart()
##         if is_restart == "off":
##             self.lineEditISUIST.setEnabled(False)
##             self.checkBoxISUIST.setEnabled(False)
        status = self.model.getRestartStatisticsStatus()
        if status == "on":
            self.checkBoxISUIST.setChecked(True)
        else:
            self.checkBoxISUIST.setChecked(False)

        nclust = self.model.getGroupOfParticlesValue()
        self.lineEditNBCLST.setText(str(nclust))

        it = self.model.getIterationStart()
        self.lineEditIDSTNT.setText(str(it))

        it = self.model.getIterSteadyStart()
        self.lineEditNSTIST.setText(str(it))

        seuil = self.model.getThresholdValue()
        self.lineEditSEUIL.setText(str(seuil))

        # volume
        status = self.model.getVolumeStatisticsStatus()
        if status == "on":
            self.groupBoxISTALA.setChecked(True)
        else:
            self.groupBoxISTALA.setChecked(False)
        self.slotISTALA()

        # boundary
        status = self.model.getBoundaryStatisticsStatus()
        if status == "on":
            self.groupBoxIENSI3.setChecked(True)
        else:
            self.groupBoxIENSI3.setChecked(False)
        self.slotIENSI3()

        self.case.undoStartGlobal()


    def _initVolumicNames(self):
        """
        Initialize names for volumic statistics.
        """
        self.modelVolumicNames = StandardItemModelVolumicNames(self.model)

        self.tableViewVolumicNames.setModel(self.modelVolumicNames)
        self.tableViewVolumicNames.setAlternatingRowColors(True)
        self.tableViewVolumicNames.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewVolumicNames.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewVolumicNames.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewVolumicNames.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewVolumicNames.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)


    def _initBoundariesNames(self):
        """
        Initialize names for volumic statistics.
        """
        self.modelBoundariesNames = StandardItemModelBoundariesNames(self.model)

        self.tableViewBoundariesNames.setModel(self.modelBoundariesNames)
        self.tableViewBoundariesNames.setAlternatingRowColors(True)
        self.tableViewBoundariesNames.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewBoundariesNames.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewBoundariesNames.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewBoundariesNames.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewBoundariesNames.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)


    @pyqtSlot()
    def slotISUIST(self):
        """
        Input ISUIST.
        """
        if self.checkBoxISUIST.isChecked():
            status = "on"
        else:
            status = "off"
        self.model.setRestartStatisticsStatus(status)


    @pyqtSlot(str)
    def slotNBCLST(self, text):
        """
        Input NBCLST.
        """
        if self.lineEditNBCLST.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setGroupOfParticlesValue(value)

    @pyqtSlot()
    def slotISTALA(self):
        """
        Input ISTALA.
        """
        if self.groupBoxISTALA.isChecked():

            self.model.setVolumeStatisticsStatus("on")
            self._initVolumicNames()

        else:

            self.model.setVolumeStatisticsStatus("off")
            if hasattr(self, "modelVolumicNames"):
                del self.modelVolumicNames


    @pyqtSlot()
    def slotIDSTNT(self):
        """
        Input IDSTNT.
        """
        if self.lineEditIDSTNT.validator().state == QValidator.Acceptable:
            text = self.lineEditIDSTNT.text()
            value = from_qvariant(text, int)
            valnds =  self.model.getIterSteadyStart()

            if value > valnds:
                self.lineEditNSTIST.setText(str(value))
                self.model.setIterSteadyStart(value)
            else:
                valndsl = from_qvariant(self.lineEditNSTIST.text(), int)
                self.model.setIterSteadyStart(valndsl)

            self.model.setIterationStart(value)


    @pyqtSlot()
    def slotNSTIST(self):
        """
        Input NSTIST.
        """
        if self.lineEditNSTIST.validator().state == QValidator.Acceptable:
            text = self.lineEditNSTIST.text()
            value = from_qvariant(text, int)
            valids =  self.model.getIterationStart()

            if value < valids:
                self.lineEditIDSTNT.setText(str(value))
                self.model.setIterationStart(value)
            else:
                validsl = from_qvariant(self.lineEditIDSTNT.text(), int)
                self.model.setIterationStart(validsl)

            self.model.setIterSteadyStart(value)


    @pyqtSlot(str)
    def slotSEUIL(self, text):
        """
        Input SEUIL.
        """
        if self.lineEditSEUIL.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setThresholdValue(value)


    @pyqtSlot()
    def slotIENSI3(self):
        """
        Input IENSI3.
        """
        if self.groupBoxIENSI3.isChecked():

            self.model.setBoundaryStatisticsStatus("on")
            self._initBoundariesNames()

        else:

            self.model.setBoundaryStatisticsStatus("off")
            if hasattr(self, "modelBoundariesNames"):
                del self.modelBoundariesNames


    @pyqtSlot(str)
    def slotNSTBOR(self, text):
        """
        Input NSTBOR.
        """
        if self.lineEditNSTBOR.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setIterationStartBoundary(value)


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
