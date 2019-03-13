# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the HeadLosses model data management.

This module contains the following classes:
- HeadLosses
- HeadLossesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, logging

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
from code_saturne.Base.QtPage import DoubleValidator, ComboModel
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.HeadLossesForm import Ui_HeadLossesForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.HeadLossesModel import HeadLossesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("HeadLossesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class to display Head Losses Zones in a QTreeView
#-------------------------------------------------------------------------------


class StandardItemModelHeadLosses(QStandardItemModel):
    def __init__(self):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"), self.tr("Zone"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))
        self.dataHeadLossesZones = []


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            return self.dataHeadLossesZones[index.row()][index.column()]
        return None

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariat(self.headers[section])
        return None


    def setData(self, index, value, role):
        self.dataChanged.emit(index, index)
        return True


    def insertItem(self, label, name, local):
        line = [label, name, local]
        self.dataHeadLossesZones.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataHeadLossesZones[row]


#-------------------------------------------------------------------------------
# Main view class
#-------------------------------------------------------------------------------

class HeadLossesView(QWidget, Ui_HeadLossesForm):

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_HeadLossesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        # Create the Page layout.

        # Model and QTreeView for Head Losses
        self.modelHeadLosses = StandardItemModelHeadLosses()
        self.treeView.setModel(self.modelHeadLosses)

        # Connections
        self.treeView.clicked[QModelIndex].connect(self.slotSelectHeadLossesZones)
        self.groupBox_3.clicked[bool].connect(self.slotTransfoMatrix)

        self.lineEdit.textChanged[str].connect(self.slotKxx)
        self.lineEdit_2.textChanged[str].connect(self.slotKyy)
        self.lineEdit_3.textChanged[str].connect(self.slotKzz)

        self.lineEdit_4.textChanged[str].connect(self.slotA11)
        self.lineEdit_5.textChanged[str].connect(self.slotA12)
        self.lineEdit_6.textChanged[str].connect(self.slotA13)
        self.lineEdit_8.textChanged[str].connect(self.slotA21)
        self.lineEdit_9.textChanged[str].connect(self.slotA22)
        self.lineEdit_7.textChanged[str].connect(self.slotA23)
        self.lineEdit_11.textChanged[str].connect(self.slotA31)
        self.lineEdit_12.textChanged[str].connect(self.slotA32)
        self.lineEdit_10.textChanged[str].connect(self.slotA33)

        # Validators

        validator = DoubleValidator(self.lineEdit, min=0.0)
        validator_2 = DoubleValidator(self.lineEdit_2, min=0.0)
        validator_3= DoubleValidator(self.lineEdit_3, min=0.0)

        validator_4 = DoubleValidator(self.lineEdit_4)
        validator_5= DoubleValidator(self.lineEdit_5)
        validator_6= DoubleValidator(self.lineEdit_6)
        validator_8 = DoubleValidator(self.lineEdit_8)
        validator_9 = DoubleValidator(self.lineEdit_9)
        validator_7 = DoubleValidator(self.lineEdit_7)
        validator_11 = DoubleValidator(self.lineEdit_11)
        validator_12 = DoubleValidator(self.lineEdit_12)
        validator_10 = DoubleValidator(self.lineEdit_10)

        # Apply validators

        self.lineEdit.setValidator(validator)
        self.lineEdit_2.setValidator(validator_2)
        self.lineEdit_3.setValidator(validator_3)

        self.lineEdit_4.setValidator(validator_4)
        self.lineEdit_5.setValidator(validator_5)
        self.lineEdit_6.setValidator(validator_6)
        self.lineEdit_8.setValidator(validator_8)
        self.lineEdit_9.setValidator(validator_9)
        self.lineEdit_7.setValidator(validator_7)
        self.lineEdit_11.setValidator(validator_11)
        self.lineEdit_12.setValidator(validator_12)
        self.lineEdit_10.setValidator(validator_10)

        # Initialize Widgets

        self.entriesNumber = 0
        d = HeadLossesModel(self.case).getNameAndLocalizationZone()
        liste=[]
        liste=list(d.items())
        t=[]
        for t in liste :
            NamLoc=t[1]
            Lab=t[0 ]
            self.modelHeadLosses.insertItem(Lab, NamLoc[0],NamLoc[1])
        self.forgetStandardWindows()

        self.case.undoStartGlobal()


    @pyqtSlot("QModelIndex")
    def slotSelectHeadLossesZones(self, index):
        model = HeadLossesModel(self.case)
        label, name, local = self.modelHeadLosses.getItem(index.row())

        if hasattr(self, "modelScalars"): del self.modelScalars
        log.debug("slotSelectHeadLossesZones label %s " % label )
        self.groupBoxDef.show()
        self.groupBox_3.show()
        kxx,kyy,kzz = model.getKCoefficients(name)
        self.lineEdit.setText(str(kxx))
        self.lineEdit_2.setText(str(kyy))
        self.lineEdit_3.setText(str(kzz))

        if model.getMatrixChoice(name,'choice') == 'on':
            self.groupBox_3.setChecked(True)
            checked = True
        else:
            self.groupBox_3.setChecked(False)
            checked = False

        self.slotTransfoMatrix(checked)


    def forgetStandardWindows(self):
        """
        For forget standard windows
        """
        self.groupBoxDef.hide()
        self.groupBox_3.hide()



    @pyqtSlot(str)
    def slotKxx(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setCoefficient(name,'kxx',value )

    @pyqtSlot(str)
    def slotKyy(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_2.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setCoefficient(name,'kyy',value )

    @pyqtSlot(str)
    def slotKzz(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_3.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setCoefficient(name,'kzz',value )

    @pyqtSlot(bool)
    def slotTransfoMatrix(self,  checked):
        self.groupBox_3.setFlat(not checked)
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)

            if checked:
                model.setMatrixChoice(name,'choice','on')
                self.groupBox_3.setChecked(True)
                self.frameTransfo.show()
                a11, a12, a13, a21, a22, a23, a31, a32, a33 = model.getMatrix(name)
            else:
                model.setMatrixChoice(name,'choice','off')
                self.frameTransfo.hide()
                a11, a12, a13, a21, a22, a23, a31, a32, a33 = 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0
                self.groupBox_3.setChecked(False)

            self.lineEdit_4.setText(str(a11))
            self.lineEdit_5.setText(str(a12))
            self.lineEdit_6.setText(str(a13))
            self.lineEdit_8.setText(str(a21))
            self.lineEdit_9.setText(str(a22))
            self.lineEdit_7.setText(str(a23))
            self.lineEdit_11.setText(str(a31))
            self.lineEdit_12.setText(str(a32))
            self.lineEdit_10.setText(str(a33))


    @pyqtSlot(str)
    def slotA11(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_4.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a11',value )

    @pyqtSlot(str)
    def slotA12(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_5.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a12',value )

    @pyqtSlot(str)
    def slotA13(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_6.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a13',value )

    @pyqtSlot(str)
    def slotA21(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_8.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a21',value )

    @pyqtSlot(str)
    def slotA22(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_9.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a22',value )

    @pyqtSlot(str)
    def slotA23(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_7.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a23',value )

    @pyqtSlot(str)
    def slotA31(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_11.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a31',value )

    @pyqtSlot(str)
    def slotA32(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_12.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a32',value )

    @pyqtSlot(str)
    def slotA33(self, text):
        cindex = self.treeView.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, name, local = self.modelHeadLosses.getItem(row)
            model = HeadLossesModel(self.case)
            if self.lineEdit_10.validator().state == QValidator.Acceptable:
                value = from_qvariant(text, float)
                model.setMatrixComposant(name,'a33',value )


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
