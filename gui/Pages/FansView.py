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
This module contains the following classes:
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

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import RegExpValidator, IntValidator, DoubleValidator
from code_saturne.Base.QtPage import ComboModel, to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.FansForm import Ui_FansForm
from code_saturne.model.FansModel import FansModel
from code_saturne.Pages.FacesSelectionView import StandardItemModelFaces

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FansView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Line edit delegate for index
#-------------------------------------------------------------------------------

class LineEditDelegateIndex(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for float
#-------------------------------------------------------------------------------

class LineEditDelegateFloat(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for mesh dimension
#-------------------------------------------------------------------------------

class MeshDimDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl):
        super(MeshDimDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        self.modelCombo.addItem(self.tr("2D"), "2")
        self.modelCombo.addItem(self.tr("3D"), "3")


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        model.setData(index, to_qvariant(value), Qt.DisplayRole)


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelFans(QStandardItemModel):

    def __init__(self, mdl=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.mdl = mdl

        self.headers = [self.tr("Fan id"),
                        self.tr("Mesh dimension"),
                        self.tr("Fan radius"),
                        self.tr("Hub radius"),
                        self.tr("Axial torque"),
                        self.tr("Blade radius")]

        self.tooltip = [self.tr("Fan id"),
                        self.tr("Mesh dimension"),
                        self.tr("Fan radius"),
                        self.tr("Hub radius"),
                        self.tr("Axial torque"),
                        self.tr("Blade radius")]

        self.setColumnCount(len(self.headers))

        self._data = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant(self.tooltip[index.column()])

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            return to_qvariant(self._data[row][col])

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        if col == 0:
            new_code = from_qvariant(value, int)
            self._data[row][col] = new_code
        elif col == 1:
            criteria = str(from_qvariant(value, to_text_string))
            self._data[row][col] = criteria
            self.mdl.setFanMeshDimension(row, criteria)
        elif col == 2:
            criteria = str(from_qvariant(value, to_text_string))
            self._data[row][col] = criteria
            self.mdl.setFanProperty(row, "fan_radius", criteria)
        elif col == 3:
            criteria = str(from_qvariant(value, to_text_string))
            self._data[row][col] = criteria
            self.mdl.setFanProperty(row, "hub_radius", criteria)
        elif col == 4:
            criteria = str(from_qvariant(value, to_text_string))
            self._data[row][col] = criteria
            self.mdl.setFanProperty(row, "axial_torque", criteria)
        elif col == 5:
            criteria = str(from_qvariant(value, to_text_string))
            self._data[row][col] = criteria
            self.mdl.setFanProperty(row, "blades_radius", criteria)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def addItem(self, idx):
        """
        Add an item in the QListView.
        """
        row = self.rowCount()
        crit = self.mdl.getFanMeshDimension(idx)
        fr   = self.mdl.getFanProperty(row, "fan_radius")
        hr   = self.mdl.getFanProperty(row, "hub_radius")
        at   = self.mdl.getFanProperty(row, "axial_torque")
        br   = self.mdl.getFanProperty(row, "blades_radius")
        fan = [idx, crit, fr, hr, at, br]
        self._data.append(fan)
        self.setRowCount(row+1)


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)
        for id in range(0, len(self.mdl.getFanList())):
            self._data[id][0] = id


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------
class FansView(QWidget, Ui_FansForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)
        Ui_FansForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = FansModel(self.case)

        # tableView fans
        self.fansModel = StandardItemModelFans(self.mdl)
        self.tableViewFans.setModel(self.fansModel)
        self.tableViewFans.setAlternatingRowColors(True)
        self.tableViewFans.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewFans.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateIdx = LineEditDelegateIndex(self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(0, delegateIdx)

        delegateMeshDim = MeshDimDelegate(self, self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(1, delegateMeshDim)

        delegateFanRadius = LineEditDelegateFloat(self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(2, delegateFanRadius)

        delegateHubRadius = LineEditDelegateFloat(self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(3, delegateHubRadius)

        delegateAxialTorque = LineEditDelegateFloat(self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(4, delegateAxialTorque)

        delegateBladeRadius = LineEditDelegateFloat(self.tableViewFans)
        self.tableViewFans.setItemDelegateForColumn(5, delegateBladeRadius)

        self.tableViewFans.resizeColumnsToContents()
        self.tableViewFans.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewFans.horizontalHeader().setResizeMode(5,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFans.horizontalHeader().setSectionResizeMode(5,QHeaderView.Stretch)

        # Validator
        validatorInX = DoubleValidator(self.lineEditInletX)
        self.lineEditInletX.setValidator(validatorInX)
        validatorInY = DoubleValidator(self.lineEditInletY)
        self.lineEditInletY.setValidator(validatorInY)
        validatorInZ = DoubleValidator(self.lineEditInletZ)
        self.lineEditInletZ.setValidator(validatorInZ)
        validatorOutX = DoubleValidator(self.lineEditOutletX)
        self.lineEditOutletX.setValidator(validatorOutX)
        validatorOutY = DoubleValidator(self.lineEditOutletY)
        self.lineEditOutletY.setValidator(validatorOutY)
        validatorOutZ = DoubleValidator(self.lineEditOutletZ)
        self.lineEditOutletZ.setValidator(validatorOutZ)
        validatorCoX = DoubleValidator(self.lineEditCoefX)
        self.lineEditCoefX.setValidator(validatorCoX)
        validatorCoY = DoubleValidator(self.lineEditCoefY)
        self.lineEditCoefY.setValidator(validatorCoY)
        validatorCoZ = DoubleValidator(self.lineEditCoefZ)
        self.lineEditCoefZ.setValidator(validatorCoZ)

        # Connections
        self.pushButtonAddSFan.clicked.connect(self.slotAddFan)
        self.pushButtonDeleteFan.clicked.connect(self.slotDeleteFan)
        self.tableViewFans.pressed[QModelIndex].connect(self.slotSelectFan)
        self.lineEditInletX.textChanged[str].connect(self.slotInletX)
        self.lineEditInletY.textChanged[str].connect(self.slotInletY)
        self.lineEditInletZ.textChanged[str].connect(self.slotInletZ)
        self.lineEditOutletX.textChanged[str].connect(self.slotOutletX)
        self.lineEditOutletY.textChanged[str].connect(self.slotOutletY)
        self.lineEditOutletZ.textChanged[str].connect(self.slotOutletZ)
        self.lineEditCoefX.textChanged[str].connect(self.slotCoefX)
        self.lineEditCoefY.textChanged[str].connect(self.slotCoefY)
        self.lineEditCoefZ.textChanged[str].connect(self.slotCoefZ)

        if self.mdl.getFanList() != None:
            for i in range(len(self.mdl.getFanList())):
                self.fansModel.addItem(i)

        self.groupBoxFanOption.hide()
        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddFan(self):
        """
        Add a fan to list
        """
        self.mdl.addFan()
        self.fansModel.addItem(len(self.mdl.getFanList()) -1)
        self.tableViewFans.clearSelection()
        self.groupBoxFanOption.hide()


    @pyqtSlot()
    def slotDeleteFan(self):
        """
        Delete the selected fan from the list
        """
        idx = self.tableViewFans.currentIndex().row()
        self.mdl.delFan(idx)
        self.fansModel.delItem(idx)
        self.tableViewFans.clearSelection()
        self.groupBoxFanOption.hide()


    @pyqtSlot("QModelIndex")
    def slotSelectFan(self, index):
        """
        Return the selected item from the list.
        """
        self.groupBoxFanOption.show()

        row = index.row()
        log.debug("slotSelectFan -> %s" % (row,))
        InX   = self.mdl.getFanProperty(row, "inlet_axis_x")
        InY   = self.mdl.getFanProperty(row, "inlet_axis_y")
        InZ   = self.mdl.getFanProperty(row, "inlet_axis_z")
        OutX  = self.mdl.getFanProperty(row, "outlet_axis_x")
        OutY  = self.mdl.getFanProperty(row, "outlet_axis_y")
        OutZ  = self.mdl.getFanProperty(row, "outlet_axis_z")
        CoefX = self.mdl.getFanProperty(row, "curve_coeffs_x")
        CoefY = self.mdl.getFanProperty(row, "curve_coeffs_y")
        CoefZ = self.mdl.getFanProperty(row, "curve_coeffs_z")
        self.lineEditInletX.setText(str(InX))
        self.lineEditInletY.setText(str(InY))
        self.lineEditInletZ.setText(str(InZ))
        self.lineEditOutletX.setText(str(OutX))
        self.lineEditOutletY.setText(str(OutY))
        self.lineEditOutletZ.setText(str(OutZ))
        self.lineEditCoefX.setText(str(CoefX))
        self.lineEditCoefY.setText(str(CoefY))
        self.lineEditCoefZ.setText(str(CoefZ))


    @pyqtSlot(str)
    def slotInletX(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditInletX.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "inlet_axis_x", val)


    @pyqtSlot(str)
    def slotInletY(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditInletY.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "inlet_axis_y", val)


    @pyqtSlot(str)
    def slotInletZ(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditInletZ.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "inlet_axis_z", val)


    @pyqtSlot(str)
    def slotOutletX(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditOutletX.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "outlet_axis_x", val)


    @pyqtSlot(str)
    def slotOutletY(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditOutletY.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "outlet_axis_y", val)


    @pyqtSlot(str)
    def slotOutletZ(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditOutletZ.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "outlet_axis_z", val)


    @pyqtSlot(str)
    def slotCoefX(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditCoefX.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "curve_coeffs_x", val)


    @pyqtSlot(str)
    def slotCoefY(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditCoefY.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "curve_coeffs_y", val)


    @pyqtSlot(str)
    def slotCoefZ(self, text):
        """
        """
        idx = self.tableViewFans.currentIndex().row()
        if self.lineEditCoefZ.validator().state == QValidator.Acceptable:
            val = from_qvariant(text, float)
            self.mdl.setFanProperty(idx, "curve_coeffs_z", val)


    def tr(self, text):
        """
        Translation
        """
        return text



#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    FanView = FanView(app)
    FanView.show()
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
