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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, RegExpValidator, DoubleValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.TurboMachineryForm import Ui_TurboMachineryForm
from code_saturne.Pages.TurboMachineryModel import TurboMachineryModel
from code_saturne.Pages.FacesSelectionView import StandardItemModelFaces

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TurboMachineryView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Line edit delegate for selection
#-------------------------------------------------------------------------------

class LineEditDelegateSelector(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = editor.text()
        model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for 'velocity' in turbomachinery table
#-------------------------------------------------------------------------------

class VelocityDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(VelocityDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        vd = DoubleValidator(editor)
        editor.setValidator(vd)
        return editor


    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelRotor(QStandardItemModel):

    def __init__(self, mdl=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.mdl = mdl

        self.headers = [self.tr("Rotation velocity (rad)"),
                        self.tr("Selection criteria")]

        self.tooltip = [self.tr("Rotation velocity"),
                        self.tr("Selection criteria string")]

        self.setColumnCount(len(self.headers))

        self._data = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant(self.tooltip[index.column()])

        if role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if index.column() in (0, 1):
                if data:
                    return to_qvariant(data)
                else:
                    return to_qvariant()

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        if col == 0:
            vel = from_qvariant(value, float)
            self._data[row][col] = vel
            self.mdl.setRotorVelocity(row, vel)
        elif col == 1:
            criteria = from_qvariant(value, to_text_string)
            self._data[row][col] = criteria
            self.mdl.setRotorCriteria(row, criteria)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def addItem(self, rotor_id):
        """
        Add an item in the QListView.
        """
        row = self.rowCount()
        vel = self.mdl.getRotorVelocity(rotor_id)
        crit = self.mdl.getRotorCriteria(rotor_id)
        rotor = [vel, crit]
        self._data.append(rotor)
        self.setRowCount(row+1)


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class TurboMachineryView(QWidget, Ui_TurboMachineryForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)
        Ui_TurboMachineryForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = TurboMachineryModel(self.case)

        # Combo model
        self.modelTurboMachineryType = ComboModel(self.comboBoxTurboMachineryType, 3, 1)
        self.modelTurboMachineryType.addItem(self.tr("None"), "off")
        self.modelTurboMachineryType.addItem(self.tr("Full transient simulation"), "transient")
        self.modelTurboMachineryType.addItem(self.tr("Transient with explicit coupling"), "transient_coupled")
        self.modelTurboMachineryType.addItem(self.tr("Frozen rotor model"), "frozen")

        # Set up validators
        self.lineEditDX.setValidator(DoubleValidator(self.lineEditDX))
        self.lineEditDY.setValidator(DoubleValidator(self.lineEditDY))
        self.lineEditDZ.setValidator(DoubleValidator(self.lineEditDZ))
        self.lineEditX1.setValidator(DoubleValidator(self.lineEditX1))
        self.lineEditY1.setValidator(DoubleValidator(self.lineEditY1))
        self.lineEditZ1.setValidator(DoubleValidator(self.lineEditZ1))

        # tableView TurboMachinery
        self.rotorModel = StandardItemModelRotor(self.mdl)
        self.tableViewTurboMachinery.setModel(self.rotorModel)
        self.tableViewTurboMachinery.resizeColumnsToContents()
        self.tableViewTurboMachinery.resizeRowsToContents()
        self.tableViewTurboMachinery.setAlternatingRowColors(True)
        self.tableViewTurboMachinery.setSelectionBehavior(QAbstractItemView.SelectRows)

        delegateVelocity = VelocityDelegate(self.tableViewTurboMachinery)
        self.tableViewTurboMachinery.setItemDelegateForColumn(0, delegateVelocity)

        delegateSelector = LineEditDelegateSelector(self.tableViewTurboMachinery)
        self.tableViewTurboMachinery.setItemDelegateForColumn(1, delegateSelector)

        self.tableViewTurboMachinery.resizeColumnsToContents()
        self.tableViewTurboMachinery.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewTurboMachinery.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewTurboMachinery.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)

        # Faces to join selection (Custom Widgets)
        model = StandardItemModelFaces(self, self.mdl, 'face_joining')
        self.widgetFacesJoin.modelFaces = model
        self.widgetFacesJoin.tableView.setModel(model)

        # Connections
        self.comboBoxTurboMachineryType.activated[str].connect(self.slotTurboModel)
        self.rotorModel.dataChanged.connect(self.dataChanged)
        self.tableViewTurboMachinery.clicked[QModelIndex].connect(self.slotChangeSelection)

        self.pushButtonAdd.clicked.connect(self.slotAddRotor)
        self.pushButtonDelete.clicked.connect(self.slotDeleteRotor)

        self.lineEditDX.textChanged[str].connect(self.slotRotationX)
        self.lineEditDY.textChanged[str].connect(self.slotRotationY)
        self.lineEditDZ.textChanged[str].connect(self.slotRotationZ)

        self.lineEditX1.textChanged[str].connect(self.slotCenterRotationX1)
        self.lineEditY1.textChanged[str].connect(self.slotCenterRotationY1)
        self.lineEditZ1.textChanged[str].connect(self.slotCenterRotationZ1)

        if self.mdl.getRotorList() != None:
            for i in range(len(self.mdl.getRotorList())):
                self.rotorModel.addItem(i)

        # Initialize widget
        self.updateView()

        self.case.undoStartGlobal()


    def __setValuesRotation(self):
        """
        Put values found in xml file as soon as mode is "rotation"
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        rx, ry, rz = self.mdl.getRotationDirection(rotor_id)
        px, py, pz = self.mdl.getRotationCenter(rotor_id)

        self.lineEditDX.setText(str(rx))
        self.lineEditDY.setText(str(ry))
        self.lineEditDZ.setText(str(rz))
        self.lineEditX1.setText(str(px))
        self.lineEditY1.setText(str(py))
        self.lineEditZ1.setText(str(pz))


    def updateView(self):
        """
        Update view
        """
        mdl = self.mdl.getTurboMachineryModel()
        self.modelTurboMachineryType.setItem(str_model = mdl)
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()

        if mdl != "off":
            if len(self.mdl.getRotorList()) == 1:
                self.pushButtonDelete.setEnabled(False)
            else:
                self.pushButtonDelete.setEnabled(True)

            self.groupBoxDefineTurboMachinery.show()
            self.groupBoxJoin.show()
            self.groupBoxRotation.hide()
            if rotor_id != -1:
                self.groupBoxRotation.show()
                self.__setValuesRotation()
            if mdl == "frozen":
                self.groupBoxJoin.setTitle("Face joining (optional)")
            elif mdl == "transient":
                self.groupBoxJoin.setTitle("Face joining")
            elif mdl == "transient_coupled":
                self.groupBoxJoin.setTitle("Face mapping")

            if mdl == "transient_coupled":
                self.widgetFacesJoin.tableView.hideColumn(1)
                self.widgetFacesJoin.tableView.hideColumn(3)
            else:
                self.widgetFacesJoin.tableView.showColumn(1)
                self.widgetFacesJoin.tableView.showColumn(3)

        else:
            self.groupBoxDefineTurboMachinery.hide()
            self.groupBoxRotation.hide()
            self.groupBoxJoin.hide()


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewTurboMachinery.resizeColumnsToContents()
        self.tableViewTurboMachinery.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewTurboMachinery.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewTurboMachinery.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)

        self.updateView()


    @pyqtSlot("QModelIndex")
    def slotChangeSelection(self, text=None):
        """
        detect change selection to update constant properties
        """
        self.updateView()


    @pyqtSlot(str)
    def slotTurboModel(self, text):
        """
        Input turbomachinery model.
        """
        for nb in range(self.rotorModel.rowCount()):
            self.rotorModel.delItem(0)

        mdl = self.modelTurboMachineryType.dicoV2M[str(text)]
        self.mdl.setTurboMachineryModel(mdl)

        if len(self.mdl.getRotorList()) > 0:
            for i in range(len(self.mdl.getRotorList())):
                self.rotorModel.addItem(i)

        self.updateView()


    @pyqtSlot()
    def slotAddRotor(self):
        """
        Add rotor
        """
        self.mdl.addRotor()
        self.rotorModel.addItem(len(self.mdl.getRotorList()) -1)
        self.tableViewTurboMachinery.clearSelection()
        self.updateView()


    @pyqtSlot()
    def slotDeleteRotor(self):
        """
        Delete the selected rotor from the list
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        self.mdl.delRotor(rotor_id)
        self.rotorModel.delItem(rotor_id)
        self.tableViewTurboMachinery.clearSelection()
        self.updateView()


    @pyqtSlot(str)
    def slotRotationX(self, text):
        """
        Periodicity rotation for X
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditDX.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_x", val)


    @pyqtSlot(str)
    def slotRotationY(self, text):
        """
        Periodicity rotation for Y
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditDY.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_y", val)


    @pyqtSlot(str)
    def slotRotationZ(self, text):
        """
        Periodicity rotation for Z
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditDZ.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_z", val)


    @pyqtSlot(str)
    def slotCenterRotationX1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditX1.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationCenter(rotor_id, "invariant_x", val)


    @pyqtSlot(str)
    def slotCenterRotationY1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditY1.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationCenter(rotor_id, "invariant_y", val)


    @pyqtSlot(str)
    def slotCenterRotationZ1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.lineEditZ1.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationCenter(rotor_id, "invariant_z", val)


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
    TurboMachineryView = TurboMachineryView(app)
    TurboMachineryView.show()
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
