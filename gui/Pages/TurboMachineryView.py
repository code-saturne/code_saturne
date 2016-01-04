# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

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

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
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
        self.tableViewTurboMachinery.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)

        # Faces to join selection (Custom Widgets)
        model = StandardItemModelFaces(self, self.mdl, 'face_joining')
        self.widgetFacesJoin.modelFaces = model
        self.widgetFacesJoin.tableView.setModel(model)

        # Connections
        self.connect(self.comboBoxTurboMachineryType, SIGNAL("activated(const QString&)"), self.slotTurboModel)
        self.connect(self.rotorModel                , SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.dataChanged)
        self.connect(self.tableViewTurboMachinery   , SIGNAL("clicked(const QModelIndex &)"), self.slotChangeSelection)

        self.connect(self.pushButtonAdd, SIGNAL("clicked()"), self.slotAddRotor)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.slotDeleteRotor)

        self.connect(self.lineEditDX, SIGNAL("textChanged(const QString &)"), self.slotRotationX)
        self.connect(self.lineEditDY, SIGNAL("textChanged(const QString &)"), self.slotRotationY)
        self.connect(self.lineEditDZ, SIGNAL("textChanged(const QString &)"), self.slotRotationZ)

        self.connect(self.lineEditX1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationX1)
        self.connect(self.lineEditY1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationY1)
        self.connect(self.lineEditZ1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationZ1)

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
                self.groupBoxJoin.setTitle("Face joining (optionnal)")
            elif mdl == "transient":
                self.groupBoxJoin.setTitle("Face joining")
        else:
            self.groupBoxDefineTurboMachinery.hide()
            self.groupBoxRotation.hide()
            self.groupBoxJoin.hide()


    @pyqtSignature("const QModelIndex &, const QModelIndex &")
    def dataChanged(self, topLeft, bottomRight):
        self.tableViewTurboMachinery.resizeColumnsToContents()
        self.tableViewTurboMachinery.resizeRowsToContents()
        self.tableViewTurboMachinery.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)

        self.updateView()


    @pyqtSignature("const QModelIndex &")
    def slotChangeSelection(self, text=None):
        """
        detect change selection to update constant properties
        """
        self.updateView()


    @pyqtSignature("const QString&")
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


    @pyqtSignature("")
    def slotAddRotor(self):
        """
        Add rotor
        """
        self.mdl.addRotor()
        self.rotorModel.addItem(len(self.mdl.getRotorList()) -1)
        self.tableViewTurboMachinery.clearSelection()
        self.updateView()


    @pyqtSignature("")
    def slotDeleteRotor(self):
        """
        Delete the selected rotor from the list
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        self.mdl.delRotor(rotor_id)
        self.rotorModel.delItem(rotor_id)
        self.tableViewTurboMachinery.clearSelection()
        self.updateView()


    @pyqtSignature("const QString&")
    def slotRotationX(self, text):
        """
        Periodicity rotation for X
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_x", val)


    @pyqtSignature("const QString&")
    def slotRotationY(self, text):
        """
        Periodicity rotation for Y
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_y", val)


    @pyqtSignature("const QString&")
    def slotRotationZ(self, text):
        """
        Periodicity rotation for Z
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationVector(rotor_id, "axis_z", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationX1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationCenter(rotor_id, "invariant_x", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationY1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setRotationCenter(rotor_id, "invariant_y", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationZ1(self, text):
        """
        Periodicity : center of rotation
        """
        rotor_id = self.tableViewTurboMachinery.currentIndex().row()
        if self.sender().validator().state == QValidator.Acceptable:
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
