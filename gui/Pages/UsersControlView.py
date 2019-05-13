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
This module defines the 'Users control' page.

This module contains the following classes:
- NameDelegate
- LocationDelegate
- StandardItemModelUsersControl
- UsersControlView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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

from code_saturne.Base.QtPage import ComboModel, RegExpValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from UsersControl import Ui_UsersControl
from code_saturne.model.UsersControlModel import UsersControlModel
from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("UsersControl")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# lineEdit delegate for the name
#-------------------------------------------------------------------------------

class NameDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_pname = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_name = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_pname = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_pname = str(editor.text())

            if new_pname in model.mdl.getUserNamesList():
                default = {}
                default['name']  = self.old_pname
                default['list']   = model.mdl.getUserNamesList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_pname = result['name']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_pname = self.old_pname

            model.setData(index, new_pname, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Location delegate
#-------------------------------------------------------------------------------

class LocationDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(LocationDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 4, 1)
        self.modelCombo.addItem(self.tr("cells"), "cells")
        self.modelCombo.addItem(self.tr("interior faces"), "internal")
        self.modelCombo.addItem(self.tr("boundary faces"), "boundary")
        self.modelCombo.addItem(self.tr("vertices"), "vertices")

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_model=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("LocationDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Dimension delegate
#-------------------------------------------------------------------------------

class UserDimensionDelegate(QItemDelegate):
    """
    Use of a combo box in the table for the user array dimension.
    """

    def __init__(self, parent):
        super(UserDimensionDelegate, self).__init__(parent)
        self.parent = parent

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 5, 1)

        self.modelCombo.addItem(self.tr("1"), "1")
        self.modelCombo.addItem(self.tr("2"), "2")
        self.modelCombo.addItem(self.tr("3"), "3")
        self.modelCombo.addItem(self.tr("6"), "6")
        self.modelCombo.addItem(self.tr("9"), "9")

        editor.installEventFilter(self)
        return editor

    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_model=string)

    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("UserDimensionDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)

    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# StandardItemModelUsersControl class
#-------------------------------------------------------------------------------

class StandardItemModelUsersControl(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Name"),
                        self.tr("Location"),
                        self.tr("Dimension")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return None

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return data
            else:
                return None

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()

        name = self._data[row][0]
        new_val = from_qvariant(value, to_text_string)

        # Name
        if col == 0:
            self.mdl.setUsersName(name, new_val)
            self._data[row][col] = new_val
        # location
        elif col == 1:
            self.mdl.setUsersLocation(name, new_val)
            self._data[row][col] = new_val
        # dimension
        elif col == 2:
            self.mdl.setUsersDim(name, new_val)
            self._data[row][col] = new_val

        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, name=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        if name == None:
            name = self.mdl.addUser("User_" + str(row + 1))
        location = self.mdl.getUsersLocation(name)
        dim      = self.mdl.getUsersDim(name)

        scalar = [name, location, dim]

        self._data.append(scalar)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        name = self._data[row][0]
        del self._data[row]
        self.mdl.deleteUser(name)
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# UsersControlView class
#-------------------------------------------------------------------------------

class UsersControlView(QWidget, Ui_UsersControl):
    """
    UsersControlView layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_UsersControl.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.mdl = UsersControlModel(self.case)

        self.tableModelUsers = StandardItemModelUsersControl(self, self.mdl)
        self.tableViewUsers.setModel(self.tableModelUsers)
        if QT_API == "PYQT4":
            self.tableViewUsers.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
            self.tableViewUsers.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
            self.tableViewUsers.horizontalHeader().setResizeMode(2,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewUsers.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
            self.tableViewUsers.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)
            self.tableViewUsers.horizontalHeader().setSectionResizeMode(2,QHeaderView.Stretch)
        self.tableViewUsers.resizeColumnsToContents()
        self.tableViewUsers.resizeRowsToContents()
        self.tableViewUsers.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewUsers.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateName     = NameDelegate(self.tableViewUsers)
        delegateLocation = LocationDelegate(self.tableViewUsers)
        delegateDim      = UserDimensionDelegate(self.tableViewUsers)

        self.tableViewUsers.setItemDelegateForColumn(0, delegateName)
        self.tableViewUsers.setItemDelegateForColumn(1, delegateLocation)
        self.tableViewUsers.setItemDelegateForColumn(2, delegateDim)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddUsers)
        self.pushButtonDelete.clicked.connect(self.slotDeleteUsers)

        # load values
        for user in self.mdl.getUserNamesList():
            self.tableModelUsers.newItem(user)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddUsers(self):
        """
        Add a scalar
        """
        self.tableViewUsers.clearSelection()
        self.tableModelUsers.newItem()


    @pyqtSlot()
    def slotDeleteUsers(self):
        """
        Delete the a scalar from the list (one by one).
        """
        row = self.tableViewUsers.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteScalar -> %s" % row)
            self.tableModelUsers.deleteItem(row)
