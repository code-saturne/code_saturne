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
This module defines the 'Users control' page.

This module contains the following classes:
- LabelDelegate
- SupportDelegate
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, RegExpValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from UsersControl import Ui_UsersControl
from UsersControlModel import UsersControlModel
from code_saturne.Base.Common import LABEL_LENGTH_MAX

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("UsersControl")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# lineEdit delegate for the label
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_plabel = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_label = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_plabel = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_plabel = str(editor.text())

            if new_plabel in model.mdl.getUserLabelList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getUserLabelList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, to_qvariant(new_plabel), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Support delegate
#-------------------------------------------------------------------------------

class SupportDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(SupportDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 3, 1)
        self.modelCombo.addItem(self.tr("cells"), "cells")
        self.modelCombo.addItem(self.tr("internal faces"), "internal")
        self.modelCombo.addItem(self.tr("boundary faces"), "boundary")

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
        log.debug("SupportDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)


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
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)

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

        self.headers = [ self.tr("Label"),
                         self.tr("Support"),
                         self.tr("Dimension")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(data)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

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
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()
        name = row + 1

        new_val = from_qvariant(value, to_text_string)
        self._data[row][col] = new_val

        # Label
        if col == 0:
            self.mdl.setUsersLabel(name, new_val)
        # support
        elif col == 1:
            self.mdl.setUsersSupport(name, new_val)
        # dimension
        elif col == 2:
            self.mdl.setUsersDim(name, new_val)

        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_users=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        name = ""
        if existing_users == None :
            name = self.mdl.addUser(row + 1)
        else :
            name = existing_users
        label   = self.mdl.getUsersLabel(name)
        support = self.mdl.getUsersSupport(name)
        dim     = self.mdl.getUsersDim(name)

        scalar = [label, support, dim]

        self._data.append(scalar)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        self.mdl.deleteUser(row + 1)
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

#        if self.case['package'].name == 'code_saturne':
#            from UsersControlModelCS import UsersControlModelCS as UsersControlMdl
#        else:
#            from UsersControlModelNCFD import UsersControlModelNCFD as UsersControlMdl

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

        delegateLabel   = LabelDelegate(self.tableViewUsers)
        delegateSupport = SupportDelegate(self.tableViewUsers)
        delegateDim     = UserDimensionDelegate(self.tableViewUsers)

        self.tableViewUsers.setItemDelegateForColumn(0, delegateLabel)
        self.tableViewUsers.setItemDelegateForColumn(1, delegateSupport)
        self.tableViewUsers.setItemDelegateForColumn(2, delegateDim)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddUsers)
        self.pushButtonDelete.clicked.connect(self.slotDeleteUsers)

        # load values
        for user in self.mdl.getUsersList():
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


