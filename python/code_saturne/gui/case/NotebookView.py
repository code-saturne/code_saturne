# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
- NotebookView
"""

#-------------------------------------------------------------------------------
# Library modules import
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

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam, GuiLabelManager
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.base.QtPage import DoubleValidator, RegExpValidator
from code_saturne.gui.base.QtPage import LabelDelegate, FloatDelegate, ComboDelegate
from code_saturne.gui.case.NotebookForm import Ui_NotebookForm
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NotebookView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# List of labels which are forbidden as notebook variables names
#-------------------------------------------------------------------------------

_forbidden_labels = GuiLabelManager().getForbidden("notebook")

#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------
class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, idx, name, value, oturns_var, editable, restart, log,
                 description):
        self.index  = idx
        self.name   = name
        self.value  = value
        self.oturns = oturns_var
        self.edit   = editable
        self.read   = restart
        self.log    = log
        self.descr  = description


    def __repr__(self):
        return ("index : %s // variable : %s // value : %s // " \
                "used with openturns : %s // editable : %s // " \
                "read at restart : %s // log : %s //" \
                % (str(self.index), str(self.name), str(self.value), \
                   str(self.oturns), str(self.edit), str(self.read), \
                   str(self.log)))

#-------------------------------------------------------------------------------
# Treeitem class
#-------------------------------------------------------------------------------
class TreeItem(object):
    '''
    a python object used to return row/column data, and keep note of
    it's parents and/or children
    '''
    def __init__(self, item, header, parentItem):
        self.item = item
        self.parentItem = parentItem
        self.header = header
        self.childItems = []


    def appendChild(self, item):
        self.childItems.append(item)


    def child(self, row):
        return self.childItems[row]


    def childCount(self):
        return len(self.childItems)


    def columnCount(self):
        return 7


    def data(self, column, role):
        if self.item is None:
            if column == 0:
                return self.header
            else:
                return None
        else:
            if column == 0 and role == Qt.DisplayRole:
                return self.item.name
            elif column == 1 and role == Qt.DisplayRole:
                return self.item.value
            elif column == 2 and role == Qt.DisplayRole:
                return self.item.oturns
            elif column == 3 and role == Qt.DisplayRole:
                return self.item.edit
            elif column == 4 and role == Qt.DisplayRole:
                return self.item.read
            elif column == 5 and role == Qt.DisplayRole:
                return self.item.log
            elif column == 6 and role == Qt.DisplayRole:
                return self.item.descr
        return None


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


    def __repr__(self):
        msg = "TreeItem: item : %s\n" % repr(self.item)

        if self.parentItem:
            msg += "parentItem : %s\n" % str(self.parentItem.header)
        else:
            msg += "parentItem : ALL\n"

        msg += "childItems : "
        for child in self.childItems:
            msg += child.header + " "
        msg += "\n"

        return msg 


#-------------------------------------------------------------------------------
# StandarItemModelOutput class
#-------------------------------------------------------------------------------

class VariableStandardItemModel(QAbstractItemModel):

    def __init__(self, parent, case, mdl):
        """
        """
        QAbstractItemModel.__init__(self)

        self.parent = parent
        self.case   = case
        self.mdl    = mdl

        self.rootItem = TreeItem(None, "ALL", None)
        self.parents = {0 : self.rootItem}

        self.disabledItem = []
        self.populateModel()


    def columnCount(self, parent = None):
        if parent and parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return 7


    def data(self, index, role):
        if not index.isValid():
            return None

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return None

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return self.tr("variable name")
            elif index.column() == 1:
                return self.tr("value")
            elif index.column() == 2:
                return self.tr("OpenTurns Variable")
            elif index.column() == 3:
                return self.tr("Editable")
            elif index.column() == 4:
                return self.tr("Read at restart")
            elif index.column() == 5:
                return self.tr("Print to default log file")
            elif index.column() == 6:
                return self.tr("Description")

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        if index.column() == 3:
            if "Yes" in self.mdl.getVariableOt(index.row()):
                return Qt.NoItemFlags
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return self.tr("Variable name")
            elif section == 1:
                return self.tr("Value")
            elif section == 2:
                return self.tr("OpenTurns\nVariable")
            elif section == 3:
                return self.tr("Editable")
            elif section == 4:
                return self.tr("Read at\nrestart")
            elif section == 5:
                return self.tr("Print to default\nlog file")
            elif section == 6:
                return self.tr("Description")
        return None


    def index(self, row, column, parent = QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        try:
            childItem = parentItem.child(row)
        except:
            childItem = None

        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()


    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        if not childItem:
            return QModelIndex()

        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)


    def rowCount(self, parent=QModelIndex()):
        if parent.column() > 0:
            return 0
        if not parent.isValid():
            p_Item = self.rootItem
        else:
            p_Item = parent.internalPointer()
        return p_Item.childCount()


    def populateModel(self):
        for idx in self.mdl.getVarList():
            parentItem = self.rootItem
            cname  = self.mdl.getVariableName(idx)
            value  = self.mdl.getVariableValue(idx)
            oturns = self.mdl.getVariableOt(idx)
            edit   = self.mdl.getVariableEditable(idx)
            read   = self.mdl.getVariableRestart(idx)
            log    = self.mdl.getVariableLog(idx)
            descr  = self.mdl.getVariableDescription(idx)
            item = item_class(idx, cname, value, oturns, edit, read, log, descr)
            new_item = TreeItem(item, cname, parentItem)
            parentItem.appendChild(new_item)


    def setData(self, index, value, role=None):
        item = index.internalPointer()

        if index.column() == 0:
            v = from_qvariant(value, to_text_string)
            item.item.name = v
            self.mdl.setVariableName(item.item.index, item.item.name)

        elif index.column() == 1:
            value = str(from_qvariant(value, float))
            item.item.value = value
            self.mdl.setVariableValue(item.item.value, idx=item.item.index)

        elif index.column() == 2:
            value = from_qvariant(value, to_text_string)
            item.item.oturns = value
            self.mdl.setVariableOt(item.item.index, item.item.oturns)
            # Input OT variables not editable and not written in log file
            if value == "Yes: Input":
                item.item.edit = "No"
                self.mdl.setVariableEditable(item.item.index, item.item.edit)
                item.item.log = "No"
                self.mdl.setVariableLog(item.item.index, item.item.log)
            # Output OT variables are editable and written in log file
            elif value == "Yes: Output":
                item.item.edit = "Yes"
                self.mdl.setVariableEditable(item.item.index, item.item.edit)
                # value already set to Yes in model
                item.item.log = "Yes"

        elif index.column() == 3:
            editable = from_qvariant(value, to_text_string)
            item.item.edit = editable
            self.mdl.setVariableEditable(item.item.index, item.item.edit)
            if editable == "Yes":
                # value already set to Yes in model
                item.item.log = "Yes"

        elif index.column() == 4:
            restart = from_qvariant(value, to_text_string)
            item.item.read = restart
            self.mdl.setVariableRestart(item.item.index, item.item.read)

        elif index.column() == 5:
            log = from_qvariant(value, to_text_string)
            item.item.log = log
            self.mdl.setVariableLog(item.item.index, item.item.log)

        elif index.column() == 6:
            description = from_qvariant(value, to_text_string)
            item.item.descr = description
            self.mdl.setVariableDescription(item.item.index, item.item.descr)
        self.dataChanged.emit(QModelIndex(), QModelIndex())

        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class NotebookView(QWidget, Ui_NotebookForm):
    """
    Class to open the Notebook Page.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NotebookForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.case.undoStopGlobal()
        self.mdl = NotebookModel(self.case)

        self.modelVar = VariableStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewNotebook.setModel(self.modelVar)
        self.treeViewNotebook.setAlternatingRowColors(True)
        self.treeViewNotebook.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewNotebook.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewNotebook.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewNotebook.expandAll()
        self.treeViewNotebook.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewNotebook.setDragEnabled(False)

        nameDelegate = LabelDelegate(self.treeViewNotebook,
                                     xml_model=self.mdl,
                                     forbidden_labels=_forbidden_labels)
        self.treeViewNotebook.setItemDelegateForColumn(0, nameDelegate)

        valDelegate = FloatDelegate(self.treeViewNotebook, self.mdl)
        self.treeViewNotebook.setItemDelegateForColumn(1, valDelegate)

        OTOptions = ["No", "Yes: Input", "Yes: Output"]
        otvalDelegate = ComboDelegate(self.treeViewNotebook,
                                      opts_list=OTOptions)
        self.treeViewNotebook.setItemDelegateForColumn(2, otvalDelegate)

        EditableOptions = ["No","Yes"]
        editableDelegate = ComboDelegate(self.treeViewNotebook,
                                         opts_list=EditableOptions)
        self.treeViewNotebook.setItemDelegateForColumn(3, editableDelegate)

        RestartOptions = ["Yes","No"]
        restartDelegate = ComboDelegate(self.treeViewNotebook,
                                        opts_list=RestartOptions)
        self.treeViewNotebook.setItemDelegateForColumn(4, restartDelegate)

        LogOptions = ["No","Yes"]
        logDelegate = ComboDelegate(self.treeViewNotebook,
                                    opts_list=LogOptions)
        self.treeViewNotebook.setItemDelegateForColumn(5, logDelegate)

        descriptionDelegate = LabelDelegate(self.treeViewNotebook)
        self.treeViewNotebook.setItemDelegateForColumn(6, descriptionDelegate)

        self.treeViewNotebook.resizeColumnToContents(0)
        self.treeViewNotebook.resizeColumnToContents(2)
        self.treeViewNotebook.resizeColumnToContents(4)
        self.treeViewNotebook.resizeColumnToContents(5)
        self.treeViewNotebook.resizeColumnToContents(6)

        # Connections
        self.toolButtonAdd.clicked.connect(self.slotAddVariable)
        self.toolButtonDelete.clicked.connect(self.slotDeleteVariable)
        self.toolButtonImport.clicked.connect(self.slotImportVariable)

        self.modelVar.dataChanged.connect(self.dataChanged)

        self.case.undoStartGlobal()


    def dataChanged(self):
        self.update()


    @pyqtSlot()
    def slotAddVariable(self):
        """
        Add one variable
        """
        self.mdl.addVariable()
        self.modelVar = VariableStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewNotebook.setModel(self.modelVar)
        self.treeViewNotebook.expandAll()


    @pyqtSlot()
    def slotDeleteVariable(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        current = self.treeViewNotebook.currentIndex()
        idx = current.row()
        self.mdl.deleteVariable(idx)
        self.modelVar = VariableStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewNotebook.setModel(self.modelVar)
        self.treeViewNotebook.expandAll()


    @pyqtSlot()
    def slotImportVariable(self):
        """
        select a csv/txt file to add and update variables
        """
        log.debug("slotImportVariable")

        data = self.case['data_path']
        title = self.tr("variable(s) file")
        filetypes = self.tr("csv file (*.csv);;txt file (*.txt);;All Files (*)")
        fle = QFileDialog.getOpenFileName(self, title, data, filetypes)[0]
        fle = str(fle)
        if not fle:
            return
        fle = os.path.abspath(fle)

        self.mdl.ImportVariableFromFile(fle)
        self.modelVar = VariableStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewNotebook.setModel(self.modelVar)
        self.treeViewNotebook.expandAll()


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
