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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import from_qvariant, to_qvariant, to_text_string
from code_saturne.Base.QtPage import DoubleValidator, RegExpValidator
from code_saturne.Pages.NotebookForm import Ui_NotebookForm
from code_saturne.Pages.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NotebookView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------
class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, idx, name, value, oturns_var, description):
        self.index  = idx
        self.name   = name
        self.value  = value
        self.oturns = oturns_var
        self.descr  = description

    def __repr__(self):
        return "variable : %s // value : %s // used with openturns : %s"\
               % (self.name, self.value, self.oturns)

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
        return 4


    def data(self, column, role):
        if self.item == None:
            if column == 0:
                return to_qvariant(self.header)
            else:
                return to_qvariant()
        else:
            if column == 0 and role == Qt.DisplayRole:
                return to_qvariant(self.item.name)
            elif column == 1 and role == Qt.DisplayRole:
                return to_qvariant(self.item.value)
            elif column == 2 and role == Qt.DisplayRole:
                return to_qvariant(self.item.oturns)
            elif column == 3 and role == Qt.DisplayRole:
                return to_qvariant(self.item.descr)
        return to_qvariant()


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None, xml_model=None):
        super(LabelDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v =  RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())
            model.setData(index, to_qvariant(p_value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None, xml_model=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        vd = DoubleValidator(editor)
        editor.setValidator(vd)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(str(value))

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class OTVariableDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """

    def __init__(self, parent = None):
        super(OTVariableDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem("No")
        editor.addItem("Yes: Input")
        editor.addItem("Yes: Output")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        #string = index.model().dataMeshes[row][col]
        string = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, to_qvariant(value))


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class DescrDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None, xml_model=None):
        super(DescrDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        rx = "[\-_A-Za-z0-9 ]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v =  RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())
            model.setData(index, to_qvariant(p_value), Qt.DisplayRole)


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
            return 4


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return to_qvariant()

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return to_qvariant(self.tr("variable name"))
            elif index.column() == 1:
                return to_qvariant(self.tr("value"))
            elif index.column() == 2:
                return to_qvariant(self.tr("OpenTurns Variable"))
            elif index.column() == 3:
                return to_qvariant(self.tr("Description"))

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        itm = index.internalPointer()
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("variable name"))
            elif section == 1:
                return to_qvariant(self.tr("value"))
            elif section == 2:
                return to_qvariant(self.tr("OpenTurns Variable"))
            elif section == 3:
                return to_qvariant(self.tr("Description"))
        return to_qvariant()


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
            descr  = self.mdl.getVariableDescription(idx)
            item = item_class(idx, cname, value, oturns, descr)
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
            self.mdl.setVariableValue(item.item.index, item.item.value)

        elif index.column() == 2:
            value = from_qvariant(value, to_text_string)
            item.item.oturns = value
            self.mdl.setVariableOt(item.item.index, item.item.oturns)

        elif index.column() == 3:
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
    Class to open the Body Forces (gravity) Page.
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

        nameDelegate = LabelDelegate(self.treeViewNotebook, self.mdl)
        self.treeViewNotebook.setItemDelegateForColumn(0, nameDelegate)

        valDelegate = ValueDelegate(self.treeViewNotebook, self.mdl)
        self.treeViewNotebook.setItemDelegateForColumn(1, valDelegate)

        otvalDelegate = OTVariableDelegate(self.treeViewNotebook)
        self.treeViewNotebook.setItemDelegateForColumn(2, otvalDelegate)

        descriptionDelegate = DescrDelegate(self.treeViewNotebook)
        self.treeViewNotebook.setItemDelegateForColumn(3, descriptionDelegate)

        self.treeViewNotebook.resizeColumnToContents(0)

        # Connections
        self.toolButtonAdd.clicked.connect(self.slotAddVariable)
        self.toolButtonDelete.clicked.connect(self.slotDeleteVariable)
        self.toolButtonImport.clicked.connect(self.slotImportVariable)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddVariable(self):
        """
        Add one variable
        """
        name = self.mdl.addVariable()
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
