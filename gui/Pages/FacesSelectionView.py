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
- LineEditDelegateVerbosity
- LineEditDelegateSelector
- StandardItemModelFaces
- FacesSelectionView
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

from code_saturne.Base.QtPage import RegExpValidator, DoubleValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Pages.FacesSelectionForm import Ui_FacesSelectionForm

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FacesSelectionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for references
#-------------------------------------------------------------------------------

class LineEditDelegateVerbosity(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator =  RegExpValidator(editor, QRegExp("^[0-9 ]*$"))
        editor.setValidator(validator)
        return editor


    def setEditorData(self, lineEdit, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        lineEdit.setText(value)


    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, to_qvariant(value), Qt.DisplayRole)

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
        validator =  RegExpValidator(editor, QRegExp("[ -~]*"))
        editor.setValidator(validator)
        return editor


    def setEditorData(self, lineEdit, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        lineEdit.setText(value)


    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for Fraction and Plane
#-------------------------------------------------------------------------------

class FractionPlaneDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(FractionPlaneDelegate, self).__init__(parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor, min=0.)
        validator.setExclusiveMin(True)
        editor.setValidator(validator)
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
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelFaces(QStandardItemModel):

    def __init__(self, parent, mdl=None, tag=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.parent = parent
        self.mdl = mdl
        self.tag = tag

        self.headers = [self.tr("Fraction"),
                        self.tr("Plane"),
                        self.tr("Verbosity"),
                        self.tr("Visualization"),
                        self.tr("Selection criteria")]

        self.tooltip = [self.tr("Relative merge tolerance (between 0 and 0.49, 0.1 by default) : maximum intersection distance"),
                        self.tr("Maximum angle between normals of coplanar faces (25 degrees by default)"),
                        self.tr("Verbosity level"),
                        self.tr("Visualization output level (0 for none)"),
                        self.tr("Selection criteria string")]

        self.setColumnCount(len(self.headers))

        self.dataFaces = []
        if tag and mdl:
            self.populateModel()


    def populateModel(self):

        # Default values
        self.default = {}
        for key in ('selector', 'fraction', 'plane', 'verbosity', 'visualization'):
            self.default[key] = self.mdl.defaultValues()[key]

        if self.tag == "face_joining":
            for j in range(self.mdl.getJoinSelectionsCount()):
                d = self.mdl.getJoinFaces(j)
                self.dataFaces.append(d)
                row = self.rowCount()
                self.setRowCount(row+1)

        elif self.tag == "face_periodicity":
            for j in range(self.mdl.getPeriodicSelectionsCount()):
                d = self.mdl.getPeriodicFaces(j)
                self.dataFaces.append(d)
                row = self.rowCount()
                self.setRowCount(row+1)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        row = index.row()
        col = index.column()

        if role == Qt.ToolTipRole:
            return to_qvariant(self.tooltip[col])

        if role == Qt.DisplayRole:
            if col == 0:
                return to_qvariant(self.dataFaces[row]['fraction'])
            elif col == 1:
                return to_qvariant(self.dataFaces[row]['plane'])
            elif col == 2:
                return to_qvariant(self.dataFaces[row]['verbosity'])
            elif col == 3:
                return to_qvariant(self.dataFaces[row]['visualization'])
            elif col == 4:
                return to_qvariant(self.dataFaces[row]['selector'])

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
            self.dataFaces[row]['fraction'] = str(from_qvariant(value, to_text_string))
        elif col == 1:
            self.dataFaces[row]['plane'] = str(from_qvariant(value, to_text_string))
        elif col == 2:
            self.dataFaces[row]['verbosity'] = str(from_qvariant(value, to_text_string))
        elif col == 3:
            self.dataFaces[row]['visualization'] = str(from_qvariant(value, to_text_string))
        elif col == 4:
            self.dataFaces[row]['selector'] = str(from_qvariant(value, to_text_string))

        if self.tag == "face_joining":
            self.mdl.replaceJoinFaces(row, self.dataFaces[row])

        elif self.tag == "face_periodicity":
            self.mdl.replacePeriodicFaces(row, self.dataFaces[row])

        log.debug("setData -> dataFaces = %s" % self.dataFaces)
        self.dataChanged.emit(index, index)
        return True


    def addItem(self):
        """
        Add an item in the QListView.
        """
        title = self.tr("Warning")

        if self.tag == "face_joining":
            self.mdl.addJoinFaces(self.default)

        elif self.tag == "face_periodicity":
            self.mdl.addPeriodicFaces(self.default)

        self.dataFaces.append(self.default.copy())
        log.debug("addItem -> dataFaces = %s" % self.dataFaces)
        row = self.rowCount()
        self.setRowCount(row+1)


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        log.debug("StandardItemModelFaces.delete row = %i %i" % (row, self.mdl.getJoinSelectionsCount()))

        if self.tag == "face_joining":
            self.mdl.deleteJoinFaces(row)
            del self.dataFaces[row]
            row = self.rowCount()
            self.setRowCount(row-1)

        elif self.tag == "face_periodicity":
            self.mdl.deletePeriodicity(row)
            del self.dataFaces[row]
            row = self.rowCount()
            self.setRowCount(row-1)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class FacesSelectionView(QWidget, Ui_FacesSelectionForm):
    """
    """
    def __init__(self, *args):
        """
        Constructor.
        """
        QWidget.__init__(self, *args)
        Ui_FacesSelectionForm.__init__(self)
        self.setupUi(self)

        self.modelFaces = None

        self.tableView.setModel(self.modelFaces)

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.tableView.horizontalHeader().setStretchLastSection(True)

        delegateFraction = FractionPlaneDelegate(self.tableView)
        self.tableView.setItemDelegateForColumn(0, delegateFraction)

        delegatePlane = FractionPlaneDelegate(self.tableView)
        self.tableView.setItemDelegateForColumn(1, delegatePlane)

        delegateVerbosity = LineEditDelegateVerbosity(self.tableView)
        self.tableView.setItemDelegateForColumn(2, delegateVerbosity)
        self.tableView.setItemDelegateForColumn(3, delegateVerbosity)

        delegateSelector = LineEditDelegateSelector(self.tableView)
        self.tableView.setItemDelegateForColumn(4, delegateSelector)

        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableView.setSelectionMode(QAbstractItemView.SingleSelection)

        # Connections

        self.pushButtonAdd.clicked.connect(self.slotAddItem)
        self.pushButtonDelete.clicked.connect(self.slotDelItem)


    def slotAddItem(self):
        """
        Create a new faces selection.
        """
        self.modelFaces.addItem()


    def slotDelItem(self):
        """
        Delete a single selected row.
        """
        for index in self.tableView.selectionModel().selectedIndexes():
            self.modelFaces.delItem(index.row())
            break


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
    FacesSelectionView = FacesSelectionView(app)
    FacesSelectionView.show()
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
