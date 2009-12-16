# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module contains the following classes:
- LineEditDelegateReferences
- LineEditDelegateGroups
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.QtPage import RegExpValidator, DoubleValidator
from Base.Toolbox import GuiParam
from FacesSelectionForm import Ui_FacesSelectionForm
from Pages.SolutionDomainModel import SolutionDomainModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FacesSelectionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for references
#-------------------------------------------------------------------------------

class LineEditDelegateReferences(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator =  RegExpValidator(editor, QRegExp("^[0-9 ]*$"))
        editor.setValidator(validator)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, lineEdit, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        lineEdit.setText(value)


    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for groups
#-------------------------------------------------------------------------------

class LineEditDelegateGroups(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator =  RegExpValidator(editor, QRegExp("[ '_A-Za-z0-9]*"))
        editor.setValidator(validator)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, lineEdit, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        lineEdit.setText(value)


    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for Fraction and Plan
#-------------------------------------------------------------------------------

class FractionPlanDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(FractionPlanDelegate, self).__init__(parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor, min=0.)
        validator.setExclusiveMin(True)
        editor.setValidator(validator)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value, ok = editor.text().toDouble()
        if editor.validator().state == QValidator.Acceptable:
            model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelFaces(QStandardItemModel):

    def __init__(self, parent, mdl=None, tag=None, perio_name=None):
        """
        """
        QStandardItemModel.__init__(self)

        self.parent = parent
        self.mdl = mdl
        self.tag = tag
        self.perio_name = perio_name

        self.headers = [self.tr("References"),
                        self.tr("Groups"),
                        self.tr("Reverse"),
                        self.tr("Fraction"),
                        self.tr("Plane"),
                        self.tr("Semi-conf")]

        self.tooltip = [self.tr("Code_Saturne Preprocessor sub-option: --color"),
                        self.tr("Code_Saturne Preprocessor sub-option: --group"),
                        self.tr("Code_Saturne Preprocessor sub-option: --invsel"),
                        self.tr("Code_Saturne Preprocessor sub-option: --fraction"),
                        self.tr("Code_Saturne Preprocessor sub-option: --plane"),
                        self.tr("Code_Saturne Preprocessor sub-option: --semi-conf")]

        if self.tag == "faces_periodic":
                self.headers.pop()
                self.tooltip.pop()

        elif self.tag == "faces_select":
            self.headers = [self.tr("References"), self.tr("Groups"), self.tr("Reverse")]
            self.tooltip = [self.tr("Code_Saturne Preprocessor sub-option: --color"),
                            self.tr("Code_Saturne Preprocessor sub-option: --group"),
                            self.tr("Code_Saturne Preprocessor sub-option: --invsel")]

        self.setColumnCount(len(self.headers))

        self.dataFaces = []
        if tag and mdl:
            self.populateModel()


    def populateModel(self):

        # Default values
        self.default = {}
        for key in ('color', 'group', 'reverse', 'fraction', 'plan', 'semiconf'):
            self.default[key] = self.mdl.defaultValues()[key]

        if self.tag == "faces_join":
            for j in range(self.mdl.getJoinSelectionsNumber()):
                d = self.mdl.getJoinFaces(j+1)
                self.dataFaces.append(d)
                row = self.rowCount()
                self.setRowCount(row+1)

        else:
            if self.tag == "faces_periodic":
                d = self.mdl.getPeriodicFaces(self.perio_name)

            elif self.tag == "faces_select":
                d = self.mdl.getSelectFaces()

            if d:
                self.dataFaces.append(d)
                row = self.rowCount()
                self.setRowCount(row+1)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        row = index.row()
        col = index.column()

        if role == Qt.ToolTipRole:
            return QVariant(self.tooltip[col])

        if role == Qt.DisplayRole:
            if col == 0:
                return QVariant(self.dataFaces[row]['color'])
            elif col == 1:
                return QVariant(self.dataFaces[row]['group'])
            elif col == 3:
                return QVariant(self.dataFaces[row]['fraction'])
            elif col == 4:
                return QVariant(self.dataFaces[row]['plan'])

        if role == Qt.CheckStateRole:
            if col == 2:
                if self.dataFaces[row]['reverse'] == "on":
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)
            elif col == 5:
                if self.dataFaces[row]['semiconf'] == "on":
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() in [2, 5]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        if col   == 0:
            self.dataFaces[row]['color'] = str(value.toString())
        elif col == 1:
            self.dataFaces[row]['group'] = str(value.toString())
        elif col == 3:
            self.dataFaces[row]['fraction'] = str(value.toString())
        elif col == 4:
            self.dataFaces[row]['plan'] = str(value.toString())
        elif col == 2:
            v, ok = value.toInt()
            if v == Qt.Checked:
                self.dataFaces[row]['reverse'] = "on"
            else:
                self.dataFaces[row]['reverse'] = "off"
        elif col == 5:
            v, ok = value.toInt()
            if v == Qt.Checked:
                self.dataFaces[row]['semiconf'] = "on"
            else:
                self.dataFaces[row]['semiconf'] = "off"

        if self.tag == "faces_join":
            self.mdl.replaceJoinFaces(row+1, self.dataFaces[row])

        elif self.tag == "faces_periodic":
            self.mdl.replacePeriodicFaces(self.perio_name, self.dataFaces[row])

        elif self.tag == "faces_select":
            self.mdl.replaceSelectFaces(self.dataFaces[row])

        log.debug("setData -> dataFaces = %s" % self.dataFaces)
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self):
        """
        Add an item in the QListView.
        """
        title = self.tr("Warning")

        if self.tag == "faces_select":
            if self.mdl.getSelectFaces() or self.rowCount() == 1:
                msg = self.tr("For interior faces selection, only a single criterion is allowed.")
                QMessageBox.information(self.parent, title, msg)
                return

        elif self.tag == "faces_periodic":
            if self.mdl.getPeriodicFaces(self.perio_name):
                msg = self.tr("For a given periodicity, only a single faces selection is allowed.")
                QMessageBox.information(self.parent, title, msg)
                return

        if self.tag == "faces_join":
            self.mdl.addJoinFaces(self.default)

        elif self.tag == "faces_periodic":
            self.mdl.addPeriodicFaces(self.perio_name, self.default)

        elif self.tag == "faces_select":
            self.mdl.addSelectFaces(self.default)

        self.dataFaces.append(self.default.copy())
        log.debug("addItem -> dataFaces = %s" % self.dataFaces)
        row = self.rowCount()
        self.setRowCount(row+1)


    def delItem(self, row):
        """
        Delete an item from the QTableView.
        """
        log.debug("StandardItemModelFaces.delete row = %i %i" % (row, self.mdl.getJoinSelectionsNumber()))

        if self.tag == "faces_join":
            self.mdl.deleteJoinFaces(str(row+1))

        elif self.tag == "faces_periodic":
            self.mdl.deletePeriodicFaces(self.perio_name)

        elif self.tag == "faces_select":
            self.mdl.deleteSelectFaces()

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
        #self.modelFaces = StandardItemModelFaces(self)
        #self.tableView.setModel(self.modelFaces)

        delegateReferences = LineEditDelegateReferences(self.tableView)
        self.tableView.setItemDelegateForColumn(0, delegateReferences)

        delegateGroups = LineEditDelegateGroups(self.tableView)
        self.tableView.setItemDelegateForColumn(1, delegateGroups)

        delegateFraction = FractionPlanDelegate(self.tableView)
        self.tableView.setItemDelegateForColumn(3, delegateFraction)

        delegatePlan = FractionPlanDelegate(self.tableView)
        self.tableView.setItemDelegateForColumn(4, delegatePlan)

##        self.tableView.resizeColumnsToContents()
##         self.tableView.resizeRowsToContents()
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableView.setSelectionMode(QAbstractItemView.SingleSelection)

        # Connections

        self.connect(self.pushButtonNew,    SIGNAL("clicked()"), self.slotAddItem)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.slotDelItem)


#    def populateModel(self, model, tag):
#        self.modelFaces.populateModel(model, tag)
#
#
#    def setPeriodicityName(self, perio_name):
#        self.modelFaces.perio_name = perio_name


    @pyqtSignature("")
    def slotAddItem(self):
        """
        Create a new faces selection.
        """
        self.modelFaces.addItem()


    @pyqtSignature("")
    def slotDelItem(self):
        """
        Delete a single selectded row.
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
