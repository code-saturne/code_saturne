# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
- OutputSurfacicVariablesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.OutputSurfacicVariablesForm import Ui_OutputSurfacicVariablesForm
from code_saturne.model.OutputControlModel import OutputControlModel
from code_saturne.model.OutputSurfacicVariablesModel import OutputSurfacicVariablesModel
from code_saturne.model.OutputSurfacicFieldsModel import OutputSurfacicFieldsModel #AZ
from code_saturne.gui.case.OutputVolumicVariablesView import LabelDelegate

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputSurfacicVariablesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModelOutput class
#-------------------------------------------------------------------------------

class StandardItemModelOutput(QStandardItemModel):

    def __init__(self, case, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.case = case
        self.mdl = mdl

        self.setColumnCount(3)
        self.dataLabel      = []
        self.dataName       = []
        self.dataPost       = []
        self.disableItem    = []
        self.populateModel()


    def populateModel(self):
        # Data initialization

        for name in self.mdl.list_name:
            row = self.rowCount()
            self.setRowCount(row + 1)

            label = self.mdl.dicoLabelName[name]
            post  = self.mdl.getPostProcessing(label)

            if not OutputControlModel(self.case).isSurfaceWriterActive():
                self.disableItem.append((row, 1))
                post = "off"

            self.dataLabel.append(label)
            self.dataName.append(name)
            self.dataPost.append(post)


    def data(self, index, role):
        if not index.isValid():
            return None

        # ToolTips BUG
        if role == Qt.ToolTipRole:
            if index.column() == 2:
                return self.tr("code_saturne keyword: ipstdv")

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 2:
                return "Post-processing"

        # Display
        if role == Qt.DisplayRole:
            row = index.row()
            if index.column() == 0:
                return self.dataLabel[row]
            elif index.column() == 1:
                return self.dataName[row]
            else:
                return None

        # CheckState
        if role == Qt.CheckStateRole:
            row = index.row()
            if index.column() == 2:
                value = self.dataPost[row]
                if value == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        if (index.row(), index.column()) in self.disableItem:
            return Qt.ItemIsSelectable
        elif index.column() == 0 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 1 :
            return  Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 2 :
            return  Qt.ItemIsEnabled | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return self.tr("Output label")
            if section == 1:
                return self.tr("Internal name")
            elif section == 2:
                return self.tr("Post-\nprocessing")
        return None


    def setData(self, index, value, role=None):
        row = index.row()
        if index.column() == 0:
            label = str(from_qvariant(value, to_text_string))
            if label == "":
                label = self.dataLabel[row]
            self.mdl.setPropertyLabel(self.dataLabel[row], label)
            self.dataLabel[row] = label

        elif index.column() == 2:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                self.dataPost[row] = "on"
            else:
                self.dataPost[row] = "off"
            if not OutputControlModel(self.case).isSurfaceWriterActive():
                self.dataPost[row] = "off"

            self.mdl.setPostProcessing(self.dataLabel[row], self.dataPost[row])

        self.dataChanged.emit(index, index)
        return True

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class OutputSurfacicVariablesView(QWidget, Ui_OutputSurfacicVariablesForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_OutputSurfacicVariablesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
            self.mdl = OutputSurfacicVariablesModel(self.case)
        else:
            self.mdl = OutputSurfacicFieldsModel(self.case)

        self.modelOutput = StandardItemModelOutput(self.case, self.mdl)
        self.tableViewOutput.setModel(self.modelOutput)
        self.tableViewOutput.setAlternatingRowColors(True)
        self.tableViewOutput.resizeColumnToContents(0)
        self.tableViewOutput.resizeRowsToContents()
        self.tableViewOutput.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewOutput.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewOutput.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewOutput.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewOutput.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        labelDelegate = LabelDelegate(self.tableViewOutput, self.mdl)
        self.tableViewOutput.setItemDelegateForColumn(0, labelDelegate)

        self.case.undoStartGlobal()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
