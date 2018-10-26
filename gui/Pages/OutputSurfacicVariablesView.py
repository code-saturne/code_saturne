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
- OutputSurfacicVariablesView
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
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.OutputSurfacicVariablesForm import Ui_OutputSurfacicVariablesForm
from code_saturne.Pages.OutputControlModel import OutputControlModel
from code_saturne.Pages.OutputSurfacicVariablesModel import OutputSurfacicVariablesModel
from code_saturne.Pages.OutputSurfacicFieldsModel import OutputSurfacicFieldsModel #AZ
from code_saturne.Pages.OutputControlModelNeptune import OutputControlModel as OutputControlModelNeptune
from code_saturne.Pages.OutputVolumicVariablesView import LabelDelegate

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

            if self.case.xmlRootNode().tagName == "Code_Saturne_GUI" :
                if OutputControlModel(self.case).getAssociatedWriterIdList("-2") == []:
                    self.disableItem.append((row, 1))
                    post = "off"
            else :
                if OutputControlModelNeptune(self.case).getAssociatedWriterIdList("-2") == []:
                    self.disableItem.append((row, 1))
                    post = "off"

            self.dataLabel.append(label)
            self.dataName.append(name)
            self.dataPost.append(post)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        # ToolTips BUG
        if role == Qt.ToolTipRole:
            if index.column() == 2:
                return to_qvariant(self.tr("Code_Saturne keyword: ipstdv"))

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 2:
                return to_qvariant("Post-processing")

        # Display
        if role == Qt.DisplayRole:
            row = index.row()
            if index.column() == 0:
                return to_qvariant(self.dataLabel[row])
            elif index.column() == 1:
                return to_qvariant(self.dataName[row])
            else:
                return to_qvariant()

        # CheckState
        if role == Qt.CheckStateRole:
            row = index.row()
            if index.column() == 2:
                value = self.dataPost[row]
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)

        return to_qvariant()


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
                return to_qvariant(self.tr("Output label"))
            if section == 1:
                return to_qvariant(self.tr("Internal name"))
            elif section == 2:
                return to_qvariant(self.tr("Post-\nprocessing"))
        return to_qvariant()


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
            if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
                if OutputControlModel(self.case).getAssociatedWriterIdList("-2") == []:
                    self.dataPost[row] = "off"
            else :
                if OutputControlModelNeptune(self.case).getAssociatedWriterIdList("-2") == []:
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


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
