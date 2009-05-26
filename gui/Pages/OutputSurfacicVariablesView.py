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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from OutputSurfacicVariablesForm import Ui_OutputSurfacicVariablesForm
import Base.QtPage as QtPage
from Pages.OutputSurfacicVariablesModel import OutputSurfacicVariablesModel
from Pages.OutputVolumicVariablesView import LabelDelegate

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

    def __init__(self, case):
        """
        """
        QStandardItemModel.__init__(self)

        self.case = case
        self.mdl = OutputSurfacicVariablesModel(self.case)

        self.setColumnCount(2)
        self.dataLabel      = []
        self.dataPost       = []
        self.disableItem    = []
        self.populateModel()


    def populateModel(self):
        # Data initialization

        for name in OutputSurfacicVariablesModel(self.case).list_name:
            row = self.rowCount()
            self.setRowCount(row + 1)

            label = self.mdl.dicoLabelName[name]
            post  = self.mdl.getPostProcessing(label)

#            if name in ('yplus', 'effort', 'all_variables'):
#                post = "on"
#                if (row,1) not in self.disableItem:
#                    self.disableItem.append((row,1))

            self.dataLabel.append(label)
            self.dataPost.append(post)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        # ToolTips BUG
        if role == Qt.ToolTipRole:
            if index.column() == 0 and index.column() > 3:
                return QVariant(self.tr("Code_Saturne keyword: NBRVAF"))
            elif index.column() == 1 and index.column() > 3:
                return QVariant(self.tr("Code_Saturne keyword: IRAYVF"))
            elif index.column() == 1 and index.column() <= 3:
                return QVariant(self.tr("Code_Saturne keyword: IPSTYP/IPSTCL/IPSTFT/IPSTFO"))

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 1:
                return QVariant("Post-processing")

        # Display
        if role == Qt.DisplayRole:
            row = index.row()
            if index.column() == 0:
                return QVariant(self.dataLabel[row])
            else:
                return QVariant()

        # CheckState
        if role == Qt.CheckStateRole:
            row = index.row()
            if index.column() == 1:
                value = self.dataPost[row]
                if value == 'on':
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        if (index.row(), index.column()) in self.disableItem:
            return Qt.ItemIsSelectable
        elif index.column() == 0 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 1 :
            return  Qt.ItemIsEnabled | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return QVariant(self.tr("Name"))
            elif section == 1:
                return QVariant(self.tr("Post-\nprocessing"))
        return QVariant()


    def setData(self, index, value, role=None):
        row = index.row()
        if index.column() == 0:
            label = str(value.toString())
            if label == "": label = self.dataLabel[row]
            self.mdl.setPropertyLabel(self.dataLabel[row], label)
            self.dataLabel[row] = label

        elif index.column() == 1:
            v, ok = value.toInt()
            if v == Qt.Checked:
                self.dataPost[row] = "on"
            else:
                self.dataPost[row] = "off"

            self.mdl.setPostProcessing(self.dataLabel[row], self.dataPost[row])

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
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
        self.mdl = OutputSurfacicVariablesModel(self.case)

        self.modelOutput = StandardItemModelOutput(self.case)
        self.tableViewOutput.setModel(self.modelOutput)
        self.tableViewOutput.setAlternatingRowColors(True)
        self.tableViewOutput.resizeColumnToContents(0)
        self.tableViewOutput.resizeRowsToContents()
        self.tableViewOutput.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewOutput.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewOutput.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewOutput.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        labelDelegate = LabelDelegate(self.tableViewOutput, self.mdl)
        self.tableViewOutput.setItemDelegateForColumn(0, labelDelegate)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------