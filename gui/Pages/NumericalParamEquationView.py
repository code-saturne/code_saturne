# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2007 EDF S.A., France
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
- StandardItemModelScheme
- StandardItemModelSolver
- NumericalParamEquationView
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
from NumericalParamEquationForm import Ui_NumericalParamEquationForm
from Pages.NumericalParamEquationModel import NumericalParamEquatModel
from Pages.TurbulenceModel import TurbulenceModel
from Pages.SteadyManagementModel import SteadyManagementModel
import Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NumericalParamEquationView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Combo box delegate for ISCHCV
#-------------------------------------------------------------------------------

class SchemeOrderDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(SchemeOrderDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem(QString("Upwind"))
        editor.addItem(QString("Centered"))
        editor.addItem(QString("SOLU"))
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"upwind": 0, "centered": 1, "solu": 2}
        row = index.row()
        string = index.model().dataScheme[row]['ischcv']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, QVariant(value))

#-------------------------------------------------------------------------------
# Line edit delegate for BLENCV
#-------------------------------------------------------------------------------

class BlendingFactorDelegate(QItemDelegate):
    def __init__(self, parent=None, xml_model=None):
        super(BlendingFactorDelegate, self).__init__(parent)
        self.parent = parent
        self.turb = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        if self.turb.getTurbulenceModel() in ('LES_Smagorinsky', 'LES_dynamique'):
            validator = QtPage.DoubleValidator(editor, min=0.95, max=1.)
        else:
            validator = QtPage.DoubleValidator(editor, min=0., max=1.)
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
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, QVariant(value))

#-------------------------------------------------------------------------------
# Delegate for Solver QTableView
#-------------------------------------------------------------------------------

class SolveurDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(SolveurDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        if index.column() == 2:
            validator = QtPage.DoubleValidator(editor, min=0., max=0.01)
            validator.setExclusiveMin(True)
        else:
            validator = QtPage.IntValidator(editor, min=1)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if index.column() == 2:
            value, ok = editor.text().toDouble()
        else:
            value, ok = editor.text().toInt()
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, QVariant(value))

#-------------------------------------------------------------------------------
# Scheme class
#-------------------------------------------------------------------------------

class StandardItemModelScheme(QStandardItemModel):

    def __init__(self, NPE):
        """
        """
        QStandardItemModel.__init__(self)
        self.NPE = NPE
        self.dataScheme = []
        # list of items to be disabled in the QTableView
        self.disabledItem = []
        self.populateModel()
        self.headers = [self.tr("Name"),
                        self.tr("Scheme"),
                        self.tr("Blending\nFactor"),
                        self.tr("Slope\nTest"),
                        self.tr("Flux\nReconstruction")]
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                if column < 3:
                    role = Qt.DisplayRole
                else:
                    role = Qt.CheckStateRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)


    def populateModel(self):
        #FIXME
        #self.dicoV2M= {"Upwind" : 'upwind', self.tr("Centered"): 'centered', "SOLU": 'solu'}
        #self.dicoM2V= {"upwind" : 'Upwind', "centered": self.tr('Centered'), "solu": 'SOLU'}
        self.dicoV2M= {"Upwind" : 'upwind', "Centered": 'centered', "SOLU": 'solu'}
        self.dicoM2V= {"upwind" : 'Upwind', "centered": 'Centered', "solu": 'SOLU'}

        for label in self.NPE.getSchemeList():
            dico           = {}
            dico['label']  = label
            dico['blencv'] = self.NPE.getBlendingFactor(label)
            dico['ischcv'] = self.NPE.getScheme(label)
            dico['isstpc'] = self.NPE.getSlopeTest(label)
            dico['ircflu'] = self.NPE.getFluxReconstruction(label)
            self.dataScheme.append(dico)
            log.debug("populateModel-> dataScheme = %s" % dico)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        row = index.row()
        dico = self.dataScheme[row]

        if role == Qt.ToolTipRole:
            if index.column() == 1:
                return QVariant(self.tr("Code_Saturne keyword: ISCHCV"))
            elif index.column() == 2:
                return QVariant(self.tr("Code_Saturne keyword: BLENCV"))
            elif index.column() == 3:
                return QVariant(self.tr("Code_Saturne keyword: ISSTPC"))
            elif index.column() == 4:
                return QVariant(self.tr("Code_Saturne keyword: IRCFLU"))

        elif role == Qt.DisplayRole:
            if index.column() == 0:
                return QVariant(dico['label'])
            elif index.column() == 1:
                return QVariant(self.dicoM2V[dico['ischcv']])
            elif index.column() == 2:
                return QVariant(dico['blencv'])
            else:
                return QVariant()

        elif role == Qt.CheckStateRole:
            if index.column() == 3:
                if dico['isstpc'] == 'on':
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)
            elif index.column() == 4:
                if dico['ircflu'] == 'on':
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsSelectable

        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 1 or index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 3 or index.column() == 4:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role=None):
        row = index.row()
        label = self.dataScheme[row]['label']

        # set ISCHCV
        if index.column() == 1:
            self.dataScheme[row]['ischcv'] = self.dicoV2M[str(value.toString())]

            if self.dataScheme[row]['ischcv'] == "upwind":
                if (row, 2) not in self.disabledItem:
                    self.disabledItem.append((row, 2))
                if (row, 3) not in self.disabledItem:
                    self.disabledItem.append((row, 3))
                self.dataScheme[row]['blencv'] = 0.0
                self.dataScheme[row]['isstpc'] = "off"
            else:
                if (row, 2) in self.disabledItem:
                    self.disabledItem.remove((row, 2))
                    self.dataScheme[row]['blencv'] = 1.0
                if (row, 3) in self.disabledItem:
                    self.disabledItem.remove((row, 3))
                    self.dataScheme[row]['isstpc'] = "on"

            self.NPE.setScheme(label, self.dataScheme[row]['ischcv'])
            self.NPE.setBlendingFactor(label, self.dataScheme[row]['blencv'])

        # set BLENCV
        elif index.column() == 2:
            if self.dataScheme[row]['ischcv'] != "upwind":
                self.dataScheme[row]['blencv'], ok = value.toDouble()
                self.NPE.setBlendingFactor(label, self.dataScheme[row]['blencv'])

        # set ISSTPC
        elif index.column() == 3:
            if self.dataScheme[row]['ischcv'] != "upwind":
                v, ok = value.toInt()
                if v == Qt.Unchecked:
                    self.dataScheme[row]['isstpc'] = "off"
                else:
                    self.dataScheme[row]['isstpc'] = "on"
                self.NPE.setSlopeTest(label, self.dataScheme[row]['isstpc'])

        # set IRCFLU
        elif index.column() == 4:
            v, ok = value.toInt()
            if v == Qt.Unchecked:
                self.dataScheme[row]['ircflu'] = "off"
            else:
                self.dataScheme[row]['ircflu'] = "on"
            self.NPE.setFluxReconstruction(label, self.dataScheme[row]['ircflu'])

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# Solver class
#-------------------------------------------------------------------------------

class StandardItemModelSolver(QStandardItemModel):
    """
    Model associated with a QTableView. 
    """
    def __init__(self, NPE, SM):
        """
        """
        QStandardItemModel.__init__(self)
        self.NPE = NPE
        self.SM = SM
        self.setColumnCount(4)
        self.dataSolver = []
        # list of items to be disabled in the view
        self.disabledItem = []
        self.populateModel()


    def populateModel(self):
        for label in self.NPE.getSolveurList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico           = {}
            dico['label']  = label
            dico['nitmax'] = self.NPE.getMaxIterNumber(label)
            dico['epsilo'] = self.NPE.getSolveurPrecision(label)
            if self.NPE.isScalar(label):
                dico['cdtvar'] = self.NPE.getScalarTimeStepFactor(label)
            else:
                dico['cdtvar'] = ""
                self.disabledItem.append((row,3))

            self.dataSolver.append(dico)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):

        if not index.isValid():
            return QVariant()

        if role == Qt.ToolTipRole:
            if index.column() == 1:
                return QVariant(self.tr("Code_Saturne keyword: NITMAX"))
            elif index.column() == 2:
                return QVariant(self.tr("Code_Saturne keyword: EPSILO"))
            elif index.column() == 3:
                return QVariant(self.tr("Code_Saturne keyword: CDTVAR"))

        elif role == Qt.DisplayRole:
            row = index.row()
            dico = self.dataSolver[row]

            if index.column() == 0:
                return QVariant(dico['label'])
            elif index.column() == 1:
                return QVariant(dico['nitmax'])
            elif index.column() == 2:
                return QVariant(dico['epsilo'])
            elif index.column() == 3:
                return QVariant(dico['cdtvar'])
            else:
                return QVariant()

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()
 
 
    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return QVariant(self.tr("Name"))
            elif section == 1:
                return QVariant(self.tr("Maximum\nIteration Number"))
            elif section == 2:
                return QVariant(self.tr("Solver\nPrecision"))
            elif section == 3:
                return QVariant(self.tr("Time Step\nFactor"))
            else:
                return QVariant()
        return QVariant()


    def setData(self, index, value, role=None):
        row = index.row()
        label = self.dataSolver[row]['label']
        
        if index.column() == 1:
            self.dataSolver[row]['nitmax'], ok = value.toInt()
            self.NPE.setMaxIterNumber(label, self.dataSolver[row]['nitmax'])

        elif index.column() == 2:
            self.dataSolver[row]['epsilo'], ok = value.toDouble()
            self.NPE.setSolveurPrecision(label, self.dataSolver[row]['epsilo'])

        elif index.column() == 3:
            self.dataSolver[row]['cdtvar'], ok = value.toDouble()
            self.NPE.setScalarTimeStepFactor(label, self.dataSolver[row]['cdtvar'])

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class NumericalParamEquationView(QWidget, Ui_NumericalParamEquationForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_NumericalParamEquationForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.NPE = NumericalParamEquatModel(self.case)
        self.SM  = SteadyManagementModel(self.case)
        self.turb = TurbulenceModel(self.case)

        # Scheme
        self.modelScheme = StandardItemModelScheme(self.NPE)
        self.tableViewScheme.setModel(self.modelScheme)
        self.tableViewScheme.setAlternatingRowColors(True)
        self.tableViewScheme.resizeColumnToContents(0)
        self.tableViewScheme.resizeRowsToContents()
        self.tableViewScheme.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewScheme.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewScheme.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewScheme.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegateISCHCV = SchemeOrderDelegate(self.tableViewScheme)
        self.tableViewScheme.setItemDelegateForColumn(1, delegateISCHCV)

        delegateBLENCV = BlendingFactorDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(2, delegateBLENCV)

        # Solveur
        self.modelSolver = StandardItemModelSolver(self.NPE, self.SM)
        self.tableViewSolver.setModel(self.modelSolver)
        self.tableViewSolver.setAlternatingRowColors(True)
        self.tableViewSolver.resizeColumnToContents(0)
        self.tableViewSolver.resizeRowsToContents()
        self.tableViewSolver.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewSolver.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewSolver.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewSolver.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        if self.SM.getSteadyFlowManagement() == 'on':
            self.tableViewSolver.setColumnHidden(3, True)

        delegate = SolveurDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegate(delegate)


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
