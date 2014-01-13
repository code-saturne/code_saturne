# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
import sys
if sys.version_info[0] == 2:
    import sip
    sip.setapi('QString', 2)

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Pages.NumericalParamEquationForm import Ui_NumericalParamEquationForm
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
        editor.addItem("Upwind")
        editor.addItem("Centered")
        editor.addItem("SOLU")
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
                model.setData(idx, value)

#-------------------------------------------------------------------------------
# Combo box delegate for IRESOL
#-------------------------------------------------------------------------------

class SolverChoiceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(SolverChoiceDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        editor.addItem("Automatic")
        editor.addItem("Conjugate gradient")
        editor.addItem("Jacobi")
        editor.addItem("BI-CGSTAB")
        editor.addItem("GMRES")
        editor.addItem("Multigrid")
        editor.installEventFilter(self)
        if index.row() > 0:
            editor.removeItem(4)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"automatic": 0, "conjugate_gradient": 1, "jacobi": 2, "bi_cgstab": 3, "gmres": 4, "multigrid": 5}
        row = index.row()
        string = index.model().dataSolver[row]['iresol']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value)

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
        if self.turb.getTurbulenceModel() in ('LES_Smagorinsky', 'LES_dynamique', 'LES_WALE'):
            validator = QtPage.DoubleValidator(editor, min=0.95, max=1.)
        else:
            validator = QtPage.DoubleValidator(editor, min=0., max=1.)
            validator.setExclusiveMin(True)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        value = str(index.model().data(index, Qt.DisplayRole))
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = float(editor.text())
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value)

#-------------------------------------------------------------------------------
# Line edit delegate for nswrsm
#-------------------------------------------------------------------------------

class RhsReconstructionDelegate(QItemDelegate):
    def __init__(self, parent=None, xml_model=None):
        super(RhsReconstructionDelegate, self).__init__(parent)
        self.parent = parent
        self.turb = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = QtPage.IntValidator(editor, min=1)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        value = str(index.model().data(index, Qt.DisplayRole))
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = float(editor.text())
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value)

#-------------------------------------------------------------------------------
# Delegate for Solver QTableView
#-------------------------------------------------------------------------------

class SolverDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(SolverDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        if index.column() == 3:
            validator = QtPage.DoubleValidator(editor, min=0., max=0.01)
            validator.setExclusiveMin(True)
        elif (index.column() == 2 or index.column() == 4):
            validator = QtPage.IntValidator(editor, min=1)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = str(index.model().data(index, Qt.DisplayRole))
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            if index.column() == 3:
                value = float(editor.text())
            elif (index.column() == 2 or index.column() == 4):
                value = int(editor.text())
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value)

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
                        self.tr("Flux\nReconstruction"),
                        self.tr("RHS Sweep\nReconstruction")]
        self.keys = ['label', 'ischcv', 'blencv', 'isstpc', 'ircflu', 'nswrsm']
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                if column == 1 or column == 2 or column == 5:
                    role = Qt.DisplayRole
                else:
                    role = Qt.CheckStateRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)


    def populateModel(self):
        self.dicoV2M= {"Upwind" : 'upwind', "Centered": 'centered', "SOLU": 'solu'}
        self.dicoM2V= {"upwind" : 'Upwind', "centered": 'Centered', "solu": 'SOLU'}

        for label in self.NPE.getSchemeList():
            dico           = {}
            dico['label']  = label
            dico['blencv'] = self.NPE.getBlendingFactor(label)
            dico['ischcv'] = self.NPE.getScheme(label)
            dico['isstpc'] = self.NPE.getSlopeTest(label)
            dico['ircflu'] = self.NPE.getFluxReconstruction(label)
            dico['nswrsm'] = self.NPE.getRhsReconstruction(label)
            self.dataScheme.append(dico)
            log.debug("populateModel-> dataScheme = %s" % dico)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return

        row = index.row()
        column = index.column()
        dico = self.dataScheme[row]
        key = self.keys[column]

        if dico[key] == None:
            return

        if role == Qt.ToolTipRole:
            if index.column() > 0:
                return self.tr("Code_Saturne keyword: " + key.upper())

        elif role == Qt.DisplayRole and not column in [3, 4]:
            if key == 'ischcv':
                return self.dicoM2V[dico[key]]
            else:
                return dico[key]

        elif role == Qt.CheckStateRole and column in [3, 4]:
            st = None
            if key in ['isstpc', 'ircflu']:
                st = dico[key]
            if st == 'on':
                return Qt.Checked
            else:
                return Qt.Unchecked

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return


    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.NoItemFlags

        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 1 or index.column() == 2 or index.column() == 5:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 3 or index.column() == 4:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        label = self.dataScheme[row]['label']

        # for Pressure, most fields are empty
        ## if column > 0 and str(value.toString()) in ['', 'None']:
        if column > 0 and str(value) in ['', 'None']:
            if (row, column) not in self.disabledItem:
                self.disabledItem.append((row, column))
            return False

        # set ISCHCV
        if column == 1:
            self.dataScheme[row]['ischcv'] = self.dicoV2M[str(value)]
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
        elif column == 2:
            if self.dataScheme[row]['ischcv'] != "upwind":
                self.dataScheme[row]['blencv'] = float(value)
                self.NPE.setBlendingFactor(label, self.dataScheme[row]['blencv'])

        # set ISSTPC
        elif column == 3:
            if self.dataScheme[row]['ischcv'] != "upwind":
                v = int(value)
                if v == Qt.Unchecked:
                    self.dataScheme[row]['isstpc'] = "off"
                else:
                    self.dataScheme[row]['isstpc'] = "on"
                self.NPE.setSlopeTest(label, self.dataScheme[row]['isstpc'])

        # set IRCFLU
        elif column == 4:
            v = int(value)
            if v == Qt.Unchecked:
                self.dataScheme[row]['ircflu'] = "off"
            else:
                self.dataScheme[row]['ircflu'] = "on"
            self.NPE.setFluxReconstruction(label, self.dataScheme[row]['ircflu'])

        # set NSWRSM
        elif column == 5:
            self.dataScheme[row]['nswrsm'] = int(value)
            self.NPE.setRhsReconstruction(label, self.dataScheme[row]['nswrsm'])

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
        self.setColumnCount(5)
        self.dataSolver = []
        # list of items to be disabled in the view
        self.disabledItem = []
        self.populateModel()


    def populateModel(self):
        self.dicoV2M= {"Multigrid": 'multigrid',"Conjugate gradient" : 'conjugate_gradient',
                       "Jacobi": 'jacobi', "BI-CGSTAB": 'bi_cgstab', "GMRES": 'gmres',
                       "Automatic": "automatic"}
        self.dicoM2V= {"multigrid" : 'Multigrid',"conjugate_gradient" : 'Conjugate gradient',
                       "jacobi": 'Jacobi', "bi_cgstab": 'BI-CGSTAB', 'gmres': "GMRES",
                       "automatic": "Automatic"}
        for label in self.NPE.getSolverList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico           = {}
            dico['label']  = label
            dico['iresol'] = self.NPE.getSolverChoice(label)
            dico['nitmax'] = self.NPE.getMaxIterNumber(label)
            dico['epsilo'] = self.NPE.getSolverPrecision(label)
            if self.NPE.isScalar(label):
                dico['cdtvar'] = self.NPE.getScalarTimeStepFactor(label)
            else:
                dico['cdtvar'] = ""
                self.disabledItem.append((row,4))

            self.dataSolver.append(dico)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):

        if not index.isValid():
            return

        if role == Qt.ToolTipRole:
            if index.column() == 1:
                return self.tr("Code_Saturne keyword: IRESOL")
            elif index.column() == 2:
                return self.tr("Code_Saturne keyword: NITMAX")
            elif index.column() == 3:
                return self.tr("Code_Saturne keyword: EPSILO")
            elif index.column() == 4:
                return self.tr("Code_Saturne keyword: CDTVAR")

        elif role == Qt.DisplayRole:
            row = index.row()
            dico = self.dataSolver[row]

            if index.column() == 0:
                return dico['label']
            elif index.column() == 1:
                return self.dicoM2V[dico['iresol']]
            elif index.column() == 2:
                return dico['nitmax']
            elif index.column() == 3:
                return dico['epsilo']
            elif index.column() == 4:
                return dico['cdtvar']
            else:
                return

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return


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
                return self.tr("Name")
            elif section == 1:
                return self.tr("Solver\nChoice")
            elif section == 2:
                return self.tr("Maximum\nIteration Number")
            elif section == 3:
                return self.tr("Solver\nPrecision")
            elif section == 4:
                return self.tr("Time Step\nFactor")
            else:
                return
        return


    def setData(self, index, value, role=None):
        row = index.row()
        label = self.dataSolver[row]['label']

        if index.column() == 1:
            self.dataSolver[row]['iresol'] = self.dicoV2M[str(value)]
            self.NPE.setSolverChoice(label, self.dataSolver[row]['iresol'])

        elif index.column() == 2:
            self.dataSolver[row]['nitmax'] = int(value)
            self.NPE.setMaxIterNumber(label, self.dataSolver[row]['nitmax'])

        elif index.column() == 3:
            self.dataSolver[row]['epsilo'] = float(value)
            self.NPE.setSolverPrecision(label, self.dataSolver[row]['epsilo'])

        elif index.column() == 4:
            self.dataSolver[row]['cdtvar'] = float(value)
            self.NPE.setScalarTimeStepFactor(label, self.dataSolver[row]['cdtvar'])

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# Line edit delegate for minimum value
#-------------------------------------------------------------------------------

class MinimumDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(MinimumDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = QtPage.DoubleValidator(editor)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        value = str(index.model().data(index, Qt.DisplayRole))
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = float(editor.text())
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    maxi = model.getData(idx)['scamax']
                    label = model.getData(idx)['label']
                    if model.checkMinMax(label, value, maxi):
                        model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for maximum value
#-------------------------------------------------------------------------------

class MaximumDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(MaximumDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = QtPage.DoubleValidator(editor)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        value = str(index.model().data(index, Qt.DisplayRole))
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = float(editor.text())
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    mini = model.getData(idx)['scamin']
                    label = model.getData(idx)['label']
                    if model.checkMinMax(label, mini, value):
                        model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelClipping(QStandardItemModel):
    """
    """
    def __init__(self, parent, NPE):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Name"),
                        self.tr("Minimal\nvalue"),
                        self.tr("Maximal\nvalue")]

        self.setColumnCount(len(self.headers))

        self.toolTipRole = [self.tr("Code_Saturne keyword: NSCAUS"),
                            self.tr("Code_Saturne keyword: SCAMIN"),
                            self.tr("Code_Saturne keyword: SCAMAX")]

        self._data = []
        self._disable = []
        self.parent = parent
        self.NPE = NPE
        self.populateModel()

    def populateModel(self):
        for label in self.NPE.getClippingList():
            row = self.rowCount()
            self.setRowCount(row + 1)
            dico             = {}
            dico['label']    = label
            dico['scamin']   = self.NPE.getMinValue(label)
            dico['scamax']   = self.NPE.getMaxValue(label)

            self._data.append(dico)
            log.debug("populateModel-> _data = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return

        row = index.row()
        col = index.column()

        if role == Qt.ToolTipRole:
            return self.toolTipRole[col]
        if role == Qt.DisplayRole:
            row = index.row()
            dico = self._data[row]
            if col == 0:
                return dico['label']
            elif col == 1:
                return dico['scamin']
            elif col == 2:
                return dico['scamax']
            else:
                return
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return


    def setData(self, index, value, role=None):
        if not index.isValid():
            return Qt.ItemIsEnabled
        row = index.row()
        label = self._data[row]['label']

        if index.column() == 1:
            self._data[row]['scamin'] = float(value)
            self.NPE.setMinValue(label, self._data[row]['scamin'])

        elif index.column() == 2:
            self._data[row]['scamax'] = float(value)
            self.NPE.setMaxValue(label, self._data[row]['scamax'])

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def checkMinMax(self, label, mini, maxi):
        """
        Verify the coherence between mini and maxi
        """
        log.debug("checkMinMax")
        OK = 1
        if mini > maxi:
            title = self.tr("Information")
            msg = self.tr("The minimal value is greater than the maximal "\
                          "value. Therefore there will be no clipping for the "\
                          "scalar named:\n\n%1").arg(label)
            QMessageBox.information(self.parent, title, msg)
            return OK

        return OK


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
        self.case.undoStopGlobal()
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

        delegateNSWRSM = RhsReconstructionDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(5, delegateNSWRSM)

        # Solver
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
            self.tableViewSolver.setColumnHidden(4, True)

        delegate = SolverDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegate(delegate)

        delegateSolverChoice = SolverChoiceDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(1, delegateSolverChoice)

        # Clipping
        self.modelClipping = StandardItemModelClipping(self, self.NPE)
        self.tableViewClipping.setModel(self.modelClipping)
        self.tableViewClipping.setAlternatingRowColors(True)
        self.tableViewClipping.resizeColumnToContents(0)
        self.tableViewClipping.resizeRowsToContents()
        self.tableViewClipping.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewClipping.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewClipping.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewClipping.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegateMin = MinimumDelegate(self.tableViewClipping)
        self.tableViewClipping.setItemDelegateForColumn(1, delegateMin)

        delegateMax = MaximumDelegate(self.tableViewClipping)
        self.tableViewClipping.setItemDelegateForColumn(2, delegateMax)

        if len(self.NPE.getClippingList()) == 0:
            self.tab_clipping.setEnabled(False)

        self.tabWidgetScheme.setCurrentIndex(self.case['current_tab'])

        self.connect(self.tabWidgetScheme, SIGNAL("currentChanged(int)"), self.slotchanged)

        self.case.undoStartGlobal()


    @pyqtSignature("int")
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


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
