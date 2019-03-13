# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, IntValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.NumericalParamEquationForm import Ui_NumericalParamEquationForm
from code_saturne.model.NumericalParamEquationModel import NumericalParamEquationModel
from code_saturne.model.TurbulenceModel import TurbulenceModel

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
        editor.addItem("Automatic")
        editor.addItem("Upwind")
        editor.addItem("Centered")
        editor.addItem("SOLU")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"automatic": 0, "upwind": 1, "centered": 2, "solu": 3}
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

        mg = index.model().dataSolver[index.row()]['mg']
        editor.addItem("Automatic")
        editor.addItem("Conjugate gradient")
        editor.addItem("Flexible conjugate gradient")
        editor.addItem("Inexact conjugate gradient")
        editor.addItem("Jacobi")
        editor.addItem("BiCGstab")
        editor.addItem("BiCGstab2")
        editor.addItem("GMRES")
        editor.addItem("Gauss Seidel")
        editor.addItem("Symmetric Gauss Seidel")
        editor.addItem("conjugate residual")
        if mg:
            editor.addItem("Multigrid, V-cycle")
            editor.addItem("Multigrid, K-cycle")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"automatic": 0,
                "conjugate_gradient": 1,
                "flexible_conjugate_gradient": 2,
                "inexact_conjugate_gradient": 3,
                "jacobi": 4,
                "bi_cgstab": 5,
                "bi_cgstab2": 6,
                "gmres": 7,
                "gauss_seidel": 8,
                "symmetric_gauss_seidel": 9,
                "PCR3": 10,
                "multigrid": 11,
                "multigrid_k_cycle": 12}
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
# Combo box delegate for preconditioning
#-------------------------------------------------------------------------------

class PreconditioningChoiceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(PreconditioningChoiceDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        editor.addItem("Automatic")
        editor.addItem("None")
        editor.addItem("Multigrid, V-cycle")
        editor.addItem("Multigrid, K-cycle")
        editor.addItem("Jacobi")
        editor.addItem("Polynomial")

        solver = index.model().dataSolver[index.row()]['iresol']
        if solver in ('multigrid', 'multigrid_k_cycle'):
            editor.model().item(1).setEnabled(False)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        string = index.model().dataSolver[row]['precond']
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

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
            validator = DoubleValidator(editor, min=0.95, max=1.)
        else:
            validator = DoubleValidator(editor, min=0., max=1.)
            # validator.setExclusiveMin(True)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
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
        validator = IntValidator(editor, min=1)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = from_qvariant(editor.text(), float)
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
            validator = DoubleValidator(editor, min=0., max=0.01)
            validator.setExclusiveMin(True)
        elif (index.column() == 2 or index.column() == 4):
            validator = IntValidator(editor, min=1)
        elif index.column() == 5:
            validator = DoubleValidator(editor, min=0.)
            validator.setExclusiveMin(True)
        editor.setValidator(validator)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            if index.column() == 3 or index.column() == 5:
                value = from_qvariant(editor.text(), float)
            elif (index.column() == 2 or index.column() == 4):
                value = from_qvariant(editor.text(), int)
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
        self.keys = ['name', 'ischcv', 'blencv', 'isstpc', 'ircflu', 'nswrsm']
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
        self.dicoV2M= {"Automatic" : 'automatic', "Upwind" : 'upwind', "Centered": 'centered', "SOLU": 'solu'}
        self.dicoM2V= {'automatic' : "Automatic", "upwind" : 'Upwind', "centered": 'Centered', "solu": 'SOLU'}

        for name in self.NPE.getSchemeList():
            dico           = {}
            dico['name']  = name
            dico['blencv'] = self.NPE.getBlendingFactor(name)
            dico['ischcv'] = self.NPE.getScheme(name)
            dico['isstpc'] = self.NPE.getSlopeTest(name)
            dico['ircflu'] = self.NPE.getFluxReconstruction(name)
            dico['nswrsm'] = self.NPE.getRhsReconstruction(name)
            self.dataScheme.append(dico)
            log.debug("populateModel-> dataScheme = %s" % dico)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        column = index.column()
        dico = self.dataScheme[row]
        key = self.keys[column]

        if dico[key] == None:
            return None

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

        return None


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
        return None


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        name = self.dataScheme[row]['name']

        # for Pressure, most fields are empty
        if column > 0 and str(from_qvariant(value, to_text_string)) in ['', 'None']:
            if (row, column) not in self.disabledItem:
                self.disabledItem.append((row, column))
            return False

        # set ISCHCV
        if column == 1:
            self.dataScheme[row]['ischcv'] = self.dicoV2M[str(from_qvariant(value, to_text_string))]
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

            self.NPE.setScheme(name, self.dataScheme[row]['ischcv'])
            self.NPE.setBlendingFactor(name, self.dataScheme[row]['blencv'])

        # set BLENCV
        elif column == 2:
            if self.dataScheme[row]['ischcv'] != "upwind":
                self.dataScheme[row]['blencv'] = from_qvariant(value, float)
                self.NPE.setBlendingFactor(name, self.dataScheme[row]['blencv'])

        # set ISSTPC
        elif column == 3:
            if self.dataScheme[row]['ischcv'] != "upwind":
                v = from_qvariant(value, int)
                if v == Qt.Unchecked:
                    self.dataScheme[row]['isstpc'] = "off"
                else:
                    self.dataScheme[row]['isstpc'] = "on"
                self.NPE.setSlopeTest(name, self.dataScheme[row]['isstpc'])

        # set IRCFLU
        elif column == 4:
            v = from_qvariant(value, int)
            if v == Qt.Unchecked:
                self.dataScheme[row]['ircflu'] = "off"
            else:
                self.dataScheme[row]['ircflu'] = "on"
            self.NPE.setFluxReconstruction(name, self.dataScheme[row]['ircflu'])

        # set NSWRSM
        elif column == 5:
            self.dataScheme[row]['nswrsm'] = from_qvariant(value, int)
            self.NPE.setRhsReconstruction(name, self.dataScheme[row]['nswrsm'])

        self.dataChanged.emit(index, index)
        return True

#-------------------------------------------------------------------------------
# Solver class
#-------------------------------------------------------------------------------

class StandardItemModelSolver(QStandardItemModel):
    """
    Model associated with a QTableView.
    """
    def __init__(self, NPE):
        """
        """
        QStandardItemModel.__init__(self)
        self.NPE = NPE
        self.setColumnCount(6)
        self.dataSolver = []
        # list of items to be disabled in the view
        self.disabledItem = []
        self.populateModel()


    def populateModel(self):
        self.dicoV2M= {"Multigrid, V-cycle"     : 'multigrid',
                       "Multigrid, K-cycle"     : 'multigrid_k_cycle',
                       "Conjugate gradient"     : 'conjugate_gradient',
                       "Flexible conjugate gradient" : 'flexible_conjugate_gradient',
                       "Inexact conjugate gradient"  : 'inexact_conjugate_gradient',
                       "Jacobi"                 : 'jacobi',
                       "BiCGstab"               : 'bi_cgstab',
                       "BiCGstab2"              : 'bi_cgstab2',
                       "GMRES"                  : 'gmres',
                       "Automatic"              : "automatic",
                       "Gauss Seidel"           : "gauss_seidel",
                       "Symmetric Gauss Seidel" : "symmetric_gauss_seidel",
                       "conjugate residual"     : "PCR3",
                       "None"                   : "none",
                       "Polynomial"             : "polynomial"}
        self.dicoM2V= {"multigrid"              : 'Multigrid, V-cycle',
                       "multigrid_k_cycle"      : 'Multigrid, K-cycle',
                       "conjugate_gradient"     : 'Conjugate gradient',
                       "inexact_conjugate_gradient"  : 'Inexact conjugate gradient',
                       "flexible_conjugate_gradient" : 'Flexible conjugate gradient',
                       "jacobi"                 : 'Jacobi',
                       "bi_cgstab"              : 'BiCGstab',
                       "bi_cgstab2"             : 'BiCGstab2',
                       'gmres'                  : "GMRES",
                       "automatic"              : "Automatic",
                       "gauss_seidel"           : "Gauss Seidel",
                       "symmetric_gauss_seidel" : "Symmetric Gauss Seidel",
                       "PCR3"                   : "conjugate residual",
                       "none"                   : "None",
                       "polynomial"             : "Polynomial"}

        for name in self.NPE.getSolverList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico            = {}
            dico['name']    = name
            dico['mg']      = self.NPE.getSolverAllowMultigrid(name)
            dico['iresol']  = self.NPE.getSolverChoice(name)
            dico['precond'] = self.NPE.getPreconditioningChoice(name)
            dico['epsilo']  = self.NPE.getSolverPrecision(name)
            dico['verbo']   = self.NPE.getVerbosity(name)
            if self.NPE.isScalar(name):
                dico['cdtvar'] = self.NPE.getScalarTimeStepFactor(name)
            else:
                dico['cdtvar'] = ""
                self.disabledItem.append((row,5))

            self.dataSolver.append(dico)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):

        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            if index.column() == 3:
                return self.tr("Field variable calculation option keyword: epsilo")
            elif index.column() == 4:
                return self.tr("Field variable calculation option keyword: iwarni")
            elif index.column() == 5:
                return self.tr("Code_Saturne keyword: CDTVAR")

        elif role == Qt.DisplayRole:
            row = index.row()
            dico = self.dataSolver[row]

            if index.column() == 0:
                return dico['name']
            elif index.column() == 1:
                return self.dicoM2V[dico['iresol']]
            elif index.column() == 2:
                return self.dicoM2V[dico['precond']]
            elif index.column() == 3:
                return dico['epsilo']
            elif index.column() == 4:
                return dico['verbo']
            elif index.column() == 5:
                return dico['cdtvar']
            else:
                return None

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


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
                return self.tr("Preconditioning\nChoice")
            elif section == 3:
                return self.tr("Solver\nPrecision")
            elif section == 4:
                return self.tr("Verbosity")
            elif section == 5:
                return self.tr("Time Step\nFactor")
            else:
                return None
        return None


    def setData(self, index, value, role=None):
        row = index.row()
        name = self.dataSolver[row]['name']

        if index.column() == 1:
            self.dataSolver[row]['iresol'] = self.dicoV2M[from_qvariant(value, to_text_string)]
            self.NPE.setSolverChoice(name, self.dataSolver[row]['iresol'])

        elif index.column() == 2:
            self.dataSolver[row]['precond'] = self.dicoV2M[from_qvariant(value, to_text_string)]
            self.NPE.setPreconditioningChoice(name, self.dataSolver[row]['precond'])

        elif index.column() == 3:
            self.dataSolver[row]['epsilo'] = from_qvariant(value, float)
            self.NPE.setSolverPrecision(name, self.dataSolver[row]['epsilo'])

        elif index.column() == 4:
            self.dataSolver[row]['verbo'] = from_qvariant(value, int)
            self.NPE.setVerbosity(name, self.dataSolver[row]['verbo'])

        elif index.column() == 5:
            self.dataSolver[row]['cdtvar'] = from_qvariant(value, float)
            self.NPE.setScalarTimeStepFactor(name, self.dataSolver[row]['cdtvar'])

        self.dataChanged.emit(index, index)
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
        v = DoubleValidator(editor)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    maxi = model.getData(idx)['scamax']
                    name = model.getData(idx)['name']
                    if model.checkMinMax(name, value, maxi):
                        model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for verbosity
#-------------------------------------------------------------------------------

class VerbosityDelegate(QItemDelegate):
    def __init__(self, parent=None, xml_model=None):
        super(VerbosityDelegate, self).__init__(parent)
        self.parent = parent
        self.turb = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0)
        validator.setExclusiveMin(False)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value)


#-------------------------------------------------------------------------------
# Line edit delegate for maximum value
#-------------------------------------------------------------------------------

class MaximumDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(MaximumDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    mini = model.getData(idx)['scamin']
                    name = model.getData(idx)['name']
                    if model.checkMinMax(name, mini, value):
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
        for name in self.NPE.getClippingList():
            row = self.rowCount()
            self.setRowCount(row + 1)
            dico             = {}
            dico['name']    = name
            dico['scamin']   = self.NPE.getMinValue(name)
            dico['scamax']   = self.NPE.getMaxValue(name)

            self._data.append(dico)
            log.debug("populateModel-> _data = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        if role == Qt.ToolTipRole:
            return self.toolTipRole[col]
        if role == Qt.DisplayRole:
            row = index.row()
            dico = self._data[row]
            if col == 0:
                return dico['name']
            elif col == 1:
                return dico['scamin']
            elif col == 2:
                return dico['scamax']
            else:
                return None
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


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
        return None


    def setData(self, index, value, role=None):
        if not index.isValid():
            return Qt.ItemIsEnabled
        row = index.row()
        name = self._data[row]['name']

        if index.column() == 1:
            self._data[row]['scamin'] = from_qvariant(value, float)
            self.NPE.setMinValue(name, self._data[row]['scamin'])

        elif index.column() == 2:
            self._data[row]['scamax'] = from_qvariant(value, float)
            self.NPE.setMaxValue(name, self._data[row]['scamax'])

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def checkMinMax(self, name, mini, maxi):
        """
        Verify the coherence between mini and maxi
        """
        log.debug("checkMinMax")
        OK = 1
        if mini > maxi:
            title = self.tr("Information")
            msg = self.tr("The minimal value is greater than the maximal "\
                          "value. Therefore there will be no clipping for the "\
                          "scalar named:\n\n%1").arg(name)
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
        self.NPE = NumericalParamEquationModel(self.case)
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
        if QT_API == "PYQT4":
            self.tableViewScheme.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewScheme.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateISCHCV = SchemeOrderDelegate(self.tableViewScheme)
        self.tableViewScheme.setItemDelegateForColumn(1, delegateISCHCV)

        delegateBLENCV = BlendingFactorDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(2, delegateBLENCV)

        delegateNSWRSM = RhsReconstructionDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(5, delegateNSWRSM)

        # Solver
        self.modelSolver = StandardItemModelSolver(self.NPE)
        self.tableViewSolver.setModel(self.modelSolver)
        self.tableViewSolver.setAlternatingRowColors(True)
        self.tableViewSolver.resizeColumnToContents(0)
        self.tableViewSolver.resizeRowsToContents()
        self.tableViewSolver.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewSolver.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewSolver.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewSolver.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSolver.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        from code_saturne.model.TimeStepModel import TimeStepModel
        idtvar = TimeStepModel(self.case).getTimePassing()
        if idtvar in [-1, 2]:
            self.tableViewSolver.setColumnHidden(5, True)

        delegate = SolverDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegate(delegate)

        delegateSolverChoice = SolverChoiceDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(1, delegateSolverChoice)

        delegatePrecondChoice = PreconditioningChoiceDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(2, delegatePrecondChoice)

        delegateVerbosity = VerbosityDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(4, delegateVerbosity)

        # Clipping
        self.modelClipping = StandardItemModelClipping(self, self.NPE)
        self.tableViewClipping.setModel(self.modelClipping)
        self.tableViewClipping.setAlternatingRowColors(True)
        self.tableViewClipping.resizeColumnToContents(0)
        self.tableViewClipping.resizeRowsToContents()
        self.tableViewClipping.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewClipping.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewClipping.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewClipping.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewClipping.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateMin = MinimumDelegate(self.tableViewClipping)
        self.tableViewClipping.setItemDelegateForColumn(1, delegateMin)

        delegateMax = MaximumDelegate(self.tableViewClipping)
        self.tableViewClipping.setItemDelegateForColumn(2, delegateMax)

        if len(self.NPE.getClippingList()) == 0:
            self.tab_clipping.setEnabled(False)

        self.tabWidgetScheme.setCurrentIndex(self.case['current_tab'])

        self.tabWidgetScheme.currentChanged[int].connect(self.slotchanged)

        self.case.undoStartGlobal()


    @pyqtSlot(int)
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
