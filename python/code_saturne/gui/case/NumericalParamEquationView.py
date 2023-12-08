# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, IntValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.NumericalParamEquationForm import Ui_NumericalParamEquationForm
from code_saturne.model.NumericalParamEquationModel import NumericalParamEquationModel
from code_saturne.model.TurbulenceModel import TurbulenceModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NumericalParamEquationView")
log.setLevel(GuiParam.DEBUG)

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
        editor.addItem("GCR")
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
                "gcr": 8,
                "gauss_seidel": 9,
                "symmetric_gauss_seidel": 10,
                "PCR3": 11,
                "multigrid": 12,
                "multigrid_k_cycle": 13}
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
        editor.addItem("Multigrid, K-cycle, HPC")
        editor.addItem("Jacobi")
        editor.addItem("Polynomial")

        solver = index.model().dataSolver[index.row()]['iresol']
        if solver in ('multigrid', 'multigrid_k_cycle', 'multigrid_k_cycle_hpc'):
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
# Combo box delegate for ISCHCV
#-------------------------------------------------------------------------------

class SchemeOrderDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, update_layout=None):
        super(SchemeOrderDelegate, self).__init__(parent)
        self.parent = parent
        self.update_layout = update_layout


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        row = index.row()
        category = index.model().dataScheme[row]['category']
        if category > 0:
            editor.addItem("Automatic")
            editor.addItem("Centered")
            editor.addItem("SOLU (centered gradient)")
            editor.addItem("SOLU (upwind gradient)")
            editor.addItem("Blending (SOLU/centered)")
            if category != 2:
                editor.addItem("NVD/TVD")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"automatic": 0, "centered": 1, "solu": 2,
                "solu_upwind_gradient": 3, "blending": 4,
                "nvd_tvd": 5}
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
        if self.update_layout != None:
            self.update_layout()


#-------------------------------------------------------------------------------
# Combo box delegate for ISSTPC
#-------------------------------------------------------------------------------

class SchemeSlopeTestDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(SchemeSlopeTestDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        row = index.row()
        category = index.model().dataScheme[row]['category']
        if category > 0:
            if index.model().dataScheme[row]['ischcv'] not in ('blending', 'nvd_tvd'):
                editor.addItem("Enabled")
            editor.addItem("Disabled")
            if category != 2:
                editor.addItem("Beta limiter")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"on": 0, "off": 1, "beta_limiter": 2}
        row = index.row()
        string = index.model().dataScheme[row]['isstpc']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value)

#-------------------------------------------------------------------------------
# Combo box delegate for TVD/NVD limiter
#-------------------------------------------------------------------------------

class SchemeNVDLimiterDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(SchemeNVDLimiterDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        row = index.row()
        category = index.model().dataScheme[row]['category']
        if category > 0:
            editor.addItem("GAMMA")
            editor.addItem("SMART")
            editor.addItem("CUBISTA")
            editor.addItem("SUPERBEE")
            editor.addItem("MUSCL")
            editor.addItem("MINMOD")
            editor.addItem("CLAM")
            editor.addItem("STOIC")
            editor.addItem("OSHER")
            editor.addItem("WASEB")
            if category == 3:
                editor.addItem("HRIC")
                editor.addItem("CICSAM")
                editor.addItem("STACS")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"gamma": 0,
                "smart": 1,
                "cubista": 2,
                "superbee": 3,
                "muscl": 4,
                "minmod": 5,
                "clam": 6,
                "stoic": 7,
                "osher": 8,
                "waseb": 9,
                "hric": 10,
                "cicsam": 11,
                "stacs": 12}
        row = index.row()
        string = index.model().dataScheme[row]['nvd_limiter']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
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
        self.nvd_count = 0
        self.populateModel()
        self.headers = [self.tr("Name"),
                        self.tr("Scheme"),
                        self.tr("Centering\nBlend"),
                        self.tr("Slope\nTest"),
                        self.tr("NVD\nLimiter"),
                        self.tr("Flux\nReconstruction"),
                        self.tr("RHS Sweep\nReconstruction")]
        self.keys = ['name', 'ischcv', 'blencv', 'isstpc', 'nvd_limiter',
                     'ircflu', 'nswrsm']
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                if column in (1, 2, 3, 4, 6):
                    role = Qt.DisplayRole
                else:
                    role = Qt.CheckStateRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)

        self.tooltips = [
            self.tr("Equation parameter: 'ischcv'\n\n"
                    "Base (second-order) convective scheme"),
            self.tr("Equation parameter: 'blencv'\n\n"
                    "Blending of chosen convective scheme with upwind scheme\n"
                    "(0: full upwind, 1: pure base convection scheme)"),
            self.tr("Equation parameter: 'isstpc'\n\n"
                    "Enable slope test to switch to an upwind convective\n"
                    "scheme under certain conditions.\n"
                    "The use of the slope test stabilizes the calculation\n"
                    "but may reduce spatial convergence order."),
            self.tr("Field keyword: 'limiter_choice'\n\n"
                    "NVD limiter choice (for NVD/TVD convective scheme)."),
            self.tr("Equation parameter: 'ircflux'\n\n"
                    "Flux reconstruction: indicate whether the convective\n"
                    "and diffusive fluxes at the faces should be reconstructed\n"
                    "at non-orthogonal mesh faces.\n"
                    "Deactivating this reconstruction can have a stabilizing\n"
                    "effect on the calculation.\n"
                    "It is sometimes useful with the k−ε model, if the mesh\n"
                    "is strongly non-orthogonal in the near-wall region,\n"
                    " where the gradients of k and ε are strong"),
            self.tr("Equation parameter: 'nswrsm'\n\n"
                    "RHS Sweep Reconstruction: number of iterations for the\n"
                    "reconstruction of the right-hand sides of the equation")
        ]

    def populateModel(self):
        self.dicoV2M = {"Automatic": 'automatic',
                        "Centered": 'centered',
                        "SOLU (centered gradient)": 'solu',
                        "SOLU (upwind gradient)": 'solu_upwind_gradient',
                        "Blending (SOLU/centered)": 'blending',
                        "NVD/TVD": 'nvd_tvd'}
        self.dicoM2V = {'automatic': "Automatic",
                        'centered': "Centered",
                        'solu': "SOLU (centered gradient)",
                        'solu_upwind_gradient': "SOLU (upwind gradient)",
                        'blending': "Blending (SOLU/centered)",
                        'nvd_tvd': "NVD/TVD"}

        self.dicoV2M_isstpc = {"Enabled": 'on',
                               "Disabled": 'off',
                               "Beta limiter": 'beta_limiter'}
        self.dicoM2V_isstpc = {'on': "Enabled",
                               'off': "Disabled",
                               'beta_limiter': "Beta limiter"}

        for v in self.NPE.getSchemeList():
            name = v[0]
            dico           = {}
            dico['name']  = name
            dico['blencv'] = self.NPE.getBlendingFactor(name)
            dico['ischcv'] = self.NPE.getScheme(name)
            if dico['ischcv'] == 'nvd_tvd':
                self.nvd_count += 1
            dico['isstpc'] = self.NPE.getSlopeTest(name)
            dico['nvd_limiter'] = self.NPE.getNVDLimiter(name)
            dico['ircflu'] = self.NPE.getFluxReconstruction(name)
            dico['nswrsm'] = self.NPE.getRhsReconstruction(name)
            dico['category'] = v[1]
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

        if dico[key] is None:
            return None

        if role == Qt.ToolTipRole:
            col = index.column()
            if col > 0 and col < 7:
                return self.tooltips[col-1]
            elif col > 0:
                return self.tr("code_saturne keyword: " + key)

        elif role == Qt.DisplayRole and not column == 5:
            if key == 'ischcv':
                return self.dicoM2V[dico[key]]
            elif key == 'isstpc':
                if self.dataScheme[row]['blencv'] > 0:
                    return self.dicoM2V_isstpc[dico[key]]
                else:
                    return ""
            elif key == 'nvd_limiter':
                v = dico[key]
                if v is not None:
                    return v.upper()
                else:
                    return v
            else:
                return dico[key]

        elif role == Qt.CheckStateRole and column == 5:
            st = None
            if key in ['ircflu']:
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

        column = index.column()
        row = index.row()
        # disabled item
        if (row, column) in self.disabledItem:
            return Qt.NoItemFlags

        if column == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif column in (1, 2, 6):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif column == 3:
            if self.dataScheme[row]['blencv'] > 0:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.NoItemFlags
        elif column == 4:
            if self.dataScheme[row]['ischcv'] == 'nvd_tvd':
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.NoItemFlags
        elif column == 5:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]

        if role == Qt.ToolTipRole:
            if section == 0:
                return self.tr("variable or equation name")
            elif section < 7:
                return self.tooltips[section-1]

        return None


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        name = self.dataScheme[row]['name']

        # for Pressure, most fields are empty
        if column > 0 and str(from_qvariant(value, to_text_string)) in ['', 'None']:
            if column not in (3, 4):
                if (row, column) not in self.disabledItem:
                    self.disabledItem.append((row, column))
            else:
                if (row, 2) in self.disabledItem:
                    self.disabledItem.append((row, column))
            return False

        # set ISCHCV
        if column == 1:
            ischcv_prev = self.dataScheme[row]['ischcv']
            ischcv = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            self.dataScheme[row]['ischcv'] = ischcv
            self.NPE.setScheme(name, ischcv)
            self.NPE.setBlendingFactor(name, self.dataScheme[row]['blencv'])
            nvd_inc = 0
            if ischcv_prev == 'nvd_tvd':
                nvd_inc -= 1
            if ischcv == 'nvd_tvd':
                nvd_inc += 1
            elif ischcv in ('blending', 'nvd_tvd'):
                if self.dataScheme[row]['isstpc'] == 'on':
                    self.dataScheme[row]['isstpc'] = 'off'
                    self.NPE.setSlopeTest(name, self.dataScheme[row]['isstpc'])

            if nvd_inc != 0:
                self.dataScheme[row]['nvd_limiter'] = self.NPE.getNVDLimiter(name)
                self.nvd_count += nvd_inc

        # set BLENCV
        elif column == 2:
            self.dataScheme[row]['blencv'] = from_qvariant(value, float)
            self.NPE.setBlendingFactor(name, self.dataScheme[row]['blencv'])

        # set ISSTPC
        elif column == 3:
            self.dataScheme[row]['isstpc'] = self.dicoV2M_isstpc[str(from_qvariant(value, to_text_string))]
            self.NPE.setSlopeTest(name, self.dataScheme[row]['isstpc'])

        # set limiter
        elif column == 4:
            self.dataScheme[row]['nvd_limiter'] = str(from_qvariant(value, to_text_string)).lower()
            self.NPE.setNVDLimiter(name, self.dataScheme[row]['nvd_limiter'])

        # set IRCFLU
        elif column == 5:
            v = from_qvariant(value, int)
            if v == Qt.Unchecked:
                self.dataScheme[row]['ircflu'] = "off"
            else:
                self.dataScheme[row]['ircflu'] = "on"
            self.NPE.setFluxReconstruction(name, self.dataScheme[row]['ircflu'])

        # set NSWRSM
        elif column == 6:
            self.dataScheme[row]['nswrsm'] = from_qvariant(value, int)
            self.NPE.setRhsReconstruction(name, self.dataScheme[row]['nswrsm'])

        self.dataChanged.emit(index, index)
        return True

#-------------------------------------------------------------------------------
# Combo box delegate for Gradient type
#-------------------------------------------------------------------------------

class GradientTypeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, boundary=False):
        super(GradientTypeDelegate, self).__init__(parent)
        self.parent = parent
        self.boundary = boundary


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        row = index.row()
        if self.boundary:
            editor.addItem("Automatic")
        else:
            editor.addItem("Global")
        editor.addItem("Green Iter")
        editor.addItem("LSQ")
        editor.addItem("LSQ Ext")
        editor.addItem("Green LSQ")
        editor.addItem("Green LSQ Ext")
        editor.addItem("Green VTX")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"automatic": 0, "global": 0,
                "green_iter": 1, "lsq": 2, "lsq_ext": 3,
                "green_lsq": 4, "green_lsq_ext": 5,
                "green_vtx": 6}
        row = index.row()
        if self.boundary:
            string = index.model().dataScheme[row]['b_gradient_r']
        else:
            string = index.model().dataScheme[row]['c_gradient_r']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value)


#-------------------------------------------------------------------------------
# Combo box delegate for Gradient type
#-------------------------------------------------------------------------------

class GradientLimiterDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None):
        super(GradientLimiterDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        row = index.row()
        editor.addItem("Disabled")
        editor.addItem("Cell gradient")
        editor.addItem("Face gradient")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"none": 0, "cell": 1, "face": 2}
        row = index.row()
        string = index.model().dataScheme[row]['imligr']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value)


#-------------------------------------------------------------------------------
# Line edit delegate for gradient epsilon
#-------------------------------------------------------------------------------

class GradientFloatDelegate(QItemDelegate):
    def __init__(self, parent=None, max_val=None):
        super(GradientFloatDelegate, self).__init__(parent)
        self.parent = parent
        self.max_val = max_val


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        if self.max_val != None:
            validator = DoubleValidator(editor, min=0., max=self.max_val)
        else:
            validator = DoubleValidator(editor, min=0.)
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
# Gradient class
#-------------------------------------------------------------------------------

class StandardItemModelGradient(QStandardItemModel):

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
                        self.tr("Volume\nGradient"),
                        self.tr("Boundary\nReconstruction"),
                        self.tr("Fixed-point\nThreshold"),
                        self.tr("Limiter\nType"),
                        self.tr("Limiter\nFactor")]
        self.keys = ['name', 'c_gradient_r', 'b_gradient_r', 'epsrgr',
                     'imligr', 'climgr']
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                role = Qt.DisplayRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)

        self.tooltips = [
            self.tr("Equation parameter: 'imrgra'\n\n"
                    "Gradient reconstruction scheme\n\n"
                    "- Green Iter: Green-Gauss with iterative face value reconstruction\n"
                    "- LSQ: least-squares on standard (face-adjacent) neighborhood\n"
                    "- LSQ Ext: least-squares on extended neighborhood\n"
                    "- Green LSQ: Green-Gauss with LSQ gradient face values\n"
                    "- Green LSQ Ext: Green-Gauss with LSQ Ext gradient face values\n"
                    "- Green VTX: Green-Gauss with vertex interpolated face values."),
            self.tr("Equation parameter: 'b_gradient_r'\n\n"
                    "Local boundary gradient for reconstruction at boundaries;\n"
                    "Using the default least-squares reconstruction\n"
                    "allows computing the gradient only at selected cells\n"
                    "for better performance\n\n."
                    "- Green Iter: Green-Gauss with iterative face value reconstruction\n"
                    "- LSQ: least-squares on standard (face-adjacent) neighborhood\n"
                    "- LSQ Ext: least-squares on extended neighborhood\n"
                    "- Green LSQ: Green-Gauss with LSQ gradient face values\n"
                    "- Green LSQ Ext: Green-Gauss with LSQ Ext gradient face values\n"
                    "- Green VTX: Green-Gauss with vertex interpolated face values."),
            self.tr("Equation parameter: 'epsrgr'\n\n"
                    "Relative precision for the iterative gradient and\n"
                    "fixed-point Neumann BC computation for least-squares gradient."),
            self.tr("Equation parameter: 'imligr'\n\n"
                    "Gradient limiter type.\n"
                    "For the default/least squares boundary reconstruction\n",
                    "the face gradient limiter is replaced\n"
                    "by the cell gradient limiter"),
            self.tr("Equation parameter: 'cmligr'\n\n"
                    "Gradient limiter factor.")
        ]


    def populateModel(self):
        self.dicoV2M = {"Automatic": 'automatic',
                        "Global": 'global',
                        "Green Iter": 'green_iter',
                        "LSQ": 'lsq',
                        "LSQ Ext": 'lsq_ext',
                        "Green LSQ": 'green_lsq',
                        "Green LSQ Ext": 'green_lsq_ext',
                        "Green VTX": 'green_vtx'}
        self.dicoM2V = {}
        for k in self.dicoV2M:
            self.dicoM2V[self.dicoV2M[k]] = k

        self.dicoV2M_imligr = {"Disabled": 'none',
                               "Cell gradient": 'cell',
                               "Face gradient": 'face'}
        self.dicoM2V_imligr = {}
        for k in self.dicoV2M_imligr:
            self.dicoM2V_imligr[self.dicoV2M_imligr[k]] = k

        for v in self.NPE.getSchemeList():
            name = v[0]
            dico = {}
            dico['name']  = name
            dico['c_gradient_r'] = self.NPE.getCellGradientType(name)
            dico['b_gradient_r'] = self.NPE.getBoundaryGradientType(name)
            dico['epsrgr'] = self.NPE.getGradientEpsilon(name)
            dico['imligr'] = self.NPE.getGradientLimiter(name)
            dico['climgr'] = self.NPE.getGradientLimitFactor(name)
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

        if dico[key] is None:
            return None

        if role == Qt.ToolTipRole:
            col = index.column()
            if col > 0 and col < 6:
                return self.tooltips[col-1]
            elif col > 0:
                return self.tr("code_saturne keyword: " + key)

        elif role == Qt.DisplayRole:
            if column in (1, 2):
                return self.dicoM2V[dico[key]]
            elif column == 4:
                return self.dicoM2V_imligr[dico[key]]
            elif column == 5:
                if self.dataScheme[row]['imligr'] != 'none':
                    return dico[key]
                else:
                    return ""
            else:
                return dico[key]

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags

        column = index.column()
        row = index.row()
        # disabled item
        if (row, column) in self.disabledItem:
            return Qt.NoItemFlags

        if column == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif column == 5:
            if self.dataScheme[row]['imligr'] != 'none':
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.NoItemFlags
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]

        if role == Qt.ToolTipRole:
            if section == 0:
                return self.tr("variable or equation name")
            elif section < 6:
                return self.tooltips[section-1]

        return None


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        name = self.dataScheme[row]['name']

        if column == 1:
            c_gradient_r = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            self.dataScheme[row]['c_gradient_r'] = c_gradient_r
            self.NPE.setCellGradientType(name, c_gradient_r)

        elif column == 2:
            b_gradient_r = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            self.dataScheme[row]['b_gradient_r'] = b_gradient_r
            self.NPE.setBoundaryGradientType(name, b_gradient_r)

        elif column == 3:
            v = from_qvariant(value, float)
            self.dataScheme[row]['epsrgr'] = v
            self.NPE.setGradientEpsilon(name, v)

        elif column == 4:
            imligr = self.dicoV2M_imligr[str(from_qvariant(value, to_text_string))]
            self.dataScheme[row]['imligr'] = imligr
            self.NPE.setGradientLimiter(name, imligr)

        elif column == 5:
            if value != "":
                v = from_qvariant(value, float)
                self.dataScheme[row]['climgr'] = v
                self.NPE.setGradientLimitFactor(name, v)

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
                       "Multigrid, K-cycle, HPC" : 'multigrid_k_cycle_hpc',
                       "Conjugate gradient"     : 'conjugate_gradient',
                       "Flexible conjugate gradient" : 'flexible_conjugate_gradient',
                       "Inexact conjugate gradient"  : 'inexact_conjugate_gradient',
                       "Jacobi"                 : 'jacobi',
                       "BiCGstab"               : 'bi_cgstab',
                       "BiCGstab2"              : 'bi_cgstab2',
                       "GMRES"                  : 'gmres',
                       "GCR"                    : 'gcr',
                       "Automatic"              : "automatic",
                       "Gauss Seidel"           : "gauss_seidel",
                       "Symmetric Gauss Seidel" : "symmetric_gauss_seidel",
                       "conjugate residual"     : "PCR3",
                       "None"                   : "none",
                       "Polynomial"             : "polynomial"}
        self.dicoM2V= {"multigrid"              : 'Multigrid, V-cycle',
                       "multigrid_k_cycle"      : 'Multigrid, K-cycle',
                       "multigrid_k_cycle_hpc"  : 'Multigrid, K-cycle, HPC',
                       "conjugate_gradient"     : 'Conjugate gradient',
                       "inexact_conjugate_gradient"  : 'Inexact conjugate gradient',
                       "flexible_conjugate_gradient" : 'Flexible conjugate gradient',
                       "jacobi"                 : 'Jacobi',
                       "bi_cgstab"              : 'BiCGstab',
                       "bi_cgstab2"             : 'BiCGstab2',
                       'gmres'                  : "GMRES",
                       'gcr'                    : "GCR",
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
                return self.tr("Equation parameter: epsilo\n\n"
                               "Convergence threshold for linear solver")
            elif index.column() == 4:
                return self.tr("Equation parameter: verbosity")
            elif index.column() == 5:
                return self.tr("code_saturne keyword: cdtvar")

        elif role == Qt.DisplayRole:
            row = index.row()
            dico = self.dataSolver[row]

            if index.column() == 0:
                return dico['name']
            elif index.column() == 1:
                return self.dicoM2V[dico['iresol']]
            elif index.column() == 2:
                if dico['iresol'] not in ('multigrid', 'jacobi',
                                          'gauss_seidel',
                                          'symmetric_gauss_seidel'):
                    return self.dicoM2V[dico['precond']]
                else:
                    return None
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

        row = index.row()
        column = index.column()
        # disable item
        if (row, column) in self.disabledItem:
            return Qt.NoItemFlags

        if column == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif column == 2:
            if self.dataSolver[row]['iresol'] not in ('multigrid', 'jacobi',
                                                      'gauss_seidel',
                                                      'symmetric_gauss_seidel'):
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.NoItemFlags
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
            return None
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

        delegateISCHCV = SchemeOrderDelegate(self.tableViewScheme, self._tableViewLayout)
        self.tableViewScheme.setItemDelegateForColumn(1, delegateISCHCV)

        delegateBLENCV = BlendingFactorDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(2, delegateBLENCV)

        delegateISSTPC = SchemeSlopeTestDelegate(self.tableViewScheme)
        self.tableViewScheme.setItemDelegateForColumn(3, delegateISSTPC)

        delegateNVDLIM = SchemeNVDLimiterDelegate(self.tableViewScheme)
        self.tableViewScheme.setItemDelegateForColumn(4, delegateNVDLIM)

        delegateNSWRSM = RhsReconstructionDelegate(self.tableViewScheme, self.turb)
        self.tableViewScheme.setItemDelegateForColumn(6, delegateNSWRSM)

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
            self.tableViewSolver.setColumnHidden(6, True)

        delegate = SolverDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegate(delegate)

        delegateSolverChoice = SolverChoiceDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(1, delegateSolverChoice)

        delegatePrecondChoice = PreconditioningChoiceDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(2, delegatePrecondChoice)

        delegateVerbosity = VerbosityDelegate(self.tableViewSolver)
        self.tableViewSolver.setItemDelegateForColumn(4, delegateVerbosity)

        # Gradient
        self.modelGradient = StandardItemModelGradient(self.NPE)
        self.tableViewGradient.setModel(self.modelGradient)
        self.tableViewGradient.setAlternatingRowColors(True)
        self.tableViewGradient.resizeColumnToContents(0)
        self.tableViewGradient.resizeRowsToContents()
        self.tableViewGradient.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewGradient.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewGradient.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewGradient.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewGradient.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateIMRGRA = GradientTypeDelegate(self.tableViewGradient)
        self.tableViewGradient.setItemDelegateForColumn(1, delegateIMRGRA)

        delegateIMRGRB = GradientTypeDelegate(self.tableViewGradient, boundary=True)
        self.tableViewGradient.setItemDelegateForColumn(2, delegateIMRGRB)

        delegateEPSRGR = GradientFloatDelegate(self.tableViewGradient, 1.0)
        self.tableViewGradient.setItemDelegateForColumn(3, delegateEPSRGR)

        delegateIMLIGR = GradientLimiterDelegate(self.tableViewGradient)
        self.tableViewGradient.setItemDelegateForColumn(4, delegateIMLIGR)

        delegateCLIMGR = GradientFloatDelegate(self.tableViewGradient)
        self.tableViewGradient.setItemDelegateForColumn(5, delegateCLIMGR)

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

        self._tableViewLayout()


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


    def _tableViewLayout(self):
        """
        Configure QTableView column number
        """
        fm = self.tableViewScheme.fontMetrics()

        if QT_API == "PYQT4":
            self.tableViewMeshes.horizontalHeader().setResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewScheme.horizontalHeader().setResizeMode(1, QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableViewScheme.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewScheme.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)

        if self.modelScheme.nvd_count == 0 :
            self.tableViewScheme.setColumnHidden(4, True)
        else:
            self.tableViewScheme.setColumnHidden(4, False)


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
