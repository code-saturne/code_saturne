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
- MatisseGeomView
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

from Base.Toolbox import GuiParam
from MatisseGeomForm import Ui_MatisseGeomForm
import Pages.MatisseTypeModel as MatisseType
from Pages.MatisseGeomModel import MatisseGeomModel
import Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseGeomView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegates
#-------------------------------------------------------------------------------

class LineEditDelegateInt(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = QtPage.IntegerValidator(editor, "validatorInt")
        editor.setValidator(validator)
        #editor.installEventFilter(self)
        return editor

    def setEditorData(self, lineEdit, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        lineEdit.setText(value)

    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, QVariant(value))


class LineEditDelegateFloat(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = QtPage.DoubleValidator(editor, "validatorFloat")
        editor.setValidator(validator)
        #editor.installEventFilter(self)
        return editor

    def setEditorData(self, lineEdit, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        lineEdit.setText(value)

    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, QVariant(value))


#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelGeom(QStandardItemModel):

    def __init__(self, case):
        """
        """
        QStandardItemModel.__init__(self)

        self.case = case
        self.model = MatisseGeomModel(self.case)
        self.model_mat_type = MatisseType.MatisseTypeModel(self.case)

        self.setColumnCount(4)
        self._initData()


    def _initData(self):

        # Int
        self.nechrg = 0
        self.nergrs = 0
        self.neclrg = 0
        self.nergch = 0
        self.neciel = 0
        self.nptran = 0
        self.nplgrs = 0
        self.nelgrs = 0
        self.nchest = 0
        self.netran = 0

        # Double
        self.jeuchr = 0.
        self.jeurcl = 0.
        self.jeuclr = 0.
        self.jeurch = 0.
        self.hbdtoi = 0.
        self.epchem = 0.
        self.epregi = 0.
        self.hconve = 0.
        self.rconve = 0.
        self.hchali = 0.
        self.hcheva = 0.
        self.hfttoi = 0.
        self.ptrres = 0.
        self.frdtra = 0.
        self.plgres = 0.
        self.epchel = 0.
        self.dmcont = 0.

        self.texts = {}
        self.texts['jeuchr'] = (6 , self.tr("Upstream space between chimney/register"), "m")
        self.texts['nechrg'] = (7 , self.tr("Number of cells between chimney/upstream register"))
        self.texts['jeurcl'] = (8 , self.tr("Space between upstream register/canisters"), "m")
        self.texts['nergrs'] = (9 , self.tr("Number of cells between upstream register/canisters"))
        self.texts['jeuclr'] = (10, self.tr("Space between canisters/downstream register"), "m")
        self.texts['neclrg'] = (11, self.tr("Number of cells between canisters/downstream register"))
        self.texts['jeurch'] = (12, self.tr("Downstrean space between register/chimney"), "m")
        self.texts['nergch'] = (13, self.tr("Number of cells between downstream register/chimney"))
        self.texts['hbdtoi'] = (18, self.tr("Height of roof edge"), "m")
        self.texts['neciel'] = (20, self.tr("Number of cells layers above canisters"))
        self.texts['epregi'] = (4 , self.tr("Upstream and downstream registers/doors thickness"), "m")
        self.texts['epchem'] = (5 , self.tr("chimneys' thickness"), "m")
        self.texts['hconve'] = (14, self.tr("Convergent height"), "m")
        self.texts['rconve'] = (15, self.tr("Convergent ratio"), "m")
        self.texts['hchali'] = (16, self.tr("Inlet chimney height"), "m")
        self.texts['hcheva'] = (17, self.tr("Outlet chimney height"), "m")
        self.texts['hfttoi'] = (19, self.tr("Roof ridge height"), "m")
        self.texts['ptrres'] = (21, self.tr("Transverse step of canisters network"), "m")
        self.texts['nptran'] = (22, self.tr("Number of tranverse steps"))
        self.texts['netran'] = (23, self.tr("Number of cells by tranverse step"))
        self.texts['frdtra'] = (24, self.tr("Transverse reduction factor model/real"), "m")
        self.texts['plgres'] = (25, self.tr("Longitudinal step of canisters network"), "m")
        self.texts['nplgrs'] = (26, self.tr("Number of longitudinal steps"))
        self.texts['nelgrs'] = (27, self.tr("Number of cells by longitudinal step"))
        self.texts['epchel'] = (28, self.tr("Cells height (storage area)"), "m")
        self.texts['nchest'] = (29, self.tr("Number of cells layers in the storage area"))
        self.texts['dmcont'] = (30, self.tr("Canister diameter"), "m")


        self.variables = [
            ['epregi', self.epregi, 'double'],
            ['epchem', self.epchem, 'double'],
            ['jeuchr', self.jeuchr, 'double'],
            ['nechrg', self.nechrg, 'int'],
            ['jeurcl', self.jeurcl, 'double'],
            ['nergrs', self.nergrs, 'int'],
            ['jeuclr', self.jeuclr, 'double'],
            ['neclrg', self.neclrg, 'int'],
            ['jeurch', self.jeurch, 'double'],
            ['nergch', self.nergch, 'int'],
            ['hconve', self.hconve, 'double'],
            ['rconve', self.rconve, 'double'],
            ['hchali', self.hchali, 'double'],
            ['hcheva', self.hcheva, 'double'],
            ['hbdtoi', self.hbdtoi, 'double'],
            ['hfttoi', self.hfttoi, 'double'],
            ['neciel', self.neciel, 'int'],
            ['ptrres', self.ptrres, 'double'],
            ['nptran', self.nptran, 'int'],
            ['netran', self.netran, 'int'],
            ['frdtra', self.frdtra, 'double'],
            ['plgres', self.plgres, 'double'],
            ['nplgrs', self.nplgrs, 'int'],
            ['nelgrs', self.nelgrs, 'int'],
            ['epchel', self.epchel, 'double'],
            ['nchest', self.nchest, 'int'],
            ['dmcont', self.dmcont, 'double']
            ]


        self.rows_disabled = []

        i = 0
        for variable in self.variables:
            if variable[2] == 'double':
                val = self.model.getMatisseGeomDoubleVar(variable[0])
            elif variable[2] == 'int':
                val = self.model.getMatisseGeomIntVar(variable[0])
            else :
                print(variable[2]+": unknown type")
                sys.exit(1)
            var = self.variables[i][1]
            var = val
            self.variables[i][1] = var

            t = self.model_mat_type.getMatisseType()
            if variable[0] in ('hfttoi', 'hbdtoi', 'neciel') :
                if (t == 'vault') or (t == 'djw'):
                    if not i in self.rows_disabled:
                        self.rows_disabled.append(i)
            elif variable[0] in ('jeuclr', 'neclrg', 'jeurch', 'nergch' ) :
                if t == 'emm' :
                    if not i in self.rows_disabled:
                        self.rows_disabled.append(i)
            else :
                pass
            i += 1

        self.setRowCount(len(self.variables))


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:
            row = index.row()
            var = self.variables[row][0]

            if index.column() == 0:
                num = self.texts[var][0]
                return QVariant(num)

            if index.column() == 1:
                txt = self.texts[var][1]
                return QVariant(txt)

            if index.column() == 2:
                val = self.variables[row][1]
                return QVariant(val)

            if index.column() == 3:
                if len(self.texts[var])>2:
                    unit = self.texts[var][2]
                else:
                    unit = ""
                return QVariant(unit)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.row() in self.rows_disabled:
            return Qt.ItemIsSelectable
        if index.column() in [0,1,3]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    #def setData(self, index, value, role):
    def setData(self, index, value):
        row = index.row()

        if index.column() == 2:
            var = self.variables[row][0]
            typ = self.variables[row][2]
            if typ == "int":
                val, ok = value.toInt()
            elif typ == "double":
                val, ok = value.toDouble()
            else:
                val = 0.
            # get attribute and set value ???
            attr = getattr(self, var)
            attr = val
            self.variables[row][1] = val
            self.model.setMatisseGeomVar(var, val)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseGeomView(QWidget, Ui_MatisseGeomForm):
    """
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseGeomForm.__init__(self)
        self.setupUi(self)

        self.case = case

        # Create the Page layout.

        self.modelGeom = StandardItemModelGeom(self.case)
        self.tableView.setModel(self.modelGeom)
        self.tableView.setAlternatingRowColors(True)
        self.tableView.resizeColumnsToContents()

##         Note: If a delegate has been assigned to both a row and a column,
##         the row delegate (i.e., this delegate) will take presedence and
##         manage the intersecting cell index.

        # First defines a float delegate (i.e. validator) for the entire column ...
        delegateFloat = LineEditDelegateFloat(self.tableView)
        self.tableView.setItemDelegateForColumn(2,delegateFloat)

        # ... then define an int delegate for certain rows!
        delegateInt = LineEditDelegateInt(self.tableView)
        self.tableView.setItemDelegateForRow(3,delegateInt)
        self.tableView.setItemDelegateForRow(5,delegateInt)
        self.tableView.setItemDelegateForRow(7,delegateInt)
        self.tableView.setItemDelegateForRow(9,delegateInt)
        self.tableView.setItemDelegateForRow(16,delegateInt)
        self.tableView.setItemDelegateForRow(18,delegateInt)
        self.tableView.setItemDelegateForRow(19,delegateInt)
        self.tableView.setItemDelegateForRow(22,delegateInt)
        self.tableView.setItemDelegateForRow(23,delegateInt)
        self.tableView.setItemDelegateForRow(25,delegateInt)


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
