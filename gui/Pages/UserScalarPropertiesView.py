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
- LabelDelegate
- TypeDelegate
- ValueDelegate
- StandardItemModelScalars
- UserScalarPropertiesView
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

from UserScalarPropertiesForm import Ui_UserScalarPropertiesForm
from DefineUserScalarsModel import DefineUserScalarsModel
from Base.Toolbox import GuiParam
from Base.Common import LABEL_LENGTH_MAX
from Base.QtPage import ComboModel, DoubleValidator, RegExpValidator

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("UserScalarPropertiesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the label
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_plabel = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_label = ""
        #editor.installEventFilter(self)
        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        self.old_plabel = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_plabel = str(editor.text())

            if new_plabel in model.mdl.getScalarLabelsList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getScalarLabelsList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)
    
                from VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, QVariant(QString(new_plabel)), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Combo box delegate for the type of coefficient 
#-------------------------------------------------------------------------------

class TypeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, stbar): 
        super(TypeDelegate, self).__init__(parent)
        self.parent   = parent
        self.stbar    = stbar


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 2, 1)
        self.modelCombo.addItem(self.tr("constant"), 'constant')
        self.modelCombo.addItem(self.tr("variable"), 'variable')
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_model=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("TypeDelegate value = %s"%value)
        if value == "variable":
            msg = self.tr("You must complete the user subroutine usphyv")
            self.stbar.showMessage(msg, 2000)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, QVariant(value), Qt.DisplayRole)


    def tr(self, text):
        return text 

#-------------------------------------------------------------------------------
# Line edit delegate for the value
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor, min=0.)
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):
    """
    """
    def __init__(self, parent,  mdl):
        """
        """
        QStandardItemModel.__init__(self)
        
        self.headers = [self.tr("Name"),
                        self.tr("Associated\nScalar"),
                        self.tr("Type of\ncoefficient"),    # constant/variable
                        self.tr("Type of\nvalue"),          # initial/reference value
                        self.tr("Value\n(m2/s)")]
        self.setColumnCount(len(self.headers))
        
        self._data = []
        self.mdl = mdl
        self.parent = parent


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        # Display Fortran message in tips !
        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            val = self._data[row][col]
            return QVariant(val)
        return QVariant()
    

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() in [0, 2, 4]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()
        scalar_label = self._data[row][1]

        # Label
        if col == 0:
            new_plabel = str(value.toString())
            self._data[row][col] = new_plabel
            self.mdl.setScalarDiffusivityLabel(scalar_label, new_plabel)

        # Coeff. type
        elif col == 2:
            typ = str(value.toString())
            self._data[row][col] = typ
            if typ == 'constant':
                self._data[row][3] = self.tr("Reference value")
            elif typ == 'variable':
                self._data[row][3] = self.tr("Initial value")

            self.mdl.setScalarDiffusivityChoice(scalar_label, typ)
            self._subroutineMessage()

        # Value
        elif col == 4:
            coeff, ok  = value.toDouble() 
            self._data[row][col] = coeff
            self.mdl.setScalarDiffusivityInitialValue(scalar_label, coeff)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, scalar):
        """
        Add an item in the table view
        """
        coeff_label = self.mdl.getScalarDiffusivityLabel(scalar)
        coeff_value = self.mdl.getScalarDiffusivityInitialValue(scalar)
        typC = self.mdl.getScalarDiffusivityChoice(scalar)
        if typC == "constant":
            typV = self.tr("Reference value")
        else:
            typV = self.tr("Initial value")

        line = [coeff_label, scalar, typC, typV, coeff_value]
        self._data.append(line)

        row = self.rowCount()
        self.setRowCount(row+1)
        
        self._subroutineMessage()


    def _subroutineMessage(self):
        self.parent.labelFortran.hide()
        for r in range(self.rowCount()):
            if self._data[r][2] == 'variable':
                self.parent.labelFortran.show()


    def getData(self, index):
        row = index.row()
        return self._data[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class UserScalarPropertiesView(QWidget, Ui_UserScalarPropertiesForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_UserScalarPropertiesForm.__init__(self)
        self.setupUi(self)

        self.case  = case
        self.stbar = stbar

        self.mdl = DefineUserScalarsModel(self.case)

        # Widget Layout

        self.modelScalars = StandardItemModelScalars(self, self.mdl)

        self.table.setModel(self.modelScalars)
        self.table.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegateLabel = LabelDelegate(self.table)
        delegateType  = TypeDelegate(self.table, self.stbar)
        delegateValue = ValueDelegate(self.table)

        self.table.setItemDelegateForColumn(0, delegateLabel)
        self.table.setItemDelegateForColumn(2, delegateType)
        self.table.setItemDelegateForColumn(4, delegateValue)

        self.connect(self.modelScalars, SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.dataChanged)

        # Initialize widgets

        for scalar in self.mdl.getUserScalarLabelsList():
            if not self.mdl.getScalarVariance(scalar):
                self.modelScalars.addItem(scalar)


    @pyqtSignature("const QModelIndex &, const QModelIndex &")
    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.table.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.table.resizeColumnToContents(col)


    def tr(self, text):
        return text 

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
