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
This module defines the 'Interfacial forces' page.

This module contains the following classes:
- FieldDelegate
- DragDelegate
- AddedMassDelegate
- LiftDelegate
- DispersionTurbulentDelegate
- WallForceDelegate
- StandardItemModelInterfacialForces
- InterfacialForcesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.Base.QtPage import ComboModel, to_qvariant, from_qvariant, to_text_string
from InterfacialForces import Ui_InterfacialForces
from code_saturne.model.InterfacialForcesModel import InterfacialForcesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialForcesView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Combo box delegate for the field choice
#-------------------------------------------------------------------------------


class FieldDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl):
        super(FieldDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        if index.column() == 0 :
            txt = index.model().getData(index)[0]
            fieldIda = self.mdl.getFieldId(txt)
            for fieldId in self.mdl.getFieldIdaList(fieldIda) :
                label = self.mdl.getLabel(fieldId[0])
                self.modelCombo.addItem(self.tr(label), label)
        else :
            row = index.row()
            fieldIda, id = self.mdl.getForceList()[row]
            label = self.mdl.getLabel(id)
            self.modelCombo.addItem(self.tr(label), label)
            for fieldId in self.mdl.getFieldIdbList(fieldIda) :
                label = self.mdl.getLabel(fieldId)
                self.modelCombo.addItem(self.tr(label), label)

        editor.setMinimumSize(editor.sizeHint())
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
        log.debug("FieldDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the drag model
#-------------------------------------------------------------------------------


class DragDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(DragDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fielda = index.model().getData(index)[0]
        fieldb = index.model().getData(index)[1]

        fieldaId = self.mdl.getFieldId(fielda)
        fieldbId = self.mdl.getFieldId(fieldb)

        for model in self.mdl.getAvailableDragModels(fieldaId, fieldbId) :
            self.modelCombo.addItem(self.tr(self.dicoM2V[model]), model)

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("DragDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the added mass model
#-------------------------------------------------------------------------------


class AddedMassDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(AddedMassDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        for model in self.mdl.getAvailableAddedMassModels() :
            self.modelCombo.addItem(self.tr(self.dicoM2V[model]), model)

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("AddedMassDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the lift
#-------------------------------------------------------------------------------


class LiftDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(LiftDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        for model in self.mdl.getAvailableLiftModels() :
            self.modelCombo.addItem(self.tr(self.dicoM2V[model]), model)

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("LiftDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the turbulent dispersion
#-------------------------------------------------------------------------------


class DispersionTurbulentDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(DispersionTurbulentDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fielda = index.model().getData(index)[0]
        fieldb = index.model().getData(index)[1]

        fieldaId = self.mdl.getFieldId(fielda)
        fieldbId = self.mdl.getFieldId(fieldb)

        for model in self.mdl.getAvailableTurbulenteDispersionModelList(fieldaId, fieldbId) :
            self.modelCombo.addItem(self.tr(self.dicoM2V[model]), model)

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("DispersionTurbulentDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Combo box delegate for the wall forces
#-------------------------------------------------------------------------------


class WallForceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(WallForceDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fielda = index.model().getData(index)[0]
        fieldb = index.model().getData(index)[1]

        fieldaId = self.mdl.getFieldId(fielda)
        fieldbId = self.mdl.getFieldId(fieldb)

        for model in self.mdl.getAvailableWallForcesModelList(fieldaId, fieldbId) :
            self.modelCombo.addItem(self.tr(self.dicoM2V[model]), model)

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("WallForceDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(self.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# StandardItemModelInterfacialForces class
#-------------------------------------------------------------------------------

class StandardItemModelInterfacialForces(QStandardItemModel):

    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Field A name"),
                         self.tr("Field B name"),
                         self.tr("Drag"),
                         self.tr("Added mass"),
                         self.tr("Lift"),
                         self.tr("Turbulent\ndispersion"),
                         self.tr("Wall forces")]

        self.setColumnCount(len(self.headers))

        self.tooltip = []
        self._data   = []
        self.mdl     = mdl
        self.parent  = parent
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(data)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        # field A choice
        if col == 0:
            oldField = self._data[row][col]
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            oldFieldId = self.mdl.getFieldId(oldField)
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            fieldaId = self.mdl.getFieldId(new_pmodel)

            self.mdl.setFielda(oldFieldId, fieldbId, fieldaId)
            self.updateItem(index)

        # field B choice
        elif col == 1:
            oldField = self._data[row][col]
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            oldFieldId = self.mdl.getFieldId(oldField)
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(new_pmodel)

            self.mdl.setFieldb(fieldaId, oldFieldId, fieldbId)
            self.updateItem(index)

        # drag
        elif col == 2:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            self.mdl.setDragModel(fieldaId, fieldbId, self.dicoV2M[new_pmodel])

        # added mass
        elif col == 3:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            self.mdl.setAddMassModel(fieldaId, fieldbId, self.dicoV2M[new_pmodel])

        # lift
        elif col == 4:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            self.mdl.setLiftModel(fieldaId, fieldbId, self.dicoV2M[new_pmodel])

        # turbulent dispersion
        elif col == 5:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            self.mdl.setTurbDispModel(fieldaId, fieldbId, self.dicoV2M[new_pmodel])

        # wall forces
        elif col == 6:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(self._data[row][1])
            self.mdl.setWallForceModel(fieldaId, fieldbId, self.dicoV2M[new_pmodel])

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, force=None):
        """
        Add/load a force.
        """
        if (len(self.mdl.getFreeCouples()) > 0 or force != None ) :
            row = self.rowCount()

            if force == None :
                couple = self.mdl.addForce()
            else :
                couple = force

            labela = self.mdl.getLabel(couple[0])
            labelb = self.mdl.getLabel(couple[1])

            drag      = self.dicoM2V[self.mdl.getDragModel(couple[0], couple[1])]
            addmass   = self.dicoM2V[self.mdl.getAddMassModel(couple[0], couple[1])]
            lift      = self.dicoM2V[self.mdl.getLiftModel(couple[0], couple[1])]
            dispturb  = self.dicoM2V[self.mdl.getTurbDispModel(couple[0], couple[1])]
            wallforce = self.dicoM2V[self.mdl.getWallForceModel(couple[0], couple[1])]

            force = [labela, labelb, drag, addmass, lift, dispturb, wallforce]

            self._data.append(force)
            self.setRowCount(row+1)
        else :
            title = self.tr("Interfacial momentum transfer")
            msg   = self.tr("No fields couple to create a new force")
            QMessageBox.information(self.parent, title, msg)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        fieldaId = self.mdl.getFieldId(self._data[row][0])
        fieldbId = self.mdl.getFieldId(self._data[row][1])
        self.mdl.deleteForce(fieldaId, fieldbId)
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def updateItem(self, index):
        """
        update item
        """
        row = index.row()
        (fieldaId, fieldbId) = self.mdl.getForceList()[row]
        self._data[row][1] = self.mdl.getLabel(fieldbId)
        self._data[row][2] = self.dicoM2V[self.mdl.getDragModel(fieldaId, fieldbId)]
        self._data[row][3] = self.dicoM2V[self.mdl.getAddMassModel(fieldaId, fieldbId)]
        self._data[row][4] = self.dicoM2V[self.mdl.getLiftModel(fieldaId, fieldbId)]
        self._data[row][5] = self.dicoM2V[self.mdl.getTurbDispModel(fieldaId, fieldbId)]
        self._data[row][6] = self.dicoM2V[self.mdl.getWallForceModel(fieldaId, fieldbId)]


#-------------------------------------------------------------------------------
#  class InterfacialForces
#-------------------------------------------------------------------------------

class InterfacialForcesView(QWidget, Ui_InterfacialForces):
    """
    Main fields layout.
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InterfacialForces.__init__(self)
        self.setupUi(self)

        self.case = case
        self.browser = tree
        self.case.undoStopGlobal()
        self.mdl = InterfacialForcesModel(self.case)

        # Dico
        self.dicoM2V= {"none"                  : 'None',
                       "LLB_model"             : 'LLB model',
                       "GTD_model"             : 'GTD model',
                       "antal"                 : 'Antal',
                       "tomiyama"              : 'Tomiyama',
                       "ishii"                 : 'Ishii',
                       "inclusions"            : 'Gobin et al.',
                       "Wen_Yu"                : 'Wen and Yu',
                       "standard"              : 'Standard',
                       "zuber"                 : 'Zuber',
                       "coef_cst"              : 'Constant coefficient',
                       "Tomiyama_SMD"          : 'Tomiyama SMD'}

        self.dicoV2M= {"None"                   : 'none',
                       "LLB model"              : 'LLB_model',
                       "GTD model"              : 'GTD_model',
                       "Antal"                  : 'antal',
                       "Tomiyama"               : 'tomiyama',
                       "Ishii"                  : 'ishii',
                       "Gobin et al."           : 'inclusions',
                       "Wen and Yu"             : 'Wen_Yu',
                       "Standard"               : 'standard',
                       "Zuber"                  : 'zuber',
                       "Constant coefficient"   : 'coef_cst',
                       "Tomiyama SMD"           : 'Tomiyama_SMD'}

        self.tableModelInterfacialForces = StandardItemModelInterfacialForces(self, self.mdl, self.dicoM2V, self.dicoV2M)
        self.tableViewInterfacialForces.setModel(self.tableModelInterfacialForces)
        self.tableViewInterfacialForces.resizeColumnsToContents()
        self.tableViewInterfacialForces.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewInterfacialForces.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInterfacialForces.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInterfacialForces.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewInterfacialForces.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInterfacialForces.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewInterfacialForces.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewInterfacialForces.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewInterfacialForces.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateFielda    = FieldDelegate(self.tableViewInterfacialForces, self.mdl)
        delegateFieldb    = FieldDelegate(self.tableViewInterfacialForces, self.mdl)
        delegateDrag      = DragDelegate(self.tableViewInterfacialForces, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateAddedMass = AddedMassDelegate(self.tableViewInterfacialForces, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateLift      = LiftDelegate(self.tableViewInterfacialForces, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateDispTurb  = DispersionTurbulentDelegate(self.tableViewInterfacialForces, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateWallForce = WallForceDelegate(self.tableViewInterfacialForces, self.mdl, self.dicoM2V, self.dicoV2M)

        self.tableViewInterfacialForces.setItemDelegateForColumn(0, delegateFielda)
        self.tableViewInterfacialForces.setItemDelegateForColumn(1, delegateFieldb)
        self.tableViewInterfacialForces.setItemDelegateForColumn(2, delegateDrag)
        self.tableViewInterfacialForces.setItemDelegateForColumn(3, delegateAddedMass)
        self.tableViewInterfacialForces.setItemDelegateForColumn(4, delegateLift)
        self.tableViewInterfacialForces.setItemDelegateForColumn(5, delegateDispTurb)
        self.tableViewInterfacialForces.setItemDelegateForColumn(6, delegateWallForce)

        # Combo models
        self.modelContinuousMomentumTransfer = ComboModel(self.comboBoxContinuousMomentumTransfer, 3, 1)
        self.modelContinuousMomentumTransfer.addItem(self.tr('none'), 'none')
        self.modelContinuousMomentumTransfer.addItem(self.tr('Large Interface Model'), 'Large_Interface_Model')
        self.modelContinuousMomentumTransfer.addItem(self.tr('Large Bubble Model'), 'Large_Bubble_Model')

        self.modelInterfacialMethod = ComboModel(self.comboBoxInterfacialMethod, 2, 1)
        self.modelInterfacialMethod.addItem(self.tr('Refined gradient method'), 'refined_gradient')
        self.modelInterfacialMethod.addItem(self.tr('Refined gradient method with 2 criteria'), 'refined_gradient_2_criteria')

        # hide/show groupBoxContinuousMomentumTransfer
        if len(self.mdl.getContinuousFieldList()) >=2:
            self.groupBoxContinuousMomentumTransfer.show()
            model = self.mdl.getContinuousCouplingModel()
            self.modelContinuousMomentumTransfer.setItem(str_model=model)
            self.updateContinuousMomemtum()
        else :
            self.groupBoxContinuousMomentumTransfer.hide()

        if len(self.mdl.getDispersedFieldList()) > 0:
            self.groupBoxDispersedMomentumTransfer.show()
        else :
            self.groupBoxDispersedMomentumTransfer.hide()

        # Connect signals to slots
        self.pushButtonAdd.clicked.connect(self.slotAddForce)
        self.pushButtonDelete.clicked.connect(self.slotDeleteForce)
        self.comboBoxContinuousMomentumTransfer.activated[str].connect(self.slotContinuousMomentumTransfer)
        self.comboBoxInterfacialMethod.activated[str].connect(self.slotInterfacialMethod)
        self.checkBoxGradPCorrection.clicked.connect(self.slotGradPCorrection)
        self.tableModelInterfacialForces.dataChanged.connect(self.dataChanged)
        self.checkBoxBubblesForLIM.clicked.connect(self.slotBubblesForLIM)
        self.checkBoxInterfaceSharpening.clicked.connect(self.slotInterfaceSharpening)
        self.checkBoxUnsharpenedCells.clicked.connect(self.slotUnsharpenedCells)
        self.checkBoxSurfaceTension.clicked.connect(self.slotSurfaceTension)

        for force in self.mdl.getForceList() :
            self.tableModelInterfacialForces.newItem(force)

        self.case.undoStartGlobal()


    def updateContinuousMomemtum(self):
        """
        update view for continuous momentum choice
        """
        model = self.mdl.getContinuousCouplingModel()

        # hide/show groupBoxSeparatePhases qui si LIM ou sep.phases
        if model != "none" :
            self.groupBoxInterfaceSharpening.show()
            isInterfaceSharpening = self.mdl.getInterfaceSharpeningStatus() == "on"
            self.checkBoxInterfaceSharpening.setChecked(isInterfaceSharpening)

            if isInterfaceSharpening:
                self.groupBoxInterfaceSharpeningOptions.show()
                isUnsharpenedCells = self.mdl.getUnsharpenedCellsStatus() == "on"
                self.checkBoxUnsharpenedCells.setChecked(isUnsharpenedCells)
                isSurfaceTension = self.mdl.getSurfaceTensionStatus() == "on"
                self.checkBoxSurfaceTension.setChecked(isSurfaceTension)

            else:
                self.groupBoxInterfaceSharpeningOptions.hide()

            if model == "Large_Interface_Model" :
                self.groupBoxBubblesForLIM.show()
                isBubblesForLIM = self.mdl.getBubblesForLIMStatus() == "on"
                self.checkBoxBubblesForLIM.setChecked(isBubblesForLIM)
                self.groupBoxSeparatePhases.hide()
            elif model == "Large_Bubble_Model" :
                self.groupBoxSeparatePhases.hide()
                self.groupBoxBubblesForLIM.hide()
            else :
                self.groupBoxBubblesForLIM.hide()
                self.groupBoxSeparatePhases.show()
                isGradCorrection = self.mdl.getGradPCorrectionStatus() == "on"
                self.checkBoxGradPCorrection.setChecked(isGradCorrection)
                if isGradCorrection :
                    self.comboBoxInterfacialMethod.setEnabled(1)
                    model = self.mdl.getGradPCorrectionModel()
                    self.modelInterfacialMethod.setItem(str_model=model)
                else :
                    self.comboBoxInterfacialMethod.setEnabled(0)
        else :
            self.groupBoxSeparatePhases.hide()
            self.groupBoxBubblesForLIM.hide()
            self.groupBoxInterfaceSharpening.hide()


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewInterfacialForces.resizeColumnsToContents()
        self.tableViewInterfacialForces.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewInterfacialForces.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewInterfacialForces.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)


    @pyqtSlot()
    def slotAddForce(self):
        """
        Add a Force.
        """
        self.tableViewInterfacialForces.clearSelection()
        self.tableModelInterfacialForces.newItem()


    @pyqtSlot()
    def slotDeleteForce(self):
        """
        Suppress a Force.
        """
        row = self.tableViewInterfacialForces.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteProfile -> %s" % row)
            self.tableModelInterfacialForces.deleteItem(row)


    @pyqtSlot(str)
    def slotContinuousMomentumTransfer(self, text):
        """
        configure momentum transfer for continuous phases
        """
        value = self.modelContinuousMomentumTransfer.dicoV2M[text]
        log.debug("slotContinuousMomentumTransfer -> %s" % value)
        self.mdl.setContinuousCouplingModel(value)

        self.updateContinuousMomemtum()


    @pyqtSlot(str)
    def slotInterfacialMethod(self, text):
        """
        configure gradP correction model for continuous phases
        """
        value = self.modelInterfacialMethod.dicoV2M[text]
        log.debug("slotInterfacialMethod -> %s" % value)
        self.mdl.setGradPCorrectionModel(value)


    @pyqtSlot(bool)
    def slotGradPCorrection(self, checked):
        """
        check box for gradP correction
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setGradPCorrectionStatus(status)

        if status == 'on' :
            self.comboBoxInterfacialMethod.setEnabled(1)
            model = self.mdl.getGradPCorrectionModel()
            self.modelInterfacialMethod.setItem(str_model=model)
        else :
            self.comboBoxInterfacialMethod.setEnabled(0)


    @pyqtSlot(bool)
    def slotBubblesForLIM(self, checked):
        """
        check box for bubbles forces for LIM
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setBubblesForLIMStatus(status)
        self.updateContinuousMomemtum()
        self.browser.configureTree(self.case)


    @pyqtSlot(bool)
    def slotInterfaceSharpening(self, checked):
        """
        check box for interface sharpening
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setInterfaceSharpeningStatus(status)
        self.updateContinuousMomemtum()


    @pyqtSlot(bool)
    def slotUnsharpenedCells(self, checked):
        """
        check box for allow unsharpened cells
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setUnsharpenedCellsStatus(status)


    @pyqtSlot(bool)
    def slotSurfaceTension(self, checked):
        """
        check box for activate surface tension
        """
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setSurfaceTensionStatus(status)


