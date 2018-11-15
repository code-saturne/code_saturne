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
This module defines the 'Interfacial enthalpy transfer' page.

This module contains the following classes:
- FieldDelegate
- StandardItemModelInterfacialEnthalpy
- InterfacialEnthalpyView
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from InterfacialEnthalpy import Ui_InterfacialEnthalpy
from InterfacialEnthalpyModel import InterfacialEnthalpyModel
from NonCondensableModel import NonCondensableModel
from InterfacialForcesModel import InterfacialForcesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InterfacialEnthalpyView")
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
            fieldIda, id = self.mdl.getEnthalpyCoupleList()[row]
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
# StandardItemModelInterfacialEnthalpy class
#-------------------------------------------------------------------------------

class StandardItemModelInterfacialEnthalpy(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Field A name"),
                         self.tr("Field B name ")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data = []
        self.mdl = mdl


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

        # field B choice
        elif col == 1:
            oldField = self._data[row][col]
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            oldFieldId = self.mdl.getFieldId(oldField)
            fieldaId = self.mdl.getFieldId(self._data[row][0])
            fieldbId = self.mdl.getFieldId(new_pmodel)

            self.mdl.setFieldb(fieldaId, oldFieldId, fieldbId)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, couplefield = None):
        """
        Add/load a couple a field for model
        """
        # TODO on bride pour l'instant a un seul couple
        # on doit controler qu'on est en eau/vap sur le couple 1/2
        #if (len(self.mdl.getFreeCouples()) > 0 or couplefield != None ) :
        if (len(self.mdl.getFreeCouples()) > 0 or couplefield != None ) and self.rowCount() == 0 :
            row = self.rowCount()

            if couplefield == None :
                couple = self.mdl.addEnthalpyCouple()
            else :
                couple = couplefield

            labela = self.mdl.getLabel(couple[0])
            labelb = self.mdl.getLabel(couple[1])

            field = [labela, labelb]

            self._data.append(field)
            self.setRowCount(row+1)
        else :
            title = self.tr("Interfacial enthalpy transfer")
            msg   = self.tr("No fields couple to create a new enthalpy transfer")
            QMessageBox.information(self.parent, title, msg)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        fieldaId = self.mdl.getFieldId(self._data[row][0])
        fieldbId = self.mdl.getFieldId(self._data[row][1])
        self.mdl.deleteEnthalpyCouple(fieldaId, fieldbId)
        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def getCouple(self, row) :
        """
        return list of field Id
        """
        fieldaId = self.mdl.getFieldId(self._data[row][0])
        fieldbId = self.mdl.getFieldId(self._data[row][1])
        return [fieldaId, fieldbId]


#-------------------------------------------------------------------------------
# InterfacialEnthalpyView class
#-------------------------------------------------------------------------------

class InterfacialEnthalpyView(QWidget, Ui_InterfacialEnthalpy):
    """
    InterfacialEnthalpyView layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InterfacialEnthalpy.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = InterfacialEnthalpyModel(self.case)

        if (len(self.mdl.getSolidFieldIdList()) > 0 or len(self.mdl.getFieldIdList()) > 2 or self.mdl.getPredefinedFlow == "None") :
            self.groupBoxSolidEnergyTransfer.show()
            # Combo box models
            self.modelSolidEnergyTransfer = ComboModel(self.comboBoxSolidEnergyTransfer, 2, 1)
            self.modelSolidEnergyTransfer.addItem(self.tr("none"), "none")
            self.modelSolidEnergyTransfer.addItem(self.tr("gas-particle"), "gas_particule")

            model = self.mdl.getSolidEnergyTransfer()
            self.modelSolidEnergyTransfer.setItem(str_model = model)
        else :
            self.groupBoxSolidEnergyTransfer.hide()

        self.tableModelLiquidGasEnergyTransfer = StandardItemModelInterfacialEnthalpy(self, self.mdl)
        self.tableViewLiquidGasEnergyTransfer.setModel(self.tableModelLiquidGasEnergyTransfer)
        self.tableViewLiquidGasEnergyTransfer.resizeColumnsToContents()
        self.tableViewLiquidGasEnergyTransfer.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)
        self.tableViewLiquidGasEnergyTransfer.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewLiquidGasEnergyTransfer.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateFielda    = FieldDelegate(self.tableViewLiquidGasEnergyTransfer, self.mdl)
        delegateFieldb    = FieldDelegate(self.tableViewLiquidGasEnergyTransfer, self.mdl)

        self.tableViewLiquidGasEnergyTransfer.setItemDelegateForColumn(0, delegateFielda)
        self.tableViewLiquidGasEnergyTransfer.setItemDelegateForColumn(1, delegateFieldb)

        # Combo models
        self.modelPonderationCoefFielda = ComboModel(self.comboBoxPonderationCoefFielda, 3, 1)
        self.modelPonderationCoefFielda.addItem(self.tr("alp1"), "alp1")
        self.modelPonderationCoefFielda.addItem(self.tr("alp2"), "alp2")
        self.modelPonderationCoefFielda.addItem(self.tr("alp1*alp2"), "alp1_alp2")

        self.modelPonderationCoefFieldb = ComboModel(self.comboBoxPonderationCoefFieldb, 3, 1)
        self.modelPonderationCoefFieldb.addItem(self.tr("alp1"), "alp1")
        self.modelPonderationCoefFieldb.addItem(self.tr("alp2"), "alp2")
        self.modelPonderationCoefFieldb.addItem(self.tr("alp1*alp2"), "alp1_alp2")

        # hide/show groupBoxLiquidGasEnergyTransfer
        if len(self.mdl.getFreeCouples()) < 1 and len(self.mdl.getEnthalpyCoupleList()) == 0 :
            self.groupBoxLiquidGasEnergyTransfer.hide()
        else :
            self.groupBoxLiquidGasEnergyTransfer.show()

            if len(NonCondensableModel(self.case).getNonCondensableLabelList()) > 0 \
            and InterfacialForcesModel(self.case).getContinuousCouplingModel() \
                == 'Large_Interface_Model':
                self.checkBoxActivatePool.show()
            else:
                self.checkBoxActivatePool.hide()

        # Validators
        validatorRelaxa = DoubleValidator(self.lineEditRelaxationTimeFielda, min = 0.0)
        validatorRelaxb = DoubleValidator(self.lineEditRelaxationTimeFieldb, min = 0.0)
        self.lineEditRelaxationTimeFielda.setValidator(validatorRelaxa)
        self.lineEditRelaxationTimeFieldb.setValidator(validatorRelaxb)

        # Connect signals to slots
        self.pushButtonAdd.clicked.connect(self.slotAddEnthalpy)
        self.pushButtonDelete.clicked.connect(self.slotDeleteEnthalpy)
        self.tableModelLiquidGasEnergyTransfer.dataChanged.connect(self.dataChanged)
        self.tableViewLiquidGasEnergyTransfer.clicked.connect(self.__slotSelectField)
        self.comboBoxFieldaModel.activated[str].connect(self.slotFieldaModel)
        self.comboBoxPonderationCoefFielda.activated[str].connect(self.slotPonderationCoefFielda)
        self.comboBoxFieldbModel.activated[str].connect(self.slotFieldbModel)
        self.comboBoxPonderationCoefFieldb.activated[str].connect(self.slotPonderationCoefFieldb)
        self.lineEditRelaxationTimeFielda.textChanged[str].connect(self.slotRelaxationTimeFielda)
        self.lineEditRelaxationTimeFieldb.textChanged[str].connect(self.slotRelaxationTimeFieldb)
        self.comboBoxSolidEnergyTransfer.activated[str].connect(self.slotSolidEnergyTransfer)
        self.checkBoxActivatePool.stateChanged.connect(self.slotPoolBoilingModel)

        for couple in self.mdl.getEnthalpyCoupleList() :
            self.tableModelLiquidGasEnergyTransfer.newItem(couple)

        self.groupBoxLiquidVaporModel.hide()

        # Initial state of Pool boiling model
        if self.mdl.getPoolBoiling() == 'on':
            self.checkBoxActivatePool.setChecked(True)

        self.case.undoStartGlobal()


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewLiquidGasEnergyTransfer.resizeColumnsToContents()
        self.tableViewLiquidGasEnergyTransfer.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewLiquidGasEnergyTransfer.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        self.__updateModel(row)


    @pyqtSlot()
    def slotAddEnthalpy(self):
        """
        Add an enthalpy interfacial force
        """
        self.tableViewLiquidGasEnergyTransfer.clearSelection()
        self.groupBoxLiquidVaporModel.hide()
        self.tableModelLiquidGasEnergyTransfer.newItem()


    @pyqtSlot()
    def slotDeleteEnthalpy(self):
        """
        Suppress an enthalpy interfacial force
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteLiquidGasEnergyTransfer -> %s" % row)
            self.tableModelLiquidGasEnergyTransfer.deleteItem(row)
        self.groupBoxLiquidVaporModel.hide()


    def __slotSelectField(self, index):
        """
        Select a Field in the QTable
        """
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(index.row())
        self.groupBoxLiquidVaporModel.show()
        self.__updateModel(index.row())


    def __updateModel(self, row) :
        """
        update if necessary
        """
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)

        self.modelFieldaModel = ComboModel(self.comboBoxFieldaModel, 1, 1)
        self.modelFieldbModel = ComboModel(self.comboBoxFieldbModel, 1, 1)

        for field in (fieldIda, fieldIdb) :
            if field == fieldIda :
                self.__currentmodel = self.modelFieldaModel
            else :
                self.__currentmodel = self.modelFieldbModel

            self.__currentmodel.addItem(self.tr("Relaxation time : alphk.rok.cpk.(Ts-Tk)/tauk"),"relaxation_time")
            self.__currentmodel.addItem(self.tr("No source term"),"no_source_term")

        # fieldIda == continuous by construction and nature(fieldIda) != nature(fieldIdb)
        if self.mdl.getFieldNature(fieldIda) == "liquid" :
            if self.mdl.getCriterion(fieldIdb) == "continuous" :
                self.modelFieldaModel.addItem(self.tr("Coste-Lavieville NURETH13 model, Wall Law Type Model"),"wall_law_type_model")
                self.modelFieldbModel.addItem(self.tr("Interfacial Sublayer Model for LI3C"),"sublayer_LI3C")
            else :
                self.modelFieldaModel.addItem(self.tr("Bulk model (Ranz-Marshall)"),"bulk")
                self.modelFieldaModel.addItem(self.tr("Flashing (Cathare)"),"flashing")
                self.modelFieldaModel.addItem(self.tr("Bubble model for liquid (Manon-Berne)"),"bubble_model_for_liquid")
                #suppression temporaire le temps que le modele soit au point
                #self.modelFieldbModel.addItem(self.tr("Bubble model for vapour"),"bubble_model_for_vapour")
                self.modelFieldbModel.addItem(self.tr("Relaxation time + subcooled gas treatment"),"relaxation_time_subcooled")
        else :
            if self.mdl.getCriterion(fieldIdb) == "continuous" :
                self.modelFieldaModel.addItem(self.tr("Interfacial Sublayer Model for LI3C"),"sublayer_LI3C")
                self.modelFieldbModel.addItem(self.tr("Coste-Lavieville NURETH13 model, Wall Law Type Model"),"wall_law_type_model")
            else :
                self.modelFieldaModel.addItem(self.tr("Bulk model (Ranz-Marshall)"),"bulk")
                if len(NonCondensableModel(self.case).getNonCondensableLabelList()) > 0:
                    self.modelFieldaModel.addItem(self.tr("Droplet model for vapour"),"droplet_model_for_vapour")
                    self.modelFieldbModel.addItem(self.tr("Droplet model for liquid"),"droplet_model_for_liquid")


        model = self.mdl.getFieldModel(fieldIda, fieldIdb, fieldIda)
        self.modelFieldaModel.setItem(str_model = model)

        if model == 'relaxation_time' :
            model = self.mdl.getPonderationCoef(fieldIda, fieldIdb, fieldIda)
            self.modelPonderationCoefFielda.setItem(str_model = model)
            value = self.mdl.getRelaxationTime(fieldIda, fieldIdb, fieldIda)
            self.lineEditRelaxationTimeFielda.setText(str(value))

            self.comboBoxPonderationCoefFielda.show()
            self.labelPonderationCoefFielda.show()
            self.lineEditRelaxationTimeFielda.show()
            self.labelRelaxationTimeFielda.show()
        else :
            self.comboBoxPonderationCoefFielda.hide()
            self.labelPonderationCoefFielda.hide()
            self.lineEditRelaxationTimeFielda.hide()
            self.labelRelaxationTimeFielda.hide()

        model = self.mdl.getFieldModel(fieldIda, fieldIdb, fieldIdb)
        self.modelFieldbModel.setItem(str_model = model)

        if model == 'relaxation_time' :
            model = self.mdl.getPonderationCoef(fieldIda, fieldIdb, fieldIdb)
            self.modelPonderationCoefFieldb.setItem(str_model = model)
            value = self.mdl.getRelaxationTime(fieldIda, fieldIdb, fieldIdb)
            self.lineEditRelaxationTimeFieldb.setText(str(value))

            self.comboBoxPonderationCoefFieldb.show()
            self.labelPonderationCoefFieldb.show()
            self.lineEditRelaxationTimeFieldb.show()
            self.labelRelaxationTimeFieldb.show()
        else :
            self.comboBoxPonderationCoefFieldb.hide()
            self.labelPonderationCoefFieldb.hide()
            self.lineEditRelaxationTimeFieldb.hide()
            self.labelRelaxationTimeFieldb.hide()


    @pyqtSlot(str)
    def slotFieldaModel(self, text):
        """
        set model for field a
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        choice = self.modelFieldaModel.dicoV2M[text]
        self.mdl.setFieldModel(fieldIda, fieldIdb, fieldIda, choice)
        self.__updateModel(row)


    @pyqtSlot(str)
    def slotPonderationCoefFielda(self, text):
        """
        set ponderation coefficient for field a
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        choice = self.modelPonderationCoefFielda.dicoV2M[text]
        self.mdl.setPonderationCoef(fieldIda, fieldIdb, fieldIda, choice)


    @pyqtSlot(str)
    def slotFieldbModel(self, text):
        """
        set model for field b
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        choice = self.modelFieldbModel.dicoV2M[text]
        self.mdl.setFieldModel(fieldIda, fieldIdb, fieldIdb, choice)
        self.__updateModel(row)


    @pyqtSlot(str)
    def slotPonderationCoefFieldb(self, text):
        """
        set ponderation coefficient for field b
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        choice = self.modelPonderationCoefFieldb.dicoV2M[text]
        self.mdl.setPonderationCoef(fieldIda, fieldIdb, fieldIdb, choice)


    @pyqtSlot(str)
    def slotRelaxationTimeFielda(self, text):
        """
        Update the relaxation time for field a
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        if self.lineEditRelaxationTimeFielda.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setRelaxationTime(fieldIda, fieldIdb, fieldIda, value)


    @pyqtSlot(str)
    def slotRelaxationTimeFieldb(self, text):
        """
        Update the relaxation time for field b
        """
        row = self.tableViewLiquidGasEnergyTransfer.currentIndex().row()
        (fieldIda, fieldIdb) = self.tableModelLiquidGasEnergyTransfer.getCouple(row)
        if self.lineEditRelaxationTimeFieldb.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setRelaxationTime(fieldIda, fieldIdb, fieldIdb, value)


    @pyqtSlot(str)
    def slotSolidEnergyTransfer(self, text):
        """
        set model for solid enthalpy transfer
        """
        choice = self.modelSolidEnergyTransfer.dicoV2M[text]
        self.mdl.setSolidEnergyTransfer(choice)


    @pyqtSlot()
    def slotPoolBoilingModel(self):
        """
        Activate or deactivate the pool boiling model
        """
        self.mdl.setPoolBoiling(state = self.checkBoxActivatePool.isChecked() )

