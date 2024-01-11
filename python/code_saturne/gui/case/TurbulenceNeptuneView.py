# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the 'Turbulence' page.

This module contains the following classes:
- TurbulenceDelegate
- CouplingDelegate
- StandardItemModelTurbulence
- TurbulenceView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from TurbulenceNeptune import Ui_Turbulence
from code_saturne.model.TurbulenceNeptuneModel import TurbulenceModel, TurbulenceModelsDescription

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TurbulenceView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Combo box delegate for the turbulence
#-------------------------------------------------------------------------------


class TurbulenceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(TurbulenceDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        field = self.mdl.mainFieldsModel.list_of_fields[index.row()]
        fieldId = field.f_id

        #TODO : move this to TurbulenceNeptuneModel
        if field.flow_type == "continuous":
            turbulence_models = TurbulenceModelsDescription.continuousTurbulenceModels
        else:
            carrier = field.carrier_id
            if self.mdl.mainFieldsModel.getPredefinedFlow() == "boiling_flow":
                turbulence_models = TurbulenceModelsDescription.bubblyFlowsTurbulenceModels
            elif self.mdl.mainFieldsModel.getPredefinedFlow() == "droplet_flow":
                turbulence_models = TurbulenceModelsDescription.dropletFlowsTurbulenceModels
            elif field.phase == "solid" or self.mdl.getTurbulenceModel(carrier) != "none":
                turbulence_models = TurbulenceModelsDescription.dispersedTurbulenceModels
            else:
                turbulence_models = ["none"]
        for turb in turbulence_models:
            self.modelCombo.addItem(self.tr(self.dicoM2V[turb]), turb)
        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("TurbulenceDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the two way coupling
#-------------------------------------------------------------------------------


class CouplingDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(CouplingDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        field = self.mdl.mainFieldsModel.list_of_fields[index.row()]
        fieldId = field.f_id

        if field.flow_type == "continuous" :
               self.modelCombo.addItem(self.tr(self.dicoM2V["none"]), "none")
               self.modelCombo.disableItem(str_model="none")
        else :
               self.modelCombo.addItem(self.tr(self.dicoM2V["none"]), "none")
               carrier = field.carrier_id

               reverseCoupling = True
               if carrier == "all" :
                  # Dispersed phase carried by all the continuous phases
                  # check if the turbulence of all the continuous phases is well defined :
                  for continuous_field in self.mdl.mainFieldsModel.getContinuousFieldList():
                      carrier = continuous_field.f_id
                      _m = self.mdl.getTurbulenceModel(carrier)
                      reverseCoupling = _m in TurbulenceModelsDescription.reverseCouplingModels
               else :
                  # Dispersed phase carried by one of the continuous phases
                  _m = self.mdl.getTurbulenceModel(carrier)
                  reverseCoupling = _m in TurbulenceModelsDescription.reverseCouplingModels

               if reverseCoupling :
                  if field.phase == "gas" :
                      # bulles
                      self.modelCombo.addItem(self.tr(self.dicoM2V["large_inclusions"]), "large_inclusions")
                  else :
                      # gouttes et solide
                      self.modelCombo.addItem(self.tr(self.dicoM2V["small_inclusions"]), "small_inclusions")

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("CouplingDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the thermal turbulent fluxes
#-------------------------------------------------------------------------------

class TurbFluxDelegate(QItemDelegate):
    """
    Use of a combobox in the table.
    """

    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(TurbFluxDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        fieldId = self.mdl.mainFieldsModel.list_of_fields[index.row()].f_id

        if self.mdl.mainFieldsModel.getFieldFromId(fieldId).enthalpy_model != 'off':
            if self.mdl.useAdvancedThermalFluxes(fieldId) == True:
                for turbFlux in TurbulenceModelsDescription.ThermalTurbFluxModels:
                    self.modelCombo.addItem(self.tr(self.dicoM2V[turbFlux]), turbFlux)

            else:
                turb_flux = TurbulenceModelsDescription.ThermalTurbFluxModels[0]
                self.modelCombo.addItem(self.tr(self.dicoM2V[turbFlux]), turbFlux)
                self.modelCombo.disableItem(index=0)
        else:
            self.modelCombo.addItem(self.tr(self.dicoM2V['none']), 'none')
            self.modelCombo.setItem(str_view='none')

        editor.setMinimumSize(editor.sizeHint())
        editor.installEventFilter(self)
        return editor

    def setEditorData(self, comboBox, index):
        col = index.column()
        string = index.model().getData(index)[col]
        self.modelCombo.setItem(str_view=string)

    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        log.debug("TurbFluxDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModelMainFields class
#-------------------------------------------------------------------------------

class StandardItemModelTurbulence(QStandardItemModel):

    def __init__(self, mdl, case, dicoM2V, dicoV2M):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("Field label"),
                         self.tr("Carrier phase"),
                         self.tr("Turbulence model"),
                         self.tr("Thermal turbulent flux"),
                         self.tr("Reverse coupling")]

        self.setColumnCount(len(self.headers))

        self.tooltip = []

        self._data = []
        self.mdl      = mdl
        self.case     = case
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return None

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return data
            else:
                return None

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):

        # NoItemsFlags is used to have a grayed out option

        field = self.mdl.mainFieldsModel.list_of_fields[index.row()]
        fieldId = field.f_id
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 2 :
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        if index.column() == 3 :
            if self.mdl.useAdvancedThermalFluxes(fieldId) == True \
                and field.enthalpy_model != 'off':
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.NoItemFlags
        elif index.column() == 1 or index.column() == 4 :
            if field.flow_type == "continuous" :
                return Qt.NoItemFlags
            else :
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()
        FieldId = self.mdl.mainFieldsModel.list_of_fields[row].f_id

        # turbulence model
        if col == 2:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            self.mdl.setTurbulenceModel(FieldId, self.dicoV2M[new_pmodel])

            # Security check for models which cannot have a formula on inlet!
            if new_pmodel in ['mixing length', 'Q2-Q12 Tchen', 'R2-R12 Tchen']:
                from code_saturne.model.LocalizationModel import LocalizationModel
                from code_saturne.model.BoundaryNeptune import Boundary

                blm = LocalizationModel('BoundaryZone', self.case)

                for zone in blm.getZones():
                    if "inlet" in zone.getNature():
                        boundary = Boundary(zone.getNature(),
                                            zone.getLabel(),
                                            self.case,
                                            FieldId)

                        tc = boundary.getTurbulenceChoice(FieldId)
                        if tc == "formula":
                            boundary.setTurbulenceChoice(FieldId, "hydraulic_diameter")

            self.updateItem()

        # Turbulent thermal fluxes (for continuous phases)
        elif col == 3:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            self.mdl.setThermalTurbulentFlux(FieldId, self.dicoV2M[new_pmodel])
            self.updateItem()

        # two way coupling
        elif col == 4:
            new_pmodel = from_qvariant(value, to_text_string)
            self._data[row][col] = new_pmodel
            self.mdl.setTwoWayCouplingModel(FieldId, self.dicoV2M[new_pmodel])
            self.updateItem()

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, fieldId):
        """
        Add/load a field in the model.
        """
        row = self.rowCount()
        field = self.mdl.mainFieldsModel.getFieldFromId(fieldId)

        label        = field.label
        carrier_id   = field.carrier_id
        carrier = self.mdl.mainFieldsModel.getFieldFromId(carrier_id)
        carrierLabel = carrier.label
        turbulence = self.dicoM2V[self.mdl.getTurbulenceModel(fieldId)]
        turb_flux  = self.dicoM2V[self.mdl.getThermalTurbulentFlux(fieldId)]
        twoway     = self.dicoM2V[self.mdl.getTwoWayCouplingModel(fieldId)]

        field = [label, carrierLabel, turbulence, turb_flux, twoway]

        self._data.append(field)
        self.setRowCount(row+1)


    def updateItem(self):
        """
        update item
        """
        for id in self.mdl.mainFieldsModel.getFieldIdList() :
            self._data[int(id)-1][2] = self.dicoM2V[self.mdl.getTurbulenceModel(id)]
            self._data[int(id)-1][3] = self.dicoM2V[self.mdl.getThermalTurbulentFlux(id)]
            self._data[int(id)-1][4] = self.dicoM2V[self.mdl.getTwoWayCouplingModel(id)]


#-------------------------------------------------------------------------------
# turbulence class
#-------------------------------------------------------------------------------

class TurbulenceView(QWidget, Ui_Turbulence):
    """
    Main fields layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_Turbulence.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = TurbulenceModel(self.case)

        # Dico
        self.dicoM2V= {"none"                        : 'none',
                       "mixing_length"               : 'mixing length',
                       "k-epsilon"                   : 'k-epsilon',
                       "rij-epsilon_ssg"             : 'Rij-epsilon SSG',
                       "rij-epsilon_ebrsm"           : 'Rij-epsilon EBRSM',
                       "k-epsilon_linear_production" : 'k-epsilon linear production',
                       "les_smagorinsky"             : 'LES (Smagorinsky)',
                       "les_wale"                    : 'LES (WALE)',
                       "q2-q12-tchen"                : 'Q2-Q12 Tchen',
                       "q2-q12"                      : 'Q2-Q12',
                       "r2-q12"                      : 'R2-Q12',
                       "r2-r12-tchen"                : 'R2-R12 Tchen',
                       "separate_phase"              : 'separate phase',
                       "separate_phase_cond"         : 'separate phase cond',
                       "small_inclusions"            : 'small inclusions',
                       "large_inclusions"            : 'large inclusions',
                       "sgdh"                        : 'SGDH',
                       "ggdh"                        : 'GGDH'}

        self.dicoV2M= {"none"                        : 'none',
                       "mixing length"               : 'mixing_length',
                       "k-epsilon"                   : 'k-epsilon',
                       "Rij-epsilon SSG"             : 'rij-epsilon_ssg',
                       "Rij-epsilon EBRSM"           : 'rij-epsilon_ebrsm',
                       "k-epsilon linear production" : 'k-epsilon_linear_production',
                       "LES (Smagorinsky)"           : 'les_smagorinsky',
                       "LES (WALE)"                  : 'les_wale',
                       "Q2-Q12 Tchen"                : 'q2-q12-tchen',
                       "Q2-Q12"                      : 'q2-q12',
                       "R2-Q12"                      : 'r2-q12',
                       "R2-R12 Tchen"                : 'r2-r12-tchen',
                       "separate phase"              : 'separate_phase',
                       "separate phase cond"         : 'separate_phase_cond',
                       "small inclusions"            : 'small_inclusions',
                       "large inclusions"            : 'large_inclusions',
                       "SGDH"                        : 'sgdh',
                       "GGDH"                        : 'ggdh'}


        # Validators
        validatorMix = DoubleValidator(self.lineEditMixingLength, min = 0.0)
        validatorMix.setExclusiveMin(False)
        self.lineEditMixingLength.setValidator(validatorMix)

        self.tableModelTurbulence = StandardItemModelTurbulence(self.mdl, self.case,
                                                                self.dicoM2V, self.dicoV2M)
        self.tableViewTurbulence.setModel(self.tableModelTurbulence)
        self.tableViewTurbulence.resizeColumnsToContents()
        self.tableViewTurbulence.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewTurbulence.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTurbulence.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTurbulence.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewTurbulence.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTurbulence.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTurbulence.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
        self.tableViewTurbulence.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewTurbulence.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateTurbulence = TurbulenceDelegate(self.tableViewTurbulence, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateTurbFlux = TurbFluxDelegate(self.tableViewTurbulence, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateCoupling = CouplingDelegate(self.tableViewTurbulence, self.mdl, self.dicoM2V, self.dicoV2M)
        self.tableViewTurbulence.setItemDelegateForColumn(2, delegateTurbulence)
        self.tableViewTurbulence.setItemDelegateForColumn(3, delegateTurbFlux)
        self.tableViewTurbulence.setItemDelegateForColumn(4, delegateCoupling)

        # Combo models
        self.modelContinuousCoupling = ComboModel(self.comboBoxContinuousCoupling, 1, 1)
        self.modelContinuousCoupling.addItem(self.tr('none'), 'none')
        if self.mdl.mainFieldsModel.getPhaseChangeTransferStatus() == "off":
            self.modelContinuousCoupling.addItem(self.tr("separate phases"), "separate_phase")
        else:
            self.modelContinuousCoupling.addItem(self.tr("separate phases + cond"), "separate_phase_cond")

        # hide groupBoxMixingLength
        self.groupBoxMixingLength.hide()

        # Connect signals to slots
        self.tableModelTurbulence.dataChanged.connect(self.dataChanged)
        self.lineEditMixingLength.textChanged[str].connect(self.slotMixingLength)
        self.tableViewTurbulence.clicked.connect(self.slotChangeSelection)
        self.comboBoxContinuousCoupling.activated[str].connect(self.slotContinuousCoupling)

        # hide/show groupBoxContinuousCoupling
        if len(self.mdl.mainFieldsModel.getContinuousFieldList()) >=2 :
            self.groupBoxContinuousCoupling.show()
            model = self.mdl.getContinuousCouplingModel()
            self.modelContinuousCoupling.setItem(str_model=model)
        else :
            self.groupBoxContinuousCoupling.hide()

        for fieldId in self.mdl.mainFieldsModel.getFieldIdList():
            self.tableModelTurbulence.newItem(fieldId)

        self.case.undoStartGlobal()


    def slotChangeSelection(self, text=None):
        """
        detect change selection to update
        """
        row = self.tableViewTurbulence.currentIndex().row()
        self.update(row)


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewTurbulence.resizeColumnsToContents()
        self.tableViewTurbulence.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewTurbulence.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewTurbulence.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        row = self.tableViewTurbulence.currentIndex().row()
        self.update(row)


    def update(self, row):
        """
        show groupBoxMixingLength if necessary
        """
        fieldId = self.mdl.mainFieldsModel.list_of_fields[row].f_id
        if fieldId != 0:
            turbModel = self.mdl.getTurbulenceModel(fieldId)
            if turbModel == "mixing_length" :
                self.groupBoxMixingLength.show()
                self.lineEditMixingLength.setText(str(self.mdl.getMixingLength(fieldId)))
            else :
                self.groupBoxMixingLength.hide()
                # If the user chose GGDH for a RSM turbulence model, we set
                # the thermal fluxes model back to SGDH for consistency
                if 'rij-epsilon' not in turbModel and \
                        self.mdl.getThermalTurbulentFlux(fieldId) == 'ggdh':
                    self.mdl.setThermalTurbulentFlux(fieldId, 'sgdh')


    @pyqtSlot(str)
    def slotMixingLength(self, text):
        """
        Update the mixing length
        """
        table_id = self.tableViewTurbulence.currentIndex().row()
        fieldId = self.mdl.mainFieldsModel.list_of_fields[table_id].f_id
        if self.lineEditMixingLength.validator().state == QValidator.Acceptable:
            mix = from_qvariant(text, float)
            self.mdl.setMixingLength(fieldId, mix)


    @pyqtSlot(str)
    def slotContinuousCoupling(self, text):
        """
        define continuous/continuous coupling modele
        """
        value = self.modelContinuousCoupling.dicoV2M[text]
        log.debug("slotContinuousCoupling -> %s" % value)
        self.mdl.setContinuousCouplingModel(value)



#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
