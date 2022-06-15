# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module define the 'ThermodynamicsField' page.
This module contains the following classes:
- MaterialsDelegate
- MethodDelegate
- StandardItemModelProperty
- ThermodynamicsFieldView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

from code_saturne.model.EosWrapper import eosWrapper

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
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel
from code_saturne.gui.base.QtPage import to_text_string
from code_saturne.gui.case.ThermodynamicsField import Ui_ThermodynamicsField
from code_saturne.model.ThermodynamicsModel import *
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.SpeciesModel import SpeciesModel
from code_saturne.model.OutputFieldsModel import OutputFieldsModel
from code_saturne.model.NonCondensableModel import NonCondensableModel
from code_saturne.model.LocalizationModel import LocalizationModel

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

from code_saturne.model.EosWrapper import eosWrapper
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermodynamicsFieldView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Combo box delegate for the material
#-------------------------------------------------------------------------------

class MaterialsDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(MaterialsDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M
        self.eos = eosWrapper()

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        self.modelCombo.addItem(self.tr(self.dicoM2V["user_material"]), 'user_material')
        fieldId= index.row() + 1
        # suppress perfect gas
        tmp = ["Argon", "Nitrogen", "Hydrogen", "Oxygen", "Helium", "Air"]
        if self.mdl.getFieldNature(fieldId) != "solid" and self.mdl.checkEOSRequirements(fieldId):
            fls = self.eos.getListOfFluids()
            for fli in fls:
                if fli not in tmp:
                    tmp.append(fli)
                    self.modelCombo.addItem(self.tr(fli), fli)
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
        log.debug("MaterialsDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.modelCombo.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the method
#-------------------------------------------------------------------------------


class MethodDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(MethodDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M
        self.eos = eosWrapper()

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fieldId= index.row() + 1
        if self.mdl.getMaterials(fieldId) == "user_material" :
            self.modelCombo.addItem(self.tr(self.dicoM2V["user_properties"]), 'user_properties')
        else :
            material = self.mdl.getMaterials(fieldId)
            fls = self.eos.getFluidMethods(material)

            # If non condensable-gases, filter the list, only Cathare and
            # Cathare2 tables are allowed with EOS
#            if len(self.parent.mdl.ncond.getNonCondensableByFieldId(fieldId)) > 0:
#                fls = [fli for fli in fls if fli in ("Cathare", "Cathare2")]

            for fli in fls:
                self.modelCombo.addItem(self.tr(fli),fli)

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
        log.debug("MethodDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.modelCombo.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the method
#-------------------------------------------------------------------------------


class ReferenceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, mdl, dicoM2V, dicoV2M):
        super(ReferenceDelegate, self).__init__(parent)
        self.parent   = parent
        self.mdl      = mdl
        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M
        self.eos = eosWrapper()

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fieldId= index.row() + 1

        material = self.mdl.getMaterials(fieldId)
        method   = self.mdl.getMethod(fieldId)
        if material == "user_material" :
            self.modelCombo.addItem(self.tr(self.dicoM2V["user_material"]), 'user_material')
        else :
            if self.eos.isActive():
                phase = self.mdl.getFieldNature(fieldId)

                if phase == "liquid":
                    ref = self.eos.getLiquidReferences(material, method)
                elif phase == "gas":
                    ref = self.eos.getVaporReferences(material, method)

                for r in ref:
                    self.modelCombo.addItem(self.tr(r), r)

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
        log.debug("ReferenceDelegate value = %s"%value)

        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, self.modelCombo.dicoM2V[value], Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandardItemModelProperty class
#-------------------------------------------------------------------------------

class StandardItemModelProperty(QStandardItemModel):

    def __init__(self, mdl, ncond, dicoM2V, dicoV2M):
        """
        """
        QStandardItemModel.__init__(self)


        self.headers = [ self.tr("Field label"),
                         self.tr("Material"),
                         self.tr("Method"),
                         self.tr("Reference")]

        self.setColumnCount(len(self.headers))

        self.tooltip = []

        self._data = []
        self.mdl = mdl
        self.ncond = ncond

        self.dicoM2V  = dicoM2V
        self.dicoV2M  = dicoV2M


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.ToolTipRole:
            return None

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if index.column() in (0, 1, 2, 3):
                if data:
                    return data
                else:
                    return None

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags
        # Lock fields with non condensable gas
        field_id = index.row() + 1
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()
        FieldId = row + 1

        # Materials
        if col == 1:
            oldMaterial = self._data[row][col]
            new_mat = value
            self._data[row][col] = new_mat
            self.mdl.setMaterials(FieldId, self.dicoV2M[new_mat])
            # refresh method to default value
            self.mdl.updateMethod(FieldId, oldMaterial)
            self._data[row][2] = self.dicoM2V[str(self.mdl.getMethod(FieldId))]

        # Method
        elif col == 2:
            new_met = value
            self._data[row][col] = new_met
            self.mdl.setMethod(FieldId, self.dicoV2M[new_met])

        # Reference
        elif col == 3:
            new_ref = value
            self._data[row][col] = new_ref
            self.mdl.setFluidReference(FieldId, new_ref)

        if row < 2:
            self.updateTable(0)
            self.updateTable(1)

        self.dataChanged.emit(index, index)
        return True


    def updateTable(self, row):
        """
        update reference. Use only for EOS
        """
        fieldId = row + 1
        self._data[row][1] = self.dicoM2V[str(self.mdl.getMaterials(fieldId))]
        self._data[row][2] = self.dicoM2V[str(self.mdl.getMethod(fieldId))]
        self._data[row][3] = self.mdl.getFluidReference(fieldId)


    def updateReference(self, row):
        """
        update reference. Use only for EOS
        """
        fieldId = row + 1
        self._data[row][3] = self.mdl.updateReference(fieldId)


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, fieldId):
        """
        load field in the model
        """
        row = self.rowCount()
        label = self.mdl.getLabel(fieldId)
        try:
            material = self.dicoM2V[self.mdl.getMaterials(fieldId)]
        except KeyError:
            material = self.dicoM2V[self.mdl.defaultValues(fieldId)["material"]]
            self.mdl.setMaterials(fieldId, material)
        try:
            method = self.dicoM2V[self.mdl.getMethod(fieldId)]
        except KeyError:
            method = self.dicoM2V[self.mdl.defaultValues(fieldId)["method"]]
            self.mdl.setMethod(fieldId, method)
        reference = self.mdl.updateReference(fieldId)

        field = [label, material, method, reference]

        self._data.append(field)
        self.setRowCount(row + 1)


    def getMethod(self, row):
        return self._data[row][2]


    def getLabel(self, row):
        return self._data[row][0]


#-------------------------------------------------------------------------------
# MainFieldsView class
#-------------------------------------------------------------------------------

class ThermodynamicsFieldView(QWidget, Ui_ThermodynamicsField):
    """
    Thermodynamics layout.
    """
    density = """# Density of air

rho = 1.293 * (273.15 / temperature);

# density for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

rho1 = 1.25051;
rho2 = 1.7832;
A = (Y1 / rho1) + (Y2 /rho2);
rho = 1.0 / A;

"""
    molecular_viscosity="""# Sutherland's Formula
# Gas             Cst    T0      mu0
# air             120    291.15  18.27e-6
# nitrogen        111    300.55  17.81e-6
# oxygen          127    292.25  20.18e-6
# carbon dioxide  240    293.15  14.8e-6
# carbon monoxide 118    288.15  17.2e-6
# hydrogen        72     293.85  8.76e-6
# ammonia         370    293.15  9.82e-6
# sulfur dioxide  416    293.65  12.54e-6
# helium          79.4   273     19e-6

CST = 120;
T0 = 291.15;
mu0 = 18.27e-6;

if ( temperature > 0 && temperature < 555) {
mu = mu0 * (T0+CST / temperature+CST) * (temperature/T0)^(3./2.);
} else {
mu = -999.0;
}

"""
    specific_heat="""# specific heat for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

Cp1 = 520.3;
Cp2 = 1040.0;
cp = Y1 * Cp1 + Y2 *Cp2;
"""
    thermal_conductivity="""# oxygen
lambda = 6.2e-5 * temperature + 8.1e-3;

# nitrogen
lambda = 6.784141e-5 * temperature + 5.564317e-3;

# hydrogen
lambda = 4.431e-4 * temperature + 5.334e-2;

"""

    temperature="""
# if total enthalpy is used you need to use the specific enthalpy insted
# of enthalpy. Specific enthalpy is computed using:
# hspec = enthalpy - 0.5*square_norm(U)

Cp = 1000
temperature = enthalpy / 1000;
"""

    def __init__(self, parent = None): #FIXME!!, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ThermodynamicsField.__init__(self)
        self.setupUi(self)

        self.case      = None
        self.mdl       = None
        self.notebook  = None
        self.zone      = None
        self.zone_name = None
        self.zone_id   = None

        self.eos = eosWrapper()

    def setup(self, case, zone_name):
        self.case = case
        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_name = zone.getLabel()
                self.zone_id   = zone.getCodeNumber()
        self.case.undoStopGlobal()

        self.mdl = ThermodynamicsModel(self.case)
        self.notebook = NotebookModel(self.case)
        self.ncond = NonCondensableModel(self.case)

        is_main_zone = (zone_name == "all_cells")
        # Dico
        self.dicoM2V= {"user_material" : 'user material',
                       "user_properties" : 'user properties'}

        self.dicoV2M= {"user material" : 'user_material',
                       "user properties" : 'user_properties'}

        eos_used = False
        for field_id in self.mdl.getFieldIdList():
            if self.mdl.checkEOSRequirements(field_id):
                eos_used = True
                break
        if eos_used:
            fls = self.eos.getListOfFluids()
            for fli in fls:
                self.dicoM2V[fli] = fli
                self.dicoV2M[fli] = fli

                for flli in self.eos.getFluidMethods(fli):
                    self.dicoM2V[flli] = flli
                    self.dicoV2M[flli] = flli

        self.list_scalars = []

        self.m_spe = SpeciesModel(self.case)
        self.m_out = OutputFieldsModel(self.case)
        self.currentFluid = 0

        label = self.m_out.getVariableLabel("none", "pressure")
        self.list_scalars.append(('pressure', label))

        self.tableModelProperties = StandardItemModelProperty(self.mdl, self.ncond, self.dicoM2V, self.dicoV2M)
        self.tableViewProperties.setModel(self.tableModelProperties)
        self.tableViewProperties.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewProperties.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateMaterials = MaterialsDelegate(self.tableViewProperties, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateMethod    = MethodDelegate(self.tableViewProperties, self.mdl, self.dicoM2V, self.dicoV2M)
        delegateReference = ReferenceDelegate(self.tableViewProperties,
                                              self.mdl,
                                              self.dicoM2V,
                                              self.dicoV2M)

        self.tableViewProperties.setItemDelegateForColumn(1, delegateMaterials)
        self.tableViewProperties.setItemDelegateForColumn(2, delegateMethod)
        self.tableViewProperties.setItemDelegateForColumn(3, delegateReference)

        # Combo models

        self.modelDensity             = ComboModel(self.comboBoxDensity, 2, 1)
        self.modelViscosity           = ComboModel(self.comboBoxViscosity, 2, 1)
        self.modelSpecificHeat        = ComboModel(self.comboBoxSpecificHeat, 2, 1)
        self.modelThermalConductivity = ComboModel(self.comboBoxThermalConductivity, 2, 1)


        self.modelDensity.addItem(self.tr('constant'), 'constant')
        self.modelDensity.addItem(self.tr('user law'), 'user_law')
        self.modelViscosity.addItem(self.tr('constant'), 'constant')
        self.modelViscosity.addItem(self.tr('user law'), 'user_law')
        self.modelSpecificHeat.addItem(self.tr('constant'), 'constant')
        self.modelSpecificHeat.addItem(self.tr('user law'), 'user_law')
        self.modelThermalConductivity.addItem(self.tr('constant'), 'constant')
        self.modelThermalConductivity.addItem(self.tr('user law'), 'user_law')

        # Validators

        validatorRho = DoubleValidator(self.lineEditDensity, min = 0.0)
        validatorMu = DoubleValidator(self.lineEditViscosity, min = 0.0)
        validatorCp = DoubleValidator(self.lineEditSpecificHeat, min = 0.0)
        validatorAl = DoubleValidator(self.lineEditThermalConductivity, min = 0.0)
        validatorEm = DoubleValidator(self.lineEditEmissivity, min = 0.0)
        validatorEc = DoubleValidator(self.lineEditElastCoef, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorMu.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)
        validatorEm.setExclusiveMin(False)
        validatorEc.setExclusiveMin(False)

        self.lineEditDensity.setValidator(validatorRho)
        self.lineEditViscosity.setValidator(validatorMu)
        self.lineEditSpecificHeat.setValidator(validatorCp)
        self.lineEditThermalConductivity.setValidator(validatorAl)
        self.lineEditEmissivity.setValidator(validatorEm)
        self.lineEditElastCoef.setValidator(validatorEc)

        # Connections

        self.lineEditDensity.textChanged[str].connect(self.slotRho)
        self.lineEditViscosity.textChanged[str].connect(self.slotMu)
        self.lineEditSpecificHeat.textChanged[str].connect(self.slotCp)
        self.lineEditThermalConductivity.textChanged[str].connect(self.slotAl)
        self.lineEditEmissivity.textChanged[str].connect(self.slotEmissivity)
        self.lineEditElastCoef.textChanged[str].connect(self.slotElastCoef)
        self.pushButtonDensity.clicked.connect(self.slotFormulaRho)
        self.pushButtonViscosity.clicked.connect(self.slotFormulaMu)
        self.pushButtonSpecificHeat.clicked.connect(self.slotFormulaCp)
        self.pushButtonThermalConductivity.clicked.connect(self.slotFormulaAl)
        self.checkBoxRadiativeTransfer.clicked.connect(self.slotRadTrans)

        self.comboBoxes = {}
        self.comboBoxes['Rho']     = self.comboBoxDensity
        self.comboBoxes['Mu']      = self.comboBoxViscosity
        self.comboBoxes['Cp']      = self.comboBoxSpecificHeat
        self.comboBoxes['Al']      = self.comboBoxThermalConductivity

        # Connect combo-boxes and disable them if not main zone
        for k in self.comboBoxes.keys():
            self.comboBoxes[k].activated[str].connect(getattr(self,"slotState"+k))
            self.comboBoxes[k].setEnabled(is_main_zone)

        self.tableModelProperties.dataChanged.connect(self.dataChanged)
        self.tableViewProperties.clicked.connect(self.slotChangeSelection)
        self.pushButtonEOS.clicked.connect(self.slotEOS)
        self.pushButtonTemperature.clicked.connect(self.slotFormulaTemperature)
        self.pushButtondRodp.clicked.connect(self.slotFormuladrodp)
        self.pushButtondRodh.clicked.connect(self.slotFormuladrodh)

        # load Field
        for fieldId in self.mdl.getFieldIdList():
            self.tableModelProperties.newItem(fieldId)

        self.tableViewProperties.resizeColumnsToContents()
        self.tableViewProperties.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewProperties.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewProperties.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)

        self.initializeWidget()

        self.case.undoStartGlobal()

    def initializeWidget(self):

        # hide groupBoxConstantProperties
        self.groupBoxConstantProperties.hide()
        if self.eos.isActive():
            self.groupBoxEOS.show()
        else:
            self.groupBoxEOS.hide()

        self.labelNonCondensableWarning.hide()
        self.labelIdenticalMaterialWarning.hide()
        if len(self.ncond.getNonCondensableLabelList()) > 0:
            self.labelNonCondensableWarning.show()
        if self.mdl.checkIdenticalMaterialsRequirements():
            self.labelIdenticalMaterialWarning.show()

    def __changeChoice(self, text, sym, tag):
        """
        Input variable state
        """

        currentFluid = self.currentFluid

        __model  = getattr(self, "model"      + sym)
        __line   = getattr(self, "lineEdit"   + sym)
        __combo  = getattr(self, "comboBox"   + sym)
        __button = getattr(self, "pushButton" + sym)

        choice = __model.dicoV2M[text]
        log.debug("__changeChoice -> %s, %s" % (text, choice))

        if choice == 'constant':
            __button.setEnabled(False)
            __line.setEnabled(True)
            __line.setText(str(self.mdl.getInitialValue(currentFluid, tag)))
            __button.setStyleSheet("background-color: None")
        elif choice == "user_law":
            __button.setEnabled(True)
            __line.setEnabled(True)
            __line.setText(str(self.mdl.getInitialValue(currentFluid, tag)))
            exp = self.mdl.getFormula(currentFluid, tag, zone=self.zone_id)
            if exp:
                __button.setStyleSheet("background-color: green")
                __button.setToolTip(exp)
            else:
                __button.setStyleSheet("background-color: red")
        else :
            __button.setEnabled(False)
            __line.setEnabled(False)
            __line.setText("")
            __button.setStyleSheet("background-color: None")

        self.mdl.setPropertyMode(currentFluid, tag, choice)


    @pyqtSlot(str)
    def slotStateRho(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        self.__changeChoice(str(text), 'Density', 'density')


    @pyqtSlot(str)
    def slotStateMu(self, text):
        """
        Method to call 'getState' with correct arguements for 'Mu'
        """
        self.__changeChoice(str(text), 'Viscosity', 'molecular_viscosity')


    @pyqtSlot(str)
    def slotStateCp(self, text):
        """
        Method to call 'getState' with correct arguements for 'Cp'
        """
        self.__changeChoice(str(text), 'SpecificHeat', 'specific_heat')


    @pyqtSlot(str)
    def slotStateAl(self, text):
        """
        Method to call 'getState' with correct arguements for 'Al'
        """
        self.__changeChoice(str(text), 'ThermalConductivity', 'thermal_conductivity')


    def slotChangeSelection(self, text=None):
        """
        detect change selection to update constant properties
        """
        if self.tableViewProperties.currentIndex().isValid():
            row = self.tableViewProperties.currentIndex().row()
            self.update(row)
        else:
            self.groupBoxConstantProperties.hide()
            self.groupBoxEOS.hide()


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewProperties.resizeColumnsToContents()
        self.tableViewProperties.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewProperties.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewProperties.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        row = self.tableViewProperties.currentIndex().row()
        self.update(row)


    def update(self, row):
        """
        show groupBoxConstantProperties if necessary
        """
        is_main_zone = (self.zone_name == "all_cells")

        method = self.tableModelProperties.getMethod(row)
        fieldId = row + 1
        self.currentFluid = fieldId
        if method == "user properties" :
            self.groupBoxEOS.hide()
            self.groupBoxConstantProperties.show()

            # Deactivate Thermal variable if no energy resolution
            mfm = MainFieldsModel(self.case)
            if mfm.getEnergyResolution(fieldId) == 'off':
                self.comboBoxSpecificHeat.setEnabled(False)
                self.lineEditSpecificHeat.setReadOnly(True)
                self.lineEditSpecificHeat.setEnabled(False)

                self.comboBoxThermalConductivity.setEnabled(False)
                self.lineEditThermalConductivity.setReadOnly(True)
                self.lineEditThermalConductivity.setEnabled(False)
            else:
                self.comboBoxSpecificHeat.setEnabled(is_main_zone)
                self.lineEditSpecificHeat.setReadOnly(False)
                self.lineEditSpecificHeat.setEnabled(is_main_zone)

                self.comboBoxThermalConductivity.setEnabled(is_main_zone)
                self.lineEditThermalConductivity.setReadOnly(False)
                self.lineEditThermalConductivity.setEnabled(is_main_zone)

            self.groupBoxCompressible.hide()
            if self.mdl.getFieldNature(fieldId) == "solid":
                self.groupBoxSolidProp.show()
                self.lineEditEmissivity.setText(str(self.mdl.getInitialValue(fieldId, 'emissivity')))
                self.lineEditElastCoef.setText(str(self.mdl.getInitialValue(fieldId, 'elasticity')))
                isRadiativeTransfer  = self.mdl.getRadiativeTransferStatus(fieldId) == "on"
                self.checkBoxRadiativeTransfer.setChecked(isRadiativeTransfer)
            else :
                self.groupBoxSolidProp.hide()

            list = [('density', 'Density'),
                    ('molecular_viscosity', 'Viscosity'),
                    ('specific_heat', 'SpecificHeat'),
                    ('thermal_conductivity', 'ThermalConductivity')]
            for tag, symbol in list :
                __line   = getattr(self, "lineEdit" + symbol)
                __button = getattr(self, "pushButton" + symbol)
                __model  = getattr(self, "model"      + symbol)

                method2 = self.mdl.getPropertyMode(fieldId, tag)

                __model.setItem(str_model=method2)

                if method2 == "constant" :
                    __button.setEnabled(False)
                    __line.setEnabled(is_main_zone)
                    __line.setText(str(self.mdl.getInitialValue(fieldId, tag)))
                    __button.setStyleSheet("background-color: None")
                elif method2 == "user_law" :
                    __button.setEnabled(True)
                    __line.setEnabled(False)
                    __line.setText(str(self.mdl.getInitialValue(fieldId, tag)))
                    exp = self.mdl.getFormula(fieldId, tag, zone=self.zone_id)
                    if exp:
                        __button.setStyleSheet("background-color: green")
                        __button.setToolTip(exp)
                    else:
                        __button.setStyleSheet("background-color: red")
                else :
                    __button.setEnabled(False)
                    __line.setEnabled(False)
                    __line.setText("")
                    __button.setStyleSheet("background-color: None")

            # Test for compressible flow
            if MainFieldsModel(self.case).getCompressibleStatus(fieldId) == "on":
                self.groupBoxCompressible.show()
                exp = self.mdl.getFormula(fieldId, 'd_rho_d_P', zone=self.zone_id)
                if exp:
                    self.pushButtondRodp.setStyleSheet("background-color: green")
                    self.pushButtondRodp.setToolTip(exp)
                else:
                    self.pushButtondRodp.setStyleSheet("background-color: red")
                exp = self.mdl.getFormula(fieldId, 'd_rho_d_h', zone=self.zone_id)
                if exp:
                    self.pushButtondRodh.setStyleSheet("background-color: green")
                    self.pushButtondRodh.setToolTip(exp)
                else:
                    self.pushButtondRodh.setStyleSheet("background-color: red")

            # Temperature / enthalpy law
            self.groupBoxTemperature.hide()
            if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
                self.groupBoxTemperature.show()

                fieldId = self.currentFluid
                label = self.m_out.getVariableLabel(str(fieldId), 'temperature')
                exp = self.mdl.getFormula(fieldId, 'temperature', zone=self.zone_id)
                if exp:
                    self.pushButtonTemperature.setStyleSheet("background-color: green")
                    self.pushButtonTemperature.setToolTip(exp)
                else:
                    self.pushButtonTemperature.setStyleSheet("background-color: red")
        else :
            self.groupBoxEOS.show()
            self.groupBoxConstantProperties.hide()


    @pyqtSlot(str)
    def slotRho(self, text):
        """
        Update the density
        """
        fieldId = self.currentFluid
        if self.lineEditDensity.validator().state == QValidator.Acceptable:
            rho = float(text)
            self.mdl.setInitialValueDensity(fieldId, rho)


    @pyqtSlot(str)
    def slotMu(self, text):
        """
        Update the molecular viscosity
        """
        fieldId = self.currentFluid
        if self.lineEditViscosity.validator().state == QValidator.Acceptable:
            mu = float(text)
            self.mdl.setInitialValueViscosity(fieldId,mu)


    @pyqtSlot(str)
    def slotCp(self, text):
        """
        Update the specific heat
        """
        fieldId = self.currentFluid
        if self.lineEditSpecificHeat.validator().state == QValidator.Acceptable:
            cp = float(text)
            self.mdl.setInitialValueHeat(fieldId,cp)


    @pyqtSlot(str)
    def slotAl(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditThermalConductivity.validator().state == QValidator.Acceptable:
            al = float(text)
            self.mdl.setInitialValueCond(fieldId,al)


    @pyqtSlot(str)
    def slotEmissivity(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditEmissivity.validator().state == QValidator.Acceptable:
            em = float(text)
            self.mdl.setInitialValueEmissivity(fieldId,em)


    @pyqtSlot(str)
    def slotElastCoef(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditElastCoef.validator().state == QValidator.Acceptable:
            ec = float(text)
            self.mdl.setInitialValueElastCoef(fieldId,ec)


    @pyqtSlot()
    def slotFormulaRho(self):
        """
        User formula for density
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_rho = self.mdl.getFormulaRhoComponents(fieldId, zone=self.zone_id)

        exa = ThermodynamicsFieldView.density

        vname = "density_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_rho,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'density', result, zone=self.zone_id)
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(exp)


    @pyqtSlot()
    def slotFormulaMu(self):
        """
        User formula for molecular viscosity
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_mu = self.mdl.getFormulaMuComponents(fieldId, zone=self.zone_id)

        exa = ThermodynamicsFieldView.molecular_viscosity

        vname = "molecular_viscosity_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_mu,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMu -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'molecular_viscosity', result, zone=self.zone_id)
            self.pushButtonViscosity.setStyleSheet("background-color: green")
            self.pushButtonViscosity.setToolTip(exp)


    @pyqtSlot()
    def slotFormulaCp(self):
        """
        User formula for specific heat
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_cp = self.mdl.getFormulaCpComponents(fieldId, zone=self.zone_id)

        exa = ThermodynamicsFieldView.specific_heat

        vname = "specific_heat_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_cp,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaCp -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'specific_heat', result, zone=self.zone_id)
            self.pushButtonSpecificHeat.setStyleSheet("background-color: green")
            self.pushButtonSpecificHeat.setToolTip(exp)


    @pyqtSlot()
    def slotFormulaAl(self):
        """
        User formula for thermal conductivity
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_al = self.mdl.getFormulaAlComponents(fieldId, zone=self.zone_id)

        exa = ThermodynamicsFieldView.thermal_conductivity

        vname = "thermal_conductivity_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_al,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaAl -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'thermal_conductivity', result, zone=self.zone_id)
            self.pushButtonThermalConductivity.setStyleSheet("background-color: green")
            self.pushButtonThermalConductivity.setToolTip(exp)

    @pyqtSlot(bool)
    def slotRadTrans(self, checked):
        """
        check box for radiative transfer
        """
        fieldId = self.currentFluid
        status = 'off'
        if checked:
            status = 'on'
        self.mdl.setRadiativeTransferStatus(fieldId, status)


    @pyqtSlot()
    def slotEOS(self):
        """
        call EOS GUI
        """
        if self.eos.isActive():
            command = None
            try:
                cfg = self.case.case['package'].config
                if cfg.libs['eos'].have:
                    eos_bin_dir = os.path.join(cfg.libs['eos'].prefix, "bin")
                    if os.path.isdir(eos_bin_dir):
                        command = os.path.join(eos_bin_dir, "eos_gui")
            except Exception:  # if case/package not available (should not happen)
                print("Warning: package configuration not available")
                pass

            if command != None:
                import subprocess
                self.runProcess = subprocess.Popen(command, shell=True)


    @pyqtSlot()
    def slotFormulaTemperature(self):
        """
        User formula for temperature as a function of enthalpy
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols = self.mdl.getFormulaTemperatureComponents(fieldId, zone=self.zone_id)

        exa = ThermodynamicsFieldView.temperature

        vname = "temperature_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTemperature -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'temperature', result, zone=self.zone_id)
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(result)


    @pyqtSlot()
    def slotFormuladrodp(self):
        """
        User formula for d(ro) / dp (compressible flow)
        """
        fieldId = self.currentFluid
        exp, req, sca, symbols = self.mdl.getFormuladrodpComponents(fieldId, zone=self.zone_id)

        exa = "d_rho_d_P = 0.;"

        vname = "d_rho_d_P_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladrodp -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'd_rho_d_P', result, zone=self.zone_id)
            self.pushButtondRodp.setStyleSheet("background-color: green")
            self.pushButtondRodp.setToolTip(result)


    @pyqtSlot()
    def slotFormuladrodh(self):
        """
        User formula for d(ro) / dh (compressible flow)
        """
        fieldId = self.currentFluid
        exp, req, sca, symbols = self.mdl.getFormuladrodhComponents(fieldId, zone=self.zone_id)

        exa = "d_rho_d_h = 0.;"

        vname = "d_rho_d_h_%s" % (str(fieldId))
        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = vname,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladrodh -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'd_rho_d_h', result, zone=self.zone_id)
            self.pushButtondRodh.setStyleSheet("background-color: green")
            self.pushButtondRodh.setToolTip(result)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

