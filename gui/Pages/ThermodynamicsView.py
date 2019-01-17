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
This module define the 'Thermodynamics' page.
This module contains the following classes:
- MaterialsDelegate
- MethodDelegate
- StandardItemModelProperty
- ThermodynamicsView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

EOS = 1
try:
   import eosAva
except:
   EOS = 0
else :
   import eosAva


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
from code_saturne.Base.QtPage import DoubleValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from Thermodynamics import Ui_Thermodynamics
from ThermodynamicsModel import *
from code_saturne.Pages.MainFieldsModel import MainFieldsModel
from code_saturne.Pages.SpeciesModel import SpeciesModel
from code_saturne.Pages.OutputFieldsModel import OutputFieldsModel
from code_saturne.Pages.NonCondensableModel import NonCondensableModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermodynamicsView")
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
        if EOS == 1 :
            self.ava = eosAva.EosAvailable()

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        self.modelCombo.addItem(self.tr(self.dicoM2V["user_material"]), 'user_material')
        fieldId= index.row() + 1
        # suppress perfect gas
        tmp = ["Argon", "Nitrogen", "Hydrogen", "Oxygen", "Helium", "Air"]
        if EOS == 1 and self.mdl.getFieldNature(fieldId) != "solid" :
            fls = self.ava.whichFluids()
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
                model.setData(idx, to_qvariant(self.modelCombo.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


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
        if EOS == 1 :
            self.ava = eosAva.EosAvailable()

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)

        fieldId= index.row() + 1
        if self.mdl.getMaterials(fieldId) == "user_material" :
            self.modelCombo.addItem(self.tr(self.dicoM2V["user_properties"]), 'user_properties')
        else :
            if EOS == 1 :
                material = self.mdl.getMaterials(fieldId)
                self.ava.setMethods(material)
                fls = self.ava.whichMethods()
                for fli in fls:
                    if fli != "Ovap" and fli != "Flica4" and fli != "StiffenedGas":
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
                model.setData(idx, to_qvariant(self.modelCombo.dicoM2V[value]), Qt.DisplayRole)


    def tr(self, text):
        return text


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
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if index.column() in (0, 1, 2, 3):
                if data:
                    return to_qvariant(data)
                else:
                    return to_qvariant()

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 or index.column() == 3:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            if len(self.ncond.getNonCondensableLabelList()) > 0 and (index.row()==0 or index.row()==1):
                return Qt.ItemIsSelectable
            else :
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


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
            new_mat = from_qvariant(value, to_text_string)
            self._data[row][col] = new_mat
            self.mdl.setMaterials(FieldId, self.dicoV2M[new_mat])
            # refresh method to default value
            self.mdl.updateMethod(FieldId, oldMaterial)
            self._data[row][2] = self.dicoM2V[str(self.mdl.getMethod(FieldId))]

        # Method
        elif col == 2:
            new_met = from_qvariant(value, to_text_string)
            self._data[row][col] = new_met
            self.mdl.setMethod(FieldId, self.dicoV2M[new_met])

        self.updateReference(row)
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
        self._data[row][3] = self.mdl.updateReference(fieldId)


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
        row       =  self.rowCount()
        label     = self.mdl.getLabel(fieldId)
        material  = self.dicoM2V[self.mdl.getMaterials(fieldId)]
        method    = self.dicoM2V[self.mdl.getMethod(fieldId)]
        reference = self.mdl.updateReference(fieldId)

        field = [label, material, method, reference]

        self._data.append(field)
        self.setRowCount(row+1)


    def getMethod(self, row):
        return self._data[row][2]


    def getLabel(self, row):
        return self._data[row][0]


#-------------------------------------------------------------------------------
# MainFieldsView class
#-------------------------------------------------------------------------------

class ThermodynamicsView(QWidget, Ui_Thermodynamics):
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

    surface_tension="""# water-air at 20°C
sigma = 0.075;

"""
    temperature="""
Cp = 1000
temperature = enthalpy / 1000;
"""

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_Thermodynamics.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.mdl = ThermodynamicsModel(self.case)
        self.notebook = NotebookModel(self.case)
        self.ncond = NonCondensableModel(self.case)

        # Dico
        self.dicoM2V= {"user_material" : 'user material',
                       "user_properties" : 'user properties'}

        self.dicoV2M= {"user material" : 'user_material',
                       "user properties" : 'user_properties'}
        if EOS == 1 :
            self.ava = eosAva.EosAvailable()
            fls = self.ava.whichFluids()
            for fli in fls:
                self.dicoM2V[fli] = fli
                self.dicoV2M[fli] = fli
                self.ava.setMethods(fli)

                flls = self.ava.whichMethods()
                for flli in flls:
                    if flli != "Ovap" and flli != "Flica4" and flli != "StiffenedGas":
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

        self.tableViewProperties.setItemDelegateForColumn(1, delegateMaterials)
        self.tableViewProperties.setItemDelegateForColumn(2, delegateMethod)

        # Combo models

        self.modelDensity             = ComboModel(self.comboBoxDensity, 2, 1)
        self.modelViscosity           = ComboModel(self.comboBoxViscosity, 2, 1)
        self.modelSpecificHeat        = ComboModel(self.comboBoxSpecificHeat, 2, 1)
        self.modelThermalConductivity = ComboModel(self.comboBoxThermalConductivity, 2, 1)
        self.modelSurfaceTension      = ComboModel(self.comboBoxSurfaceTension, 2, 1)

        self.modelHsat                = ComboModel(self.comboBoxHsat, 2, 1)
        self.modeldHsatdp             = ComboModel(self.comboBoxdHsatdp, 2, 1)

        self.modelDensity.addItem(self.tr('constant'), 'constant')
        self.modelDensity.addItem(self.tr('user law'), 'user_law')
        self.modelViscosity.addItem(self.tr('constant'), 'constant')
        self.modelViscosity.addItem(self.tr('user law'), 'user_law')
        self.modelSpecificHeat.addItem(self.tr('constant'), 'constant')
        self.modelSpecificHeat.addItem(self.tr('user law'), 'user_law')
        self.modelThermalConductivity.addItem(self.tr('constant'), 'constant')
        self.modelThermalConductivity.addItem(self.tr('user law'), 'user_law')
        self.modelSurfaceTension.addItem(self.tr('constant'), 'constant')
        self.modelSurfaceTension.addItem(self.tr('user law'), 'user_law')
        self.modelHsat.addItem(self.tr('Liquid'), 'Liquid')
        self.modelHsat.addItem(self.tr('Gas'), 'Gas')
        self.modeldHsatdp.addItem(self.tr('Liquid'), 'Liquid')
        self.modeldHsatdp.addItem(self.tr('Gas'), 'Gas')

        # Validators

        validatorRho = DoubleValidator(self.lineEditDensity, min = 0.0)
        validatorMu = DoubleValidator(self.lineEditViscosity, min = 0.0)
        validatorCp = DoubleValidator(self.lineEditSpecificHeat, min = 0.0)
        validatorAl = DoubleValidator(self.lineEditThermalConductivity, min = 0.0)
        validatorSt = DoubleValidator(self.lineEditSurfaceTension, min = 0.0)
        validatorEm = DoubleValidator(self.lineEditEmissivity, min = 0.0)
        validatorEc = DoubleValidator(self.lineEditElastCoef, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorMu.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)
        validatorSt.setExclusiveMin(True)
        validatorEm.setExclusiveMin(False)
        validatorEc.setExclusiveMin(False)

        self.lineEditDensity.setValidator(validatorRho)
        self.lineEditViscosity.setValidator(validatorMu)
        self.lineEditSpecificHeat.setValidator(validatorCp)
        self.lineEditThermalConductivity.setValidator(validatorAl)
        self.lineEditSurfaceTension.setValidator(validatorSt)
        self.lineEditEmissivity.setValidator(validatorEm)
        self.lineEditElastCoef.setValidator(validatorEc)

        # Connections

        self.lineEditDensity.textChanged[str].connect(self.slotRho)
        self.lineEditViscosity.textChanged[str].connect(self.slotMu)
        self.lineEditSpecificHeat.textChanged[str].connect(self.slotCp)
        self.lineEditThermalConductivity.textChanged[str].connect(self.slotAl)
        self.lineEditSurfaceTension.textChanged[str].connect(self.slotSt)
        self.lineEditEmissivity.textChanged[str].connect(self.slotEmissivity)
        self.lineEditElastCoef.textChanged[str].connect(self.slotElastCoef)
        self.pushButtonDensity.clicked.connect(self.slotFormulaRho)
        self.pushButtonViscosity.clicked.connect(self.slotFormulaMu)
        self.pushButtonSpecificHeat.clicked.connect(self.slotFormulaCp)
        self.pushButtonThermalConductivity.clicked.connect(self.slotFormulaAl)
        self.pushButtonSurfaceTension.clicked.connect(self.slotFormulaSt)
        self.checkBoxRadiativeTransfer.clicked.connect(self.slotRadTrans)
        self.comboBoxDensity.activated[str].connect(self.slotStateRho)
        self.comboBoxViscosity.activated[str].connect(self.slotStateMu)
        self.comboBoxSpecificHeat.activated[str].connect(self.slotStateCp)
        self.comboBoxThermalConductivity.activated[str].connect(self.slotStateAl)
        self.comboBoxSurfaceTension.activated[str].connect(self.slotStateSt)
        self.tableModelProperties.dataChanged.connect(self.dataChanged)
        self.tableViewProperties.clicked.connect(self.slotChangeSelection)
        self.pushButtonEOS.clicked.connect(self.slotEOS)
        self.pushButtonTemperature.clicked.connect(self.slotFormulaTemperature)
        self.pushButtondRodp.clicked.connect(self.slotFormuladrodp)
        self.pushButtondRodh.clicked.connect(self.slotFormuladrodh)
        self.comboBoxHsat.activated[str].connect(self.slotStateHsat)
        self.comboBoxdHsatdp.activated[str].connect(self.slotStatedHsatdp)
        self.pushButtonHsat.clicked.connect(self.slotFormulaHsat)
        self.pushButtondHsatdp.clicked.connect(self.slotFormuladHsatdp)
        self.pushButtonTsat.clicked.connect(self.slotFormulaTsat)
        self.pushButtondTsatdp.clicked.connect(self.slotFormuladTsatdp)
        self.pushButtonHlat.clicked.connect(self.slotFormulaHlat)

        # load Field
        for fieldId in self.mdl.getFieldIdList():
            self.tableModelProperties.newItem(fieldId)

        self.tableViewProperties.resizeColumnsToContents()
        self.tableViewProperties.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewProperties.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewProperties.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        # hide groupBoxConstantProperties
        self.groupBoxConstantProperties.hide()
        if EOS == 1 :
            self.groupBoxEOS.show()
        else :
            self.groupBoxEOS.hide()

        self.groupBoxEauvap.hide()

        # hide or not surface tension
        self.__updateSurfTension()

        self.case.undoStartGlobal()


    def __updateSurfTension(self):
        """
        surface tension only if we don't use eauvap and EOS
        """
        self.groupBoxNoFieldProperties.hide()

        need = 1
        # if we use EOS for one gas phas we have surface tension
        # calculated by EOS component
        for field in self.mdl.getGasPhaseList():
            if self.mdl.getMaterials(field) != "user_material":
                need = 0

        if need:
            self.groupBoxNoFieldProperties.show()
            tag = 'surface_tension'
            currentFluid = 'none'
            sym = 'SurfaceTension'

            __model  = getattr(self, "model"      + sym)
            __line   = getattr(self, "lineEdit"   + sym)
            __button = getattr(self, "pushButton" + sym)

            choice = self.mdl.getPropertyMode(currentFluid, tag)
            __model.setItem(str_model=choice)

            if choice == 'constant':
                __button.setEnabled(False)
                __line.setEnabled(True)
                __line.setText(str(self.mdl.getInitialValue(currentFluid, tag)))
                __button.setStyleSheet("background-color: None")
            elif choice == "user_law":
                __button.setEnabled(True)
                __line.setEnabled(True)
                __line.setText(str(self.mdl.getInitialValue(currentFluid, tag)))
                exp = self.mdl.getFormula('none', 'surface_tension')
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


    def __changeChoice(self, text, sym, tag):
        """
        Input variable state
        """
        if tag == 'surface_tension':
            currentFluid = 'none'
        else:
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
            exp = self.mdl.getFormula(currentFluid, tag)
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


    @pyqtSlot(str)
    def slotStateSt(self, text):
        """
        Method to call 'getState' with correct arguements for 'St'
        """
        self.__changeChoice(str(text), 'SurfaceTension', 'surface_tension')


    def slotChangeSelection(self, text=None):
        """
        detect change selection to update constant properties
        """
        row = self.tableViewProperties.currentIndex().row()
        self.update(row)


    def dataChanged(self, topLeft, bottomRight):
        self.tableViewProperties.resizeColumnsToContents()
        self.tableViewProperties.resizeRowsToContents()
        if QT_API == "PYQT4":
            self.tableViewProperties.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewProperties.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        row = self.tableViewProperties.currentIndex().row()
        self.update(row)


    @pyqtSlot(str)
    def slotStateHsat(self, text):
        """
        Method to change phas for Hsat
        """
        choice = self.comboBoxHsat.currentText()
        if choice == "Liquid":
            name = "SaturationEnthalpyLiquid"
        else:
            name = "SaturationEnthalpyGas"
        label = self.m_out.getVariableLabel('none', name)
        exp = self.mdl.getFormula('none', name)
        if exp:
            self.pushButtonHsat.setStyleSheet("background-color: green")
            self.pushButtonHsat.setToolTip(exp)
        else:
            self.pushButtonHsat.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def slotStatedHsatdp(self, text):
        """
        Method to change phas for Hsat
        """
        choice = self.comboBoxdHsatdp.currentText()
        if choice == "Liquid":
            name = "d_Hsat_d_P_Liquid"
        else:
            name = "d_Hsat_d_P_Gas"
        label = self.m_out.getVariableLabel('none', name)
        exp = self.mdl.getFormula('none', name)
        if exp:
            self.pushButtondHsatdp.setStyleSheet("background-color: green")
            self.pushButtondHsatdp.setToolTip(exp)
        else:
            self.pushButtondHsatdp.setStyleSheet("background-color: red")


    def update(self, row):
        """
        show groupBoxConstantProperties if necessary
        """
        method = self.tableModelProperties.getMethod(row)
        fieldId = row + 1
        self.currentFluid = fieldId
        self.groupBoxEauvap.hide()
        if method == "user properties" :
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
                self.comboBoxSpecificHeat.setEnabled(True)
                self.lineEditSpecificHeat.setReadOnly(False)
                self.lineEditSpecificHeat.setEnabled(True)

                self.comboBoxThermalConductivity.setEnabled(True)
                self.lineEditThermalConductivity.setReadOnly(False)
                self.lineEditThermalConductivity.setEnabled(True)

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
                    __line.setEnabled(True)
                    __line.setText(str(self.mdl.getInitialValue(fieldId, tag)))
                    __button.setStyleSheet("background-color: None")
                elif method2 == "user_law" :
                    __button.setEnabled(True)
                    __line.setEnabled(False)
                    __line.setText(str(self.mdl.getInitialValue(fieldId, tag)))
                    exp = self.mdl.getFormula(fieldId, tag)
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
                exp = self.mdl.getFormula(fieldId, 'd_rho_d_P')
                if exp:
                    self.pushButtondRodp.setStyleSheet("background-color: green")
                    self.pushButtondRodp.setToolTip(exp)
                else:
                    self.pushButtondRodp.setStyleSheet("background-color: red")
                exp = self.mdl.getFormula(fieldId, 'd_rho_d_h')
                if exp:
                    self.pushButtondRodh.setStyleSheet("background-color: green")
                    self.pushButtondRodh.setToolTip(exp)
                else:
                    self.pushButtondRodh.setStyleSheet("background-color: red")

            # Test for water / steam
            if (MainFieldsModel(self.case).getPredefinedFlow() != "None" and \
                MainFieldsModel(self.case).getPredefinedFlow() != "particles_flow" and \
                row < 2):

                self.groupBoxEauvap.show()

                # Tsat
                label = self.m_out.getVariableLabel('none', 'SaturationTemperature')
                exp = self.mdl.getFormula('none', 'SaturationTemperature')
                if exp:
                    self.pushButtonTsat.setStyleSheet("background-color: green")
                    self.pushButtonTsat.setToolTip(exp)
                else:
                    self.pushButtonTsat.setStyleSheet("background-color: red")

                # Tsatdp
                label = self.m_out.getVariableLabel('none', 'd_Tsat_d_P')
                exp = self.mdl.getFormula('none', 'd_Tsat_d_P')
                if exp:
                    self.pushButtondTsatdp.setStyleSheet("background-color: green")
                    self.pushButtondTsatdp.setToolTip(exp)
                else:
                    self.pushButtondTsatdp.setStyleSheet("background-color: red")

                # Hlat
                label = self.m_out.getVariableLabel('none', 'LatentHeat')
                exp = self.mdl.getFormula('none', 'LatentHeat')
                if exp:
                    self.pushButtonHlat.setStyleSheet("background-color: green")
                    self.pushButtonHlat.setToolTip(exp)
                else:
                    self.pushButtonHlat.setStyleSheet("background-color: red")

                # Hsat
                choice = self.comboBoxHsat.currentText()
                if choice == "Liquid":
                    name = "SaturationEnthalpyLiquid"
                else:
                    name = "SaturationEnthalpyGas"
                label = self.m_out.getVariableLabel('none', name)
                exp = self.mdl.getFormula('none', name)
                if exp:
                    self.pushButtonHsat.setStyleSheet("background-color: green")
                    self.pushButtonHsat.setToolTip(exp)
                else:
                    self.pushButtonHsat.setStyleSheet("background-color: red")

                # Hsatdp
                choice = self.comboBoxdHsatdp.currentText()
                if choice == "Liquid":
                    name = "d_Hsat_d_P_Liquid"
                else:
                    name = "d_Hsat_d_P_Gas"
                label = self.m_out.getVariableLabel('none', name)
                exp = self.mdl.getFormula('none', name)
                if exp:
                    self.pushButtondHsatdp.setStyleSheet("background-color: green")
                    self.pushButtondHsatdp.setToolTip(exp)
                else:
                    self.pushButtondHsatdp.setStyleSheet("background-color: red")

                # TODO gerer les mots cle liquid / gas

            # Temperature / enthalpy law
            self.groupBoxTemperature.hide()
            if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
                self.groupBoxTemperature.show()

                fieldId = self.currentFluid
                label = self.m_out.getVariableLabel(str(fieldId), 'temperature')
                exp = self.mdl.getFormula(fieldId, 'temperature')
                if exp:
                    self.pushButtonTemperature.setStyleSheet("background-color: green")
                    self.pushButtonTemperature.setToolTip(exp)
                else:
                    self.pushButtonTemperature.setStyleSheet("background-color: red")
        else :
            self.groupBoxConstantProperties.hide()

        self.__updateSurfTension()


    @pyqtSlot(str)
    def slotRho(self, text):
        """
        Update the density
        """
        fieldId = self.currentFluid
        if self.lineEditDensity.validator().state == QValidator.Acceptable:
            rho = from_qvariant(text, float)
            self.mdl.setInitialValueDensity(fieldId, rho)


    @pyqtSlot(str)
    def slotMu(self, text):
        """
        Update the molecular viscosity
        """
        fieldId = self.currentFluid
        if self.lineEditViscosity.validator().state == QValidator.Acceptable:
            mu = from_qvariant(text, float)
            self.mdl.setInitialValueViscosity(fieldId,mu)


    @pyqtSlot(str)
    def slotCp(self, text):
        """
        Update the specific heat
        """
        fieldId = self.currentFluid
        if self.lineEditSpecificHeat.validator().state == QValidator.Acceptable:
            cp = from_qvariant(text, float)
            self.mdl.setInitialValueHeat(fieldId,cp)


    @pyqtSlot(str)
    def slotAl(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditThermalConductivity.validator().state == QValidator.Acceptable:
            al = from_qvariant(text, float)
            self.mdl.setInitialValueCond(fieldId,al)


    @pyqtSlot(str)
    def slotSt(self, text):
        """
        Update the surface tension
        """
        fieldId = self.currentFluid
        if self.lineEditSurfaceTension.validator().state == QValidator.Acceptable:
            st = from_qvariant(text, float)
            self.mdl.setInitialValueTens(st)


    @pyqtSlot(str)
    def slotEmissivity(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditEmissivity.validator().state == QValidator.Acceptable:
            em = from_qvariant(text, float)
            self.mdl.setInitialValueEmissivity(fieldId,em)


    @pyqtSlot(str)
    def slotElastCoef(self, text):
        """
        Update the thermal conductivity
        """
        fieldId = self.currentFluid
        if self.lineEditElastCoef.validator().state == QValidator.Acceptable:
            ec = from_qvariant(text, float)
            self.mdl.setInitialValueElastCoef(fieldId,ec)


    def getFormulaComponents(self, fieldId, tag):
        """
        Get the formula components for a given tag
        """

        if tag == 'density':
            return self.getFormulaRhoComponents(fieldId)

        elif tag == 'molecular_viscosity':
            return self.getFormulaMuComponents(fieldId)

        elif tag == 'specific_heat':
            return self.getFormulaCpComponents(fieldId)

        elif tag == 'thermal_conductivity':
            return self.getFormulaAlComponents(fieldId)

        else:
            msg = 'Formula is not available for field %s_%s in MEG' % (tag,str(fieldId))
            raise Exception(msg)


    def getFormulaRhoComponents(self, fieldId):
        """
        User formula for density
        """
        exp = self.mdl.getFormula(fieldId, 'density')
        if not exp:
            exp = "rho = 1.8;"
        req = [('rho', 'Density')]
        exa = ThermodynamicsView.density

        symbols_rho = []
        for s in self.list_scalars:
           symbols_rho.append(s)

        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols_rho.append((label, "enthalpy_"+str(fieldId)))
        rho0_value = self.mdl.getInitialValue(fieldId, 'density')
        symbols_rho.append(('rho0', 'Density (reference value) = '+str(rho0_value)))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols_rho.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols_rho.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols_rho, exa


    @pyqtSlot()
    def slotFormulaRho(self):
        """
        User formula for density
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_rho, exa = self.getFormulaRhoComponents(fieldId)

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_rho,
                                examples   = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'density', result)
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(exp)


    def getFormulaMuComponents(self, fieldId):
        """
        User formula for molecular viscosity
        """
        exp = self.mdl.getFormula(fieldId, 'molecular_viscosity')
        if not exp:
            exp = "mu = 4.56e-05;"
        req = [('mu', 'Molecular Viscosity')]
        exa = ThermodynamicsView.molecular_viscosity

        symbols_mu = []
        for s in self.list_scalars:
           symbols_mu.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols_mu.append((label, 'enthalpy_'+str(fieldId)))
        mu0_val = self.mdl.getInitialValue(fieldId, 'molecular_viscosity')
        symbols_mu.append(('mu0', 'Viscosity (reference value) = '+str(mu0_val)))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols_mu.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols_mu.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols_mu, exa


    @pyqtSlot()
    def slotFormulaMu(self):
        """
        User formula for molecular viscosity
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_mu, exa = self.getFormulaMuComponents(fieldId)

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_mu,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMu -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'molecular_viscosity', result)
            self.pushButtonViscosity.setStyleSheet("background-color: green")
            self.pushButtonViscosity.setToolTip(exp)


    def getFormulaCpComponents(self, fieldId):
        """
        User formula for specific heat
        """
        exp = self.mdl.getFormula(fieldId, 'specific_heat')

        if not exp:
            exp = "cp = 4000.;"
        req = [('cp', 'Specific heat')]
        exa = ThermodynamicsView.specific_heat

        symbols_cp = []
        for s in self.list_scalars:
           symbols_cp.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols_cp.append((label, "enthalpy_"+str(fieldId)))
        cp0_val = self.mdl.getInitialValue(fieldId, "specific_heat")
        symbols_cp.append(('cp0', 'Specific heat (reference value) = '+str(cp0_val)))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols_cp.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols_cp.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols_cp, exa


    @pyqtSlot()
    def slotFormulaCp(self):
        """
        User formula for specific heat
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_cp, exa = self.getFormulaCpComponents(fieldId)

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_cp,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaCp -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'specific_heat', result)
            self.pushButtonSpecificHeat.setStyleSheet("background-color: green")
            self.pushButtonSpecificHeat.setToolTip(exp)


    def getFormulaAlComponents(self, fieldId):
        """
        User formula for thermal conductivity
        """
        exp = self.mdl.getFormula(fieldId, 'thermal_conductivity')
        if not exp:
            exp = "lambda = 1.e-5;"
        req = [('lambda', 'Thermal conductivity')]
        exa = ThermodynamicsView.thermal_conductivity

        symbols_al = []
        for s in self.list_scalars:
           symbols_al.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols_al.append((label, 'enthalpy_'+str(fieldId)))
        l0_val = self.mdl.getInitialValue(fieldId, 'thermal_conductivity')
        symbols_al.append(('lambda0', 'Thermal conductivity (reference value) = '+str(l0_val)))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols_al.append((s, s))

        for (nme, val) in self.notebook.getNotebookList():
            symbols_al.append((nme, 'value (notebook) = ' + str(val)))

        return exp, req, self.list_scalars, symbols_al, exa


    @pyqtSlot()
    def slotFormulaAl(self):
        """
        User formula for thermal conductivity
        """
        fieldId = self.currentFluid

        exp, req, sca, symbols_al, exa = self.getFormulaAlComponents(fieldId)

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_al,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaAl -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'thermal_conductivity', result)
            self.pushButtonThermalConductivity.setStyleSheet("background-color: green")
            self.pushButtonThermalConductivity.setToolTip(exp)


    @pyqtSlot()
    def slotFormulaSt(self):
        """
        User formula for surface tension
        """
        exp = self.mdl.getFormula('none', 'surface_tension')
        if not exp:
            exp = "sigma = 0.075;"
        req = [('sigma', 'Surface Tension')]
        exa = ThermodynamicsView.surface_tension

        symbols_st = []
        for s in self.list_scalars:
           symbols_st.append(s)
        for fieldId in self.mdl.getFieldIdList():
            if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
                label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
                symbols_st.append((label, self.tr("Field variable")))
        symbols_st.append(('sigma0', 'Surface tension (reference value)'))

        for s in self.m_spe.getScalarNameList():
              symbols_st.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_st,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaSt -> %s" % str(result))
            self.mdl.setFormula('none', 'surface_tension', result)
            self.pushButtonSurfaceTension.setStyleSheet("background-color: green")
            self.pushButtonSurfaceTension.setToolTip(exp)


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
        if EOS == 1:
            import subprocess
            import cs_config

            cfg = cs_config.config()
            if cfg.libs['eos'].have == "yes":
                eos_bin_dir = os.path.join(cfg.libs['eos'].prefix, "bin")
                if os.path.isdir(eos_bin_dir):
                    command = os.path.join(eos_bin_dir, "eos_gui")
                    self.runProcess = subprocess.Popen(command, shell=True)


    @pyqtSlot()
    def slotFormulaTemperature(self):
        """
        User formula for temperature as a function of enthalpy
        """
        fieldId = self.currentFluid
        label = self.m_out.getVariableLabel(str(fieldId), 'temperature')
        exp = self.mdl.getFormula(fieldId, 'temperature')
        if not exp:
            exp = label + " = 273.15;"
        req = [(label, 'temperature')]
        exa = ThermodynamicsView.temperature

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, self.tr("Field variable")))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTemperature -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'temperature', result)
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(result)


    @pyqtSlot()
    def slotFormuladrodp(self):
        """
        User formula for d(ro) / dp (compressible flow)
        """
        fieldId = self.currentFluid
        exp = self.mdl.getFormula(fieldId, 'd_rho_d_P')
        if not exp:
            exp = "d_rho_d_P = 0.;"
        req = [('d_rho_d_P', 'Partial derivative of density with respect to pressure')]
        exa = "d_rho_d_P = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, self.tr("Field variable")))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladrodp -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'd_rho_d_P', result)
            self.pushButtondRodp.setStyleSheet("background-color: green")
            self.pushButtondRodp.setToolTip(result)


    @pyqtSlot()
    def slotFormuladrodh(self):
        """
        User formula for d(ro) / dh (compressible flow)
        """
        fieldId = self.currentFluid
        exp = self.mdl.getFormula(fieldId, 'd_rho_d_h')
        if not exp:
            exp = "d_rho_d_h = 0.;"
        req = [('d_rho_d_h', 'Partial derivative of density with respect to enthalpy')]
        exa = "d_rho_d_h = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)
        if MainFieldsModel(self.case).getEnergyResolution(fieldId) == "on":
            label = self.m_out.getVariableLabel(str(fieldId), "enthalpy")
            symbols.append((label, self.tr("Field variable")))

        for s in self.m_spe.getScalarByFieldId(fieldId):
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladrodh -> %s" % str(result))
            self.mdl.setFormula(str(fieldId), 'd_rho_d_h', result)
            self.pushButtondRodh.setStyleSheet("background-color: green")
            self.pushButtondRodh.setToolTip(result)


    @pyqtSlot()
    def slotFormulaHsat(self):
        """
        User formula for enthalpy of saturation (water-steam flow)
        """
        choice = self.comboBoxHsat.currentText()
        if choice == "Liquid":
            name = "SaturationEnthalpyLiquid"
        else:
            name = "SaturationEnthalpyGas"

        label = self.m_out.getVariableLabel('none', name)
        exp = self.mdl.getFormula('none', name)
        if not exp:
            exp = label + " = 0.;"
        req = [(label, 'enthalpy of saturation')]
        exa = label + " = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaHsat -> %s" % str(result))
            self.mdl.setFormula('none', name, result)
            self.pushButtonHsat.setStyleSheet("background-color: green")
            self.pushButtonHsat.setToolTip(result)


    @pyqtSlot()
    def slotFormuladHsatdp(self):
        """
        User formula for d(Hsat) / dp (water-steam flow)
        """
        choice = self.comboBoxdHsatdp.currentText()
        if choice == "Liquid":
            name = "d_Hsat_d_P_Liquid"
        else:
            name = "d_Hsat_d_P_Gas"

        label = self.m_out.getVariableLabel('none', name)
        exp = self.mdl.getFormula('none', name)
        if not exp:
            exp = label + " = 0.;"
        req = [(label, 'Partial derivative of enthalpy of saturation with respect to pressure')]
        exa = label + " = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladHsatdp -> %s" % str(result))
            self.mdl.setFormula('none', name, result)
            self.pushButtondHsatdp.setStyleSheet("background-color: green")
            self.pushButtondHsatdp.setToolTip(result)


    @pyqtSlot()
    def slotFormulaTsat(self):
        """
        User formula for temperature of saturation (water-steam flow)
        """
        label = self.m_out.getVariableLabel('none', 'SaturationTemperature')
        exp = self.mdl.getFormula('none', 'SaturationTemperature')
        if not exp:
            exp = label + " = 273.15;"
        req = [(label, 'temperature of saturation')]
        exa = label + " = 273.15;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTsat -> %s" % str(result))
            self.mdl.setFormula('none', 'SaturationTemperature', result)
            self.pushButtonTsat.setStyleSheet("background-color: green")
            self.pushButtonTsat.setToolTip(result)


    @pyqtSlot()
    def slotFormuladTsatdp(self):
        """
        User formula for d(Tsat) / dp (water-steam flow)
        """
        label = self.m_out.getVariableLabel('none', 'd_Tsat_d_P')
        exp = self.mdl.getFormula('none', 'd_Tsat_d_P')
        if not exp:
            exp = label + " = 0.;"
        req = [(label , 'Partial derivative of temperature of saturation with respect to pressure')]
        exa = label + " = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormuladTsatdp -> %s" % str(result))
            self.mdl.setFormula('none', 'd_Tsat_d_P', result)
            self.pushButtondTsatdp.setStyleSheet("background-color: green")
            self.pushButtondTsatdp.setToolTip(result)


    @pyqtSlot()
    def slotFormulaHlat(self):
        """
        User formula for latent heat (water-steam flow)
        """
        label = self.m_out.getVariableLabel('none', 'LatentHeat')
        exp = self.mdl.getFormula('none', 'LatentHeat')
        if not exp:
            exp = label + " = 0.;"
        req = [(label, 'latent heat')]
        exa = label + " = 0.;"

        symbols = []
        for s in self.list_scalars:
           symbols.append(s)

        for s in self.m_spe.getScalarNameList():
            symbols.append((s, self.tr("Additional species")))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaHlat -> %s" % str(result))
            self.mdl.setFormula('none', 'LatentHeat', result)
            self.pushButtonHlat.setStyleSheet("background-color: green")
            self.pushButtonHlat.setToolTip(result)


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

