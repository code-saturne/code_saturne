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
This module contains the following classes and function:
- NameDelegate
- ChemicalFormulaDelegate
- ValueDelegate
- StandardItemModelSpecies
- GasCombustionView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging, os

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.base.QtPage import ComboModel, RegExpValidator
from code_saturne.gui.base.QtPage import DoubleValidator, from_qvariant
from code_saturne.gui.base.QtPage import to_text_string
from code_saturne.gui.base.QtPage import IntValidator
from code_saturne.gui.case.GasCombustionForm import Ui_GasCombustionForm
from code_saturne.model.GasCombustionModel import GasCombustionModel
from code_saturne.model.GasCombustionModel import ThermochemistryData

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("GasCombustionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the species label
#-------------------------------------------------------------------------------

class NameDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_pname = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_pname = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_pname = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_pname = str(editor.text())

            if new_pname in model.mdl.getSpeciesNamesList():
                default = {}
                default['name']  = self.old_pname
                default['list']   = model.mdl.getSpeciesNamesList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_pname = result['name']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_pname = self.old_pname

            model.setData(index, new_pname, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for the chemical formula
#-------------------------------------------------------------------------------

class ChemicalFormulaDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_pname = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_pname = ""
        rx = "[chonslCHONSL()][CHONSLchonsl()0-9]{0," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_pname = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_pname = str(editor.text())

            if new_pname in model.mdl.getSpeciesNamesList():
                default = {}
                default['name']  = self.old_pname
                default['list']   = model.mdl.getSpeciesNamesList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_pname = result['name']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_pname = self.old_pname

            model.setData(index, new_pname, Qt.DisplayRole)

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
                    model.setData(idx, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelSpecies(QStandardItemModel):
    """
    """
    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Species"),
                        self.tr("Chemical Formula"),
                        self.tr("Nb of moles (Fuel)"),
                        self.tr("Nb of moles (Oxidiser)"),
                        self.tr("Nb of moles (Product)"),
                        self.tr("Absorption Coeff")]

        self.setColumnCount(len(self.headers))

        self._data = []
        self.parent = parent
        self.mdl  = mdl


    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        if role == Qt.DisplayRole:
            return self._data[row][col]

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() != 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsSelectable


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

        # Label (nothing to do for the species label)
        if col == 0:
            pass

        # Chemical formula
        elif col == 1:
            ChemicalFormula = str(from_qvariant(value, to_text_string))
            self._data[row][col] = ChemicalFormula
            label = self._data[row][0]
            self.mdl.setSpeciesChemicalFormula(label, ChemicalFormula)

        # Fuel Composition
        elif col == 2:
            CompFuel = str(from_qvariant(value, to_text_string))
            self._data[row][col] = CompFuel
            label = self._data[row][0]
            self.mdl.setCompFuel(label, CompFuel)

        # Oxi Composition
        elif col == 3:
            CompOxi = str(from_qvariant(value, to_text_string))
            self._data[row][col] = CompOxi
            label = self._data[row][0]
            self.mdl.setCompOxi(label, CompOxi)

        # Product Composition
        elif col == 4:
            CompProd = str(from_qvariant(value, to_text_string))
            self._data[row][col] = CompProd
            label = self._data[row][0]
            self.mdl.setCompProd(label, CompProd)

        # Coeff absorption
        elif col == 5:
            CoeffAbsorp = str(from_qvariant(value, to_text_string))
            self._data[row][col] = CoeffAbsorp
            label = self._data[row][0]
            self.mdl.setCoeffAbsorp(label, CoeffAbsorp)

        self.dataChanged.emit(index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_name=None):
        """
        Add an item in the table view
        """
        row = self.rowCount()

        label = self.mdl.addSpecies(existing_name)

        ChemicalFormula = self.mdl.getSpeciesChemicalFormula(label)
        CompFuel = self.mdl.getCompFuel(label)
        CompOxi = self.mdl.getCompOxi(label)
        CompProd = self.mdl.getCompProd(label)
        CoeffAbsorp = self.mdl.getCoeffAbsorp(label)

        species = [label, ChemicalFormula, CompFuel, CompOxi, CompProd, CoeffAbsorp]

        self.setRowCount(row+1)
        self._data.append(species)


    def getItem(self, row):
        """
        Return the values for an item.
        """
        [label, ChemicalFormula, CompFuel, CompOxi, CompProd, CoeffAbsorp] = self._data[row]
        return label, ChemicalFormula, CompFuel, CompOxi, CompProd, CoeffAbsorp


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        log.debug("deleteItem row = %i " % row)

        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class GasCombustionView(QWidget, Ui_GasCombustionForm):
    """
    Class to open the Gas Combustion option Page.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_GasCombustionForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = GasCombustionModel(self.case)
        self.thermodata = ThermochemistryData(self.case)

        # Model for table View
        self.modelSpecies = StandardItemModelSpecies(self, self.thermodata)
        self.tableViewSpecies.setModel(self.modelSpecies)

        # Delegates
        delegateLabel            = NameDelegate(self.tableViewSpecies)
        delegateChemicalFormula  = ChemicalFormulaDelegate(self.tableViewSpecies)
        delegateCompFuel  = ValueDelegate(self.tableViewSpecies)
        delegateCompOxi   = ValueDelegate(self.tableViewSpecies)
        delegateCompProd  = ValueDelegate(self.tableViewSpecies)
        delegateCoeffAbsorp  = ValueDelegate(self.tableViewSpecies)

        self.tableViewSpecies.setItemDelegateForColumn(0, delegateLabel)
        self.tableViewSpecies.setItemDelegateForColumn(1, delegateChemicalFormula)
        self.tableViewSpecies.setItemDelegateForColumn(2, delegateCompFuel)
        self.tableViewSpecies.setItemDelegateForColumn(3, delegateCompOxi)
        self.tableViewSpecies.setItemDelegateForColumn(4, delegateCompProd)
        self.tableViewSpecies.setItemDelegateForColumn(5, delegateCoeffAbsorp)

        # tableView
        if QT_API == "PYQT4":
            self.tableViewSpecies.horizontalHeader().setResizeMode(QHeaderView.Stretch)
            self.tableViewSpecies.horizontalHeader().setResizeMode(0,QHeaderView.ResizeToContents)
            self.tableViewSpecies.horizontalHeader().setResizeMode(1,QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableViewSpecies.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.tableViewSpecies.horizontalHeader().setSectionResizeMode(0,QHeaderView.ResizeToContents)
            self.tableViewSpecies.horizontalHeader().setSectionResizeMode(1,QHeaderView.ResizeToContents)

        # Set models and number of elements for combo boxes
        self.modelGasCombustionOption = ComboModel(self.comboBoxGasCombustionOption,1,1)

        # Combo models to choose the mode to create the Janaf File

        self.userModeForChemicalReaction   = ComboModel(self.comboBoxUserChoice, 2, 1)

        self.userModeForChemicalReaction.addItem("Automatic definition", 'auto')
        self.userModeForChemicalReaction.addItem("Defined by user", 'user')

        # Connections
        self.comboBoxUserChoice.activated[str].connect(self.slotUserChoice)
        self.comboBoxGasCombustionOption.activated[str].connect(self.slotGasCombustionOption)
        self.pushButtonThermochemistryData.pressed.connect(self.__slotSearchThermochemistryData)
        self.radioButtonCreateJanafFile.clicked.connect(self.slotCreateJanafFile)
        self.lineEditNbPointsTabu.textChanged[str].connect(self.slotNbPointsTabu)
        self.lineEditMaximumTemp.textChanged[str].connect(self.slotMaximumTemp)
        self.lineEditMinimumTemp.textChanged[str].connect(self.slotMinimumTemp)
        self.pushButtonAddSpecies.clicked.connect(self.slotAddSpecies)
        self.pushButtonDeleteSpecies.clicked.connect(self.slotDeleteSpecies)
        self.pushButtonGenerateJanafFile.clicked.connect(self.slotGenerateJanafFile)
        self.lineEditFuel.textChanged[str].connect(self.slotFuel)
        self.lineEditO2.textChanged[str].connect(self.slotVolPropO2)
        self.lineEditN2.textChanged[str].connect(self.slotVolPropN2)
        self.lineEditCOyield.textChanged[str].connect(self.slotCOyield)
        self.lineEditCSyield.textChanged[str].connect(self.slotCSyield)
        self.modelSpecies.dataChanged.connect(self.dataChanged)

        # Validators
        validatorNbPointsTabu = IntValidator(self.lineEditNbPointsTabu, min=1)
        validatorMaximumTemp  = DoubleValidator(self.lineEditMaximumTemp, min=273.0)
        validatorMinimumTemp  = DoubleValidator(self.lineEditMinimumTemp, min=273.0)
        rx = "[chonlCHONL()][CHONLchonl()0-9]{0," + str(LABEL_LENGTH_MAX-1) + "}"
        validatorFuel         = RegExpValidator(self.lineEditFuel,QRegExp(rx))
        validatorO2  = DoubleValidator(self.lineEditO2, min=1e-12, max=1.0)
        validatorN2  = DoubleValidator(self.lineEditN2, min=0.0, max=1.0)
        validatorCOyield  = DoubleValidator(self.lineEditCOyield, min=0.0)
        validatorCSyield  = DoubleValidator(self.lineEditCSyield, min=0.0)

        self.lineEditNbPointsTabu.setValidator(validatorNbPointsTabu)
        self.lineEditMaximumTemp.setValidator(validatorMaximumTemp)
        self.lineEditMinimumTemp.setValidator(validatorMinimumTemp)
        self.lineEditFuel.setValidator(validatorFuel)
        self.lineEditO2.setValidator(validatorO2)
        self.lineEditN2.setValidator(validatorN2)
        self.lineEditCOyield.setValidator(validatorCOyield)
        self.lineEditCSyield.setValidator(validatorCSyield)

        NbPointsTabu = self.thermodata.getNbPointsTabu()
        MaximumTemp  = self.thermodata.getMaximumTemp()
        MinimumTemp  = self.thermodata.getMinimumTemp()
        Option_UserMode = self.thermodata.getUserModeForChemicalReaction()
        ChemicalFormulaFuel = self.thermodata.getChemicalFormulaFuel()
        VolPropO2 = self.thermodata.getVolPropO2()
        VolPropN2 = self.thermodata.getVolPropN2()
        COyield = self.thermodata.getCOyield()
        CSyield = self.thermodata.getCSyield()

        self.lineEditNbPointsTabu.setText(str(NbPointsTabu))
        self.lineEditMaximumTemp.setText(str(MaximumTemp))
        self.lineEditMinimumTemp.setText(str(MinimumTemp))
        self.lineEditFuel.setText(str(ChemicalFormulaFuel))
        self.lineEditO2.setText(str(VolPropO2))
        self.lineEditN2.setText(str(VolPropN2))
        self.lineEditCOyield.setText(str(COyield))
        self.lineEditCSyield.setText(str(CSyield))

        # Initialize Widgets

        self.tableViewSpecies.reset()
        self.modelSpecies = StandardItemModelSpecies(self, self.thermodata)
        self.tableViewSpecies.setModel(self.modelSpecies)

        model = self.mdl.getGasCombustionModel()

        if model == 'd3p':
            self.modelGasCombustionOption.addItem(self.tr("adiabatic model"), "adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("non adiabatic model"), "extended")
        elif model == 'ebu':
            self.modelGasCombustionOption.addItem(self.tr("reference Spalding model"), "spalding")
            self.modelGasCombustionOption.addItem(self.tr("extended model with enthalpy source term"), "enthalpy_st")
            self.modelGasCombustionOption.addItem(self.tr("extended model with mixture fraction transport"), "mixture_st")
            self.modelGasCombustionOption.addItem(self.tr("extended model with enthalpy and mixture fraction transport"), "enthalpy_mixture_st")
        elif model == 'lwp':
            self.modelGasCombustionOption.addItem(self.tr("reference two-peak model with adiabatic condition"), "2-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference two-peak model with enthalpy source term"), "2-peak_enthalpy")
            self.modelGasCombustionOption.addItem(self.tr("reference three-peak model with adiabatic condition"), "3-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference three-peak model with enthalpy source term"), "3-peak_enthalpy")
            self.modelGasCombustionOption.addItem(self.tr("reference four-peak model with adiabatic condition"), "4-peak_adiabatic")
            self.modelGasCombustionOption.addItem(self.tr("reference four-peak model with enthalpy source term"), "4-peak_enthalpy")

        option = self.mdl.getGasCombustionOption()
        self.modelGasCombustionOption.setItem(str_model= option)

        name = self.mdl.getThermoChemistryDataFileName()
        if name != None and name != "":
            self.labelThermochemistryFile.setText(str(name))
            self.pushButtonThermochemistryData.setStyleSheet("background-color: green")
        else:
            self.pushButtonThermochemistryData.setStyleSheet("background-color: red")

        self.radioButtonCreateJanafFile.hide()

        if self.thermodata.getCreateThermoDataFile() == 'on':
            self.radioButtonCreateJanafFile.setChecked(True)
            self.userModeForChemicalReaction.setItem(str_model= Option_UserMode)
            for label in self.thermodata.getSpeciesNamesList():
                self.modelSpecies.newItem(label)
        else:
            self.radioButtonCreateJanafFile.setChecked(False)

        # for the moment the option to create Janaf file in the GUI is only available with d3p
        if model == 'd3p':
            self.radioButtonCreateJanafFile.show()
            self.groupBoxTabulation.show()
            self.groupBoxChemicalReaction.show()
            self.groupBoxDefinedByUser.show()
            self.groupBoxGenerateDataFile.show()
            self.groupBoxAutomatic.show()

        self.slotCreateJanafFile()

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotGasCombustionOption(self, text):
        """
        Private slot.
        Binding method for gas combustion models.
        """
        option = self.modelGasCombustionOption.dicoV2M[str(text)]
        self.mdl.setGasCombustionOption(option)


    @pyqtSlot()
    def __slotSearchThermochemistryData(self):
        """
        Select a properties file of data for electric arc
        """
        data = self.case['data_path']
        if not data:
            data = "."
        title = self.tr("Thermochemistry file of data.")
        filetypes = self.tr("Thermochemistry (*dp_*);;All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, data, filetypes)[0]
        file = str(file)
        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(data):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the DATA directory")
            QMessageBox.information(self, title, msg)
        else:
            self.labelThermochemistryFile.setText(str(file))
            self.mdl.setThermoChemistryDataFileName(file)
            self.pushButtonThermochemistryData.setStyleSheet("background-color: green")

    @pyqtSlot()
    def slotCreateJanafFile(self):
        """
        Determine if the Thermochemistry file is created with the GUI.
        """
        if self.radioButtonCreateJanafFile.isChecked():
            self.thermodata.setCreateThermoDataFile("on")
            self.groupBoxTabulation.show()
            self.groupBoxChemicalReaction.show()
            self.lineEditNbPointsTabu.setText(str(self.thermodata.getNbPointsTabu()))
            self.lineEditMaximumTemp.setText(str(self.thermodata.getMaximumTemp()))
            self.lineEditMinimumTemp.setText(str(self.thermodata.getMinimumTemp()))
            self.slotUserChoice()
            return
        else:
            self.thermodata.setCreateThermoDataFile("off")
            self.groupBoxTabulation.hide()
            self.groupBoxChemicalReaction.hide()
            self.groupBoxDefinedByUser.hide()
            self.groupBoxGenerateDataFile.hide()
            self.groupBoxAutomatic.hide()

    @pyqtSlot(str)
    def slotUserChoice(self):
        """
        """
        model = self.userModeForChemicalReaction.dicoV2M[str(self.comboBoxUserChoice.currentText())]
        self.thermodata.setUserModeForChemicalReaction(model)
        self.groupBoxGenerateDataFile.show()
        if model == 'auto':
            self.groupBoxAutomatic.show()
            self.groupBoxDefinedByUser.hide()
            #update of the xml from the GUI
            self.slotFuel(self.lineEditFuel.text())
            self.slotVolPropO2(self.lineEditO2.text())
            self.slotVolPropN2(self.lineEditN2.text())
            self.slotCOyield(self.lineEditCOyield.text())
            self.slotCSyield(self.lineEditCSyield.text())
        elif model == 'user':
            self.groupBoxDefinedByUser.show()
            self.groupBoxAutomatic.hide()
            #update of the xml from the tableView
            row_tab = self.modelSpecies.rowCount()
            for row in range(row_tab):
                data = self.modelSpecies.getItem(row)
                self.thermodata.updateSpecies(data)

    @pyqtSlot(str)
    def slotNbPointsTabu(self, text):
        """
        Input Number of points for the tabulation (ENTH-TEMP)
        """
        if self.lineEditNbPointsTabu.validator().state == QValidator.Acceptable:
            NbPointsTabu = from_qvariant(text, int)
            self.thermodata.setNbPointsTabu(NbPointsTabu)

    @pyqtSlot(str)
    def slotMaximumTemp(self, text):
        """
        Input Maximum temperature for the tabulation (ENTH-TEMP)
        """
        if self.lineEditMaximumTemp.validator().state == QValidator.Acceptable:
            MaximumTemp = from_qvariant(text, float)
            self.thermodata.setMaximumTemp(MaximumTemp)

    @pyqtSlot(str)
    def slotMinimumTemp(self, text):
        """
        Input Minimum temperature for the tabulation (ENTH-TEMP)
        """
        if self.lineEditMinimumTemp.validator().state == QValidator.Acceptable:
            MinimumTemp = from_qvariant(text, float)
            self.thermodata.setMinimumTemp(MinimumTemp)

    @pyqtSlot(str)
    def slotFuel(self, text):
        """
        Input the chemical formula for the Fuel
        """
        if self.lineEditFuel.validator().state == QValidator.Acceptable:
            ChemicalFormula = from_qvariant(text, to_text_string)
            self.thermodata.setChemicalFormulaFuel(ChemicalFormula)

    @pyqtSlot(str)
    def slotVolPropO2(self, text):
        """
        Input volume proportion for O2
        """
        if self.lineEditO2.validator().state == QValidator.Acceptable:
            VolPropO2 = from_qvariant(text, float)
            self.thermodata.setVolPropO2(VolPropO2)
            self.lineEditN2.setText(str(round(1.0 - VolPropO2,12)))

    @pyqtSlot(str)
    def slotVolPropN2(self, text):
        """
        Input volume proportion for N2
        """
        if self.lineEditN2.validator().state == QValidator.Acceptable:
            VolPropN2 = from_qvariant(text, float)
            self.thermodata.setVolPropN2(VolPropN2)
            self.lineEditO2.setText(str(round(1.0 - VolPropN2,12)))

    @pyqtSlot(str)
    def slotCOyield(self, text):
        """
        Input the CO yield
        """
        if self.lineEditCOyield.validator().state == QValidator.Acceptable:
            COyield = from_qvariant(text, float)
            self.thermodata.setCOyield(COyield)

    @pyqtSlot(str)
    def slotCSyield(self, text):
        """
        Input the CS yield
        """
        if self.lineEditCSyield.validator().state == QValidator.Acceptable:
            CSyield = from_qvariant(text, float)
            self.thermodata.setCSyield(CSyield)

    @pyqtSlot()
    def slotAddSpecies(self):
        """
        Add a new item in the table when the 'Create' button is pushed.
        """
        self.tableViewSpecies.clearSelection()
        self.modelSpecies.newItem()

    @pyqtSlot()
    def slotDeleteSpecies(self):
        """
        Just delete the current selected entries from the table and
        of course from the XML file.
        """
        lst = []
        for index in self.tableViewSpecies.selectionModel().selectedRows():
            row = index.row()
            lst.append(row)

        lst.sort()
        lst.reverse()

        for row in lst:
            label = self.modelSpecies.getItem(row)[0]
            self.thermodata.deleteSpecies(label)
            self.modelSpecies.deleteItem(row)

        self.tableViewSpecies.clearSelection()

    @pyqtSlot()
    def slotGenerateJanafFile(self):
        """
        Generate the Thermochemistry file.
        """

        data = self.case['data_path']
        if not data:
            data = "."
        filename = "dp_ThermochemistryFromGui"
        file_path = os.path.join(data, filename)

        self.thermodata.WriteThermochemistryDataFile(file_path)

        self.mdl.setThermoChemistryDataFileName(filename)

        if self.thermodata.Error_GUI :
            self.labelThermochemistryFile.setText(str("Error in : "+filename))
            self.pushButtonGenerateJanafFile.setStyleSheet("background-color: red")
            self.pushButtonThermochemistryData.setStyleSheet("background-color: orange")
        else :
            self.labelThermochemistryFile.setText(str(filename))
            self.pushButtonGenerateJanafFile.setStyleSheet("")
            self.pushButtonThermochemistryData.setStyleSheet("background-color: green")


    @pyqtSlot("QModelIndex, QModelIndex")
    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewSpecies.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewSpecies.resizeColumnToContents(col)


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
