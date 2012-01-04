# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
This module contains the following classes:
- StandardItemModelCoals
- StandardItemModelClasses
- CoalCombustionView
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

from CoalCombustionForm import Ui_CoalCombustionForm

from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, DoubleValidator

import Pages.CoalThermoChemistry as CoalThermoChemistry
from Pages.Boundary import Boundary
from Pages.LocalizationModel import LocalizationModel
from Pages.CoalCombustionModel import CoalCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CoalCombustionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Delegate for diameter
#-------------------------------------------------------------------------------

class DiameterDelegate(QItemDelegate):
    def __init__(self, parent):
        super(DiameterDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor, min=0.)
        v.setExclusiveMin()
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Delegate for oxydant composition
#-------------------------------------------------------------------------------

class OxydantDelegate(QItemDelegate):
    def __init__(self, parent):
        super(OxydantDelegate, self).__init__(parent)
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
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            model.setData(index, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel for Coals
#-------------------------------------------------------------------------------

class StandardItemModelCoals(QStandardItemModel):
    def __init__(self):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Number of coal types")]
        self.setColumnCount(len(self.headers))
        self.dataCoals = []


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataCoals[index.row()])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, label):
        """
        Add a row in the table.
        """
        self.dataCoals.append("Coal " + str(label))
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataCoals[row]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataCoals[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def deleteAll(self):
        """
        Delete all the rows in the model
        """
        self.dataCoals = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# StandarItemModel for Coal Classes
#-------------------------------------------------------------------------------

class StandardItemModelClasses(QStandardItemModel):
    def __init__(self, coalThermoChModel):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Number of classes"),
                        self.tr("Initial diameter")]
        self.setColumnCount(len(self.headers))
        self.dataClasses = []
        self.coalNumber = None
        self.coalThermoChModel = coalThermoChModel


    def setCoalNumber(self, coalNumber):
        self.coalNumber = coalNumber


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataClasses[index.row()][index.column()])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
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
        if col == 1:
            newDiameter, ok = value.toDouble()
            self.dataClasses[row][col] = newDiameter

            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.updateInitDiameterClasses(row+1, newDiameter)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, num, diameter):
        """
        Add a row in the table.
        """
        label = "Class " + str(num)
        item = [label, diameter]
        self.dataClasses.append(item)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataClasses[row]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataClasses[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def deleteAll(self):
        """
        Delete all the rows in the model
        """
        self.dataClasses = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# StandarItemModel for Oxydant
#-------------------------------------------------------------------------------

class StandardItemModelOxydant(QStandardItemModel):
    def __init__(self, coalThermoChModel):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Oxydant\nnumber"),
                        self.tr("     O2      "),
                        self.tr("     N2      "),
                        self.tr("     H2O     "),
                        self.tr("     CO2     ")]
        self.setColumnCount(len(self.headers))
        self.dataClasses = []
        self.coalThermoChModel = coalThermoChModel


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataClasses[index.row()][index.column()])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() != 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
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
        v, ok = value.toDouble()
        self.dataClasses[row][col] = v

        oxy = self.coalThermoChModel.getOxydants().getOxydant(row+1)
        if col == 1:
            oxy.setO2(v)
        elif col == 2:
            oxy.setN2(v)
        elif col == 3:
            oxy.setH2O(v)
        elif col == 4:
            oxy.setCO2(v)
        self.coalThermoChModel.getOxydants().updateOxydant(row+1, oxy)
        self.coalThermoChModel.save()

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, num, oxy):
        """
        Add a row in the table.
        """
        label = str(num)
        item = [label, oxy.getO2(), oxy.getN2(), oxy.getH2O(),oxy.getCO2()]
        self.dataClasses.append(item)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataClasses[row]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataClasses[row]
        row = self.rowCount()
        self.setRowCount(row-1)
        self.coalThermoChModel.getOxydants().deleteOxydant(row+1)
        self.coalThermoChModel.save()

    def deleteAll(self):
        """
        Delete all the rows in the model
        """
        self.dataClasses = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class CoalCombustionView(QWidget, Ui_CoalCombustionForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_CoalCombustionForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.stbar = stbar

        self.model = CoalCombustionModel(self.case)
        self.coalThermoChModel = CoalThermoChemistry.CoalThermoChemistryModel("dp_FCP", self.case)

        self.coalNumber = 0
        ###
        ### Attention voir comment creer un charbon quand il n'y a pas de fichier dp_FCP !
        ####

#        self.selectedEntryNb1p1h = None

        # widgets layout.

        # Models
        # ------
        self.modelCoals = StandardItemModelCoals()
        self.treeViewCoals.setModel(self.modelCoals)

        self.modelClasses = StandardItemModelClasses(self.coalThermoChModel)
        self.treeViewClasses.setModel(self.modelClasses)
        self.treeViewClasses.resizeColumnToContents(0)
        self.treeViewClasses.resizeColumnToContents(1)

        self.modelOxydants = StandardItemModelOxydant(self.coalThermoChModel)
        self.tableViewOxydants.setModel(self.modelOxydants)
        self.tableViewOxydants.resizeColumnsToContents()
        self.tableViewOxydants.resizeRowsToContents()

        # set Coal number in modelClasses
        coalNumber = self.treeViewCoals.currentIndex().row()
        log.debug("coalNumber row = %i " % coalNumber)
        self.modelClasses.setCoalNumber(coalNumber)

        delegateDiameter = DiameterDelegate(self.treeViewClasses)
        self.treeViewClasses.setItemDelegateForColumn(1, delegateDiameter)
        delegateOxydant = OxydantDelegate(self.tableViewOxydants)
        self.tableViewOxydants.setItemDelegate(delegateOxydant)

        # Combo box
        # ---------
        self.modelPCI = ComboModel(self.comboBoxPCIList,2,1)
        self.modelPCI.addItem(self.tr("on dry"), "sec")
        self.modelPCI.addItem(self.tr("on pure"), "pur")

        self.modelReactTypeO2 = ComboModel(self.comboBoxReactO2,2,1)
        self.modelReactTypeO2.addItem(self.tr("0.5"), "0.5")
        self.modelReactTypeO2.addItem(self.tr("1"), "1")

        self.modelReactTypeCO2 = ComboModel(self.comboBoxReactCO2,2,1)
        self.modelReactTypeCO2.addItem(self.tr("0.5"), "0.5")
        self.modelReactTypeCO2.addItem(self.tr("1"), "1")

        # Connections
        # -----------
        self.connect(self.treeViewCoals, SIGNAL("clicked(const QModelIndex &)"), self.slotSelectCoal)
        self.connect(self.pushButtonAddCoal,     SIGNAL("clicked()"), self.slotCreateCoal)
        self.connect(self.pushButtonDeleteCoal,  SIGNAL("clicked()"), self.slotDeleteCoal)
        self.connect(self.pushButtonAddClass,    SIGNAL("clicked()"), self.slotCreateClass)
        self.connect(self.pushButtonDeleteClass, SIGNAL("clicked()"), self.slotDeleteClass)
        self.connect(self.pushButtonAddOxydant,    SIGNAL("clicked()"), self.slotCreateOxydant)
        self.connect(self.pushButtonDeleteOxydant, SIGNAL("clicked()"), self.slotDeleteOxydant)

        self.connect(self.lineEditC, SIGNAL("textChanged(const QString &)"), self.slotCComposition)
        self.connect(self.lineEditH, SIGNAL("textChanged(const QString &)"), self.slotHComposition)
        self.connect(self.lineEditO, SIGNAL("textChanged(const QString &)"), self.slotOComposition)
        self.connect(self.lineEditPCI,      SIGNAL("textChanged(const QString &)"), self.slotPCI)
        self.connect(self.lineEditCp,       SIGNAL("textChanged(const QString &)"), self.slotThermalCapacity)
        self.connect(self.lineEditDensity,  SIGNAL("textChanged(const QString &)"), self.slotDensity)
        self.connect(self.comboBoxPCIList,  SIGNAL("activated(const QString&)"), self.slotPCIType)

        self.connect(self.lineEditCokeC,   SIGNAL("textChanged(const QString &)"), self.slotCCompositionCoke)
        self.connect(self.lineEditCokeH,   SIGNAL("textChanged(const QString &)"), self.slotHCompositionCoke)
        self.connect(self.lineEditCokeO,   SIGNAL("textChanged(const QString &)"), self.slotOCompositionCoke)
        self.connect(self.lineEditCokePCI, SIGNAL("textChanged(const QString &)"), self.slotPCICoke)

        self.connect(self.lineEditAshesRatio,    SIGNAL("textChanged(const QString &)"), self.slotAshesRatio)
        self.connect(self.lineEditAshesEnthalpy, SIGNAL("textChanged(const QString &)"), self.slotAshesFormingEnthalpy)
        self.connect(self.lineEditAshesCp,       SIGNAL("textChanged(const QString &)"), self.slotAshesThermalCapacity)
        self.connect(self.lineEditHumidity, SIGNAL("textChanged(const QString &)"), self.slotHumidity)

        self.connect(self.lineEditCoefY1, SIGNAL("textChanged(const QString &)"), self.slotY1CH)
        self.connect(self.lineEditCoefY2, SIGNAL("textChanged(const QString &)"), self.slotY2CH)
        self.connect(self.lineEditCoefA1, SIGNAL("textChanged(const QString &)"), self.slotA1CH)
        self.connect(self.lineEditCoefA2, SIGNAL("textChanged(const QString &)"), self.slotA2CH)
        self.connect(self.lineEditCoefE1, SIGNAL("textChanged(const QString &)"), self.slotE1CH)
        self.connect(self.lineEditCoefE2, SIGNAL("textChanged(const QString &)"), self.slotE2CH)
        self.connect(self.checkBoxY1, SIGNAL("clicked()"), self.slotIY1)
        self.connect(self.checkBoxY2, SIGNAL("clicked()"), self.slotIY2)

        self.connect(self.lineEditConstO2,  SIGNAL("textChanged(const QString &)"), self.slotPreExpoCstO2)
        self.connect(self.lineEditEnergyO2, SIGNAL("textChanged(const QString &)"), self.slotActivEnergyO2)
        self.connect(self.comboBoxReactO2,  SIGNAL("activated(const QString&)"), self.slotReactTypeO2)
        self.connect(self.lineEditConstCO2,  SIGNAL("textChanged(const QString &)"), self.slotPreExpoCstCO2)
        self.connect(self.lineEditEnergyCO2, SIGNAL("textChanged(const QString &)"), self.slotActivEnergyCO2)
        self.connect(self.comboBoxReactCO2,  SIGNAL("activated(const QString&)"), self.slotReactTypeCO2)

        # Validators
        # ----------
        validatorC = DoubleValidator(self.lineEditC, min=0., max=100.)
        validatorH = DoubleValidator(self.lineEditH, min=0., max=100.)
        validatorO = DoubleValidator(self.lineEditO, min=0., max=100.)
        validatorPCI = DoubleValidator(self.lineEditPCI, min=0.)
        validatorCp = DoubleValidator(self.lineEditCp, min=0.)
        validatorDensity = DoubleValidator(self.lineEditDensity, min=0.)
        validatorHumidity = DoubleValidator(self.lineEditHumidity, min=0., max=100.)

        validatorCCoke = DoubleValidator(self.lineEditCokeC, min=0., max=100.)
        validatorHCoke = DoubleValidator(self.lineEditCokeH, min=0., max=100.)
        validatorOCoke = DoubleValidator(self.lineEditCokeO, min=0., max=100.)
        validatorPCICoke = DoubleValidator(self.lineEditCokePCI, min=0.)

        validatorAshesRatio = DoubleValidator(self.lineEditAshesRatio, min=0., max=100.)
        validatorAshesEnthalpy = DoubleValidator(self.lineEditAshesEnthalpy, min=0.)
        validatorAshesCp = DoubleValidator(self.lineEditAshesCp, min=0.)

        validatorY1 = DoubleValidator(self.lineEditCoefY1, min=0.)
        validatorY2 = DoubleValidator(self.lineEditCoefY2, min=0.)
        validatorA1 = DoubleValidator(self.lineEditCoefA1, min=0.)
        validatorA2 = DoubleValidator(self.lineEditCoefA2, min=0.)
        validatorE1 = DoubleValidator(self.lineEditCoefE1, min=0.)
        validatorE2 = DoubleValidator(self.lineEditCoefE2, min=0.)

        validatorConstO2 = DoubleValidator(self.lineEditConstO2, min=0.)
        validatorEnergyO2 = DoubleValidator(self.lineEditEnergyO2, min=0.)
        validatorConstCO2 = DoubleValidator(self.lineEditConstCO2, min=0.)
        validatorEnergyCO2 = DoubleValidator(self.lineEditEnergyCO2, min=0.)

        self.lineEditC.setValidator(validatorC)
        self.lineEditH.setValidator(validatorH)
        self.lineEditO.setValidator(validatorO)
        self.lineEditPCI.setValidator(validatorPCI)
        self.lineEditCp.setValidator(validatorCp)
        self.lineEditDensity.setValidator(validatorDensity)

        self.lineEditCokeC.setValidator(validatorCCoke)
        self.lineEditCokeH.setValidator(validatorHCoke)
        self.lineEditCokeO.setValidator(validatorOCoke)
        self.lineEditCokePCI.setValidator(validatorPCICoke)

        self.lineEditAshesRatio.setValidator(validatorAshesRatio)
        self.lineEditAshesEnthalpy.setValidator(validatorAshesEnthalpy)
        self.lineEditAshesCp.setValidator(validatorAshesCp)
        self.lineEditHumidity.setValidator(validatorHumidity)

        self.lineEditCoefY1.setValidator(validatorY1)
        self.lineEditCoefY2.setValidator(validatorY2)
        self.lineEditCoefA1.setValidator(validatorA1)
        self.lineEditCoefA2.setValidator(validatorA2)
        self.lineEditCoefE1.setValidator(validatorE1)
        self.lineEditCoefE2.setValidator(validatorE2)

        self.lineEditConstO2.setValidator(validatorConstO2)
        self.lineEditEnergyO2.setValidator(validatorEnergyO2)
        self.lineEditConstCO2.setValidator(validatorConstCO2)
        self.lineEditEnergyCO2.setValidator(validatorEnergyCO2)

        # Initialize widgets

        coals = self.coalThermoChModel.getCoals()
        CoalsNumber = coals.getNumber()
        for coalNumber in range(0, CoalsNumber):
            self.modelCoals.addItem(coalNumber + 1)

        if coals.getNumber() >= 3:
            self.pushButtonAddCoal.setDisabled(True)
        if coals.getNumber() <= 1:
            self.pushButtonDeleteCoal.setDisabled(True)

        self.slotSelectCoal()

        if self.model.getCoalCombustionModel() == 'coal_homo':
            self.lineEditHumidity.setText(QString(str(0.)))
            self.labelHumidity.setDisabled(True)
            self.lineEditHumidity.setDisabled(True)
            self.labelUnitHumidity.setDisabled(True)

        num = self.coalThermoChModel.getOxydants().getNumber()
        for index in range(0, num):
            oxy = self.coalThermoChModel.getOxydants().getOxydant(index+1)
            self.modelOxydants.addItem(index+1, oxy)

        # Update buttons
        self.pushButtonAddOxydant.setEnabled(True)
        self.pushButtonDeleteOxydant.setEnabled(True)
        if num <= 1:
            self.pushButtonDeleteOxydant.setDisabled(True)
        if num >= 3:
            self.pushButtonAddOxydant.setDisabled(True)


    @pyqtSignature("const QModelIndex &")
    def slotSelectCoal(self, text=None):
        """
        Display values for the current coal selected in the view.
        """
        #self.treeViewCoals.setCurrentIndex(index)
        row = self.treeViewCoals.currentIndex().row()
        log.debug("selectCoal row = %i "%row)
        if row == -1:
            row = 0

        coals = self.coalThermoChModel.getCoals()
        number = row + 1
        self.coalNumber = int(number)
        coal = coals.getCoal(self.coalNumber)

        # set Coal number in modelClasses
        self.modelClasses.setCoalNumber(self.coalNumber)

        # Classes
        ClassesNumber  = coal.getClassesNumber()
        initDiameters  = coal.getInitDiameterClasses()

        self.modelClasses.deleteAll()
        for number in range(0, ClassesNumber):
            self.modelClasses.addItem(number + 1, initDiameters[number])
#        self.coalThermoChModel.save()

        log.debug("selectCoal getClassesNumber = %i " % coal.getClassesNumber())
        self.pushButtonDeleteClass.setDisabled(False)
        if coal.getClassesNumber() >= 10:
            self.pushButtonAddClass.setDisabled(True)
        if coal.getClassesNumber() <= 1:
            self.pushButtonDeleteClass.setDisabled(True)

        # General (composition)
        self.lineEditC.setText(QString(str(coal.getCDryComposition())))
        self.lineEditH.setText(QString(str(coal.getHDryComposition())))
        self.lineEditO.setText(QString(str(coal.getODryComposition())))
        self.lineEditPCI.setText(QString(str(coal.getPCIValue())))
        self.lineEditCp.setText(QString(str(coal.getThermalCapacity())))
        self.lineEditDensity.setText(QString(str(coal.getDensity())))
        self.lineEditHumidity.setText(QString(str(coal.getHumidity())))
        dry = coal.getDry()
        if dry == 0:
            key = "pur"
        elif dry == 1:
            key = "sec"
        self.modelPCI.setItem(str_model=key)

        # Coke
        self.lineEditCokeC.setText(QString(str(coal.getCDryCompositionCoke())))
        self.lineEditCokeH.setText(QString(str(coal.getHDryCompositionCoke())))
        self.lineEditCokeO.setText(QString(str(coal.getODryCompositionCoke())))
        self.lineEditCokePCI.setText(QString(str(coal.getPCICokeValue())))

        # Ashes
        self.lineEditAshesRatio.setText(QString(str(coal.getAshesRatio())))
        self.lineEditAshesEnthalpy.setText(QString(str(coal.getAshesFormingEnthalpy())))
        self.lineEditAshesCp.setText(QString(str(coal.getAshesThermalCapacity())))

        # Devolatilisation
        stat = self.coalThermoChModel.getCoals().getCoal(self.coalNumber).getIY1CH()
        if stat == 0 :
            bstat = "on"
            self.checkBoxY1.setChecked(True)
        else:
            bstat = "off"
            self.checkBoxY1.setChecked(False)
        self.ActiveWidgetIY1(bstat)

        self.lineEditCoefY1.setText(QString(str(coal.getY1CH())))

        stat = self.coalThermoChModel.getCoals().getCoal(self.coalNumber).getIY2CH()
        if stat == 0 :
            bstat = "on"
            self.checkBoxY2.setChecked(True)
        else:
            bstat = "off"
            self.checkBoxY2.setChecked(False)
        self.ActiveWidgetIY2(bstat)

        self.lineEditCoefY2.setText(QString(str(coal.getY2CH())))

        self.lineEditCoefA1.setText(QString(str(coal.getA1CH())))
        self.lineEditCoefA2.setText(QString(str(coal.getA2CH())))
        self.lineEditCoefE1.setText(QString(str(coal.getE1CH())))
        self.lineEditCoefE2.setText(QString(str(coal.getE2CH())))

        # Combustion heterogene
        self.lineEditConstO2.setText(QString(str(coal.getAHETCH_O2())))
        self.lineEditEnergyO2.setText(QString(str(coal.getEHETCH_O2())))
        self.lineEditConstCO2.setText(QString(str(coal.getAHETCH_CO2())))
        self.lineEditEnergyCO2.setText(QString(str(coal.getEHETCH_CO2())))

        if coal.getIOCHET_O2() == 0:
            key = "0.5"
        else:
            key = "1"
        self.modelReactTypeO2.setItem(str_model=key)

        if coal.getIOCHET_CO2() == 0:
            key = "0.5"
        else:
            key = "1"
        self.modelReactTypeCO2.setItem(str_model=key)


    @pyqtSignature("")
    def slotCreateCoal(self):
        """ create a new coal"""
        diameter = self.model.defaultValues()['diameter']

        # update model
        coal = CoalThermoChemistry.Coal()
        self.coalThermoChModel.getCoals().addCoal(coal)
        CoalsNumber = self.coalThermoChModel.getCoals().getNumber()
        coal = self.coalThermoChModel.getCoals().getCoal(CoalsNumber)
        coal.updateInitDiameterClasses(1, diameter)

        # Init
        self.modelCoals.deleteAll()
#        CoalsNumber = self.coalThermoChModel.getCoals().getNumber()
        for number in range(0, CoalsNumber):
            self.modelCoals.addItem(number + 1)
            log.debug("slotCreateCoal number + 1 = %i " % (number + 1))
        self.coalThermoChModel.save()

        # update Properties and scalars
        self.model.createCoalModelScalarsAndProperties(self.coalThermoChModel)

        # update Buttons
        if self.coalThermoChModel.getCoals().getNumber() >= 3:
            self.pushButtonAddCoal.setDisabled(True)
        self.pushButtonDeleteCoal.setEnabled(True)

        # Coal created is selected
        self.slotSelectCoal()


    @pyqtSignature("")
    def slotDeleteCoal(self):
        """ cancel a coal"""
        row = self.treeViewCoals.currentIndex().row()
        log.debug("slotDeleteCoal row = %i "%row)

        if row == -1:
            return

        number = row + 1

        # update Scalars and properties
        self.model.deleteCoalModelScalarsAndProperties(self.coalThermoChModel, number - 1)

        # update boundary conditions (1/2)
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            if nature == "inlet":
                bc = Boundary("coal_inlet", label, self.case)
                bc.deleteCoalFlow(number-1, self.coalThermoChModel.getCoals().getNumber())

        self.coalThermoChModel.getCoals().deleteCoal(number)

        # Init
        self.modelCoals.deleteAll()
        CoalsNumber = self.coalThermoChModel.getCoals().getNumber()
        for number in range(0, CoalsNumber):
            self.modelCoals.addItem(number + 1)
            log.debug("slotDeleteCoal number + 1 = %i " % (number + 1))
        self.coalThermoChModel.save()

        # Update buttons
        self.pushButtonAddCoal.setEnabled(True)
        if self.coalThermoChModel.getCoals().getNumber() <= 1:
            self.pushButtonDeleteCoal.setDisabled(True)
        else:
            self.pushButtonDeleteCoal.setEnabled(True)

        # First coal is selected
        self.slotSelectCoal()


    @pyqtSignature("")
    def slotCreateClass(self):
        """Create a new class"""
        diameter = self.model.defaultValues()['diameter']

        coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
        coal.addInitDiameterClasses(diameter)
        self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)

        # Init
        ClassesNumber = coal.getClassesNumber()
        initDiameters = coal.getInitDiameterClasses()
        self.modelClasses.addItem(ClassesNumber, initDiameters[ClassesNumber -1])
        log.debug("slotCreateClass number + 1 = %i " % (ClassesNumber))
        self.coalThermoChModel.save()
        self.model.createClassModelScalarsAndProperties(self.coalThermoChModel, self.coalNumber)

        # FIXME: bug ici
        # Update boundary conditions
        log.debug("slotCreateClass: number of classes: %i " % coal.getClassesNumber())
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            if zone.getNature() == "inlet":
                b = Boundary("coal_inlet", zone.getLabel(), self.case)
                #b.getCoalRatios(self.coalNumber-1)
                b.updateCoalRatios(self.coalNumber-1)

        # Update buttons
        self.pushButtonDeleteClass.setEnabled(True)
        if self.coalThermoChModel.getCoals().getCoal(self.coalNumber).getClassesNumber() >= 10:
            self.pushButtonAddClass.setDisabled(True)


    @pyqtSignature("")
    def slotDeleteClass(self):
        """ cancel a class diameter"""
        row = self.treeViewClasses.currentIndex().row()
        log.debug("slotDeleteClass  number = %i " % row)
        if row == -1:
            return

        number = row + 1

        self.model.deleteClassModelProperties(self.coalThermoChModel, self.coalNumber, number - 1)
        self.model.deleteClassModelScalars(self.coalThermoChModel, self.coalNumber, number - 1)

        coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
        coal.cancelInitDiameterClasses(number)
        self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()

        # Init
        self.modelClasses.deleteAll()
        ClassesNumber = coal.getClassesNumber()
        initDiameters = coal.getInitDiameterClasses()
        for number in range(0, ClassesNumber):
            self.modelClasses.addItem(number+1, initDiameters[number])
            log.debug("slotDeleteClass number + 1 = %i " % (number+1))

        # Update boundary conditions
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            if zone.getNature() == "inlet":
                bc = Boundary("coal_inlet", zone.getLabel(), self.case)
                bc.updateCoalRatios(self.coalNumber-1)

        # Update buttons
        self.pushButtonAddClass.setEnabled(True)
        if self.coalThermoChModel.getCoals().getCoal(self.coalNumber).getClassesNumber() <= 1:
            self.pushButtonDeleteClass.setDisabled(True)


    @pyqtSignature("")
    def slotCreateOxydant(self):
        """Create a new oxydant"""
        new_ox = CoalThermoChemistry.Oxydant()
        self.coalThermoChModel.getOxydants().addOxydant(new_ox)
        num = self.coalThermoChModel.getOxydants().getNumber()

        self.modelOxydants.addItem(num, new_ox)
        log.debug("slotCreateOxydant number = %i " % num)
        self.coalThermoChModel.save()
        #self.model.createOxydantModelScalarsAndProperties(self.coalThermoChModel, num)

        # Update buttons
        self.pushButtonDeleteOxydant.setEnabled(True)
        if self.coalThermoChModel.getOxydants().getNumber() >= 3:
            self.pushButtonAddOxydant.setDisabled(True)


    @pyqtSignature("")
    def slotDeleteOxydant(self):
        """ delete an oxydant"""
        row = self.tableViewOxydants.currentIndex().row()
        log.debug("slotDeleteOxydants number = %i " % row)
        if row == -1:
            return

        number = row + 1
        self.coalThermoChModel.getOxydants().deleteOxydant(number)
        self.coalThermoChModel.save()

        # Update boundary conditions
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            if nature == "inlet":
                bc = Boundary("coal_inlet", label, self.case)
                oxy_max = bc.getOxydantNumber()
                if oxy_max >= number:
                    bc.setOxydantNumber(oxy_max-1)

        self.modelOxydants.deleteAll()
        for number in range(0, self.coalThermoChModel.getOxydants().getNumber()):
            oxy = self.coalThermoChModel.getOxydants().getOxydant(number+1)
            self.modelOxydants.addItem(number+1, oxy)

        # Update buttons
        self.pushButtonAddOxydant.setEnabled(True)
        if self.coalThermoChModel.getOxydants().getNumber() <= 1:
            self.pushButtonDeleteOxydant.setDisabled(True)


    @pyqtSignature("const QString&")
    def slotCComposition(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setCDryComposition(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotHComposition(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setHDryComposition(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotOComposition(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setODryComposition(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotPCI(self, text):
        PCI, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setPCIValue(PCI)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotThermalCapacity(self, text):
        Cp, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setThermalCapacity(Cp)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotDensity(self, text):
        density, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setDensity(density)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotHumidity(self, text):
        humidity, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setHumidity(humidity)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotPCIType(self, text):
        """
        Return the selected type of boundary condition. Whatever the
        language (i18n), type is always in ('wall','inlet','outlet','symmetry').
        """
        key = self.modelPCI.dicoV2M[str(text)]
        coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
        if key == "sec" :
            coal.setDry(1)
        elif key == "pur":
            coal.setDry(0)
        self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()
        return key


    @pyqtSignature("const QString&")
    def slotCCompositionCoke(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setCDryCompositionCoke(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotHCompositionCoke(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setHDryCompositionCoke(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotOCompositionCoke(self, text):
        composition, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setODryCompositionCoke(composition)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotPCICoke(self, text):
        PCICoke, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setPCICokeValue(PCICoke)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotAshesRatio(self, text):
        ashesRatio, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setAshesRatio(ashesRatio)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def slotAshesFormingEnthalpy(self, text):
        ashesFormingEnthalpy, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setAshesFormingEnthalpy(ashesFormingEnthalpy)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotAshesThermalCapacity(self, text):
        ashesThermalCapacity, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setAshesThermalCapacity(ashesThermalCapacity)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotY1CH(self, text):
        Y1CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setY1CH(Y1CH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotY2CH(self, text):
        Y2CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setY2CH(Y2CH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotA1CH(self, text):
        A1CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setA1CH(A1CH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotA2CH(self, text):
        A2CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
                coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
                coal.setA2CH(A2CH)
                self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
                self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotE1CH(self, text):
        E1CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setE1CH(E1CH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotE2CH(self, text):
        E2CH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setE2CH(E2CH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("")
    def slotIY1(self):
        if self.checkBoxY1.isChecked():
            stat = "on"
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setIY1CH(0)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        else :
            stat = "off"
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setIY1CH(1)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()
        self.ActiveWidgetIY1(stat)


    @pyqtSignature("")
    def slotIY2(self):
        if self.checkBoxY2.isChecked():
            stat = "on"
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setIY2CH(0)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        else :
            stat = "off"
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setIY2CH(1)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()
        self.ActiveWidgetIY2(stat)


    def ActiveWidgetIY1(self, stat) :
        if stat == 'on':
            self.labelCoefY1.setDisabled(True)
            self.lineEditCoefY1.setDisabled(True)
        else :
            self.labelCoefY1.setEnabled(True)
            self.lineEditCoefY1.setEnabled(True)


    def ActiveWidgetIY2(self, stat) :
        if stat == 'on':
            self.labelCoefY2.setDisabled(True)
            self.lineEditCoefY2.setDisabled(True)
        else :
            self.labelCoefY2.setEnabled(True)
            self.lineEditCoefY2.setEnabled(True)


    @pyqtSignature("const QString&")
    def slotPreExpoCstO2(self, text):
        AHETCH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setAHETCH_O2(AHETCH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotActivEnergyO2(self, text):
        EHETCH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setEHETCH_O2(EHETCH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotReactTypeO2(self, text):
        key = self.modelReactTypeO2.dicoV2M[str(text)]
        coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
        if key == "0.5" :
            coal.setIOCHET_O2(0)
        else:
            coal.setIOCHET_O2(1)
        self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()
        return key


    @pyqtSignature("const QString&")
    def slotPreExpoCstCO2(self, text):
        AHETCH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setAHETCH_CO2(AHETCH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotActivEnergyCO2(self, text):
        EHETCH, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
            coal.setEHETCH_CO2(EHETCH)
            self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
            self.coalThermoChModel.save()


    @pyqtSignature("const QString&")
    def slotReactTypeCO2(self, text):
        key = self.modelReactTypeCO2.dicoV2M[str(text)]
        coal = self.coalThermoChModel.getCoals().getCoal(self.coalNumber)
        if key == "0.5" :
            coal.setIOCHET_CO2(0)
        else:
            coal.setIOCHET_CO2(1)
        self.coalThermoChModel.getCoals().updateCoal(self.coalNumber, coal)
        self.coalThermoChModel.save()
        return key


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
