# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
- StandardItemModelOxidant
- StandardItemModelRefusal
- CoalCombustionView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

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

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, RegExpValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string

from code_saturne.Pages.CoalCombustionForm import Ui_CoalCombustionForm
from code_saturne.model.Boundary import Boundary
from code_saturne.model.CoalCombustionModel import CoalCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CoalCombustionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the label fel
#-------------------------------------------------------------------------------

class LabelFuelDelegate(QItemDelegate):
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
        rx = "[_A-Za-z0-9 \(\)]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_plabel = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_plabel = str(editor.text())

            if new_plabel in model.mdl.getLabelIdList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getLabelIdList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, str(new_plabel), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Combo box delegate for the fuel type
#-------------------------------------------------------------------------------

class TypeFuelDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(TypeFuelDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem("biomass")
        editor.addItem("coal")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"biomass": 0, "coal": 1}
        row = index.row()
        string = index.model().dataCoals[row]['type']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


    def paint(self, painter, option, index):
        row = index.row()
        fueltype = index.model().dataCoals[row]['type']
        isValid = fueltype != None and fueltype != ''

        if isValid:
            QItemDelegate.paint(self, painter, option, index)
        else:
            painter.save()
            # set background color
            if option.state & QStyle.State_Selected:
                painter.setBrush(QBrush(Qt.darkRed))
            else:
                painter.setBrush(QBrush(Qt.red))
            # set text color
            painter.setPen(QPen(Qt.NoPen))
            painter.drawRect(option.rect)
            painter.setPen(QPen(Qt.black))
            value = index.data(Qt.DisplayRole)
            if value.isValid():
                text = from_qvariant(value, to_text_string)
                painter.drawText(option.rect, Qt.AlignLeft, text)
            painter.restore()

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
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Delegate for refusal
#-------------------------------------------------------------------------------

class RefusalDelegate(QItemDelegate):
    def __init__(self, parent):
        super(RefusalDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor, min=0.)
        v.setExclusiveMin()
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Delegate for oxidant composition
#-------------------------------------------------------------------------------

class OxidantDelegate(QItemDelegate):
    def __init__(self, parent):
        super(OxidantDelegate, self).__init__(parent)
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
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel for Coals
#-------------------------------------------------------------------------------

class StandardItemModelCoals(QStandardItemModel):
    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Name"), self.tr("Type")]
        self.setColumnCount(len(self.headers))
        self.mdl = mdl
        self.dataCoals = []
        self.defaultItem = []
        self.populateModel()


    def populateModel(self):
        self.dicoV2M= {"biomass": 'biomass',
                       "coal" : 'coal'}
        self.dicoM2V= {"biomass" : 'biomass',
                       "coal" : 'coal'}
        for id in self.mdl.getFuelIdList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico  = {}
            dico['name'] = self.mdl.getFuelLabel(id)
            dico['type'] = self.mdl.getFuelType(id)

            self.dataCoals.append(dico)
            if int(id) < 0:
                self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            dico = self.dataCoals[row]

            if index.column() == 0:
                return dico['name']
            elif index.column() == 1:
                return self.dicoM2V[dico['type']]
            else:
                return None

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            old_plabel = self.dataCoals[row]['name']
            new_plabel = str(from_qvariant(value, to_text_string))
            self.dataCoals[row]['name'] = new_plabel
            self.mdl.setFuelLabel(row + 1, new_plabel)

        elif col == 1:
            self.dataCoals[row]['type'] = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            self.mdl.setFuelType(row + 1, self.dataCoals[row]['type'])

        self.dataChanged.emit(index, index)
        return True


    def addItem(self, name = None, fuel_type = None):
        """
        Add a row in the table.
        """
        dico = {}
        if (name != None and fuel_type != None):
            dico['name'] = name
            dico['type'] = fuel_type
        else:
            self.mdl.createCoal()
            number = self.mdl.getCoalNumber()
            dico['name'] = self.mdl.getFuelLabel(number)
            dico['type'] = self.mdl.getFuelType(number)
        self.dataCoals.append(dico)

        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        """
        Returns the name of the fuel file.
        """
        return self.dataCoals[row]


    def deleteItem(self, row):
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
    def __init__(self, model, fuel):
        """
        """
        QStandardItemModel.__init__(self)
        self.model = model
        self.fuel = fuel
        diameter_type = self.model.getDiameterType(self.fuel)

        if diameter_type == 'automatic' :
            self.headers = [self.tr("class number"),
                            self.tr("Initial diameter (m)")]
        elif diameter_type == 'rosin-rammler_law':
            self.headers = [self.tr("class number"),
                            self.tr("Mass percent")]

        self.setColumnCount(len(self.headers))
        self.dataClasses = []


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            return self.dataClasses[index.row()][index.column()]
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 1:
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
        ClassId = row + 1

        if col == 1:
            newDiameter = from_qvariant(value, float)
            self.dataClasses[row][col] = newDiameter

            diameter_type = self.model.getDiameterType(self.fuel)
            if diameter_type == 'automatic' :
                self.model.setDiameter(self.fuel, ClassId, newDiameter)
            elif diameter_type == 'rosin-rammler_law':
                self.model.setMassPercent(self.fuel, ClassId, newDiameter)

        self.dataChanged.emit(index, index)
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
# StandarItemModel for Oxidant
#-------------------------------------------------------------------------------

class StandardItemModelOxidant(QStandardItemModel):
    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Oxidant\nnumber"),
                        self.tr("     O2      "),
                        self.tr("     N2      "),
                        self.tr("     H2O     "),
                        self.tr("     CO2     ")]
        self.setColumnCount(len(self.headers))
        self.dataClasses = []
        self.model = model


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            return self.dataClasses[index.row()][index.column()]
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() != 0:
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
        v = from_qvariant(value, float)
        self.dataClasses[row][col] = v
        oxId = row + 1

        if col == 1:
            self.model.setElementComposition(oxId, "O2", v)
        elif col == 2:
            self.model.setElementComposition(oxId, "N2", v)
        elif col == 3:
            self.model.setElementComposition(oxId, "H2O", v)
        elif col == 4:
            self.model.setElementComposition(oxId, "CO2", v)

        self.dataChanged.emit(index, index)
        return True


    def addItem(self, num):
        """
        Add a row in the table.
        """
        label = str(num)
        O2  = self.model.getElementComposition(num, "O2")
        N2  = self.model.getElementComposition(num, "N2")
        H2O = self.model.getElementComposition(num, "H2O")
        CO2 = self.model.getElementComposition(num, "CO2")
        item = [label, O2, N2, H2O, CO2]
        self.dataClasses.append(item)
        row = self.rowCount()
        self.setRowCount(row + 1)


    def getItem(self, row):
        return self.dataClasses[row]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataClasses[row]
        row = self.rowCount()
        self.setRowCount(row-1)
        self.model.deleteOxidant(row+1)

    def deleteAll(self):
        """
        Delete all the rows in the model
        """
        self.dataClasses = []
        self.setRowCount(0)


#-------------------------------------------------------------------------------
# StandarItemModel for Refusal
#-------------------------------------------------------------------------------

class StandardItemModelRefusal(QStandardItemModel):
    def __init__(self, model, fuel):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Refusal"),
                        self.tr("diameter (m)"),
                        self.tr("value")]
        self.setColumnCount(len(self.headers))
        self.dataClasses = []
        self.model = model
        self.fuel = fuel


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            return self.dataClasses[index.row()][index.column()]
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() != 0:
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
        v = from_qvariant(value, float)
        self.dataClasses[row][col] = v

        if col == 1:
            self.model.setRefusalDiameter(self.fuel, row+1, v)
        elif col == 2:
            self.model.setRefusalValue(self.fuel, row+1, v)

        self.dataChanged.emit(index, index)
        return True


    def addItem(self, num, item):
        """
        Add a row in the table.
        """
        label = str(num)
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
        self.model.deleteRefusal(self.fuel, row+1)

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
        self.case.undoStopGlobal()

        self.stbar = stbar

        self.model = CoalCombustionModel(self.case)

        # widgets layout.
        self.fuel = 1

        # Models
        # ------
        self.modelCoals = StandardItemModelCoals(self.model)
        self.treeViewCoals.setModel(self.modelCoals)
        delegate_label_fuel = LabelFuelDelegate(self.treeViewCoals)
        delegate_type       = TypeFuelDelegate(self.treeViewCoals, self.model)
        self.treeViewCoals.setItemDelegateForColumn(0, delegate_label_fuel)
        self.treeViewCoals.setItemDelegateForColumn(1, delegate_type)

        self.modelClasses = StandardItemModelClasses(self.model, self.fuel)
        self.treeViewClasses.setModel(self.modelClasses)
        self.treeViewClasses.resizeColumnToContents(0)
        self.treeViewClasses.resizeColumnToContents(1)

        self.modelRefusal = StandardItemModelRefusal(self.model, self.fuel)
        self.treeViewRefusal.setModel(self.modelRefusal)
        self.treeViewRefusal.resizeColumnToContents(0)
        self.treeViewRefusal.resizeColumnToContents(1)
        self.treeViewRefusal.resizeColumnToContents(2)

        self.modelOxidants = StandardItemModelOxidant(self.model)
        self.tableViewOxidants.setModel(self.modelOxidants)
        self.tableViewOxidants.resizeColumnsToContents()
        self.tableViewOxidants.resizeRowsToContents()

        delegateDiameter = DiameterDelegate(self.treeViewClasses)
        self.treeViewClasses.setItemDelegateForColumn(1, delegateDiameter)
        delegateRefusal = RefusalDelegate(self.treeViewRefusal)
        self.treeViewRefusal.setItemDelegate(delegateRefusal)
        delegateOxidant = OxidantDelegate(self.tableViewOxidants)
        self.tableViewOxidants.setItemDelegate(delegateOxidant)

        # Combo box
        # ---------

        self.modelKineticModel = ComboModel(self.comboBoxKineticModel,2,1)
        self.modelKineticModel.addItem(self.tr("unused"),       "unused")
        self.modelKineticModel.addItem(self.tr("CO2 mass fraction transport"),
                                       "co2_ym_transport")

        self.modelPCI = ComboModel(self.comboBoxPCIList,3,1)
        self.modelPCI.addItem(self.tr("LHV"), "LHV")
        self.modelPCI.addItem(self.tr("HHV"), "HHV")
        self.modelPCI.addItem(self.tr("IGT correlation"), "IGT_correlation")

        self.modelPCIType = ComboModel(self.comboBoxPCIType,3,1)
        self.modelPCIType.addItem(self.tr("dry basis"),    "dry_basis")
        self.modelPCIType.addItem(self.tr("dry ash free"), "dry_ash_free")
        self.modelPCIType.addItem(self.tr("as received"),  "as_received")

        self.modelY1Y2 = ComboModel(self.comboBoxY1Y2,3,1)
        self.modelY1Y2.addItem(self.tr("user define"),       "user_define")
        self.modelY1Y2.addItem(self.tr("automatic CHONS"),   "automatic_CHONS")
        self.modelY1Y2.addItem(self.tr("automatic formula"), "automatic_formula")

        self.modelDiameter = ComboModel(self.comboBoxDiameter,2,1)
        self.modelDiameter.addItem(self.tr("user define"),       "automatic")
        self.modelDiameter.addItem(self.tr("Rosin-Rammler law"), "rosin-rammler_law")

        self.modelReactTypeO2 = ComboModel(self.comboBoxReactO2,2,1)
        self.modelReactTypeO2.addItem(self.tr("0.5"), "0.5")
        self.modelReactTypeO2.addItem(self.tr("1"),   "1")

        self.modelReactTypeCO2 = ComboModel(self.comboBoxReactCO2,2,1)
        self.modelReactTypeCO2.addItem(self.tr("0.5"), "0.5")
        self.modelReactTypeCO2.addItem(self.tr("1"),   "1")

        self.modelReactTypeH2O = ComboModel(self.comboBoxReactH2O,2,1)
        self.modelReactTypeH2O.addItem(self.tr("0.5"), "0.5")
        self.modelReactTypeH2O.addItem(self.tr("1"),   "1")

        self.modelOxidantType = ComboModel(self.comboBoxOxidant,2,1)
        self.modelOxidantType.addItem(self.tr("volumic percentage"), "volumic_percent")
        self.modelOxidantType.addItem(self.tr("molar"),              "molar")

        self.modelReburning = ComboModel(self.comboBoxReburning,3,1)
        self.modelReburning.addItem(self.tr("unused"), "unused")
        self.modelReburning.addItem(self.tr("Model of Chen et al."), "chen")
        self.modelReburning.addItem(self.tr("Model of Dimitriou et al."),
                                    "dimitriou")

        # Connections
        # -----------
        self.comboBoxKineticModel.activated[str].connect(self.slotKineticModel)
        self.treeViewCoals.clicked[QModelIndex].connect(self.slotSelectCoal)
        self.pushButtonAddCoal.clicked.connect(self.slotCreateCoal)
        self.pushButtonDeleteCoal.clicked.connect(self.slotDeleteCoal)
        self.pushButtonAddClass.clicked.connect(self.slotCreateClass)
        self.pushButtonDeleteClass.clicked.connect(self.slotDeleteClass)
        self.pushButtonAddRefusal.clicked.connect(self.slotCreateRefusal)
        self.pushButtonDeleteRefusal.clicked.connect(self.slotDeleteRefusal)
        self.comboBoxDiameter.activated[str].connect(self.slotDiameterType)
        self.pushButtonAddOxidant.clicked.connect(self.slotCreateOxidant)
        self.pushButtonDeleteOxidant.clicked.connect(self.slotDeleteOxidant)

        self.lineEditC.textChanged[str].connect(self.slotCComposition)
        self.lineEditH.textChanged[str].connect(self.slotHComposition)
        self.lineEditO.textChanged[str].connect(self.slotOComposition)
        self.lineEditN.textChanged[str].connect(self.slotNComposition)
        self.lineEditS.textChanged[str].connect(self.slotSComposition)
        self.lineEditPCI.textChanged[str].connect(self.slotPCI)
        self.lineEditVolatileMatter.textChanged[str].connect(self.slotVolatileMatter)
        self.lineEditMoisture.textChanged[str].connect(self.slotMoisture)
        self.lineEditCp.textChanged[str].connect(self.slotThermalCapacity)
        self.lineEditThermalCond.textChanged[str].connect(self.slotThermalConductivity)
        self.lineEditDensity.textChanged[str].connect(self.slotDensity)
        self.comboBoxPCIList.activated[str].connect(self.slotPCIChoice)
        self.comboBoxPCIType.activated[str].connect(self.slotPCIType)
        self.lineEditCCoke.textChanged[str].connect(self.slotCCompositionCoke)
        self.lineEditHCoke.textChanged[str].connect(self.slotHCompositionCoke)
        self.lineEditOCoke.textChanged[str].connect(self.slotOCompositionCoke)
        self.lineEditNCoke.textChanged[str].connect(self.slotNCompositionCoke)
        self.lineEditSCoke.textChanged[str].connect(self.slotSCompositionCoke)

        self.lineEditAshesRatio.textChanged[str].connect(self.slotAshesRatio)
        self.lineEditAshesEnthalpy.textChanged[str].connect(self.slotAshesFormingEnthalpy)
        self.lineEditAshesCp.textChanged[str].connect(self.slotAshesThermalCapacity)

        self.comboBoxY1Y2.activated[str].connect(self.slotY1Y2)
        self.lineEditCoefY1.textChanged[str].connect(self.slotY1CH)
        self.lineEditCoefY2.textChanged[str].connect(self.slotY2CH)
        self.lineEditCoefA1.textChanged[str].connect(self.slotA1CH)
        self.lineEditCoefA2.textChanged[str].connect(self.slotA2CH)
        self.lineEditCoefE1.textChanged[str].connect(self.slotE1CH)
        self.lineEditCoefE2.textChanged[str].connect(self.slotE2CH)

        self.lineEditConstO2.textChanged[str].connect(self.slotPreExpoCstO2)
        self.lineEditEnergyO2.textChanged[str].connect(self.slotActivEnergyO2)
        self.comboBoxReactO2.activated[str].connect(self.slotReactTypeO2)
        self.lineEditConstCO2.textChanged[str].connect(self.slotPreExpoCstCO2)
        self.lineEditEnergyCO2.textChanged[str].connect(self.slotActivEnergyCO2)
        self.comboBoxReactCO2.activated[str].connect(self.slotReactTypeCO2)
        self.lineEditConstH2O.textChanged[str].connect(self.slotPreExpoCstH2O)
        self.lineEditEnergyH2O.textChanged[str].connect(self.slotActivEnergyH2O)
        self.comboBoxReactH2O.activated[str].connect(self.slotReactTypeH2O)
        self.comboBoxOxidant.activated[str].connect(self.slotOxidantType)
        self.comboBoxReburning.activated[str].connect(self.slotReburning)

        self.lineEditQPR.textChanged[str].connect(self.slotQPR)
        self.lineEditNitrogenConcentration.textChanged[str].connect(self.slotNitrogenConcentration)
        self.lineEditKobayashi1.textChanged[str].connect(self.slotKobayashi1)
        self.lineEditKobayashi2.textChanged[str].connect(self.slotKobayashi2)
        self.lineEditNitrogenLowTemp.textChanged[str].connect(self.slotNLowTemp)
        self.lineEditNitrogenHighTemp.textChanged[str].connect(self.slotNHighTemp)
        self.lineEditHCNChar.textChanged[str].connect(self.slotHCNChar)

        self.checkBoxNOxFormation.clicked[bool].connect(self.slotNOxFormation)
        self.checkBoxNOxFormationFeature.clicked[bool].connect(self.slotNOxFeature)
        self.checkBoxCO2Kinetics.clicked[bool].connect(self.slotCO2Kinetics)
        self.checkBoxH2OKinetics.clicked[bool].connect(self.slotH2OKinetics)

        self.tabWidget.currentChanged[int].connect(self.slotchanged)

        # Validators
        # ----------
        validatorC   = DoubleValidator(self.lineEditC, min=0., max=100.)
        validatorH   = DoubleValidator(self.lineEditH, min=0., max=100.)
        validatorO   = DoubleValidator(self.lineEditO, min=0., max=100.)
        validatorN   = DoubleValidator(self.lineEditN, min=0., max=100.)
        validatorS   = DoubleValidator(self.lineEditS, min=0., max=100.)
        validatorPCI = DoubleValidator(self.lineEditPCI, min=0.)
        validatorCp  = DoubleValidator(self.lineEditCp, min=0.)
        validatorla  = DoubleValidator(self.lineEditThermalCond, min=0.)
        validatorDensity = DoubleValidator(self.lineEditDensity, min=0.)
        validatorMoisture = DoubleValidator(self.lineEditMoisture, min=0., max=100.)
        validatorVolatileMatter = DoubleValidator(self.lineEditVolatileMatter, min=0., max=100.)
        validatorCCoke = DoubleValidator(self.lineEditCCoke, min=0., max=100.)
        validatorHCoke = DoubleValidator(self.lineEditHCoke, min=0., max=100.)
        validatorOCoke = DoubleValidator(self.lineEditOCoke, min=0., max=100.)
        validatorNCoke = DoubleValidator(self.lineEditNCoke, min=0., max=100.)
        validatorSCoke = DoubleValidator(self.lineEditSCoke, min=0., max=100.)

        validatorAshesRatio = DoubleValidator(self.lineEditAshesRatio, min=0., max=100.)
        validatorAshesEnthalpy = DoubleValidator(self.lineEditAshesEnthalpy, min=0.)
        validatorAshesCp = DoubleValidator(self.lineEditAshesCp, min=0.)

        validatorY1 = DoubleValidator(self.lineEditCoefY1, min=0.)
        validatorY2 = DoubleValidator(self.lineEditCoefY2, min=0.)
        validatorA1 = DoubleValidator(self.lineEditCoefA1, min=0.)
        validatorA2 = DoubleValidator(self.lineEditCoefA2, min=0.)
        validatorE1 = DoubleValidator(self.lineEditCoefE1, min=0.)
        validatorE2 = DoubleValidator(self.lineEditCoefE2, min=0.)

        validatorConstO2   = DoubleValidator(self.lineEditConstO2, min=0.)
        validatorEnergyO2  = DoubleValidator(self.lineEditEnergyO2, min=0.)
        validatorConstCO2  = DoubleValidator(self.lineEditConstCO2, min=0.)
        validatorEnergyCO2 = DoubleValidator(self.lineEditEnergyCO2, min=0.)
        validatorConstH2O  = DoubleValidator(self.lineEditConstH2O, min=0.)
        validatorEnergyH2O = DoubleValidator(self.lineEditEnergyH2O, min=0.)

        validatorQPR = DoubleValidator(self.lineEditQPR, min=0.)
        validatorNitrogenConcentration = DoubleValidator(self.lineEditNitrogenConcentration, min=0.)
        validatorKobayashi1 = DoubleValidator(self.lineEditKobayashi1, min=0., max = 1.)
        validatorKobayashi2 = DoubleValidator(self.lineEditKobayashi2, min=0., max = 1.)
        validatorNitrogenLowTemp = DoubleValidator(self.lineEditNitrogenLowTemp, min=0.)
        validatorNitrogenHighTemp = DoubleValidator(self.lineEditNitrogenHighTemp, min=0.)
        validatorHCNChar = DoubleValidator(self.lineEditHCNChar, min=0.)

        self.lineEditC.setValidator(validatorC)
        self.lineEditH.setValidator(validatorH)
        self.lineEditO.setValidator(validatorO)
        self.lineEditN.setValidator(validatorN)
        self.lineEditS.setValidator(validatorS)
        self.lineEditPCI.setValidator(validatorPCI)
        self.lineEditCp.setValidator(validatorCp)
        self.lineEditThermalCond.setValidator(validatorla)
        self.lineEditDensity.setValidator(validatorDensity)
        self.lineEditCCoke.setValidator(validatorCCoke)
        self.lineEditHCoke.setValidator(validatorHCoke)
        self.lineEditOCoke.setValidator(validatorOCoke)
        self.lineEditNCoke.setValidator(validatorNCoke)
        self.lineEditSCoke.setValidator(validatorSCoke)

        self.lineEditAshesRatio.setValidator(validatorAshesRatio)
        self.lineEditAshesEnthalpy.setValidator(validatorAshesEnthalpy)
        self.lineEditAshesCp.setValidator(validatorAshesCp)
        self.lineEditMoisture.setValidator(validatorMoisture)
        self.lineEditVolatileMatter.setValidator(validatorVolatileMatter)

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
        self.lineEditConstH2O.setValidator(validatorConstH2O)
        self.lineEditEnergyH2O.setValidator(validatorEnergyH2O)

        self.lineEditQPR.setValidator(validatorQPR)
        self.lineEditNitrogenConcentration.setValidator(validatorNitrogenConcentration)
        self.lineEditKobayashi1.setValidator(validatorKobayashi1)
        self.lineEditKobayashi2.setValidator(validatorKobayashi2)
        self.lineEditNitrogenLowTemp.setValidator(validatorNitrogenLowTemp)
        self.lineEditNitrogenHighTemp.setValidator(validatorNitrogenHighTemp)
        self.lineEditHCNChar.setValidator(validatorHCNChar)

        # Initialize widgets
        self.initializeView()

        num = self.model.getOxidantNumber()
        for index in range(0, num):
            self.modelOxidants.addItem(index + 1)

        # Update buttons
        self._updateCoalButton()
        self._updateOxidantButton()

        self.tabWidget.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


    def _updateCoalButton(self):
        """
        control solid fuel number between 1 and 3
        """
        CoalNumber = self.model.getCoalNumber()
        self.pushButtonDeleteCoal.setEnabled(True)
        self.pushButtonAddCoal.setEnabled(True)
        if CoalNumber >= 5:
            self.pushButtonAddCoal.setDisabled(True)
        elif CoalNumber <= 1:
            self.pushButtonDeleteCoal.setDisabled(True)


    def _updateClassButton(self):
        """
        control class number between 1 and 10 for a define solid fuel
        """
        ClassNumber = self.model.getClassNumber(self.fuel)

        self.pushButtonDeleteClass.setEnabled(True)
        self.pushButtonAddClass.setEnabled(True)
        if ClassNumber >= 10:
            self.pushButtonAddClass.setDisabled(True)
        elif ClassNumber <= 1:
            self.pushButtonDeleteClass.setDisabled(True)

        diameter_type = self.model.getDiameterType(self.fuel)

        if diameter_type == 'rosin-rammler_law':
            self._updateRefusalButton()


    def _updateRefusalButton(self):
        """
        control refusal number between 1 and number of class for a define solid fuel
        """
        ClassNumber   = self.model.getClassNumber(self.fuel)
        RefusalNumber = self.model.getRefusalNumber(self.fuel)

        self.pushButtonDeleteRefusal.setEnabled(True)
        self.pushButtonAddRefusal.setEnabled(True)
        if RefusalNumber >= ClassNumber:
            self.pushButtonAddRefusal.setDisabled(True)
        elif RefusalNumber <= 1:
            self.pushButtonDeleteRefusal.setDisabled(True)


    def _updateOxidantButton(self):
        """
        control oxidant number between 1 and 3
        """
        OxidantNumber = self.model.getOxidantNumber()

        self.pushButtonAddOxidant.setEnabled(True)
        self.pushButtonDeleteOxidant.setEnabled(True)
        if OxidantNumber >= 3:
            self.pushButtonAddOxidant.setDisabled(True)
        elif OxidantNumber <= 1:
            self.pushButtonDeleteOxidant.setDisabled(True)


    def initializeKineticModel(self):
        """
        initialize view with kinetic model type choice
        """
        key = self.model.getKineticModel()
        self.modelKineticModel.setItem(str_model=key)


    def initializeDiameter(self):
        """
        initialize view with diameter type choice
        """
        self.modelClasses.deleteAll()

        key = self.model.getDiameterType(self.fuel)
        self.modelDiameter.setItem(str_model=key)

        ClassesNumber = self.model.getClassNumber(self.fuel)

        if key == 'automatic':
            self.treeViewRefusal.hide()
            self.pushButtonDeleteRefusal.hide()
            self.pushButtonAddRefusal.hide()

            for number in range(0, ClassesNumber):
                diam  = self.model.getDiameter(self.fuel, number+1)
                self.modelClasses.addItem(number + 1, diam)
        else:
            self.treeViewRefusal.show()
            self.modelRefusal.deleteAll()
            RefusalNumber = self.model.getRefusalNumber(self.fuel)

            for number in range(0, ClassesNumber):
                diam  = self.model.getMassPercent(self.fuel, number+1)
                self.modelClasses.addItem(number + 1, diam)

            for number in range(0, RefusalNumber):
                refusal = self.model.getRefusal(self.fuel, number+1)
                self.modelRefusal.addItem(number+1, refusal)
                log.debug("slotDeleteRefusal number + 1 = %i " % (number+1))
            self._updateRefusalButton()
            self.pushButtonDeleteRefusal.show()
            self.pushButtonAddRefusal.show()


    def initializeNOxView(self):
        """
        initialize NOx tabview for a define solid fuel
        """
        self.labelReburning.hide()
        self.comboBoxReburning.hide()
        self.checkBoxNOxFormationFeature.hide()
        if self.model.getNOxFormationStatus() == 'on':
            self.checkBoxNOxFormation.setChecked(True)
            self.groupBoxNOxFormation.show()
            self.checkBoxNOxFormationFeature.show()
            self.lineEditQPR.setText(str(self.model.getNOxFormationParameter(self.fuel, 'nitrogen_fraction')))
            self.lineEditNitrogenConcentration.setText \
                (str(self.model.getNOxFormationParameter(self.fuel, 'nitrogen_concentration')))
            self.lineEditKobayashi1.setText(str(self.model.getHCNParameter \
                (self.fuel, "HCN_NH3_partitionning_reaction_1")))
            self.lineEditKobayashi2.setText(str(self.model.getHCNParameter \
                (self.fuel, "HCN_NH3_partitionning_reaction_2")))
            self.lineEditNitrogenLowTemp.setText(str(self.model.getNOxFormationParameter \
                (self.fuel, "nitrogen_in_char_at_low_temperatures")))
            self.lineEditNitrogenHighTemp.setText(str(self.model.getNOxFormationParameter \
                (self.fuel, "nitrogen_in_char_at_high_temperatures")))
            self.lineEditHCNChar.setText(str(self.model.getNOxFormationParameter \
                (self.fuel, "percentage_HCN_char_combustion")))
            if self.model.getNOxFormationFeature(self.fuel) == 'on':
                self.checkBoxNOxFormationFeature.setChecked(True)
                self.labelReburning.show()
                self.comboBoxReburning.show()
                mdl = self.model.getReburning(self.fuel)
                self.modelReburning.setItem(str_model=mdl)
            else:
                self.checkBoxNOxFormationFeature.setChecked(False)

        else:
            self.groupBoxNOxFormation.hide()


    def initializeKineticsView(self):
        """
        initialize kinetic tabview for a define solid fuel
        """
        if self.model.getCO2KineticsStatus() == 'on':
            self.checkBoxCO2Kinetics.setChecked(True)
            self.groupBoxParametersCO2.show()
            self.lineEditConstCO2.setText(str(self.model.getPreExponentialConstant(self.fuel, "CO2")))
            self.lineEditEnergyCO2.setText(str(self.model.getEnergyOfActivation(self.fuel, "CO2")))

            key = self.model.getOrderOfReaction(self.fuel, "CO2")
            self.modelReactTypeCO2.setItem(str_model=key)
            if key =='1':
                self.labelUnitConstCO2.setText('kg/m<sup>2</sup>/s/atm')
            elif key =='0.5':
                self.labelUnitConstCO2.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')
        else:
            self.checkBoxCO2Kinetics.setChecked(False)
            self.groupBoxParametersCO2.hide()

        if self.model.getH2OKineticsStatus() == 'on':
            self.checkBoxH2OKinetics.setChecked(True)
            self.groupBoxParametersH2O.show()
            self.lineEditConstH2O.setText(str(self.model.getPreExponentialConstant(self.fuel, "H2O")))
            self.lineEditEnergyH2O.setText(str(self.model.getEnergyOfActivation(self.fuel, "H2O")))

            key = self.model.getOrderOfReaction(self.fuel, "H2O")
            self.modelReactTypeH2O.setItem(str_model=key)
            if key =='1':
                self.labelUnitConstH2O.setText('kg/m<sup>2</sup>/s/atm')
            elif key =='0.5':
                self.labelUnitConstH2O.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')
        else:
            self.checkBoxH2OKinetics.setChecked(False)
            self.groupBoxParametersH2O.hide()


    def initializeView(self):
        """
        initialize view for a define solid fuel
        """
        self.modelClasses = StandardItemModelClasses(self.model, self.fuel)
        self.treeViewClasses.setModel(self.modelClasses)
        self.treeViewClasses.resizeColumnToContents(0)
        self.treeViewClasses.resizeColumnToContents(1)

        self.initializeKineticModel()
        self.initializeDiameter()
        self.initializeKineticsView()
        self.initializeNOxView()
        self._updateClassButton()

        # General (composition)
        self.lineEditC.setText(str(self.model.getComposition(self.fuel, "C")))
        self.lineEditH.setText(str(self.model.getComposition(self.fuel, "H")))
        self.lineEditO.setText(str(self.model.getComposition(self.fuel, "O")))
        self.lineEditN.setText(str(self.model.getComposition(self.fuel, "N")))
        self.lineEditS.setText(str(self.model.getComposition(self.fuel, "S")))
        self.lineEditCCoke.setText(str(self.model.getCokeComposition(self.fuel, "C")))
        self.lineEditHCoke.setText(str(self.model.getCokeComposition(self.fuel, "H")))
        self.lineEditOCoke.setText(str(self.model.getCokeComposition(self.fuel, "O")))
        self.lineEditNCoke.setText(str(self.model.getCokeComposition(self.fuel, "N")))
        self.lineEditSCoke.setText(str(self.model.getCokeComposition(self.fuel, "S")))
        self.lineEditPCI.setText(str(self.model.getPCIValue(self.fuel)))
        self.lineEditCp.setText(str(self.model.getProperty(self.fuel, "specific_heat_average")))

        if self.model.getCoalCombustionModel() == 'homogeneous_fuel_moisture':
            if self.model.node_lagr != 'off':
                self.labelThermalCond.show()
                self.labelUnitThermalCond.show()
                self.lineEditThermalCond.show()
                self.lineEditThermalCond.setText(str(self.model.getProperty(self.fuel, "thermal_conductivity")))
            else:
                self.labelThermalCond.hide()
                self.labelUnitThermalCond.hide()
                self.lineEditThermalCond.hide()

        self.lineEditDensity.setText(str(self.model.getProperty(self.fuel, "density")))
        self.lineEditMoisture.setText(str(self.model.getProperty(self.fuel, "moisture")))
        self.lineEditVolatileMatter.setText(str(self.model.getProperty(self.fuel, "volatile_matter")))

        PCIChoice = self.model.getPCIChoice(self.fuel)
        self.modelPCI.setItem(str_model=PCIChoice)
        if PCIChoice == 'IGT_correlation':
            self.lineEditPCI.hide()
            self.comboBoxPCIType.hide()
            self.labelUnitPCI.hide()
        else:
            self.lineEditPCI.show()
            self.comboBoxPCIType.show()
            self.labelUnitPCI.show()
            PCIType = self.model.getPCIType(self.fuel)
            self.modelPCIType.setItem(str_model=PCIType)
            self.lineEditPCI.setText(str(self.model.getPCIValue(self.fuel)))

        # Ashes
        self.lineEditAshesRatio.setText(str(self.model.getProperty(self.fuel, "rate_of_ashes_on_mass")))
        self.lineEditAshesEnthalpy.setText(str(self.model.getProperty(self.fuel, "ashes_enthalpy")))
        self.lineEditAshesCp.setText(str(self.model.getProperty(self.fuel, "ashes_thermal_capacity")))

        # Devolatilisation
        Y1Y2Choice = self.model.getY1Y2(self.fuel)
        self.modelY1Y2.setItem(str_model=Y1Y2Choice)
        if Y1Y2Choice == 'automatic_CHONS':
            self.frameY1Y2.hide()
        else:
            self.frameY1Y2.show()
        self.lineEditCoefY1.setText(str(self.model.getY1StoichiometricCoefficient(self.fuel)))
        self.lineEditCoefY2.setText(str(self.model.getY2StoichiometricCoefficient(self.fuel)))

        A1 = self.model.getDevolatilisationParameter(self.fuel, "A1_pre-exponential_factor")
        A2 = self.model.getDevolatilisationParameter(self.fuel, "A2_pre-exponential_factor")
        E1 = self.model.getDevolatilisationParameter(self.fuel, "E1_energy_of_activation")
        E2 = self.model.getDevolatilisationParameter(self.fuel, "E2_energy_of_activation")
        self.lineEditCoefA1.setText(str(A1))
        self.lineEditCoefA2.setText(str(A2))
        self.lineEditCoefE1.setText(str(E1))
        self.lineEditCoefE2.setText(str(E2))

        # Combustion heterogene
        self.lineEditConstO2.setText(str(self.model.getPreExponentialConstant(self.fuel, "O2")))
        self.lineEditEnergyO2.setText(str(self.model.getEnergyOfActivation(self.fuel, "O2")))

        key = self.model.getOrderOfReaction(self.fuel, "O2")
        self.modelReactTypeO2.setItem(str_model=key)
        if key =='1':
            self.labelUnitConstO2.setText('kg/m<sup>2</sup>/s/atm')
        elif key =='0.5':
            self.labelUnitConstO2.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')

        key = self.model.getOxidantType()
        self.modelOxidantType.setItem(str_model=key)

        if self.model.getCoalCombustionModel() == 'homogeneous_fuel':
            moisture = self.model.getProperty(self.fuel, "moisture")
            self.lineEditMoisture.setText(str(moisture))
            self.labelMoisture.setDisabled(True)
            self.lineEditMoisture.setDisabled(True)


    @pyqtSlot(str)
    def slotKineticModel(self, text):
        """
        Change the diameter type
        """
        key = self.modelKineticModel.dicoV2M[str(text)]
        self.model.setKineticModel(key)


    @pyqtSlot("QModelIndex")
    def slotSelectCoal(self, text=None):
        """
        Display values for the current coal selected in the view.
        """
        row = self.treeViewCoals.currentIndex().row()
        log.debug("selectCoal row = %i "%row)

        self.fuel = row + 1

        self.initializeView()


    @pyqtSlot()
    def slotCreateCoal(self):
        """ create a new coal"""
        # Init
        self.treeViewCoals.clearSelection()
        self.modelCoals.addItem()

        self.initializeView()

        # update Properties and scalars
        self.model.createCoalModelScalarsAndProperties()

        # update Buttons
        self._updateCoalButton()


    @pyqtSlot()
    def slotDeleteCoal(self):
        """ cancel a coal"""
        row = self.treeViewCoals.currentIndex().row()
        log.debug("slotDeleteCoal row = %i "%row)

        if row == -1:
            return

        number = row + 1

        self.model.deleteSolidFuel(number)

        # suppress item
        self.modelCoals.deleteItem(row)

        # First coal is selected
        row = 0
        self.fuel = row + 1

        self.initializeView()

        # Update buttons
        self._updateCoalButton()


    @pyqtSlot()
    def slotCreateClass(self):
        """Create a new class"""
        self.model.createClass(self.fuel)

        # Init
        ClassNumber = self.model.getClassNumber(self.fuel)
        diameter_type = self.model.getDiameterType(self.fuel)
        if diameter_type == 'automatic':
            diam = self.model.getDiameter(self.fuel, ClassNumber)
        elif diameter_type == 'rosin-rammler_law':
            diam = self.model.getMassPercent(self.fuel, ClassNumber)
        self.modelClasses.addItem(ClassNumber, diam)

        log.debug("slotCreateClass number + 1 = %i " % (ClassNumber))

        self.model.createClassModelScalarsAndProperties(self.fuel)

        # Update buttons
        self._updateClassButton()


    @pyqtSlot()
    def slotDeleteClass(self):
        """ cancel a class diameter"""
        row = self.treeViewClasses.currentIndex().row()
        log.debug("slotDeleteClass  number = %i " % row)
        if row == -1:
            return

        number = row + 1

        self.model.deleteClass(self.fuel, number)

        # Init
        self.initializeDiameter()

        # Update buttons
        self._updateClassButton()


    @pyqtSlot()
    def slotCreateRefusal(self):
        """Create a new refusal"""
        diameter = self.model.defaultValues()['diameter']

        self.model.createRefusal(self.fuel)

        # Init
        RefusalNumber = self.model.getRefusalNumber(self.fuel)
        refusal = self.model.getRefusal(self.fuel, RefusalNumber)
        self.modelRefusal.addItem(RefusalNumber, refusal)
        log.debug("slotCreateRefusal number + 1 = %i " % (RefusalNumber))

        # Update buttons
        self._updateRefusalButton()


    @pyqtSlot()
    def slotDeleteRefusal(self):
        """ cancel a refusal"""
        row = self.treeViewRefusal.currentIndex().row()
        log.debug("slotDeleteRefusal  number = %i " % row)
        if row == -1:
            return

        number = row + 1

        self.model.deleteRefusal(self.fuel, number)

        # Init
        self.modelRefusal.deleteAll()
        RefusalNumber = self.model.getRefusalNumber(self.fuel)
        for number in range(0, RefusalNumber):
            refusal = self.model.getRefusal(self.fuel, number+1)
            self.modelRefusal.addItem(number+1, refusal)
            log.debug("slotDeleteRefusal number + 1 = %i " % (number+1))

        # Update buttons
        self._updateRefusalButton()


    @pyqtSlot()
    def slotCreateOxidant(self):
        """Create a new oxidant"""
        self.model.createOxidant()
        num = self.model.getOxidantNumber()
        self.modelOxidants.addItem(str(num))

        log.debug("slotCreateOxidant number = %i " % num)

        # Update buttons
        self._updateOxidantButton()


    @pyqtSlot()
    def slotDeleteOxidant(self):
        """ delete an oxidant"""
        row = self.tableViewOxidants.currentIndex().row()
        log.debug("slotDeleteOxidants number = %i " % row)
        if row == -1:
            return

        number = row + 1
        self.model.deleteOxidant(number)

        self.modelOxidants.deleteAll()
        for number in range(0, self.model.getOxidantNumber()):
            self.modelOxidants.addItem(number+1)

        # Update buttons
        self._updateOxidantButton()


    @pyqtSlot(str)
    def slotCComposition(self, text):
        """
        Change the C composition
        """
        if self.lineEditC.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setComposition(self.fuel, "C", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotHComposition(self, text):
        """
        Change the H composition
        """
        if self.lineEditH.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setComposition(self.fuel, "H", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotOComposition(self, text):
        """
        Change the O composition
        """
        if self.lineEditO.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setComposition(self.fuel, "O", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotNComposition(self, text):
        """
        Change the N composition
        """
        if self.lineEditN.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setComposition(self.fuel, "N", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotSComposition(self, text):
        """
        Change the S composition
        """
        if self.lineEditS.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setComposition(self.fuel, "S", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotPCI(self, text):
        """
        Change the PCI value
        """
        if self.lineEditPCI.validator().state == QValidator.Acceptable:
            PCI = from_qvariant(text, float)
            self.model.setPCIValue(self.fuel, PCI)


    @pyqtSlot(str)
    def slotPCIType(self, text):
        """
        Change the PCI type
        """
        key = self.modelPCIType.dicoV2M[str(text)]
        self.model.setPCIType(self.fuel, key)


    @pyqtSlot(str)
    def slotPCIChoice(self, text):
        """
        Change the PCI choice
        """
        key = self.modelPCI.dicoV2M[str(text)]
        self.model.setPCIChoice(self.fuel, key)
        if key == 'IGT_correlation':
            self.lineEditPCI.hide()
            self.comboBoxPCIType.hide()
            self.labelUnitPCI.hide()
        else:
            self.lineEditPCI.show()
            self.comboBoxPCIType.show()
            self.labelUnitPCI.show()
            PCIType = self.model.getPCIType(self.fuel)
            self.modelPCIType.setItem(str_model=PCIType)
            self.lineEditPCI.setText(str(self.model.getPCIValue(self.fuel)))


    @pyqtSlot(str)
    def slotCCompositionCoke(self, text):
        """
        Change the C composition for coke
        """
        if self.lineEditCCoke.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setCokeComposition(self.fuel, "C", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotHCompositionCoke(self, text):
        """
        Change the H composition for coke
        """
        if self.lineEditHCoke.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setCokeComposition(self.fuel, "H", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotOCompositionCoke(self, text):
        """
        Change the O composition for coke
        """
        if self.lineEditOCoke.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setCokeComposition(self.fuel, "O", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotNCompositionCoke(self, text):
        """
        Change the N composition for coke
        """
        if self.lineEditNCoke.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setCokeComposition(self.fuel, "N", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotSCompositionCoke(self, text):
        """
        Change the S composition for coke
        """
        if self.lineEditSCoke.validator().state == QValidator.Acceptable:
            composition = from_qvariant(text, float)
            self.model.setCokeComposition(self.fuel, "S", composition)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotDiameterType(self, text):
        """
        Change the diameter type
        """
        key = self.modelDiameter.dicoV2M[str(text)]
        self.model.setDiameterType(self.fuel, key)

        self.modelClasses = StandardItemModelClasses(self.model, self.fuel)
        self.treeViewClasses.setModel(self.modelClasses)
        self.treeViewClasses.resizeColumnToContents(0)
        self.treeViewClasses.resizeColumnToContents(1)

        self.initializeDiameter()
        self._updateClassButton()


    @pyqtSlot(str)
    def slotVolatileMatter(self, text):
        """
        Change the volatile matter
        """
        if self.lineEditVolatileMatter.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "volatile_matter", value)


    @pyqtSlot(str)
    def slotThermalCapacity(self, text):
        """
        Change the thermal capacity
        """
        if self.lineEditCp.validator().state == QValidator.Acceptable:
            Cp = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "specific_heat_average", Cp)


    @pyqtSlot(str)
    def slotThermalConductivity(self, text):
        """
        Change the thermal conductivity
        """
        if self.lineEditThermalCond.validator().state == QValidator.Acceptable:
            lam = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "thermal_conductivity", lam)


    @pyqtSlot(str)
    def slotDensity(self, text):
        """
        Change the density
        """
        if self.lineEditDensity.validator().state == QValidator.Acceptable:
            density = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "density", density)


    @pyqtSlot(str)
    def slotMoisture(self, text):
        """
        Change the moisture
        """
        if self.lineEditMoisture.validator().state == QValidator.Acceptable:
            moisture = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "moisture", moisture)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotAshesRatio(self, text):
        """
        Change the ashes ratio
        """
        if self.lineEditAshesRatio.validator().state == QValidator.Acceptable:
            ashesRatio = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "rate_of_ashes_on_mass", ashesRatio)
        else:
            msg = self.tr("This value must be between 0 and 100.")
            self.stbar.showMessage(msg, 2000)


    @pyqtSlot(str)
    def slotAshesFormingEnthalpy(self, text):
        """
        Change the ashes forming enthalpy
        """
        if self.lineEditAshesEnthalpy.validator().state == QValidator.Acceptable:
            ashesFormingEnthalpy = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "ashes_enthalpy", ashesFormingEnthalpy)


    @pyqtSlot(str)
    def slotAshesThermalCapacity(self, text):
        """
        Change the ashes thermal capacity
        """
        if self.lineEditAshesCp.validator().state == QValidator.Acceptable:
            ashesThermalCapacity = from_qvariant(text, float)
            self.model.setProperty(self.fuel, "ashes_thermal_capacity", ashesThermalCapacity)


    @pyqtSlot(str)
    def slotY1CH(self, text):
        """
        Change the Y1 stoichiometric coefficient
        """
        if self.lineEditCoefY1.validator().state == QValidator.Acceptable:
            Y1CH = from_qvariant(text, float)
            self.model.setY1StoichiometricCoefficient(self.fuel, Y1CH)


    @pyqtSlot(str)
    def slotY2CH(self, text):
        """
        Change the Y2 stoichiometric coefficient
        """
        if self.lineEditCoefY2.validator().state == QValidator.Acceptable:
            Y2CH = from_qvariant(text, float)
            self.model.setY2StoichiometricCoefficient(self.fuel, Y2CH)


    @pyqtSlot(str)
    def slotY1Y2(self, text):
        """
        Change the Y1Y2 type
        """
        key = self.modelY1Y2.dicoV2M[str(text)]
        self.model.setY1Y2(self.fuel, key)
        if key == 'automatic_CHONS':
            self.frameY1Y2.hide()
        else:
            self.frameY1Y2.show()
            self.lineEditCoefY1.setText(str(self.model.getY1StoichiometricCoefficient(self.fuel)))
            self.lineEditCoefY2.setText(str(self.model.getY2StoichiometricCoefficient(self.fuel)))


    @pyqtSlot(str)
    def slotA1CH(self, text):
        """
        Change the pre exponential factor A1
        """
        if self.lineEditCoefA1.validator().state == QValidator.Acceptable:
            A1CH = from_qvariant(text, float)
            self.model.setDevolatilisationParameter(self.fuel, "A1_pre-exponential_factor", A1CH)


    @pyqtSlot(str)
    def slotA2CH(self, text):
        """
        Change the pre exponentiel factor A2
        """
        if self.lineEditCoefA2.validator().state == QValidator.Acceptable:
            A2CH = from_qvariant(text, float)
            self.model.setDevolatilisationParameter(self.fuel, "A2_pre-exponential_factor", A2CH)


    @pyqtSlot(str)
    def slotE1CH(self, text):
        """
        Change the energy of activation E1
        """
        if self.lineEditCoefE1.validator().state == QValidator.Acceptable:
            E1CH = from_qvariant(text, float)
            self.model.setDevolatilisationParameter(self.fuel, "E1_energy_of_activation", E1CH)


    @pyqtSlot(str)
    def slotE2CH(self, text):
        """
        Change the Energy of activation E2
        """
        if self.lineEditCoefE2.validator().state == QValidator.Acceptable:
            E2CH = from_qvariant(text, float)
            self.model.setDevolatilisationParameter(self.fuel, "E2_energy_of_activation", E2CH)


    @pyqtSlot(str)
    def slotPreExpoCstO2(self, text):
        """
        Change the pre exponential constant for O2
        """
        if self.lineEditConstO2.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setPreExponentialConstant(self.fuel, "O2", value)


    @pyqtSlot(str)
    def slotActivEnergyO2(self, text):
        """
        Change the energy of activation for O2
        """
        if self.lineEditEnergyO2.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEnergyOfActivation(self.fuel, "O2", value)


    @pyqtSlot(str)
    def slotReactTypeO2(self, text):
        """
        Change the order of reaction of O2
        """
        key = self.modelReactTypeO2.dicoV2M[str(text)]
        self.model.setOrderOfReaction(self.fuel, "O2", key)
        if text =='1':
            self.labelUnitConstO2.setText('kg/m<sup>2</sup>/s/atm')
        elif text =='0.5':
            self.labelUnitConstO2.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')


    @pyqtSlot(str)
    def slotPreExpoCstCO2(self, text):
        """
        Change the preexponential constant for CO2
        """
        if self.lineEditConstCO2.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setPreExponentialConstant(self.fuel, "CO2", value)


    @pyqtSlot(str)
    def slotActivEnergyCO2(self, text):
        """
        Change the energy of activation for CO2
        """
        if self.lineEditEnergyCO2.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEnergyOfActivation(self.fuel, "CO2", value)


    @pyqtSlot(str)
    def slotReactTypeCO2(self, text):
        """
        Change the order of reaction for CO2
        """
        key = self.modelReactTypeCO2.dicoV2M[str(text)]
        self.model.setOrderOfReaction(self.fuel, "CO2", key)
        if text =='1':
            self.labelUnitConstCO2.setText('kg/m<sup>2</sup>/s/atm')
        elif text =='0.5':
            self.labelUnitConstCO2.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')


    @pyqtSlot(str)
    def slotPreExpoCstH2O(self, text):
        """
        Change the pre exponential constant for H2O
        """
        if self.lineEditConstH2O.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setPreExponentialConstant(self.fuel, "H2O", value)


    @pyqtSlot(str)
    def slotActivEnergyH2O(self, text):
        """
        Change the energy of activation for H2O
        """
        if self.lineEditEnergyH2O.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEnergyOfActivation(self.fuel, "H2O", value)


    @pyqtSlot(str)
    def slotReactTypeH2O(self, text):
        """
        Change the order of reaction
        """
        key = self.modelReactTypeH2O.dicoV2M[str(text)]
        self.model.setOrderOfReaction(self.fuel, "H2O", key)
        if text =='1':
            self.labelUnitConstH2O.setText('kg/m<sup>2</sup>/s/atm')
        elif text =='0.5':
            self.labelUnitConstH2O.setText('kg/m<sup>2</sup>/s/atm<sup>1/2</sup>')


    @pyqtSlot(str)
    def slotQPR(self, text):
        """
        Change the nitrogen fraction
        """
        if self.lineEditQPR.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setNOxFormationParameter(self.fuel, 'nitrogen_fraction', value)


    @pyqtSlot(str)
    def slotNitrogenConcentration(self, text):
        """
        Change the nitrogen concentration
        """
        if self.lineEditNitrogenConcentration.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setNOxFormationParameter(self.fuel, 'nitrogen_concentration', value)


    @pyqtSlot(str)
    def slotKobayashi1(self, text):
        """
        Change the nitrogen partition reaction of reaction 1
        """
        if self.lineEditKobayashi1.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setHCNParameter(self.fuel, "HCN_NH3_partitionning_reaction_1", value)


    @pyqtSlot(str)
    def slotKobayashi2(self, text):
        """
        Change the Nitrogen partition reaction of reaction 2
        """
        if self.lineEditKobayashi2.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setHCNParameter(self.fuel, "HCN_NH3_partitionning_reaction_2", value)


    @pyqtSlot(str)
    def slotNLowTemp(self, text):
        """
        Change the nitrogen in char at low temperatures
        """
        if self.lineEditNitrogenLowTemp.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setNOxFormationParameter(self.fuel, 'nitrogen_in_char_at_low_temperatures', value)


    @pyqtSlot(str)
    def slotNHighTemp(self, text):
        """
        Change the nitrogen in char at  temperatures
        """
        if self.lineEditNitrogenHighTemp.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setNOxFormationParameter(self.fuel, 'nitrogen_in_char_at_high_temperatures', value)


    @pyqtSlot(str)
    def slotHCNChar(self, text):
        """
        Change the nitrogen percentage in char combustion
        """
        if self.lineEditHCNChar.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setNOxFormationParameter(self.fuel, 'percentage_HCN_char_combustion', value)


    @pyqtSlot(str)
    def slotOxidantType(self, text):
        """
        Change the oxidant type
        """
        key = self.modelOxidantType.dicoV2M[str(text)]
        self.model.setOxidantType(key)


    @pyqtSlot(str)
    def slotReburning(self, text):
        """
        Change the reburning type
        """
        key = self.modelReburning.dicoV2M[str(text)]
        self.model.setReburning(self.fuel, key)


    @pyqtSlot(bool)
    def slotNOxFormation(self, checked):
        """
        check box for NOx formation
        """
        status = 'off'
        if checked:
            status = 'on'
        self.model.setNOxFormationStatus(status)
        self.initializeNOxView()


    @pyqtSlot(bool)
    def slotNOxFeature(self, checked):
        """
        check box for NOx formation
        """
        status = 'off'
        if checked:
            status = 'on'
        self.model.setNOxFormationFeature(self.fuel, status)
        self.initializeNOxView()


    @pyqtSlot(bool)
    def slotCO2Kinetics(self, checked):
        """
        check box for CO2 kinetics
        """
        status = 'off'
        if checked:
            status = 'on'
        self.model.setCO2KineticsStatus(status)
        self.initializeKineticsView()


    @pyqtSlot(bool)
    def slotH2OKinetics(self, checked):
        """
        check box for H2O kinetics
        """
        status = 'off'
        if checked:
            status = 'on'
        self.model.setH2OKineticsStatus(status)
        self.initializeKineticsView()


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
