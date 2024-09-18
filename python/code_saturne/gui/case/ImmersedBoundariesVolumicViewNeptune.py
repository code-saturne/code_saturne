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
This module defines the Immersed volume zones view data management.

This module contains the following classes and function:
- ImmersedBoundariesVolumicViewNeptune
"""

#-------------------------------------------------------------------------------
# Standard modules
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
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.ImmersedBoundariesVolumicFormNeptune import Ui_ImmersedBoundariesVolumicFormNeptune
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesVolumicViewNeptune")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QComboBox delegate for the moving type : Set motion or computed from fluid forces
#-------------------------------------------------------------------------------

class MovingTypeDelegate(QItemDelegate):
    """
    Use of a combobox to set the fsi moving type
    """

    def __init__(self, parent, mdl):
        super(MovingTypeDelegate, self).__init__(parent)
        self.parent  = parent
        self.mdl     = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ["fixed", "imposed", "computed"]:
            editor.addItem(itm)

        self.updateComboBoxItems(editor, index)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().data_vol[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)


    def updateComboBoxItems(self, comboBox, index):
        row = index.row()

        is_fsi = self.mdl.getObjectFSI(row+1)

        if is_fsi == 'on':
            for i in range(comboBox.count()):
                comboBox.model().item(i).setEnabled(True)
        else:
            for i in range(comboBox.count()):
                item = comboBox.model().item(i)
                if item.text() == "computed":
                    item.setEnabled(False)
                else:
                    item.setEnabled(True)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemVolume(QStandardItemModel):

    def __init__(self, model, case, is_thermal, tree):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Object name"),
                        self.tr("Fluid Stucture Interaction"),
                        self.tr("Moving"),
                        self.tr("Conjugate heat transfer"),
                        self.tr("Initialization"),
                        self.tr("Physical properties"),
                        self.tr("thermal source term")]

        self.tooltip = [self.tr("Name of the solid object"),
                        self.tr("Solve interaction between the fluid and the structure"),
                        self.tr("Type of motion interaction with the flow"),
                        self.tr("Conjugate heat transfer between the fluid and the structure"),
                        self.tr("Initialization of the solid object"),
                        self.tr("Physical properties for the solid object"),
                        self.tr("Thermal source term")]

        self.setColumnCount(len(self.headers))
        self.data_vol = []
        self.__model = model
        self.case = case
        self.is_thermal = is_thermal
        self.browser = tree

    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        # Tooltips
        if role == Qt.ToolTipRole:
            return self.tooltip[col]

        elif role == Qt.CheckStateRole:
            if col == 1 or col == 3 or col == 4 or col == 5 or col == 6:
                data_vol = self.data_vol[row][col]
                if data_vol == "on":
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        # Display
        elif role == Qt.DisplayRole:
            return self.data_vol[row][col]

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        if col == 0:
            return Qt.ItemIsSelectable

        elif (col == 1):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable

        elif (col == 2):
            is_fsi = self.data_vol[row][1]
            if is_fsi == 'on':
                # The moving type is always 'computed' when solving FSI
                return Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

        elif col == 3:
            if self.is_thermal == False:
                # CHT can't be activated without thermal
                return Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable

        elif (col == 4):
            return Qt.ItemIsSelectable

        elif (col == 5):
            return Qt.ItemIsSelectable

        elif (col == 6):
            return Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        num = row + 1

        # Check if thermal (only if energy resolution is activated)
        if (self.is_thermal == False):
            self.data_vol[row][3] = "off"

        if col == 1:
            state = from_qvariant(value, int)

            if state == Qt.Unchecked:
                self.data_vol[row][col] = "off"
            else:
                self.data_vol[row][col] = "on"

            self.__model.setObjectFSI(num, self.data_vol[row][1])

            is_fsi = self.data_vol[row][1]
            is_moving = self.data_vol[row][2]
            is_CHT = self.data_vol[row][3]

            if is_fsi == 'on':
                self.data_vol[row][2] = "computed"
                self.data_vol[row][5] = "on"
                self.__model.setObjectMoving(num, self.data_vol[row][2])
                self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])
            else:
                if (is_moving == "computed"): # computed is only for FSI
                    self.data_vol[row][2] = "imposed"
                    self.__model.setObjectMoving(num, self.data_vol[row][2])

                if (is_moving in ['fixed'] and is_CHT == 'off'):
                    self.data_vol[row][5] = "off"
                    self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])
                else:
                    self.data_vol[row][5] = "on"
                    self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])

            if (is_CHT == "on"): # or (is_fsi == "on" and self.is_thermal == True)):
                self.data_vol[row][4] = "on"
                self.__model.setObjectInit(num, self.data_vol[row][4])
            else:
                self.data_vol[row][4] = "off"
                self.__model.setObjectInit(num, self.data_vol[row][4])


        elif col == 2:
            self.data_vol[row][col] = str(from_qvariant(value, to_text_string))
            self.__model.setObjectMoving(num, self.data_vol[row][col])

            is_moving = self.data_vol[row][col]
            is_fsi = self.data_vol[row][1]
            is_CHT = self.data_vol[row][3]

            if (is_moving == "imposed"):
                self.data_vol[row][5] = "on"
                self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])

            elif (is_moving in ["computed","fixed"]):
                if (is_CHT == 'off'):
                    self.data_vol[row][5] = "off"
                    self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])


        elif col == 3: #CHT
            state = from_qvariant(value, int)

            if (self.is_thermal == False):
                self.data_vol[row][col] = "off"
            else:
                if state == Qt.Unchecked:
                    self.data_vol[row][col] = "off"
                else:
                    self.data_vol[row][col] = "on"

            self.__model.setObjectCHT(num, self.data_vol[row][col])
            self.browser.configureTree(self.case)

            is_fsi = self.data_vol[row][1]
            is_moving = self.data_vol[row][2]
            is_CHT = self.data_vol[row][3]

            if (is_CHT == "on"):
                self.data_vol[row][4] = "on"
                self.data_vol[row][5] = "on"
                self.data_vol[row][6] = "on"
                self.__model.setObjectInit(num, self.data_vol[row][4])
                self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])
                self.__model.setObjectThermalSourceTerm(num, self.data_vol[row][6])
            else:

                self.data_vol[row][6] = "off"
                self.__model.setObjectThermalSourceTerm(num, self.data_vol[row][6])

                if (is_fsi == 'off'):
                    if (is_moving == 'computed'):
                        self.data_vol[row][2] == "imposed"
                        self.__model.setObjectMoving(num, self.data_vol[row][col])

                    if (is_moving == 'fixed'):
                        self.data_vol[row][4] = "off"
                        self.data_vol[row][5] = "off"
                        self.__model.setObjectInit(num, self.data_vol[row][4])
                        self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])
                    elif (is_moving == 'imposed'):
                        self.data_vol[row][4] = "off"
                        self.data_vol[row][5] = "on"
                        self.__model.setObjectInit(num, self.data_vol[row][4])
                        self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])
                else: # FSI and off-CHT
                    self.data_vol[row][5] = "on"
                    self.__model.setObjectPhysicalProperties(num, self.data_vol[row][5])

                    if (self.is_thermal == True):
                        self.data_vol[row][4] = "on"
                        self.__model.setObjectInit(num, self.data_vol[row][4])
                    else:
                        self.data_vol[row][4] = "off"
                        self.__model.setObjectInit(num, self.data_vol[row][4])


        elif col == 4:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.data_vol[row][col] = "off"
            else:
                self.data_vol[row][col] = "on"

            self.__model.setObjectInit(num, self.data_vol[row][col])

        elif col == 5:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.data_vol[row][col] = "off"
            else:
                self.data_vol[row][col] = "on"

            self.__model.setObjectPhysicalProperties(num, self.data_vol[row][col])

        elif col == 6:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.data_vol[row][col] = "off"
            else:
                self.data_vol[row][col] = "on"

            self.__model.setObjectThermalSourceTerm(num, self.data_vol[row][col])

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getData(self, index):
        row = index.row()
        return self.data_vol[row]

    def addItem(self, object_name, object_is_fsi, object_moving,
                object_is_cht, object_is_init, object_is_phy, object_st):
        """
        Add a row in the table.
        """
        self.data_vol.append([object_name, object_is_fsi, object_moving,
                              object_is_cht, object_is_init,
                              object_is_phy, object_st])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.data_vol[row]
        row = self.rowCount()
        self.setRowCount(row-1)

    def getItem(self, row):
        return self.data_vol[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesVolumicViewNeptune(QWidget, Ui_ImmersedBoundariesVolumicFormNeptune):
    """
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesVolumicFormNeptune.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.ibm = ImmersedBoundariesModel(self.case)
        self.mfm = MainFieldsModel(self.case)
        self.current_obj = None
        self.browser = tree

        # Check if thermal (only if energy resolution is activated)
        self.is_thermal = False
        for fId in self.mfm.getFieldIdList():
            if self.mfm.getFieldFromId(fId).enthalpy_model != 'off':
                self.is_thermal = True

        # Models
        self.model_vol = StandardItemVolume(self.ibm, self.case,
                                            self.is_thermal, self.browser)
        self.tableViewIBMVolumezone.setModel(self.model_vol)

        for obj in range(1,self.ibm.getNumberOfObjects()+1):

            if (self.ibm.getObjectMoving(obj) == "computed" and
                self.ibm.getObjectFSI(obj) == "off"):
                self.ibm.setObjectMoving(obj, "imposed")

            #initialization
            if (self.is_thermal == False):
                self.ibm.setObjectCHT(obj, "off")
                self.ibm.setObjectThermalSourceTerm(obj, "off")
                self.ibm.setObjectInit(obj, "off")
                self.ibm.setObjectBoundaryEnergyMode(obj, "off")

                if (self.ibm.getObjectMoving(obj) == "imposed" or
                    self.ibm.getObjectFSI(obj) == "on"):

                    self.ibm.setObjectPhysicalProperties(obj, "on")
                else:
                    self.ibm.setObjectPhysicalProperties(obj, "off")
            else:
                if (self.ibm.getObjectCHT(obj) == "on"):
                    self.ibm.setObjectInit(obj, "on") # for temp
                    self.ibm.setObjectPhysicalProperties(obj, "on")
                    self.ibm.setObjectBoundaryEnergyMode(obj, "off")
                    self.ibm.setObjectThermalSourceTerm(obj, "on")
                else:
                    self.ibm.setObjectInit(obj, "off")
                    self.ibm.setObjectThermalSourceTerm(obj, "off")

                    if (self.ibm.getObjectFSI(obj) == "on"):
                        self.ibm.setObjectPhysicalProperties(obj, "on")
                        self.ibm.setObjectBoundaryEnergyMode(obj, "off")
                    else:
                        self.ibm.setObjectInit(obj, "off")
                        self.ibm.setObjectPhysicalProperties(obj, "off")
                        if (self.ibm.getObjectMoving(obj) == "imposed"):
                            self.ibm.setObjectPhysicalProperties(obj, "on")
                        self.ibm.setObjectBoundaryEnergyMode(obj, "temperature")

            self.model_vol.addItem(self.ibm.getObjectName(obj),
                                   self.ibm.getObjectFSI(obj),
                                   self.ibm.getObjectMoving(obj),
                                   self.ibm.getObjectCHT(obj),
                                   self.ibm.getObjectInit(obj),
                                   self.ibm.getObjectPhysicalProperties(obj),
                                   self.ibm.getObjectThermalSourceTerm(obj))

        if QT_API == "PYQT4":
            self.tableViewIBMVolumezone.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewIBMVolumezone.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            #self.tableViewIBMVolumezone.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewIBMVolumezone.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewIBMVolumezone.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            #self.tableViewIBMVolumezone.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        delegateMovingType = MovingTypeDelegate(self.tableViewIBMVolumezone, self.ibm)
        self.tableViewIBMVolumezone.setItemDelegateForColumn(2, delegateMovingType)

        self.model_vol.dataChanged.connect(self.dataChanged)

        self.tableViewIBMVolumezone.clicked[QModelIndex].connect(self.slotChangedSelection)

        self.updatePageView()

        self.case.undoStartGlobal()


    @pyqtSlot("QModelIndex")
    def slotChangedSelection(self, index):
        """
        detect change in selection and update view
        """
        row = self.tableViewIBMVolumezone.currentIndex().row()
        self.current_obj = row + 1
        self.updatePageView()

    def dataChanged(self, topLeft, bottomRight):
        self.updatePageView()

    def updatePageView(self):

        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            self.groupBoxVolumeZone.hide()
            return

        self.groupBoxVolumeZone.show()

        current_obj = self.tableViewIBMVolumezone.currentIndex().row() + 1

        # Desactivate boundary condition for CHT computation
        if (self.ibm.getObjectCHT(current_obj) == "on"):
            self.ibm.setObjectBoundaryEnergyMode(current_obj, "off")

        self.browser.configureTree(self.case)



#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
