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
This module defines the Immersed boundaries view data management.

This module contains the following classes and function:
- SyrthesVerbosityDelegate
- ProjectionAxisDelegate
- SelectionCriteriaDelegate
- StandardItemModelSyrthes
- ConjugateHeatTransferView
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
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.ImmersedBoundariesNeptune import Ui_ImmersedBoundariesNeptune
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesViewNeptune")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QLineEdit delegate to attach a label to the FSI object
#-------------------------------------------------------------------------------

class FSIObjectNameDelegate(QItemDelegate):

    def __init__(self, parent = None):
        super(FSIObjectNameDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                                   to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "":
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QComboBox delegate for the FSI type : Set motion or computed from fluid forces
#-------------------------------------------------------------------------------

class FSIMovingDelegate(QItemDelegate):
    """
    USe of a comboBox to set the moving attribute of the object
    """
    def __init__(self, parent):
        super(FSIMovingDelegate, self).__init__(parent)
        self.parent  = parent

    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ['non_moving', 'moving']:
            editor.addItem(itm)

        editor.installEventFilter(self)
        return editor

    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataFSI[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QComboBox delegate for the FSI type : Set motion or computed from fluid forces
#-------------------------------------------------------------------------------

class FSITypeDelegate(QItemDelegate):
    """
    Use of a combobox to set the fsi interaction type
    """

    def __init__(self, parent, mdl):
        super(FSITypeDelegate, self).__init__(parent)
        self.parent  = parent
        self.mdl     = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ["off", "imposed", "computed"]:
            editor.addItem(itm)

        if self.mdl.getObjectMotion(index.row()+1) == 'moving':
            editor.model().item(0).setEnabled(False)
        else:
            editor.model().item(0).setEnabled(True)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataFSI[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for validation of Syrthes verbosity or visualization
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# QComboBox delegate for Axis Projection in Conjugate Heat Transfer table
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelFSI(QStandardItemModel):

    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Object name"),
                        self.tr("Object motion"),
                        self.tr("Interaction type")]
        self.tooltip = [self.tr("Name of solid object"),
                        self.tr("Is the object moving or not"),
                        self.tr("Type of interaction with the flow")]

        self.setColumnCount(len(self.headers))
        self.dataFSI = []
        self.__model = model


    def data(self, index, role):
        if not index.isValid():
            return None

        # Tooltips
        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]

        # Display
        if role == Qt.DisplayRole:
            return self.dataFSI[index.row()][index.column()]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        if index.column() in (0, 1):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 2:
            if self.__model.getObjectMotion(index.row()+1) == 'non_moving':
                return Qt.NoItemFlags
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return

        row = index.row()
        col = index.column()

        self.dataFSI[row][col] = str(from_qvariant(value, to_text_string))

        num = row + 1
        self.__model.setObjectName(num, self.dataFSI[row][0])
        self.__model.setObjectMotion(num, self.dataFSI[row][1])
        if self.dataFSI[row][1] == 'non_moving':
            self.dataFSI[row][2] = 'off'
        else:
            if self.dataFSI[row][2] == 'off':
                self.dataFSI[row][2] = 'imposed'

        self.__model.setObjectInteraction(num, self.dataFSI[row][2])

#        self.dataChanged.emit(index, index)

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getData(self, index):
        row = index.row()
        return self.dataFSI[row]

    def addItem(self, object_name, motion_type, interaction_type):
        """
        Add a row in the table.
        """
        self.dataFSI.append([object_name, motion_type, interaction_type])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataFSI[row]
        row = self.rowCount()
        self.setRowCount(row-1)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesViewNeptune(QWidget, Ui_ImmersedBoundariesNeptune):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesNeptune.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.ibm = ImmersedBoundariesModel(self.case)

        self.current_obj = None

        # Models
        self.modelFSI = StandardItemModelFSI(self.ibm)
        self.tableViewFSI.setModel(self.modelFSI)

        if QT_API == "PYQT4":
            self.tableViewFSI.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setResizeMode(4, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFSI.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)

        self.modelFSI.dataChanged.connect(self.dataChanged)

        delegateObjectLabel  = FSIObjectNameDelegate(self.tableViewFSI)
        self.tableViewFSI.setItemDelegateForColumn(0, delegateObjectLabel)

        delegateObjectMotion = FSIMovingDelegate(self.tableViewFSI)
        self.tableViewFSI.setItemDelegateForColumn(1, delegateObjectMotion)

        delegateObjectType   = FSITypeDelegate(self.tableViewFSI, self.ibm)
        self.tableViewFSI.setItemDelegateForColumn(2, delegateObjectType)

        self.checkBoxActivate.stateChanged.connect(self.slotCheckActivate)

        self.tableViewFSI.clicked.connect(self.slotChangedSelection)


        for ind in ['Explicit', 'MEDCoupling']:
            eval('self.radioButton'+ind+'.toggled.connect(self.slotRadioButton)')

        # Connections
        self.pushButtonAddFSI.clicked.connect(self.slotAddFSI)
        self.pushButtonDeleteFSI.clicked.connect(self.slotDeleteFSI)

        # Check for MEDCoupling presence
        import cs_config
        cfg = cs_config.config()
        self.has_medcoupling = cfg.libs['medcoupling'].have == 'yes'
        # deactivated for the moment
        self.has_medcoupling = False
        self.radioButtonMEDCoupling.setEnabled(self.has_medcoupling)
        if self.ibm.getMethod() == 'medcoupling' and self.has_medcoupling == False:
            self.setMethod('explicit')

        # Show/hide widgets on start
        if self.ibm.getOnOff() == 'off':
            self.groupBoxMethod.hide()
            self.groupBoxObjects.hide()
            self.groupBoxObjProperties.hide()
            self.groupBoxExplicit.hide()
            self.groupBoxMEDCoupling.hide()
        else:
            self.groupBoxMethod.show()
            self.groupBoxObjects.show()
            if self.ibm.getMethod() == 'explicit':
                self.groupBoxExplicit.show()
            else:
                self.groupBoxMEDCoupling.show()

        self.updatePageView()

        self.case.undoStartGlobal()


    def slotChangedSelection(self, text=None):
        """
        detect change in selection and update view
        """
        row = self.tableViewFSI.currentIndex().row()
        self.current_obj = row + 1
        self.updatePageView()


    @pyqtSlot(int)
    def slotCheckActivate(self, val):

        # Set the method state
        if val == 0:
            self.ibm.setOnOff('off')
        else:
            self.ibm.setOnOff('on')

        # Update the view if needed
        self.updatePageView()

    def dataChanged(self, topLeft, bottomRight):
#        self.tableViewFSI.resizeColumnsToContents()
#        self.tableViewFSI.resizeRowsToContents()

#        if QT_API == "PYQT4":
#            self.tableViewFSI.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
#        elif QT_API == "PYQT5":
#            self.tableViewFSI.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)

        self.updatePageView()

    def updatePageView(self):

        if self.ibm.getOnOff() == 'off':
            self.groupBoxMethod.hide()
            self.groupBoxObjects.hide()
            self.groupBoxExplicit.hide()
            self.groupBoxMEDCoupling.hide()

        else:
            self.groupBoxMethod.show()
            self.groupBoxObjects.show()

            # If motion is set to 'fixed', interaction is switched back to off
            for obj in range(self.ibm.getNumberOfFSIObjects()):
                if self.ibm.getObjectMotion(obj+1) == 'non_moving':
                    self.ibm.setObjectInteraction(obj+1, 'off')

            # Which button to show for the solid definition
            if self.current_obj:
                if self.ibm.getMethod() == 'explicit':
                    self.groupBoxExplicit.show()
                    self.groupBoxMEDCoupling.hide()
                    self.radioButtonExplicit.setChecked(True)
                elif self.ibm.getMethod() == 'medcoupling':
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.show()
                    self.radioButtonMEDCoupling.setChecked(True)
                else:
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.hide()

                if self.ibm.getObjectMotion(self.current_obj) == 'moving':
                    self.groupBoxObjProperties.show()
                else:
                    self.groupBoxObjProperties.hide()


    @pyqtSlot()
    def slotRadioButton(self):

        for ind in ['Explicit', 'MEDCoupling']:

            radioButton = eval('self.radioButton'+ind)
            if radioButton.isChecked():
                self.ibm.setMethod(ind.lower())

        self.updatePageView()


    @pyqtSlot()
    def slotAddFSI(self):

        name        = '_'.join([self.ibm.defaultValues()['fsi_object_name'],
                                str(self.ibm.getNumberOfFSIObjects()+1)])
        motion      = self.ibm.defaultValues()['fsi_moving']
        interaction = self.ibm.defaultValues()['fsi_interaction']

        num = self.ibm.addFSIObject(name, motion, interaction)
        self.modelFSI.addItem(name, motion, interaction)


    @pyqtSlot()
    def slotDeleteFSI(self):
        row = self.tableViewFSI.currentIndex().row()
        log.debug("slotDeleteFSI -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing object")
            QMessageBox.information(self, title, msg)
        else:
            self.modelFSI.deleteRow(row)
            self.ibm.deleteFSIObject(row+1)

    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
