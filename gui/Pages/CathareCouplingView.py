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
This module defines the cathare coupling view data management.

This module contains the following classes and function:
- CathareEltDelegate
- NcfdBcDelegate
- SelectionCriteriaDelegate
- StandardItemModelCathare
- CathareCouplingView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging
import os

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *
import code_saturne.Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.CathareCouplingForm import Ui_CathareCouplingForm
from code_saturne.model.CathareCouplingModel import CathareCouplingModel
from code_saturne.model.LocalizationModelNeptune import LocalizationModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CathareCouplingView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QLineEdit delegate for the activation of the Cathare coupling
#-------------------------------------------------------------------------------

class CathareEltDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(CathareEltDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "" :
            model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QComboBox delegate for Axis Projection in Conjugate Heat Transfer table
#-------------------------------------------------------------------------------

class NcfdBcDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent = None, case=None):
        super(NcfdBcDelegate, self).__init__(parent)
        self.parent = parent
        self.case = case


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        editor.addItem("off")

        if self.case:
            d = LocalizationModel('BoundaryZone', self.case)

            for zone in d.getZones():
                editor.addItem(zone.getLabel())

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataCathare[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for location
#-------------------------------------------------------------------------------

class SelectionCriteriaDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(SelectionCriteriaDelegate, self).__init__(parent)
        self.parent = parent
        self.__model = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)

        validator =  RegExpValidator(editor, QRegExp("[ -~]*"))
        editor.setValidator(validator)

        # Autocompletion for selection criteria!
        comp_list = ['all[]']
        comp_list += ['plane[a, b, c, d, epsilon]',
                      'plane[a, b, c, d, inside]',
                      'plane[a, b, c, d, outside]',
                      'plane[n_x, n_y, n_z, x0, y0, z0, epsilon]',
                      'plane[n_x, n_y, n_z, x0, y0, z0, inside]',
                      'plane[n_x, n_y, n_z, x0, y0, z0, outside]',
                      'box[xmin, ymin, zmin, xmax, ymax, zmax]',
                      'box[x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3]',
                      'cylinder[x0, y0, z0, x1, y1, z1, radius]',
                      'sphere[x_c, y_c, z_c, radius]']

        completer = QCompleter()
        editor.setCompleter(completer)
        model = QStringListModel()
        completer.setModel(model)
        model.setStringList(comp_list)

        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "" :
            model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelCathare(QStandardItemModel):

    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(5)

        self.headers = [self.tr("Cathare Element"),
                        self.tr("First cell"),
                        self.tr("Last cell"),
                        self.tr("Neptune BC"),
                        self.tr("Neptune 1D volume")]

        self.tooltip = [self.tr("Name of the coupled Cathare element"),
                        self.tr("First cell before the coupling frontier"),
                        self.tr("First cell after the coupling frontier"),
                        self.tr("NEPTUNE_CFD coupled boundary condition"),
                        self.tr("NEPTUNE_CFD volume equivalent to Cathare's 'Last cell'")]

        self.setColumnCount(len(self.headers))
        self.dataCathare = []
        self.__model = model


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]
        if role == Qt.DisplayRole:
            return self.dataCathare[index.row()][index.column()]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return

        row = index.row()
        if index.column() in (0, 3, 4):
            self.dataCathare[row][index.column()] = str(from_qvariant(value, to_text_string))
        else:
            self.dataCathare[row][index.column()] = from_qvariant(value, int)

        num = row + 1
        self.__model.setCathareEltName(num, self.dataCathare[row][0])
        self.__model.setCathareFCell(num, self.dataCathare[row][1])
        self.__model.setCathareLCell(num, self.dataCathare[row][2])

        bc_zone = None
        bc_type = None
        d = LocalizationModel('BoundaryZone', self.__model.case)
        for zone in d.getZones():
            if zone.getLabel() == self.dataCathare[row][3]:
                bc_zone = zone.getLocalization()
                bc_type = zone.getNature()
                break
        self.__model.setNeptuneBc(num, self.dataCathare[row][3])
        self.__model.setNeptune2dZone(num,bc_zone)
        self.__model.setNeptuneBcType(num, bc_type)

        self.__model.setNeptune1dZone(num, self.dataCathare[row][4])

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def addItem(self, CathareElt, CFirstCell, CLastCell, NeptuneBC, Neptune1dZone):
        """
        Add a row in the table.
        """
        self.dataCathare.append([CathareElt,
                                 CFirstCell,
                                 CLastCell,
                                 NeptuneBC,
                                 Neptune1dZone])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataCathare[row]
        row = self.rowCount()
        self.setRowCount(row-1)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class CathareCouplingView(QWidget, Ui_CathareCouplingForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_CathareCouplingForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.case.undoStopGlobal()

        self.__model = CathareCouplingModel(self.case)

        # Main combo box
        self.activateCathareCpl = QtPage.ComboModel(self.comboBoxActiveCpl,2,1)
        self.activateCathareCpl.addItem(self.tr("No coupling"),       "off")
        self.activateCathareCpl.addItem(self.tr("Activate coupling"), "on")


        # Models
        self.modelCathare = StandardItemModelCathare(self.__model)
        self.tableViewCathare.setModel(self.modelCathare)

        if QT_API == "PYQT4":
            self.tableViewCathare.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCathare.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCathare.horizontalHeader().setResizeMode(4, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewCathare.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCathare.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewCathare.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)

        delegateCathareElt = CathareEltDelegate(self.tableViewCathare)
        self.tableViewCathare.setItemDelegateForColumn(0, delegateCathareElt)

        delegateCathareFCell = CathareEltDelegate(self.tableViewCathare)
        self.tableViewCathare.setItemDelegateForColumn(1, delegateCathareFCell)

        delegateCathareLCell = CathareEltDelegate(self.tableViewCathare)
        self.tableViewCathare.setItemDelegateForColumn(2, delegateCathareLCell)

        delegateNeptuneBc = NcfdBcDelegate(self.tableViewCathare, case=case)
        self.tableViewCathare.setItemDelegateForColumn(3, delegateNeptuneBc)

        delegateNeptune1dZone = SelectionCriteriaDelegate(self.tableViewCathare,
                                                          self.__model)
        self.tableViewCathare.setItemDelegateForColumn(4, delegateNeptune1dZone)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddCathare)
        self.pushButtonDelete.clicked.connect(self.slotDeleteCathare)

        self.comboBoxActiveCpl.activated[str].connect(self.slotActivateCpl)

        self.radioButtonOnePhase.clicked.connect(self.slotOnePhase)
        self.radioButtonAllPhases.clicked.connect(self.slotAllPhases)

        self.lineEditCathareFile.textChanged[str].connect(self.slotCathareFile)
        self.toolButtonCathareFile.pressed.connect(self.searchCathareJDD)
        self.lineEditCplName.textChanged[str].connect(self.slotCplName)

        self.lineEditCplTime.textChanged[str].connect(self.slotCplTime)
        self.lineEditCathareInitTime.textChanged[str].connect(self.slotCathareTime)

        self.lineEditCathareInstance.textChanged[str].connect(self.slotCathareInstanceName)
        self.toolButtonCathareInstance.pressed.connect(self.searchCathareDir)
        self.lineEditNeptuneInstance.textChanged[str].connect(self.slotNeptuneInstanceName)
        self.toolButtonNeptuneInstance.pressed.connect(self.searchNeptuneDir)

        # Insert list of Cathare couplings for view
        for c in self.__model.getCathareCouplingList():
            [cathare_elt, cathare_first_cell, cathare_last_cell,
             neptune_bc, neptune_1d_zone] = c

            self.modelCathare.addItem(cathare_elt,
                                      cathare_first_cell, cathare_last_cell,
                                      neptune_bc, neptune_1d_zone)


        # ------------------------------------------
        # Activate the coupling parameters if needed
        if self.__model.getCathareActivationStatus() != 0:
            self.activateCathareCpl.setItem(str_model="on")


        if self.__getActivationState() == 'off':
            self.groupBoxCathareCpls.hide()
            self.groupBoxCplParameters.hide()
        else:
            self.groupBoxCathareCpls.show()
            self.groupBoxCplParameters.show()

        # ------------------------------------------
        if self.__model.getNphases() == 0 or self.__model.getNphases() == 1:
            self.radioButtonOnePhase.setChecked(True)
            self.radioButtonAllPhases.setChecked(False)
        else:
            self.radioButtonOnePhase.setChecked(False)
            self.radioButtonAllPhases.setChecked(True)

        # ------------------------------------------
        if self.__getActivationState() == 'on':
            self.lineEditCathareFile.setText(str(self.__model.getCathareFile()))
            self.lineEditCplName.setText(str(self.__model.getCplName()))

            self.lineEditCplTime.setText(str(self.__model.getCplTime()))
            self.lineEditCathareInitTime.setText(str(self.__model.getCathareTime()))

            self.lineEditCathareInstance.setText(str(self.__model.getCathareInstanceName()))
            self.lineEditNeptuneInstance.setText(str(self.__model.getNeptuneInstanceName()))

        # ------------------------------------------
        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotOnePhase(self):
        """
        Set the number of coupled phases to one
        """

        self.__model.setNphases(1)

    @pyqtSlot()
    def slotAllPhases(self):


        np = self.__model.getNumberOfFluids()

        self.__model.setNphases(np)

    @pyqtSlot()
    def slotAddCathare(self):
        """
        Set in view label and variables to see on profile
        """
        cathare_elt     = self.__model.defaultValues()['cathare_elt']
        cathare_fcell   = self.__model.defaultValues()['cathare_first_cell']
        cathare_lcell   = self.__model.defaultValues()['cathare_last_cell']
        neptune_bc      = self.__model.defaultValues()['neptune_bc']
        neptune_1d_zone = self.__model.defaultValues()['neptune_1d_zone']

        num = self.__model.addCathareCoupling(cathare_elt,
                                              cathare_fcell, cathare_lcell,
                                              neptune_bc, neptune_1d_zone)

        self.modelCathare.addItem(cathare_elt,
                                  cathare_fcell, cathare_lcell,
                                  neptune_bc, neptune_1d_zone)


    @pyqtSlot()
    def slotDeleteCathare(self):
        """
        Delete the profile from the list (one by one).
        """
        row = self.tableViewCathare.currentIndex().row()
        log.debug("slotDeleteCathare -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing coupling")
            QMessageBox.information(self, title, msg)
        else:
            self.modelCathare.deleteRow(row)
            self.__model.deleteCathareCoupling(row+1)


    def __getActivationState(self):

        combo = self.comboBoxActiveCpl
        dico  = self.activateCathareCpl.dicoV2M

        cpl_state = dico[str(combo.currentText())]

        return cpl_state

    @pyqtSlot()
    def slotActivateCpl(self):
        """
        Activate or Deactivate the NEPTUNE_CFD/CATHARE coupling
        """

        cpl_state = self.__getActivationState()

        if cpl_state == "off":
            self.groupBoxCathareCpls.hide()
            self.groupBoxCplParameters.hide()
            self.__model.setApiType(0)
        else:
            self.groupBoxCathareCpls.show()
            self.groupBoxCplParameters.show()
            self.__model.setApiType(1)

            self.lineEditCathareFile.setText(str(self.__model.getCathareFile()))
            self.lineEditCplName.setText(str(self.__model.getCplName()))

            self.lineEditCplTime.setText(str(self.__model.getCplTime()))
            self.lineEditCathareInitTime.setText(str(self.__model.getCathareTime()))

    @pyqtSlot()
    def slotCathareFile(self):

        value = str(self.lineEditCathareFile.text())
        self.__model.setCathareFile(value)

    @pyqtSlot()
    def slotCplName(self):

        value = str(self.lineEditCplName.text())
        self.__model.setCplName(value)

    @pyqtSlot()
    def slotCplTime(self):

        value = float(self.lineEditCplTime.text())
        self.__model.setCplTime(value)


    @pyqtSlot()
    def slotCathareTime(self):

        value = float(self.lineEditCathareInitTime.text())
        self.__model.setCathareTime(value)


    @pyqtSlot()
    def slotCathareInstanceName(self):

        value = str(self.lineEditCathareInstance.text())
        self.__model.setCathareInstanceName(value)

    @pyqtSlot()
    def slotNeptuneInstanceName(self):

        value = str(self.lineEditNeptuneInstance.text())
        self.__model.setNeptuneInstanceName(value)

    def searchCathareDir(self):
        """
        Open a File Dialog in order to search the case directory.
        """
        title    = self.tr("Select NEPTUNE_CFD instance for coupling")
        default  = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setFileMode(QFileDialog.Directory)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()
            instance_name = os.path.split(str(s[0]))[-1]

            self.lineEditCathareInstance.setText(instance_name)
            self.__model.setCathareInstanceName(instance_name)

    def searchNeptuneDir(self):
        """
        Open a File Dialog in order to search the case directory.
        """
        title    = self.tr("Select CATHARE2 instance for coupling")
        default  = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setFileMode(QFileDialog.Directory)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()
            instance_name = os.path.split(str(s[0]))[-1]

            self.lineEditNeptuneInstance.setText(instance_name)
            self.__model.setNeptuneInstanceName(instance_name)

    def searchCathareJDD(self):
        """
        Open a File Dialog in order to search the case directory.
        """
        title    = self.tr("Select CATHARE2 data file for coupling")
        default  = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setNameFilter(self.tr("*.dat"))

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()
            file_name = os.path.split(str(s[0]))[-1]

            self.lineEditCathareFile.setText(file_name)
            self.__model.setCathareFile(file_name)

    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
