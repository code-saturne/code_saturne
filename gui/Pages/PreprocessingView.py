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
This module contains the following classes and function:
- StandardItemModelThinWall
- StandardItemModelExtrude
- PreprocessingView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, logging
try:
    import ConfigParser
    configparser = ConfigParser
except Exception:
    import configparser

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, RegExpValidator, IntValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.PreprocessingForm import Ui_PreprocessingForm
from code_saturne.model.SolutionDomainModel import RelOrAbsPath, MeshModel, SolutionDomainModel
from code_saturne.Pages.SolutionDomainView import MeshNumberDelegate
from code_saturne.Pages.FacesSelectionView import StandardItemModelFaces

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PreprocessingView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for selection
#-------------------------------------------------------------------------------

class LineEditDelegateSelector(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator =  RegExpValidator(editor, QRegExp("[ -~]*"))
        editor.setValidator(validator)
        return editor


    def setEditorData(self, lineEdit, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        lineEdit.setText(value)


    def setModelData(self, lineEdit, model, index):
        value = lineEdit.text()
        model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Line edit delegate for float (thickness and reason)
#-------------------------------------------------------------------------------

class FloatDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(FloatDelegate, self).__init__(parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor, min=0.)
        validator.setExclusiveMin(True)
        editor.setValidator(validator)
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
# Line edit delegate for integer
#-------------------------------------------------------------------------------

class IntDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(IntDelegate, self).__init__(parent)


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0)
        validator.setExclusiveMin(True)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandardItemModelThinWall class
#-------------------------------------------------------------------------------

class StandardItemModelThinWall(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("zone id"),
                         self.tr("selector")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl


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
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 :
            return Qt.ItemIsSelectable
        else:
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
        name = row

        if col == 1:
            new_sup = from_qvariant(value, to_text_string)
            self._data[row][col] = new_sup
            self.mdl.replaceThinWall(name, new_sup)

        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_tw=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        name = ""
        if existing_tw == None :
            self.mdl.addThinWall()
            name = row
        else:
            name = existing_tw

        support = self.mdl.getThinWall(name)

        thin = [str(name), support]

        self._data.append(thin)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        self.mdl.deleteThinWall(row)
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# StandardItemModelExtrude class
#-------------------------------------------------------------------------------

class StandardItemModelExtrude(QStandardItemModel):

    def __init__(self, parent, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [ self.tr("zone id"),
                         self.tr("n layers"),
                         self.tr("thickness"),
                         self.tr("expansion factor"),
                         self.tr("selector")]

        self.setColumnCount(len(self.headers))
        self.parent = parent

        self.tooltip = []

        self._data  = []
        self.mdl    = mdl


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
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 :
            return Qt.ItemIsSelectable
        else:
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
        name = row

        if col == 1:
            new_sup = from_qvariant(value, int)
            self._data[row][col] = new_sup
            self.mdl.setExtrudeLayer(name, new_sup)
        elif col == 2:
            new_sup = from_qvariant(value, float)
            self._data[row][col] = new_sup
            self.mdl.setExtrudeThickness(name, new_sup)
        elif col == 3:
            new_sup = from_qvariant(value, float)
            self._data[row][col] = new_sup
            self.mdl.setExtrudeReason(name, new_sup)
        elif col == 4:
            new_sup = from_qvariant(value, to_text_string)
            self._data[row][col] = new_sup
            self.mdl.setExtrudeSelector(name, new_sup)

        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_tw=None):
        """
        Add/load a scalar in the model.
        """
        row = self.rowCount()

        name = ""
        if existing_tw == None :
            self.mdl.addExtrude()
            name = row
        else:
            name = existing_tw

        nlayer  = self.mdl.getExtrudeLayer(name)
        thick   = self.mdl.getExtrudeThickness(name)
        reason  = self.mdl.getExtrudeReason(name)
        support = self.mdl.getExtrudeSelector(name)

        thin = [str(name), nlayer, thick, reason, support]

        self._data.append(thin)
        self.setRowCount(row + 1)


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        del self._data[row]
        self.mdl.deleteExtrude(row)
        row = self.rowCount()
        self.setRowCount(row - 1)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class PreprocessingView(QWidget, Ui_PreprocessingForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_PreprocessingForm.__init__(self)
        self.setupUi(self)

        self.stbar = stbar
        self.case = case
        self.case.undoStopGlobal()
        self.mdl = SolutionDomainModel(self.case)

        # 2) Meshe preprocessing layout

        # 2.2) Connections

        self.groupBoxWarp.clicked[bool].connect(self.slotFacesCutting)
        self.lineEditWarp.textChanged[str].connect(self.slotWarpParam)
        self.groupBoxMeshSmooth.clicked[bool].connect(self.slotMeshSmooth)
        self.lineEditMeshSmooth.textChanged[str].connect(self.slotMeshSmoothParam)

        # 2.3) Set up validators
        validatorWarp = DoubleValidator(self.lineEditWarp, min=0.0)
        self.lineEditWarp.setValidator(validatorWarp)
        validatorSmooth = DoubleValidator(self.lineEditMeshSmooth, min=0.0, max=90.0)
        self.lineEditMeshSmooth.setValidator(validatorSmooth)

        # 2.4) Faces to join selection (Custom Widgets)

        model = StandardItemModelFaces(self, self.mdl, 'face_joining')
        self.widgetFacesJoin.modelFaces = model
        self.widgetFacesJoin.tableView.setModel(model)

        # 2.5) Thin wall

        self.tableModelThinWall = StandardItemModelThinWall(self, self.mdl)
        self.tableViewThinWall.setModel(self.tableModelThinWall)
        if QT_API == "PYQT4":
            self.tableViewThinWall.horizontalHeader().setResizeMode(0,QHeaderView.Stretch)
            self.tableViewThinWall.horizontalHeader().setResizeMode(1,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewThinWall.horizontalHeader().setSectionResizeMode(0,QHeaderView.Stretch)
            self.tableViewThinWall.horizontalHeader().setSectionResizeMode(1,QHeaderView.Stretch)
        self.tableViewThinWall.resizeColumnsToContents()
        self.tableViewThinWall.resizeRowsToContents()
        self.tableViewThinWall.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewThinWall.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateLabel   = MeshNumberDelegate(self.tableViewThinWall)
        delegateSupport = LineEditDelegateSelector(self.tableViewThinWall)

        self.tableViewThinWall.setItemDelegateForColumn(0, delegateLabel)
        self.tableViewThinWall.setItemDelegateForColumn(1, delegateSupport)

        # Connections
        self.pushButtonAddThinWall.clicked.connect(self.slotAddThinWall)
        self.pushButtonDeleteThinWall.clicked.connect(self.slotDeleteThinWall)

        # load values
        for tw in range(self.mdl.getThinWallSelectionsCount()):
            self.tableModelThinWall.newItem(tw)

        # 2.5) Extrude

        self.tableModelExtrude = StandardItemModelExtrude(self, self.mdl)
        self.tableViewExtrude.setModel(self.tableModelExtrude)
        if QT_API == "PYQT4":
            self.tableViewExtrude.horizontalHeader().setResizeMode(4,QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewExtrude.horizontalHeader().setSectionResizeMode(4,QHeaderView.Stretch)
        self.tableViewExtrude.resizeColumnsToContents()
        self.tableViewExtrude.resizeRowsToContents()
        self.tableViewExtrude.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewExtrude.setSelectionMode(QAbstractItemView.SingleSelection)

        delegateLabel   = MeshNumberDelegate(self.tableViewExtrude)
        delegateLayer   = IntDelegate(self.tableViewExtrude)
        delegateSupport = LineEditDelegateSelector(self.tableViewExtrude)
        delegateFloat = FloatDelegate(self.tableViewExtrude)

        self.tableViewExtrude.setItemDelegateForColumn(0, delegateLabel)   # Id
        self.tableViewExtrude.setItemDelegateForColumn(1, delegateLayer)   # nlayers
        self.tableViewExtrude.setItemDelegateForColumn(2, delegateFloat)   # thickness
        self.tableViewExtrude.setItemDelegateForColumn(3, delegateFloat)   # reason
        self.tableViewExtrude.setItemDelegateForColumn(4, delegateSupport) # criteria

        # Connections
        self.pushButtonAddExtrude.clicked.connect(self.slotAddExtrude)
        self.pushButtonDeleteExtrude.clicked.connect(self.slotDeleteExtrude)

        # load values
        for tw in range(self.mdl.getExtrudeSelectionsCount()):
            self.tableModelExtrude.newItem(tw)

        # 3) Periodicities

        self.perio_mode = ""

        # Model for periodicities

        model = StandardItemModelFaces(self, self.mdl, 'face_periodicity')
        self.widgetFacesPerio.modelFaces = model
        self.widgetFacesPerio.tableView.setModel(model)

        # Combo model for type of periodicity
        self.modelComboPeriod = ComboModel(self.comboBoxPeriodicity, 3, 1)
        self.modelComboPeriod.addItem(self.tr("Periodicity by translation"), "translation")
        self.modelComboPeriod.addItem(self.tr("Periodicity by rotation (defined by angle and direction)"), "rotation")
        self.modelComboPeriod.addItem(self.tr("Composite periodicity (defined by matrix)"), "mixed")

        # Display
        self.groupBoxMode.hide()
        self.groupBoxTranslation.hide()
        self.groupBoxRotation.hide()
        self.groupBoxMixed.hide()

        # Set up validators

        # 4)
        self.lineEditTX.setValidator(DoubleValidator(self.lineEditTX))
        self.lineEditTY.setValidator(DoubleValidator(self.lineEditTY))
        self.lineEditTZ.setValidator(DoubleValidator(self.lineEditTZ))
        self.lineEditAngle.setValidator(DoubleValidator(self.lineEditAngle))
        self.lineEditDX.setValidator(DoubleValidator(self.lineEditDX))
        self.lineEditDY.setValidator(DoubleValidator(self.lineEditDY))
        self.lineEditDZ.setValidator(DoubleValidator(self.lineEditDZ))
        self.lineEditX1.setValidator(DoubleValidator(self.lineEditX1))
        self.lineEditY1.setValidator(DoubleValidator(self.lineEditY1))
        self.lineEditZ1.setValidator(DoubleValidator(self.lineEditZ1))
        self.lineEditM11.setValidator(DoubleValidator(self.lineEditM11))
        self.lineEditM12.setValidator(DoubleValidator(self.lineEditM12))
        self.lineEditM13.setValidator(DoubleValidator(self.lineEditM13))
        self.lineEditM14.setValidator(DoubleValidator(self.lineEditM14))
        self.lineEditM21.setValidator(DoubleValidator(self.lineEditM21))
        self.lineEditM22.setValidator(DoubleValidator(self.lineEditM22))
        self.lineEditM23.setValidator(DoubleValidator(self.lineEditM23))
        self.lineEditM24.setValidator(DoubleValidator(self.lineEditM24))
        self.lineEditM31.setValidator(DoubleValidator(self.lineEditM31))
        self.lineEditM32.setValidator(DoubleValidator(self.lineEditM32))
        self.lineEditM33.setValidator(DoubleValidator(self.lineEditM33))
        self.lineEditM34.setValidator(DoubleValidator(self.lineEditM34))

        # Connections

        selectionModel = self.widgetFacesPerio.tableView.selectionModel()
        selectionModel.selectionChanged.connect(self.slotUpdatePeriodicity)

        self.widgetFacesPerio.pushButtonDelete.clicked.connect(self.slotDeletePeriodicity)

        self.comboBoxPeriodicity.activated[str].connect(self.slotPeriodicityMode)

        self.lineEditTX.textChanged[str].connect(self.slotTranslationX)
        self.lineEditTY.textChanged[str].connect(self.slotTranslationY)
        self.lineEditTZ.textChanged[str].connect(self.slotTranslationZ)

        self.lineEditAngle.textChanged[str].connect(self.slotAngleRotation)

        self.lineEditDX.textChanged[str].connect(self.slotRotationX)
        self.lineEditDY.textChanged[str].connect(self.slotRotationY)
        self.lineEditDZ.textChanged[str].connect(self.slotRotationZ)

        self.lineEditX1.textChanged[str].connect(self.slotCenterRotationX1)
        self.lineEditY1.textChanged[str].connect(self.slotCenterRotationY1)
        self.lineEditZ1.textChanged[str].connect(self.slotCenterRotationZ1)

        self.lineEditM11.textChanged[str].connect(self.slotMatrix11)
        self.lineEditM12.textChanged[str].connect(self.slotMatrix12)
        self.lineEditM13.textChanged[str].connect(self.slotMatrix13)
        self.lineEditM14.textChanged[str].connect(self.slotMatrix14)
        self.lineEditM21.textChanged[str].connect(self.slotMatrix21)
        self.lineEditM22.textChanged[str].connect(self.slotMatrix22)
        self.lineEditM23.textChanged[str].connect(self.slotMatrix23)
        self.lineEditM24.textChanged[str].connect(self.slotMatrix24)
        self.lineEditM31.textChanged[str].connect(self.slotMatrix31)
        self.lineEditM32.textChanged[str].connect(self.slotMatrix32)
        self.lineEditM33.textChanged[str].connect(self.slotMatrix33)
        self.lineEditM34.textChanged[str].connect(self.slotMatrix34)
        self.tabWidget.currentChanged[int].connect(self.slotchanged)

        # Warped faces cutting

        if self.mdl.getCutStatus() == 'on':
            self.groupBoxWarp.setChecked(True)
            self.slotFacesCutting(True)
        else:
            self.groupBoxWarp.setChecked(False)
            self.slotFacesCutting(False)

        v = self.mdl.getCutAngle()
        self.warp = v
        self.lineEditWarp.setText(str(self.warp))


        # Mesh Smoothing

        if self.mdl.getSmoothingStatus() == 'on':
            self.groupBoxMeshSmooth.setChecked(True)
            self.slotMeshSmooth(True)
        else:
            self.groupBoxMeshSmooth.setChecked(False)
            self.slotMeshSmooth(False)

        v = self.mdl.getSmoothAngle()
        self.smooth = v
        self.lineEditMeshSmooth.setText(str(self.smooth))

        # tab Widget
        self.tabWidget.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


    @pyqtSlot(bool)
    def slotFacesCutting(self, checked):
        """
        Private slot.

        Do we cut any warp faces ?

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox warp parameters
        """
        self.groupBoxWarp.setFlat(not checked)
        if checked:
            self.mdl.setCutStatus("on")
            self.frameWarp.show()
        else:
            self.mdl.setCutStatus("off")
            self.frameWarp.hide()


    @pyqtSlot(str)
    def slotWarpParam(self, text):
        """
        Private slot.

        @type text: C{QString}
        @param text: max angle of warped faces
        """
        if self.lineEditWarp.validator().state == QValidator.Acceptable:
            var = float(text)
            self.mdl.setCutAngle(var)


    @pyqtSlot(bool)
    def slotMeshSmooth(self, checked):
        """
        Private slot.

        Do we use mesh smoothing ?

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox mesh smooth parameters
        """
        self.groupBoxMeshSmooth.setFlat(not checked)
        if checked:
            self.mdl.setSmoothingStatus("on")
            self.frameSmooth.show()
        else:
            self.mdl.setSmoothingStatus("off")
            self.frameSmooth.hide()


    @pyqtSlot(str)
    def slotMeshSmoothParam(self, text):
        """
        Private slot.

        @type text: C{QString}
        @param text: angle for mesh smoothing
        """
        if self.lineEditMeshSmooth.validator().state == QValidator.Acceptable:
            var = float(text)
            self.mdl.setSmoothAngle(var)


    @pyqtSlot()
    def slotDeletePeriodicity(self):
        """
        Delete a periodicity from the list.
        """

        log.debug("slotDeletePeriodicity  = %s " % self.perio_mode)

        sel_rows = self.widgetFacesPerio.tableView.selectionModel().selectedIndexes()

        perio_id = None

        for row in sel_rows:

            perio_id = row.row()
            perio_mode = self.mdl.getPeriodicityMode(perio_id)

            self.perio_id = perio_id
            self.perio_mode = perio_mode
            self.modelComboPeriod.setItem(str_model=perio_mode)
            txt = str(self.comboBoxPeriodicity.currentText())
            self.slotPeriodicityMode(txt)

        if perio_id == None:

            self.groupBoxMode.hide()
            self.groupBoxTranslation.hide()
            self.groupBoxRotation.hide()
            self.groupBoxMixed.hide()


    def __setValuesTranslation(self, perio_id):
        """
        Put values found in xml file as soon as mode is "translation"
        """
        dx, dy, dz = self.mdl.getTranslationDirection(perio_id)

        self.lineEditTX.setText(str(dx))
        self.lineEditTY.setText(str(dy))
        self.lineEditTZ.setText(str(dz))


    def __setValuesRotation(self, perio_id):
        """
        Put values found in xml file as soon as mode is "rotation"
        """
        angle = self.mdl.getRotationAngle(perio_id)
        rx, ry, rz = self.mdl.getRotationDirection(perio_id)
        px, py, pz = self.mdl.getRotationCenter(perio_id)

        self.lineEditAngle.setText(str((angle)))
        self.lineEditDX.setText(str(rx))
        self.lineEditDY.setText(str(ry))
        self.lineEditDZ.setText(str(rz))
        self.lineEditX1.setText(str(px))
        self.lineEditY1.setText(str(py))
        self.lineEditZ1.setText(str(pz))


    def __setValuesMixed(self, perio_id):
        """
        Put values found in xml file as soon as mode is "rotation"2
        """
        m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34 = self.mdl.getTransformationMatrix(perio_id)

        self.lineEditM11.setText(str(m11))
        self.lineEditM12.setText(str(m12))
        self.lineEditM13.setText(str(m13))
        self.lineEditM14.setText(str(m14))
        self.lineEditM21.setText(str(m21))
        self.lineEditM22.setText(str(m22))
        self.lineEditM23.setText(str(m23))
        self.lineEditM24.setText(str(m24))
        self.lineEditM31.setText(str(m31))
        self.lineEditM32.setText(str(m32))
        self.lineEditM33.setText(str(m33))
        self.lineEditM34.setText(str(m34))


    def __setValuesPeriodicTransformation(self, perio, mode):
        """
        Put values found in xml file as soon as mode of
        transformation is choosen
        """
        log.debug("__setValuesPeriodicTransformation perio mode = %s %s "% (perio, mode))
        if mode == "translation" :
            self.__setValuesTranslation(perio)
        if mode == "rotation":
            self.__setValuesRotation(perio)
        if mode == "mixed":
            self.__setValuesMixed(perio)


    @pyqtSlot("QItemSelection")
    def slotUpdatePeriodicity(self, current):
        """
        This slot updates the display for the periodicity selected
        in the table view.
        """
        index = self.widgetFacesPerio.tableView.currentIndex()
        log.debug("slotUpdatePeriodicity index.row() = %i " % index.row())

        self.groupBoxMode.show()

        perio_id = index.row()
        if perio_id < 0:
            return

        perio_mode = self.mdl.getPeriodicityMode(perio_id)
        self.perio_id = perio_id
        self.perio_mode = perio_mode

        self.modelComboPeriod.setItem(str_model=perio_mode)
        txt = str(self.comboBoxPeriodicity.currentText())
        self.slotPeriodicityMode(txt)

    @pyqtSlot(str)
    def slotPeriodicityMode(self, text):
        """
        Do we have a periodicity ?
        """

        self.perio_mode = self.modelComboPeriod.dicoV2M[str(text)]

        log.debug("slotPeriodicityMode  = %s " % self.perio_mode)

        self.groupBoxTranslation.hide()
        self.groupBoxRotation.hide()
        self.groupBoxMixed.hide()

        if self.perio_mode == "":
            self.modelComboPeriod(str_model='translation')

        if self.perio_mode == "translation":
            self.groupBoxTranslation.show()
            self.groupBoxRotation.hide()
            self.groupBoxMixed.hide()

        elif self.perio_mode == "rotation":
            self.groupBoxTranslation.hide()
            self.groupBoxRotation.show()
            self.groupBoxMixed.hide()

        elif self.perio_mode == "mixed":
            self.groupBoxTranslation.hide()
            self.groupBoxRotation.hide()
            self.groupBoxMixed.show()

        sel_rows = self.widgetFacesPerio.tableView.selectionModel().selectedIndexes()

        for row in sel_rows:

            perio_id = row.row()

            self.mdl.updatePeriodicityMode(perio_id, self.perio_mode)

            if self.perio_mode == "":
                self.__setValuesTranslation(perio_id)

            if self.perio_mode == "translation":
                self.__setValuesTranslation(perio_id)

            elif self.perio_mode == "rotation":
                self.__setValuesRotation(perio_id)

            elif self.perio_mode == "mixed":
                self.__setValuesMixed(perio_id)

            self.__setValuesPeriodicTransformation(self.perio_id, self.perio_mode)


    @pyqtSlot(str)
    def slotTranslationX(self, text):
        """
        Periodicity translation for X
        """
        if self.perio_mode != "rotation" or self.perio_mode != "mixed":
            if self.lineEditTX.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setTranslationDirection(self.perio_id, 'translation_x', val)


    @pyqtSlot(str)
    def slotTranslationY(self, text):
        """
        Periodicity translation for Y
        """
        if self.perio_mode != "rotation" or self.perio_mode != "mixed":
            if self.lineEditTY.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setTranslationDirection(self.perio_id, 'translation_y', val)


    @pyqtSlot(str)
    def slotTranslationZ(self, text):
        """
        Periodicity translation for Z
        """
        if self.perio_mode != "rotation" or self.perio_mode != "mixed":
            if self.lineEditTZ.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setTranslationDirection(self.perio_id, 'translation_z', val)


    @pyqtSlot(str)
    def slotAngleRotation(self, text):
        """
        Periodicity rotation angle
        """
        if self.perio_mode == "rotation":
            if self.lineEditAngle.validator().state == QValidator.Acceptable:
                angle = float(text)
                self.mdl.setRotationAngle(self.perio_id, angle)


    @pyqtSlot(str)
    def slotRotationX(self, text):
        """
        Periodicity rotation for X
        """
        if self.perio_mode == "rotation":
            if self.lineEditDX.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationVector(self.perio_id, "axis_x", val)


    @pyqtSlot(str)
    def slotRotationY(self, text):
        """
        Periodicity rotation for Y
        """
        if self.perio_mode == "rotation":
            if self.lineEditDY.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationVector(self.perio_id, "axis_y", val)


    @pyqtSlot(str)
    def slotRotationZ(self, text):
        """
        Periodicity rotation for Z
        """
        if self.perio_mode == "rotation":
            if self.lineEditDZ.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationVector(self.perio_id, "axis_z", val)


    @pyqtSlot(str)
    def slotCenterRotationX1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.lineEditX1.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_x", val)


    @pyqtSlot(str)
    def slotCenterRotationY1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.lineEditY1.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_y", val)


    @pyqtSlot(str)
    def slotCenterRotationZ1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.lineEditZ1.validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_z", val)


    # Methods for matrix components
    @pyqtSlot(str)
    def slotMatrix11(self, text):
        self.__cmdTransformationMatrix("matrix_11", text)


    @pyqtSlot(str)
    def slotMatrix12(self, text):
        self.__cmdTransformationMatrix("matrix_12", text)


    @pyqtSlot(str)
    def slotMatrix13(self, text):
        self.__cmdTransformationMatrix("matrix_13", text)


    @pyqtSlot(str)
    def slotMatrix14(self, text):
        self.__cmdTransformationMatrix("matrix_14", text)


    @pyqtSlot(str)
    def slotMatrix21(self, text):
        self.__cmdTransformationMatrix("matrix_21", text)


    @pyqtSlot(str)
    def slotMatrix22(self, text):
        self.__cmdTransformationMatrix("matrix_22", text)


    @pyqtSlot(str)
    def slotMatrix23(self, text):
        self.__cmdTransformationMatrix("matrix_23", text)


    @pyqtSlot(str)
    def slotMatrix24(self, text):
        self.__cmdTransformationMatrix("matrix_24", text)


    @pyqtSlot(str)
    def slotMatrix31(self, text):
        self.__cmdTransformationMatrix("matrix_31", text)


    @pyqtSlot(str)
    def slotMatrix32(self, text):
        self.__cmdTransformationMatrix("matrix_32", text)


    @pyqtSlot(str)
    def slotMatrix33(self, text):
        self.__cmdTransformationMatrix("matrix_33", text)


    @pyqtSlot(str)
    def slotMatrix34(self, text):
        self.__cmdTransformationMatrix("matrix_34", text)


    def __cmdTransformationMatrix(self, pos, text):
        """
        Periodicity translation
        """
        if self.perio_mode == "mixed":
            if self.sender().validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setTransformationMatrix(self.perio_id, pos, val)


    @pyqtSlot(str)
    def slotCenterRotationX2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.sender().validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_x", val)


    @pyqtSlot(str)
    def slotCenterRotationY2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.sender().validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_y", val)


    @pyqtSlot(str)
    def slotCenterRotationZ2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            if self.sender().validator().state == QValidator.Acceptable:
                val = float(text)
                self.mdl.setRotationCenter(self.perio_id, "invariant_z", val)


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


    @pyqtSlot()
    def slotAddThinWall(self):
        """
        Add a thin wall
        """
        self.tableViewThinWall.clearSelection()
        self.tableModelThinWall.newItem()


    @pyqtSlot()
    def slotDeleteThinWall(self):
        """
        Delete the a thin wall from the list (one by one).
        """
        row = self.tableViewThinWall.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteThinWall -> %s" % row)
            self.tableModelThinWall.deleteItem(row)


    @pyqtSlot()
    def slotAddExtrude(self):
        """
        Add a thin wall
        """
        self.tableViewExtrude.clearSelection()
        self.tableModelExtrude.newItem()


    @pyqtSlot()
    def slotDeleteExtrude(self):
        """
        Delete the a thin wall from the list (one by one).
        """
        row = self.tableViewExtrude.currentIndex().row()
        if row >= 0 :
            log.debug("slotDeleteExtrude -> %s" % row)
            self.tableModelExtrude.deleteItem(row)


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
