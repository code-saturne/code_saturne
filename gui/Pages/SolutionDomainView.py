# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module contains the following classes and function:
- SpinBoxDelegate
- StandardItemModelMeshes
- StandardItemModelPeriod
- SolutionDomainMeshFormatDialogView
- SolutionDomainView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from SolutionDomainForm import Ui_SolutionDomainForm
from SolutionDomainMeshFormatDialogForm import Ui_SolutionDomainMeshFormatDialogForm
from Pages.SolutionDomainModel import SolutionDomainModel
from Pages.SolutionDomainModel import MeshModel
from FacesSelectionView import StandardItemModelFaces
from Base.QtPage import ComboModel, DoubleValidator

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SolutionDomainView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Spin box delegate for 'Number' in Meshes table
#-------------------------------------------------------------------------------

class MeshNumberDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(MeshNumberDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QSpinBox(parent)
        editor.setMinimum(1)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, spinBox, index):
        value, ok = index.model().data(index, Qt.DisplayRole).toInt()
        spinBox.setValue(value)


    def setModelData(self, spinBox, model, index):
        spinBox.interpretText()
        value = spinBox.value()
        model.setData(index, QVariant(value))
#        selectionModel = self.parent.selectionModel()
#        for idx in selectionModel.selectedIndexes():
#            if idx.column() == index.column():
#                model.setData(idx, QVariant(value))


#-------------------------------------------------------------------------------
# QComboBox delegate for 'Group cells' and 'Group Faces' in Meshes table
#-------------------------------------------------------------------------------


class GroupDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent = None):
        super(GroupDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem(QString("off"))
        editor.addItem(QString("section"))
        editor.addItem(QString("zone"))
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataMeshes[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, QVariant(value))


#-------------------------------------------------------------------------------
# StandarItemModelMeshes class
#-------------------------------------------------------------------------------

class StandardItemModelMeshes(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = mdl
        self.dataMeshes = []
        # list of items to be disabled in the QTableView
        self.disabledItem = []

        self.populateModel()

        self.headers = [self.tr("Meshes"),
                        self.tr("Format"),
                        self.tr("Number"),
                        self.tr("Add groups of faces"),
                        self.tr("Add groups of cells")]

        self.tooltip = [self.tr("Code_Saturne Preprocessor option: --mesh"),
                        self.tr("Code_Saturne Preprocessor sub-option: --format"),
                        self.tr("Code_Saturne Preprocessor sub-option: --num"),
                        self.tr("Code_Saturne Preprocessor sub-option: --grp-fac"),
                        self.tr("Code_Saturne Preprocessor sub-option: --grp-cel")]

        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                role = Qt.DisplayRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)


    def populateModel(self):

        for mesh in self.mdl.getMeshList():
            format = self.mdl.getMeshFormat(mesh)
            list   = []
            list.append(mesh)
            list.append(MeshModel().getMeshFormatDescription(format))
            if format == 'med':
                num = self.mdl.getMeshNumber(mesh)
                if not num:
                    num = 1
                list.append(num)
            else:
                list.append("")
            if format == 'cgns':
                list.append(self.mdl.getMeshGroupFaces(mesh))
                list.append(self.mdl.getMeshGroupCells(mesh))
            else:
                list.append("")
                list.append("")
            self.dataMeshes.append(list)
            log.debug("populateModel-> dataMeshes = %s" % list)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.ToolTipRole:
            return QVariant(self.tooltip[index.column()])

        elif role == Qt.DisplayRole:
            d = self.dataMeshes[index.row()][index.column()]
            if d:
                return QVariant(d)
            else:
                return QVariant()

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        self.__disableData(index.row())

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        if index.column() in [0, 1]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() in [2, 3, 4]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()

        mesh = self.dataMeshes[row][0]

        if col == 2:
            v, ok = value.toInt()
            if ok:
                self.dataMeshes[row][2] = v
                self.mdl.setMeshNumber(mesh, v)

        elif col == 3:
            v = str(value.toString())
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshGroupFaces(mesh, v)

        elif col == 4:
            v = str(value.toString())
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshGroupCells(mesh, v)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addRow(self, name, format):
        """
        Add a row in the table.
        """
        txt_format = MeshModel().getMeshFormatDescription(format)

        if format == 'med':
            item = [name, txt_format, 1, "", ""]
        elif format == 'cgns':
            item = [name, txt_format, "", "off", "off"]
        else:
            item = [name, txt_format, "", "", ""]

        self.dataMeshes.append(item)

        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataMeshes[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def __disableData(self, row):
        if self.mdl.getMeshFormat(self.dataMeshes[row][0]) != "med":
            if (row, 2) not in self.disabledItem:
                self.disabledItem.append((row, 2))
        else:
            if (row, 2) in self.disabledItem:
                self.disabledItem.remove((row, 2))

        if self.mdl.getMeshFormat(self.dataMeshes[row][0]) != "cgns":
            if (row, 3) not in self.disabledItem:
                self.disabledItem.append((row, 3))
            if (row, 4) not in self.disabledItem:
                self.disabledItem.append((row, 4))
        else:
            if (row, 3) in self.disabledItem:
                self.disabledItem.remove((row, 3))
            if (row, 4) in self.disabledItem:
                self.disabledItem.remove((row, 4))


    def getName(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataMeshes[row][0]

#-------------------------------------------------------------------------------
# StandarItemModelPeriod class
#-------------------------------------------------------------------------------

class StandardItemModelPeriod(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = mdl
        self.headers = [ self.tr("Number"), self.tr("Name")]
        self.tooltip = [self.tr("Code_Saturne Preprocessor option: --perio"),
                        self.tr("Code_Saturne Preprocessor option: --perio")]

        self.setColumnCount(len(self.headers))
        self.dataPeriod = []

        self.populateModel()


    def populateModel(self):
        l = self.mdl.getPeriodicityListName()
        for i in range(self.mdl.getPeriodicityNumber()):
            list   = []
            list.append(i+1)
            list.append(l[i])
            self.dataPeriod.append(list)
            log.debug("populateModel-> dataPeriod = %s" % list)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataPeriod[index.row()][index.column()])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        if index.column() == 1:
            row = index.row()
            new_name = str(value.toString())
            if new_name:
                self.mdl.changePeriodicityName(self.dataPeriod[row][1], new_name)
                self.dataPeriod[row][1] = new_name

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addRow(self):
        """
        Add a row in the model.
        """
        row = self.rowCount()
        self.mdl.addPeriodicity("perio" + str(row+1))
        item = [row+1, "perio" + str(row+1)]
        self.dataPeriod.append(item)
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model.
        """
        self.mdl.deletePeriodicity(self.dataPeriod[row][1])
        del self.dataPeriod[row]
        row = self.rowCount()
        self.setRowCount(row-1)

#-------------------------------------------------------------------------------
# Popup window: Mesh Format dialog window
#-------------------------------------------------------------------------------

class SolutionDomainMeshFormatDialogView(QDialog, Ui_SolutionDomainMeshFormatDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)
        Ui_SolutionDomainMeshFormatDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Mesh format selection"))
        self.default = default
        self.result  = self.default

        # Combo models and items

        self.dico_txt_format = {}
        liste = MeshModel().getBuildFormatList()
        for (num, v, txt) in liste:
            self.dico_txt_format[txt] = v
            self.comboBox.addItem(QString(txt))

        # Activate the default format: "NOPO Simail (.des)"

        self.comboBox.setCurrentIndex(2)


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        txt_format = str(self.comboBox.currentText())
        self.result = self.dico_txt_format[txt_format]
        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class SolutionDomainView(QWidget, Ui_SolutionDomainForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_SolutionDomainForm.__init__(self)
        self.setupUi(self)

        self.stbar = stbar
        self.case = case
        self.mdl = SolutionDomainModel(self.case)

        # 1) MESHES SELECTION LAYOUT

        # 1.1) Model for meshes table
        self.modelMeshes = StandardItemModelMeshes(self.mdl)
        self.tableViewMeshes.setModel(self.modelMeshes)
        self.tableViewMeshes.resizeColumnsToContents()
        self.tableViewMeshes.resizeRowsToContents()

        delegateNumber = MeshNumberDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(2, delegateNumber)

        delegateGroupFaces = GroupDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(3, delegateGroupFaces)

        delegateGroupCells = GroupDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(4, delegateGroupCells)

        # 1.2) Connections

        self.connect(self.pushButtonAddMesh, SIGNAL("clicked()"), self.slotSearchMesh)
        self.connect(self.pushButtonDeleteMesh, SIGNAL("clicked()"), self.slotDeleteMesh)
        self.connect(self.groupBoxJoin, SIGNAL("clicked(bool)"), self.slotJoinMeshes)
        self.connect(self.groupBoxWarp, SIGNAL("clicked(bool)"), self.slotFacesCutting)
        self.connect(self.lineEditWarp, SIGNAL("textChanged(const QString &)"), self.slotWarpParam)
        self.connect(self.groupBoxOrient, SIGNAL("clicked(bool)"), self.slotOrientation)

        # 1.3) Set up validators
        validatorWarp = DoubleValidator(self.lineEditWarp, min=0.0)
        self.lineEditWarp.setValidator(validatorWarp)

        # 1.4) Faces to join selection (Custom Widgets)

        model = StandardItemModelFaces(self, self.mdl, 'faces_join')
        #self.widgetFacesJoin.tableView.reset()
        self.widgetFacesJoin.modelFaces = model
        self.widgetFacesJoin.tableView.setModel(model)

        # 2) PERIODICITIES

        self.perio_name = ""
        self.perio_mode = ""

        # Model for periodicities
        self.modelPeriod = StandardItemModelPeriod(self.mdl)
        self.treeViewPeriod.setModel(self.modelPeriod)
        self.treeViewPeriod.setSelectionBehavior(QAbstractItemView.SelectRows)

        # Combo model for type of periodicity
        self.modelComboPeriod = ComboModel(self.comboBoxPeriodicity, 5, 1)
        self.modelComboPeriod.addItem(self.tr("Periodicity of translation"), "translation")
        self.modelComboPeriod.addItem(self.tr("Periodicity of rotation (defined by angle and direction)"), "rotation1")
        self.modelComboPeriod.addItem(self.tr("Periodicity of rotation (defined by matrix)"), "rotation2")
        self.modelComboPeriod.addItem(self.tr("Translation + rotation (defined by angle and direction)"), "tr+rota1")
        self.modelComboPeriod.addItem(self.tr("Translation + rotation (defined by matrix)"), "tr+rota2")

        # Display
        self.groupBoxMode.hide()
        self.groupBoxFaces.hide()
        self.groupBoxTranslation.hide()
        self.groupBoxRotation1.hide()
        self.groupBoxRotation2.hide()

        # Set up validators

        # 2)
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
        self.lineEditM21.setValidator(DoubleValidator(self.lineEditM21))
        self.lineEditM22.setValidator(DoubleValidator(self.lineEditM22))
        self.lineEditM23.setValidator(DoubleValidator(self.lineEditM23))
        self.lineEditM31.setValidator(DoubleValidator(self.lineEditM31))
        self.lineEditM32.setValidator(DoubleValidator(self.lineEditM32))
        self.lineEditM33.setValidator(DoubleValidator(self.lineEditM33))
        self.lineEditX2.setValidator(DoubleValidator(self.lineEditX2))
        self.lineEditY2.setValidator(DoubleValidator(self.lineEditY2))
        self.lineEditZ2.setValidator(DoubleValidator(self.lineEditZ2))

        # Connections

        self.connect(self.treeViewPeriod, SIGNAL("clicked(const QModelIndex &)"), self.slotUpdatePeriodicity)
        self.connect(self.modelPeriod, SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.slotUpdateByModel)
        self.connect(self.pushButtonAddPeriod,    SIGNAL("clicked()"), self.slotAddPeriodicity)
        self.connect(self.pushButtonDeletePeriod, SIGNAL("clicked()"), self.slotDeletePeriodicity)

        self.connect(self.comboBoxPeriodicity, SIGNAL("activated(const QString&)"), self.slotPeriodicityMode)

        self.connect(self.lineEditTX, SIGNAL("textChanged(const QString &)"), self.slotTranslationX)
        self.connect(self.lineEditTY, SIGNAL("textChanged(const QString &)"), self.slotTranslationY)
        self.connect(self.lineEditTZ, SIGNAL("textChanged(const QString &)"), self.slotTranslationZ)

        self.connect(self.lineEditAngle, SIGNAL("textChanged(const QString &)"), self.slotAngleRotation)

        self.connect(self.lineEditDX, SIGNAL("textChanged(const QString &)"), self.slotRotationX)
        self.connect(self.lineEditDY, SIGNAL("textChanged(const QString &)"), self.slotRotationY)
        self.connect(self.lineEditDZ, SIGNAL("textChanged(const QString &)"), self.slotRotationZ)

        self.connect(self.lineEditX1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationX1)
        self.connect(self.lineEditY1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationY1)
        self.connect(self.lineEditZ1, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationZ1)

        self.connect(self.lineEditM11, SIGNAL("textChanged(const QString &)"), self.slotMatrix11)
        self.connect(self.lineEditM12, SIGNAL("textChanged(const QString &)"), self.slotMatrix12)
        self.connect(self.lineEditM13, SIGNAL("textChanged(const QString &)"), self.slotMatrix13)
        self.connect(self.lineEditM21, SIGNAL("textChanged(const QString &)"), self.slotMatrix21)
        self.connect(self.lineEditM22, SIGNAL("textChanged(const QString &)"), self.slotMatrix22)
        self.connect(self.lineEditM23, SIGNAL("textChanged(const QString &)"), self.slotMatrix22)
        self.connect(self.lineEditM31, SIGNAL("textChanged(const QString &)"), self.slotMatrix31)
        self.connect(self.lineEditM32, SIGNAL("textChanged(const QString &)"), self.slotMatrix32)
        self.connect(self.lineEditM33, SIGNAL("textChanged(const QString &)"), self.slotMatrix33)

        self.connect(self.lineEditX2, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationX2)
        self.connect(self.lineEditY2, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationY2)
        self.connect(self.lineEditZ2, SIGNAL("textChanged(const QString &)"), self.slotCenterRotationZ2)

        # 3) INITIALIZE MESHES LIST

        # 3.1) Meshes list

        msg = ""
        nameList = self.mdl.getMeshList()
        log.debug("__init__ -> nameList = %s " % nameList)

        if nameList:
            for mesh in nameList:
                if mesh not in os.listdir(self.case['mesh_path']):
                    msg = msg  + mesh + '\n'

            if msg:
                msg =  msg + '\n'
                title = self.tr("WARNING")
                msg2  = self.tr("The following meshes are not in the meshes "  \
                                "directory given in the 'Identity and paths' " \
                                "section:\n\n" +
                                msg +
                                "Verify existence and location of the mesh " \
                                "files, and the 'Identity and Paths' section." )
                QMessageBox.warning(self, title, msg2)

#            for mesh in nameList:
#                format = self.mdl.getMeshFormat(mesh)
#                txt_format = self.dico_format[format]
#                self.modelMeshes.addRow(mesh, txt_format)

        elif os.listdir(self.case['mesh_path']):
            for mesh in os.listdir(self.case['mesh_path']):
                if MeshModel().getMeshExtension(mesh):
                    format = self._meshFormatFromFileName(mesh)
                    self.modelMeshes.addRow(mesh, format)
                    self.mdl.addMesh(mesh, format)

        self._tableViewLayout()


        # 3.2) Join parameters

        if self.mdl.getJoinMeshesStatus() == 'on':
            self.groupBoxJoin.setChecked(True)
            self.slotJoinMeshes(True)
        else:
            self.groupBoxJoin.setChecked(False)
            self.slotJoinMeshes(False)

        # 3.3) Warped faces cutting

        if self.mdl.getCutStatus() == 'on':
            self.groupBoxWarp.setChecked(True)
            self.slotFacesCutting(True)
        else:
            self.groupBoxWarp.setChecked(False)
            self.slotFacesCutting(False)

        v = self.mdl.getCutAngle()
        self.warp = v
        self.lineEditWarp.setText(str(self.warp))

        # 3.4) Reorientation

        if self.mdl.getOrientation() == 'on':
            self.groupBoxOrient.setChecked(True)
        else:
            self.groupBoxOrient.setChecked(False)


    def setLink(self, fullfile):
        """
        Tab1: set link between asked file and case's directory 'MESH'
        and keep dependency between file in a directory if destruction
        """
        name = os.path.basename(fullfile)
        new_full = self.case['mesh_path'] + "/" + name
        os.symlink(fullfile, new_full)


    def _meshFormatFromFileName(self, mesh):
        """
        Tab1: return the format and the listbox name fo the given elementNode.
        Put the format in attribute of the elementNode if necessary.
        """
        # Search the format of the mesh
        format = MeshModel().getMeshFormat(mesh)

        # If the format is not found ask to the user
        if not format:
            default = ""
            dialog = SolutionDomainMeshFormatDialogView(self, default)
            if dialog.exec_():
                format = dialog.get_result()

        return format


    def _tableViewLayout(self):
        """
        Configure QTableView column number
        """
        format_list = []
        for mesh in self.mdl.getMeshList():
            format_list.append(self.mdl.getMeshFormat(mesh))

        self.tableViewMeshes.resizeColumnsToContents()
        self.tableViewMeshes.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)

        if "med" in format_list:
            self.tableViewMeshes.setColumnHidden(2, False)
            self.tableViewMeshes.horizontalHeader().setResizeMode(2, QHeaderView.ResizeToContents)
        else:
            self.tableViewMeshes.setColumnHidden(2, True)

        if "cgns" in format_list:
            self.tableViewMeshes.setColumnHidden(3, False)
            self.tableViewMeshes.setColumnHidden(4, False)
            self.tableViewMeshes.horizontalHeader().setResizeMode(3, QHeaderView.ResizeToContents)
            self.tableViewMeshes.horizontalHeader().setResizeMode(4, QHeaderView.ResizeToContents)
        else:
            self.tableViewMeshes.setColumnHidden(3, True)
            self.tableViewMeshes.setColumnHidden(4, True)

        self.tableViewMeshes.horizontalHeader().setResizeMode(0, QHeaderView.Stretch)
        self.tableViewMeshes.horizontalHeader().setResizeMode(1, QHeaderView.Stretch)


    def _addMeshInList(self,  m):
        """
        Tab1: add input new meshes in the listbox and case.
        """
        input_meshes = string.split(m)
        if not input_meshes: return

        # 1) Verify if the all new meshes are the mesh_path directory.

        for mesh in input_meshes:

            if mesh not in os.listdir(self.case['mesh_path']):
                title = self.tr("Warning")
                msg   = self.tr("The new mesh entry is not in meshes directory " \
                                "given in the 'Identity and paths' section."  \
                                "One link is created in this directory.\n\n"  \
                                "Verify existence and location of the meshes " \
                                "files, and the 'Identity and Pathes' section." )
                QMessageBox.information(self, title, msg)

        # 2) Verify that the new mesh is not allready in the case

        for mesh in input_meshes:

            if mesh in self.mdl.getMeshList():
                title = self.tr("Warning")
                msg   = self.tr("Warning, the following input is already " \
                                "uploaded in the list:\n\n" + mesh)
                QMessageBox.information(self, title, msg)

            else:

                # 3) Update View and model

                format = self._meshFormatFromFileName(mesh)
                if format:
                    self.modelMeshes.addRow(mesh, format)
                    self.mdl.addMesh(mesh, format)

        self._tableViewLayout()


    @pyqtSignature("")
    def slotSearchMesh(self):
        msg = self.tr("Select a mesh file (the file must be " \
                      "in the directory %s)."%self.case['mesh_path'])
        self.stbar.showMessage(msg, 2000)

        title = self.tr("New mesh")
        if self.case['mesh_path'] and os.path.isdir(self.case['mesh_path']):
            path = self.case['mesh_path']
        else:
            path = os.getcwd()

        filetypes = ""
        for Format in MeshModel().getFileFormatList():
            filetypes += "%s (%s);;"%(Format[0], Format[1])

        filt = "All files (*)"
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes, filt)
        file_name = str(file_name)

        if file_name:
            if os.path.dirname(file_name) != self.case['mesh_path']:
                self.setLink(file_name)
                msg = self.tr("Select a mesh file (the file must be " \
                              "in the directory %s)."%self.case['mesh_path'])
                self.stbar.showMessage(msg, 2000)
            m = os.path.basename(file_name)
            self._addMeshInList(m)


    @pyqtSignature("")
    def slotDeleteMesh(self):
        """
        Delete the selected mesh from the list
        """
        # 1) Is there a mesh to delete ?

        selectionModel = self.tableViewMeshes.selectionModel()
        for index in selectionModel.selectedRows():
            row = index.row()
            mesh = self.modelMeshes.getName(row)

            # 2) Delete mesh from view and from model

            self.modelMeshes.deleteRow(row)
            self.mdl.delMesh(mesh)

        self._tableViewLayout()


    @pyqtSignature("bool")
    def slotJoinMeshes(self, checked):
        """
        Private slot.

        Do we join any meshes ?

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox join parameters
        """
        self.groupBoxJoin.setFlat(not checked)
        if checked:
            self.mdl.setJoinMeshesStatus("on")
            self.frameJoin.show()
        else:
            self.mdl.setJoinMeshesStatus("off")
            self.frameJoin.hide()


    @pyqtSignature("bool")
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


    @pyqtSignature("bool")
    def slotOrientation(self, checked):
        """
        Private slot.

        Do we cut any warp faces ?

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox warp parameters
        """
        #self.groupBoxOrient.setFlat(not checked)
        if checked:
            self.mdl.setOrientation ("on")
        else:
            self.mdl.setOrientation("off")


    @pyqtSignature("const QString&")
    def slotWarpParam(self, text):
        """
        Private slot.

        Input '--cwf' parameters.

        @type text: C{QString}
        @param text: max angle of warped faces
        """
        var, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setCutAngle(var)


    @pyqtSignature("")
    def slotAddPeriodicity(self):
        """
        Add a periodicity in the list.
        """
        self.modelPeriod.addRow()


    @pyqtSignature("")
    def slotDeletePeriodicity(self):
        """
        Delete a periodicity from the list.
        """
        self.groupBoxMode.hide()
        self.groupBoxFaces.hide()
        self.groupBoxTranslation.hide()
        self.groupBoxRotation1.hide()
        self.groupBoxRotation2.hide()

        for index in self.treeViewPeriod.selectionModel().selectedIndexes():
            self.modelPeriod.deleteRow(index.row())
            break


    def __setValuesTranslation(self, perio_name):
        """
        Put values found in xml file as soon as mode is "translation"
        """
        dx, dy, dz = self.mdl.getTranslationDirection(perio_name)

        self.lineEditTX.setText(QString(dx))
        self.lineEditTY.setText(QString(dy))
        self.lineEditTZ.setText(QString(dz))


    def __setValuesRotation1(self, perio_name):
        """
        Put values found in xml file as soon as mode is "rotation1"
        """
        angle = self.mdl.getRotationAngle(perio_name)
        rx, ry, rz = self.mdl.getRotationDirection(perio_name)
        px, py, pz = self.mdl.getRotationCenter(perio_name)

        self.lineEditAngle.setText(QString((angle)))
        self.lineEditDX.setText(QString(rx))
        self.lineEditDY.setText(QString(ry))
        self.lineEditDZ.setText(QString(rz))
        self.lineEditX1.setText(QString(px))
        self.lineEditY1.setText(QString(py))
        self.lineEditZ1.setText(QString(pz))


    def __setValuesRotation2(self, perio_name):
        """
        Put values found in xml file as soon as mode is "rotation"2
        """
        m11,m12,m13,m21,m22,m23,m31,m32,m33 = self.mdl.getRotationMatrix(perio_name)
        px, py, pz = self.mdl.getRotationCenter(perio_name)

        self.lineEditM11.setText(QString(m11))
        self.lineEditM12.setText(QString(m12))
        self.lineEditM13.setText(QString(m13))
        self.lineEditM21.setText(QString(m21))
        self.lineEditM22.setText(QString(m22))
        self.lineEditM23.setText(QString(m23))
        self.lineEditM31.setText(QString(m31))
        self.lineEditM32.setText(QString(m32))
        self.lineEditM33.setText(QString(m33))
        self.lineEditX2.setText(QString(px))
        self.lineEditY2.setText(QString(py))
        self.lineEditZ2.setText(QString(pz))


    def __setValuesPeriodicTransformation(self, perio, mode):
        """
        Put values found in xml file as soon as mode of
        transformation is choosen
        """
        log.debug("__setValuesPeriodicTransformation perio mode = %s %s "% (perio, mode))
        if mode == "translation" :
            self.__setValuesTranslation(perio)
        if mode == "rotation1":
            self.__setValuesRotation1(perio)
        if mode == "rotation2":
            self.__setValuesRotation2(perio)
        if mode == "tr+rota1":
            self.__setValuesTranslation(perio)
            self.__setValuesRotation1(perio)
        if mode == "tr+rota2":
            self.__setValuesTranslation(perio)
            self.__setValuesRotation2(perio)


    @pyqtSignature("const QModelIndex&, const QModelIndex&")
    def slotUpdateByModel(self, index1, index2):
        """
        This slot update the display for the periodicity modified
        in the tree view.
        """
        log.debug("slotUpdateByModel index.row() = %i " % index1.row())
        self.slotUpdatePeriodicity(index1)


    @pyqtSignature("const QModelIndex&")
    def slotUpdatePeriodicity(self, index):
        """
        This slot update the display for the periodicity selected
        in the tree view.
        """
        log.debug("slotUpdatePeriodicity index.row() = %i " % index.row())

        self.groupBoxMode.show()
        self.groupBoxFaces.show()

        perio_name = self.mdl.getPeriodicityListName()[index.row()]
        perio_mode = self.mdl.getPeriodicityMode(perio_name)
        self.perio_name = perio_name
        self.perio_mode = perio_mode

        self.modelComboPeriod.setItem(str_model=perio_mode)
        txt = str(self.comboBoxPeriodicity.currentText())
        self.slotPeriodicityMode(txt)

        # Update FaceSelection
#        self.widgetFacesPeriodic.setPeriodicityName(perio_name)
#        self.widgetFacesPeriodic.populateModel(self.mdl, 'faces_periodic')
        model = StandardItemModelFaces(self, self.mdl, 'faces_periodic', perio_name)
        self.widgetFacesPeriodic.modelFaces = model
        self.widgetFacesPeriodic.tableView.setModel(model)

    @pyqtSignature("const QString&")
    def slotPeriodicityMode(self, text):
        """
        Do we have a periodicity ?
        """
        self.perio_mode = self.modelComboPeriod.dicoV2M[str(text)]
        #perio_number = self.treeViewPeriod.currentIndex().row()
        #self.perio_number = str(perio_number+1)

        log.debug("slotPeriodicityMode  = %s " % self.perio_mode)

        self.mdl.updatePeriodicityMode(self.perio_name, self.perio_mode)

        # exit if no periodicity
#        if self.perio_number == "0":
#            log.debug("slotPeriodicityMode -> no periodicity ")
#            return

        self.groupBoxTranslation.hide()
        self.groupBoxRotation1.hide()
        self.groupBoxRotation2.hide()

        if self.perio_mode == "":
            self.modelComboPeriod(str_model='translation')
            self.__setValuesTranslation(self.perio_name)

        if self.perio_mode == "translation":
            self.groupBoxTranslation.show()
            self.groupBoxRotation1.hide()
            self.groupBoxRotation2.hide()
            self.__setValuesTranslation(self.perio_name)

        elif self.perio_mode == "rotation1":
            self.groupBoxTranslation.hide()
            self.groupBoxRotation1.show()
            self.groupBoxRotation2.hide()
            self.__setValuesRotation1(self.perio_name)

        elif self.perio_mode == "rotation2":
            self.groupBoxTranslation.hide()
            self.groupBoxRotation1.hide()
            self.groupBoxRotation2.show()
            self.__setValuesRotation2(self.perio_name)

        elif self.perio_mode == "tr+rota1":
            self.groupBoxTranslation.show()
            self.groupBoxRotation1.show()
            self.groupBoxRotation2.hide()
            self.__setValuesTranslation(self.perio_name)
            self.__setValuesRotation1(self.perio_name)

        elif self.perio_mode == "tr+rota2":
            self.groupBoxTranslation.show()
            self.groupBoxRotation1.hide()
            self.groupBoxRotation2.show()
            self.__setValuesTranslation(self.perio_name)
            self.__setValuesRotation2(self.perio_name)

#        self.mdl.addPeriodicity(self.perio_name, self.perio_mode)
        self.__setValuesPeriodicTransformation(self.perio_name, self.perio_mode)


    @pyqtSignature("const QString&")
    def slotTranslationX(self, text):
        """
        Periodicity translation for X
        """
        if self.perio_mode != "rotation1" or self.perio_mode != "rotation2":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setTranslationDirection(self.perio_name, 'translation_x', val)


    @pyqtSignature("const QString&")
    def slotTranslationY(self, text):
        """
        Periodicity translation for Y
        """
        if self.perio_mode != "rotation1" or self.perio_mode != "rotation2":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setTranslationDirection(self.perio_name, 'translation_y', val)


    @pyqtSignature("const QString&")
    def slotTranslationZ(self, text):
        """
        Periodicity translation for Z
        """
        if self.perio_mode != "rotation1" or self.perio_mode != "rotation2":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setTranslationDirection(self.perio_name, 'translation_z', val)


    @pyqtSignature("const QString&")
    def slotAngleRotation(self, text):
        """
        Periodicity rotation angle
        """
        if self.perio_mode == "rotation1" or self.perio_mode == "tr+rota1":
            angle, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationAngle(self.perio_name, angle)


    @pyqtSignature("const QString&")
    def slotRotationX(self, text):
        """
        Periodicity rotation for X
        """
        if self.perio_mode == "rotation1" or self.perio_mode =="tr+rota1":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationVector(self.perio_name, "rotation_x", val)


    @pyqtSignature("const QString&")
    def slotRotationY(self, text):
        """
        Periodicity rotation for Y
        """
        if self.perio_mode == "rotation1" or self.perio_mode =="tr+rota1":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationVector(self.perio_name, "rotation_y", val)


    @pyqtSignature("const QString&")
    def slotRotationZ(self, text):
        """
        Periodicity rotation for Z
        """
        if self.perio_mode == "rotation1" or self.perio_mode =="tr+rota1":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationVector(self.perio_name, "rotation_z", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationX1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_x", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationY1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_y", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationZ1(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_z", val)


    # Methods for matrix components
    @pyqtSignature("const QString&")
    def slotMatrix11(self, text):
        self.__cmdMatrixRotation("rotation_matrix_11", text)


    @pyqtSignature("const QString&")
    def slotMatrix12(self, text):
        self.__cmdMatrixRotation("rotation_matrix_12", text)


    @pyqtSignature("const QString&")
    def slotMatrix13(self, text):
        self.__cmdMatrixRotation("rotation_matrix_13", text)


    @pyqtSignature("const QString&")
    def slotMatrix21(self, text):
        self.__cmdMatrixRotation("rotation_matrix_21", text)


    @pyqtSignature("const QString&")
    def slotMatrix22(self, text):
        self.__cmdMatrixRotation("rotation_matrix_22", text)


    @pyqtSignature("const QString&")
    def slotMatrix23(self, text):
        self.__cmdMatrixRotation("rotation_matrix_23", text)


    @pyqtSignature("const QString&")
    def slotMatrix31(self, text):
        self.__cmdMatrixRotation("rotation_matrix_31", text)


    @pyqtSignature("const QString&")
    def slotMatrix32(self, text):
        self.__cmdMatrixRotation("rotation_matrix_32", text)


    @pyqtSignature("const QString&")
    def slotMatrix33(self, text):
        self.__cmdMatrixRotation("rotation_matrix_33", text)


    def __cmdMatrixRotation(self, pos, text):
        """
        Periodicity translation
        """
        if self.perio_mode == "rotation2" or self.perio_mode == "tr+rota2":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationMatrix(self.perio_name, pos, val)


    @pyqtSignature("const QString&")
    def slotCenterRotationX2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_x", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationY2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_y", val)


    @pyqtSignature("const QString&")
    def slotCenterRotationZ2(self, text):
        """
        Periodicity : center of rotation
        """
        if self.perio_mode != "translation":
            val, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                self.mdl.setRotationCenter(self.perio_name, "rotation_center_z", val)


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
