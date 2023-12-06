# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
- StandardItemModelMeshes
- SolutionDomainView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, logging
import configparser

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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator, RegExpValidator, IntValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.SolutionDomainForm import Ui_SolutionDomainForm
from code_saturne.model.SolutionDomainModel import RelOrAbsPath, MeshModel, SolutionDomainModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SolutionDomainView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Label delegate for 'Name' in Meshes table
#-------------------------------------------------------------------------------

class MeshNameDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(MeshNameDelegate, self).__init__(parent)
        self.parent = parent

    def paint(self, painter, option, index):
        row = index.row()
        isValid = index.model().dataMeshes[row][8]

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
            if value != None:
                text = from_qvariant(value, to_text_string)
                painter.drawText(option.rect, Qt.AlignLeft, text)
            painter.restore()


#-------------------------------------------------------------------------------
# QComboBox delegate for Format in Meshes table
#-------------------------------------------------------------------------------

class MeshFormatDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent = None, updateLayout = None):
        super(MeshFormatDelegate, self).__init__(parent)
        self.parent = parent
        self.updateLayout = updateLayout
        self.lst = MeshModel().getBuildFormatList()
        # Compute width based on longest possible string and font metrics
        fm = self.parent.fontMetrics()
        self.textSize = fm.size(Qt.TextSingleLine, 'I-deas universal')
        self.textSize.setHeight(1)
        for i in range(len(self.lst)):
            w = fm.size(Qt.TextSingleLine, str(self.lst[i][1])).width()
            if w > self.textSize.width():
                self.textSize.setWidth(w)


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for i in range(len(self.lst)):
            fmt = self.lst[i]
            editor.addItem(str(fmt[1] + fmt[2]))
        return editor


    def setEditorData(self, comboBox, index):
        key = index.model().dataMeshes[index.row()][1]
        string = ''
        for i in range(len(self.lst)):
            if key == self.lst[i][0]:
                comboBox.setCurrentIndex(i)


    def setModelData(self, comboBox, model, index):
        value = str(comboBox.currentText())
        key = ''
        for i in range(len(self.lst)):
            if value == self.lst[i][1] + self.lst[i][2]:
                key = self.lst[i][0]
        model.setData(index, key)
        if self.updateLayout != None:
            self.updateLayout()


    def sizeHint(self, option, index):
        return self.textSize


    def paint(self, painter, option, index):
        row = index.row()
        format = index.model().dataMeshes[row][1]
        isValid = format != None and format != ''

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
            if value != None:
                if value.isValid():
                    text = from_qvariant(value, to_text_string)
                    painter.drawText(option.rect, Qt.AlignLeft, text)
            painter.restore()


#-------------------------------------------------------------------------------
# Line edit delegate for 'Number' in Meshes table
#-------------------------------------------------------------------------------

class MeshNumberDelegate(QItemDelegate):
    def __init__(self, parent = None):
        super(MeshNumberDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        vd = RegExpValidator(editor, QRegExp("[0-9- ]*"))
        editor.setValidator(vd)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, lineEdit, index):
        string = index.model().dataMeshes[index.row()][2]
        lineEdit.setText(string)


    def setModelData(self, lineEdit, model, index):
        value = str(lineEdit.text()).strip()
        model.setData(index, value)


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
        editor.addItem("off")
        editor.addItem("section")
        editor.addItem("zone")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataMeshes[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value)

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
# Line edit delegate for float (thickness and reason)
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

        lst = MeshModel().getBuildFormatList()
        self.formatDict = {'':''}
        for i in range(len(lst)):
            self.formatDict[lst[i][0]] = lst[i][1]

        self.populateModel()


        self.headers = [self.tr("File name"),
                        self.tr("Format"),
                        self.tr("Numbers"),
                        self.tr("Reorient"),
                        self.tr("Fix"),
                        self.tr("Add face groups"),
                        self.tr("Add cell groups"),
                        self.tr("Path")]

        self.tooltip = [self.tr("Base mesh file name"),
                        self.tr("Mesh format"),
                        self.tr("Mesh number if multiple meshes are present in the same file"),
                        self.tr("Try to reorient badly oriented cells"),
                        self.tr("Discard cells with negative volumes or incorrect topology"),
                        self.tr("Add groups to faces"),
                        self.tr("Add groups to cells"),
                        self.tr("File path")]

        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                role = Qt.DisplayRole
                index = self.index(row, column)
                value = self.data(index, role)
                if column != 1:
                    self.setData(index, value)
                else:
                    self.setData(index, self.dataMeshes[row][1])


    def populateModel(self):

        for mesh in self.mdl.getMeshList():
            format = self.mdl.getMeshFormat(mesh)
            lst   = []
            lst.append(mesh[0])
            lst.append(format)
            if format in ['cgns', 'med', 'ensight']:
                num = self.mdl.getMeshNumbers(mesh)
                if not num:
                    num = ''
                lst.append(num)
            else:
                lst.append("")
            lst.append(self.mdl.getMeshReorient(mesh))
            lst.append(self.mdl.getMeshDiscardBadCells(mesh))
            if format == 'cgns':
                lst.append(self.mdl.getMeshGroupFaces(mesh))
                lst.append(self.mdl.getMeshGroupCells(mesh))
            else:
                lst.append("")
                lst.append("")
            lst.append(mesh[1])
            lst.append(self.__isMeshPathValid(mesh))

            self.dataMeshes.append(lst)
            log.debug("populateModel-> dataMeshes = %s" % lst)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return None

        col = index.column()

        if role == Qt.ToolTipRole:
            return self.tooltip[col]

        elif role == Qt.DisplayRole:
            d = self.dataMeshes[index.row()][col]
            if d:
                if col == 1:
                    return self.formatDict[d]
                elif col in (3, 4):
                    return None
                else:
                    return d
            else:
                return None

        elif role == Qt.TextAlignmentRole:
            if col == 7:
                return Qt.AlignLeft
            else:
                return Qt.AlignCenter

        elif role == Qt.CheckStateRole:
            if col in (3, 4):
                d = self.dataMeshes[index.row()][col]
                if d == True:
                    return Qt.Checked
                else:
                    return Qt.Unchecked
            else:
                return None

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        self.__disableData(index.row())

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        if index.column() in (0, 7):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() in (1, 2, 5, 6):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() in (3, 4):
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()


        mesh = (self.dataMeshes[row][0], self.dataMeshes[row][7])

        if col == 1:
            v = from_qvariant(value, to_text_string)
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshFormat(mesh, v)

        elif col == 2:
            v = from_qvariant(value, to_text_string)
            self.dataMeshes[row][col] = v
            if value:
                self.mdl.setMeshNumbers(mesh, v)

        elif col == 3 and role == Qt.CheckStateRole:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.dataMeshes[row][col] = False
            else:
                self.dataMeshes[row][col] = True
            self.mdl.setMeshReorient(mesh, self.dataMeshes[row][col])

        elif col == 4 and role == Qt.CheckStateRole:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.dataMeshes[row][col] = False
            else:
                self.dataMeshes[row][col] = True
            self.mdl.setMeshDiscardBadCells(mesh, self.dataMeshes[row][col])

        elif col == 5:
            if value:
                v = from_qvariant(value, to_text_string)
            else:
                v = ''
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshGroupFaces(mesh, v)

        elif col == 6:
            if value:
                v = from_qvariant(value, to_text_string)
            else:
                v = ''
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshGroupCells(mesh, v)

        self.dataChanged.emit(index, index)
        return True


    def addRow(self, mesh, format):
        """
        Add a row in the table.
        """
        if format == 'med' or format == 'ensight':
            item = [mesh[0], format, "1", "", "", "", "", mesh[1]]
        elif format == 'cgns':
            item = [mesh[0], format, "1", "", "", "off", "off", mesh[1]]
        else:
            item = [mesh[0], format, "", "", "", "", "", mesh[1]]
        item.append(self.__isMeshPathValid(mesh))

        self.dataMeshes.append(item)

        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        self.setRowCount(self.rowCount() - 1)
        del self.dataMeshes[row]


    def __disableData(self, row):

        mesh = (self.dataMeshes[row][0], self.dataMeshes[row][7])
        self.dataMeshes[row][8] = self.__isMeshPathValid(mesh)

        if not self.dataMeshes[row][1] in ["cgns", "med", "ensight"]:
            if (row, 2) not in self.disabledItem:
                self.disabledItem.append((row, 2))
        else:
            if (row, 2) in self.disabledItem:
                self.disabledItem.remove((row, 2))

        if self.dataMeshes[row][1] != "cgns":
            if (row, 5) not in self.disabledItem:
                self.disabledItem.append((row, 5))
            if (row, 6) not in self.disabledItem:
                self.disabledItem.append((row, 6))
        else:
            if (row, 5) in self.disabledItem:
                self.disabledItem.remove((row, 5))
            if (row, 6) in self.disabledItem:
                self.disabledItem.remove((row, 6))


    def __isMeshPathValid(self, mesh):
        """
        Public method. Check if a mesh named "mesh" matches an existing file
        """
        if mesh[1] != None:
            path = os.path.join(mesh[1], mesh[0])
        else:
            path = mesh[0]
        if not os.path.isabs(path) and self.mdl.case['mesh_path'] != None:
            path = os.path.join(self.mdl.case['mesh_path'], path)

        if not (os.path.isfile(path) and os.path.isabs(path)):
            return False

        return True


#-------------------------------------------------------------------------------
# File dialog to select either file or directory
#-------------------------------------------------------------------------------

class MeshInputDialog(QFileDialog):

    def __init__(self,
                 parent = None,
                 search_dirs = []):

        if len(search_dirs) == 0:
            directory = ""
        else:
            directory = str(search_dirs[0])

        try:
            QFileDialog.__init__(self,
                                 parent = parent,
                                 directory = directory)
        except (AttributeError, TypeError):
            QFileDialog.__init__(self)  # for older PyQt versions

        # Self.tr is only available once the parent class __init__ has been called,
        # so we may now set the caption, filter, and selection label

        caption =  self.tr("Select input mesh file or directory")
        self.setWindowTitle(caption)

        self.name_filter = str(self.tr("Imported or preprocessed meshes (mesh_input mesh_output *.csm)"))
        self.setNameFilter(self.name_filter)

        self.select_label = str(self.tr("Select"))
        self.setLabelText(QFileDialog.Accept, self.select_label)

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog
        if hasattr(self, 'setOptions'):
            self.setOptions(options)
        search_urls = []
        for d in search_dirs:
            search_urls.append(QUrl.fromLocalFile(d))
        self.setSidebarUrls(search_urls)

        self.currentChanged[str].connect(self._selectChange)

        self.setFileMode(QFileDialog.ExistingFile)


    def _selectChange(self, spath):
        mode = QFileDialog.ExistingFile
        path = str(spath)
        if os.path.basename(path) == 'mesh_input':
            if os.path.isdir(path):
                mode = QFileDialog.Directory
        self.setFileMode(mode)
        self.setNameFilter(self.name_filter)
        self.setLabelText(QFileDialog.Accept, self.select_label)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class SolutionDomainView(QWidget, Ui_SolutionDomainForm):
    """
    """
    def __init__(self, parent, case, stbar, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_SolutionDomainForm.__init__(self)
        self.setupUi(self)

        self.root = parent
        self.stbar = stbar
        self.case = case
        self.browser = tree

        self.case.undoStopGlobal()
        self.mdl = SolutionDomainModel(self.case)

        # 0) Mesh Input

        self.mesh_origin = self.mdl.getMeshOrigin()

        self.mesh_input = self.mdl.getMeshInput()

        if self.mesh_origin == "mesh_import":
            self.radioButtonImport.setChecked(True)
            self.radioButtonExists.setChecked(False)
            self.radioButtonCartesianMesh.setChecked(False)

            self.frameMeshInput.hide()
            self.frameMeshCartesian.hide()
            self.frameMeshImport.show()

        elif self.mesh_origin == "mesh_input":
            self.radioButtonImport.setChecked(False)
            self.radioButtonExists.setChecked(True)
            self.radioButtonCartesianMesh.setChecked(False)

            self.frameMeshImport.hide()
            self.frameMeshCartesian.hide()
            self.frameMeshInput.show()
            self.lineEditMeshInput.setText(self.mesh_input)

        elif self.mesh_origin == "mesh_cartesian":
            self.radioButtonImport.setChecked(False)
            self.radioButtonExists.setChecked(False)
            self.radioButtonCartesianMesh.setChecked(True)

            self.frameMeshInput.hide()
            self.frameMeshImport.hide()
            self.frameMeshCartesian.show()


        self.radioButtonImport.clicked.connect(self.slotSetImportMesh)
        self.radioButtonExists.clicked.connect(self.slotSetInputMesh)
        self.radioButtonCartesianMesh.clicked.connect(self.slotSetCartesianMesh)
        self.toolButtonMeshInput.pressed.connect(self.selectInputMesh)
        self.lineEditMeshInput.textChanged[str].connect(self.modifyInputMesh)

        self.cartParams = {"x_ncells":self.lineEditXNcells,
                           "y_ncells":self.lineEditYNcells,
                           "z_ncells":self.lineEditZNcells,
                           "x_min":self.lineEditXmin,
                           "x_max":self.lineEditXmax,
                           "y_min":self.lineEditYmin,
                           "y_max":self.lineEditYmax,
                           "z_min":self.lineEditZmin,
                           "z_max":self.lineEditZmax,
                           "x_prog":self.lineEditXprog,
                           "y_prog":self.lineEditYprog,
                           "z_prog":self.lineEditZprog}


        # We use lambda to connect all 9 lineEdits to the same function.
        # val corresponds to the string, and "key=k" is necessary for the connect
        # method to evaluate 'k' during creation. Otherwise last value of 'k' is
        # used for all.
        for k in self.cartParams.keys():
            self.cartParams[k].textChanged[str].connect(lambda val,key=k:
                    self.slotSetCartesianParam(val, key))
            # set validator
            if k.split("_")[1] == "ncells":
                _v = IntValidator(self.cartParams[k], min=1)
                self.cartParams[k].setValidator(_v)
            else:
                _v = DoubleValidator(self.cartParams[k])
                self.cartParams[k].setValidator(_v)

        self.cartParams["x_law"] = ComboModel(self.comboBoxXlaw, 3, 1)
        self.cartParams["y_law"] = ComboModel(self.comboBoxYlaw, 3, 1)
        self.cartParams["z_law"] = ComboModel(self.comboBoxZlaw, 3, 1)
        for k in ("x_law", "y_law", "z_law"):
            self.cartParams[k].addItem(self.tr("constant"), "constant")
            self.cartParams[k].addItem(self.tr("geometric"), "geometric")
            self.cartParams[k].addItem(self.tr("parabolic"), "parabolic")

        self.comboBoxXlaw.activated[str].connect(lambda val:
                self.slotSetCartesianParam(val, "x_law"))
        self.comboBoxYlaw.activated[str].connect(lambda val:
                self.slotSetCartesianParam(val, "y_law"))
        self.comboBoxZlaw.activated[str].connect(lambda val:
                self.slotSetCartesianParam(val, "z_law"))

        # Set initial values
        self.mesh_origin = self.mdl.getMeshOrigin()

        if self.mesh_origin == "mesh_cartesian":
            self.slotSetCartesianMesh()

        # 1) Meshes directory

        self.mesh_dirs = [None]

        d = self.mdl.getMeshDir()
        study_path = os.path.split(self.case['case_path'])[0]

        if d is None:
            d = os.path.join(study_path, 'MESH')
        elif not os.path.abspath(d):
            d = os.path.join(self.case['case_path'], d)

        if d != None:
            if os.path.isdir(d):
                self.lineEditMeshDir.setText(RelOrAbsPath(d, self.case['case_path']))
                self.mesh_dirs[0] = d

        self.case['mesh_path'] = self.mesh_dirs[0]

        package = self.case['package']

        # User and global mesh directories

        for config_file in [package.get_user_configfile(),
                            package.get_global_configfile()]:
            cfg = configparser.ConfigParser()
            cfg.read(config_file)
            if cfg.has_option('run', 'meshpath'):
                cfg_mesh_dirs = cfg.get('run', 'meshpath').split(':')
                for d in cfg_mesh_dirs:
                    self.mesh_dirs.append(d)

        del(package)

        # 2) Meshes selection layout

        # 2.1) Model for meshes table
        self.modelMeshes = StandardItemModelMeshes(self.mdl)
        self.tableViewMeshes.setModel(self.modelMeshes)
        self.tableViewMeshes.resizeColumnsToContents()
        self.tableViewMeshes.resizeRowsToContents()

        delegateName = MeshNameDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(0, delegateName)

        delegateFormat = MeshFormatDelegate(self.tableViewMeshes, self._tableViewLayout)
        self.tableViewMeshes.setItemDelegateForColumn(1, delegateFormat)

        delegateNumber = MeshNumberDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(2, delegateNumber)

        delegateGroupFaces = GroupDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(5, delegateGroupFaces)

        delegateGroupCells = GroupDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(6, delegateGroupCells)

        self.groupBoxMeshes.resizeEvent = self.MeshesResizeEvent

        # 2.2) Connections

        self.pushButtonAddMesh.clicked.connect(self.slotSearchMesh)
        self.pushButtonDeleteMesh.clicked.connect(self.slotDeleteMesh)

        # 3) Initialize meshes list

        # 3.1) Meshes default directory

        self.toolButtonMeshDir.pressed.connect(self.searchDir)
        self.toolButtonMeshDirClear.pressed.connect(self.clearDir)

        # 3.2) Meshes list

        msg = ""
        nameList = self.mdl.getMeshList()
        log.debug("__init__ -> nameList = %s " % nameList)

        if nameList:
            for i in range(len(nameList)):
                mesh = nameList[i]
                if mesh[1] != None:
                    path = os.path.join(mesh[1], mesh[0])
                else:
                    path = mesh[0]
                if not os.path.isabs(path) and self.case['mesh_path'] != None:
                    path = os.path.join(self.case['mesh_path'], path)
                if not (os.path.isfile(path) and os.path.isabs(path)):
                    msg = msg  + path + '\n'

        if msg != "":
            msg =  msg + '\n'
            title = self.tr("WARNING")
            msg2  = self.tr("The following mesh files are not present or\n" +
                            "in the meshes directory search path:\n\n" +
                            msg +
                            "Verify existence and location of the mesh files,\n" +
                            "and the 'Mesh Directory' section." )
            QMessageBox.warning(self, title, msg2)

        self._tableViewLayout()

        # Combomodels

        self.modelArg_cs_verif = ComboModel(self.comboBoxRunType, 4, 1)

        self.modelArg_cs_verif.addItem(self.tr("Import mesh only"), 'none')
        self.modelArg_cs_verif.addItem(self.tr("Mesh preprocessing only"),
                                       'mesh preprocess')
        self.modelArg_cs_verif.addItem(self.tr("Mesh quality criteria only"),
                                       'mesh quality')
        self.modelArg_cs_verif.addItem(self.tr("Standard Computation"), 'standard')

        self.modelArg_cs_verif.setItem(str_model=self.case['run_type'])

        if self.mesh_input != None:
            self.modelArg_cs_verif.disableItem(str_model='none')
        else:
            self.modelArg_cs_verif.enableItem(str_model='none')

        self.comboBoxRunType.activated[str].connect(self.slotArgRunType)

        # Checkboxes

        self.checkBoxMeshRestart.setChecked(not self.mdl.getPreprocessOnRestart())
        self.checkBoxMeshRestart.clicked.connect(self.slotMeshRestart)

        self.checkBoxMeshSave.setChecked(self.mdl.getMeshSaveOnModify())
        self.checkBoxMeshSave.clicked.connect(self.slotMeshSaveOnModify)

        # Undo/redo

        self.case.undoStartGlobal()

        # Update tree

        self.browser.configureTree(self.case)


    def MeshesResizeEvent(self, event):
        QWidget.resizeEvent(self, event)
        self._tableViewLayout()


    def searchDir(self):
        """
        Open a File Dialog in order to search the case directory.
        """
        title    = self.tr("Select input mesh directory")
        default  = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []
        for i in range(0, len(self.mesh_dirs)):
            if self.mesh_dirs[i] != None:
                l_mesh_dirs.append(QUrl.fromLocalFile(self.mesh_dirs[i]))

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setFileMode(QFileDialog.Directory)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()

            dir_name = str(s[0])
            dir_name = os.path.abspath(dir_name)

            self.lineEditMeshDir.setText(RelOrAbsPath(dir_name,
                                                 self.case['case_path']))
            self.mdl.setMeshDir(dir_name)
            self.mesh_dirs[0] = dir_name

            mesh_list = self.mdl.getMeshList()
            for i in range(len(mesh_list)):
                self.modelMeshes.dataMeshes[i][7] = mesh_list[i][1]


    def clearDir(self):
        """
        Clear the case directory.
        """
        self.lineEditMeshDir.setText('')
        if self.case['mesh_path'] != None:
            self.mesh_dirs[0] = None
            self.mdl.setMeshDir(None)

            mesh_list = self.mdl.getMeshList()
            for i in range(len(mesh_list)):
                self.modelMeshes.dataMeshes[i][7] = mesh_list[i][1]


    def selectMeshFiles(self):
        """
        Open a File Dialog in order to select mesh files.
        """
        mesh_files = []

        title    = self.tr("Select input mesh file(s)")

        default = self.mesh_dirs[0]
        if default is None:
            default  = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []
        for i in range(0, len(self.mesh_dirs)):
            if self.mesh_dirs[i] != None:
                l_mesh_dirs.append(QUrl.fromLocalFile(self.mesh_dirs[i]))

        filetypes = ""
        for Format in MeshModel().getFileFormatList():
            filetypes += "%s (%s);;"%(Format[0], Format[1])

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)
        dialog.setNameFilter(filetypes)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        if dialog.exec_() == 1:
            s = dialog.selectedFiles()
            count = len(s)
            for i in range(count):
                el = str(s[0])
                s = s[1:]
                mesh_files.append(el)

        return mesh_files


    def _tableViewLayout(self):
        """
        Configure QTableView column number
        """
        fm = self.tableViewMeshes.fontMetrics()
        last_width = 0
        n_groups = 0
        n_num = 0
        for row in self.modelMeshes.dataMeshes:
            if row[1] == 'cgns':
                n_groups += 1
                n_num += 1
            elif row[1] in ['med', 'ensight']:
                n_num += 1
            if row[7] != None:
                cmp_width = fm.size(Qt.TextSingleLine, str(row[7])).width()
                if cmp_width > last_width:
                    last_width = cmp_width

        if QT_API == "PYQT4":
            self.tableViewMeshes.horizontalHeader().setResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewMeshes.horizontalHeader().setResizeMode(1, QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableViewMeshes.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewMeshes.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        if n_num == 0:
            self.tableViewMeshes.setColumnHidden(2, True)
        else:
            self.tableViewMeshes.setColumnHidden(2, False)
            if QT_API == "PYQT4":
                self.tableViewMeshes.horizontalHeader().setResizeMode(2, QHeaderView.ResizeToContents)
            elif QT_API == "PYQT5":
                self.tableViewMeshes.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)

        if n_groups == 0:
            self.tableViewMeshes.setColumnHidden(5, True)
            self.tableViewMeshes.setColumnHidden(6, True)
        else:
            self.tableViewMeshes.setColumnHidden(5, False)
            self.tableViewMeshes.setColumnHidden(6, False)
            if QT_API == "PYQT4":
                self.tableViewMeshes.horizontalHeader().setResizeMode(5, QHeaderView.ResizeToContents)
                self.tableViewMeshes.horizontalHeader().setResizeMode(6, QHeaderView.ResizeToContents)
            elif QT_API == "PYQT5":
                self.tableViewMeshes.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeToContents)
                self.tableViewMeshes.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeToContents)

        # We have a bug under KDE (but not Gnome) if the line below is not commented.
        # self.tableViewMeshes.horizontalHeader().setSectionResizeMode(7, QHeaderView.ResizeToContents)

        self.tableViewMeshes.horizontalHeader().setStretchLastSection(True)

        if last_width > self.tableViewMeshes.columnWidth(7):
            self.tableViewMeshes.horizontalHeader().setStretchLastSection(False)
            self.tableViewMeshes.resizeColumnsToContents()


    def _addMeshInList(self,  file_name):
        """
        Tab1: add input new meshes in the listbox and case.
        """
        (d, m) = os.path.split(file_name)

        index = -1
        if self.mesh_dirs[0] != None:
            index = d.find(self.mesh_dirs[0])
            if index == 0:
                d = d[len(self.mesh_dirs[0])+1:]

        if d == '':
            d = None

        mesh = (m, d)

        # 1) Verify that the new mesh is not already in the case

        if mesh in self.mdl.getMeshList():
            title = self.tr("Warning")
            msg   = self.tr("Warning, the following input is already " \
                                "uploaded in the list:\n\n" + file_name)
            QMessageBox.information(self, title, msg)

        else:

            # 2) Update View and model

            format = MeshModel().getMeshFormat(mesh[0])
            self.mdl.addMesh(mesh)
            self.modelMeshes.addRow(mesh, format)

        self._tableViewLayout()


    @pyqtSlot()
    def selectInputMesh(self):

        # Open a File Dialog in order to search the mesh_input file or directory.

        search_dirs = []
        study_dir = os.path.split(self.case['case_path'])[0]
        for d in [os.path.join(study_dir, 'RESU_COUPLING'),
                  os.path.join(self.case['case_path'], 'RESU'),
                  study_dir]:
            if os.path.isdir(d):
                search_dirs.append(d)

        dialog = MeshInputDialog(search_dirs = search_dirs)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()

            mi = str(s[0])
            mi = os.path.abspath(mi)
            mi = RelOrAbsPath(mi, self.case['case_path'])

            self.modifyInputMesh(mi)

    @pyqtSlot(str)
    def modifyInputMesh(self, text):
        """
        Modify the mesh_input/mesh_output value
        """

        self.lineEditMeshInput.setText(text)
        self.mdl.setMeshInput(text)
        self.mesh_input = text


    @pyqtSlot()
    def slotSetInputMesh(self):

        self.radioButtonImport.setChecked(False)
        self.radioButtonExists.setChecked(True)
        self.radioButtonCartesianMesh.setChecked(False)

        self.frameMeshInput.show()
        self.frameMeshImport.hide()
        self.frameMeshCartesian.hide()

        if self.mesh_input is None:
            self.selectInputMesh()

        # If Dialog was canceled and no previous mesh_input was selected
        if self.mesh_input is None:
            self.slotSetImportMesh()
        else:
            self.setMeshOriginChoice("mesh_input")


    @pyqtSlot()
    def slotSetImportMesh(self):

        self.radioButtonImport.setChecked(True)
        self.radioButtonExists.setChecked(False)
        self.radioButtonCartesianMesh.setChecked(False)

        self.frameMeshInput.hide()
        self.frameMeshImport.show()
        self.frameMeshCartesian.hide()

        if self.mesh_input:
            self.lineEditMeshInput.setText("")
            self.mesh_input = None
            self.mdl.setMeshInput(self.mesh_input)

        self.setMeshOriginChoice("mesh_import")


    @pyqtSlot()
    def slotSetCartesianMesh(self):
        self.radioButtonImport.setChecked(False)
        self.radioButtonExists.setChecked(False)
        self.radioButtonCartesianMesh.setChecked(True)

        self.frameMeshInput.hide()
        self.frameMeshImport.hide()
        self.frameMeshCartesian.show()

        if self.mesh_input:
            self.lineEditMeshInput.setText("")
            self.mesh_input = None
            self.mdl.setMeshInput(self.mesh_input)

        self.setMeshOriginChoice("mesh_cartesian")

        for k in self.cartParams.keys():
            d = k.split("_")[0] + "_direction"
            p = k.split("_")[1]
            val = self.mdl.getCartesianParam(d, p)
            self.slotSetCartesianParam(val, k)


    @pyqtSlot(str)
    def slotSetCartesianParam(self, text, name):

        val = text
        if name in ("x_law", "y_law", "z_law"):
            self.cartParams[name].setItem(str_model=val)
            l0 = name[0]
            if val == 'constant':
                self.cartParams[l0+'_prog'].hide()
                eval('self.label'+l0.upper()+'prog.hide()')
            else:
                self.cartParams[l0+'_prog'].show()
                eval('self.label'+l0.upper()+'prog.show()')
        else:
            self.cartParams[name].setText(val)

        d = name.split("_")[0] + "_direction"
        p = name.split("_")[1]
        self.mdl.setCartesianParam(d, p,val)


    @pyqtSlot()
    def slotSearchMesh(self):
        msg = self.tr("Select a mesh file.")
        self.stbar.showMessage(msg, 2000)

        file_names = self.selectMeshFiles()

        for file_name in file_names:
            self._addMeshInList(file_name)


    @pyqtSlot()
    def slotDeleteMesh(self):
        """
        Delete the selected mesh from the list
        """
        # 1) Is there a mesh to delete ?

        selectionModel = self.tableViewMeshes.selectionModel()
        for index in selectionModel.selectedRows():
            row = index.row()
            mesh = (self.modelMeshes.dataMeshes[row][0],
                    self.modelMeshes.dataMeshes[row][7])

            # 2) Delete mesh from view and from model

            self.modelMeshes.deleteRow(row)
            self.mdl.delMesh(mesh)

        self._tableViewLayout()


    @pyqtSlot()
    def slotMeshRestart(self):
        """
        Input if preprocessing should be done on restart or not
        """
        self.mdl.setPreprocessOnRestart(not self.checkBoxMeshRestart.isChecked())


    @pyqtSlot(str)
    def slotArgRunType(self, text):
        """
        Input run type option.
        """
        run_type = self.modelArg_cs_verif.dicoV2M[str(text)]
        self.mdl.setRunType(run_type)
        self.root.initCase()
        self.browser.configureTree(self.case)


    @pyqtSlot()
    def slotMeshSaveOnModify(self):
        """
        Input if mesh should be saved if modified or not
        """
        self.mdl.setMeshSaveOnModify(self.checkBoxMeshSave.isChecked())


    def setMeshOriginChoice(self, choice):
        """
        Set mesh origin choice.
        """
        self.mesh_origin = choice
        self.mdl.setMeshOrigin(choice)

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
