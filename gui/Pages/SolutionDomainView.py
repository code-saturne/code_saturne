# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
- StandardItemModelThinWall
- StandardItemModelExtrude
- SolutionDomainView
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, RegExpValidator, IntValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.SolutionDomainForm import Ui_SolutionDomainForm
from code_saturne.Pages.SolutionDomainModel import RelOrAbsPath, MeshModel, SolutionDomainModel
from code_saturne.Pages.FacesSelectionView import StandardItemModelFaces

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
        isValid = index.model().dataMeshes[row][7]

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
        model.setData(index, to_qvariant(key))
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
        model.setData(index, to_qvariant(value))


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
        model.setData(index, to_qvariant(value))

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
        model.setData(index, to_qvariant(value), Qt.DisplayRole)


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
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


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
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


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
                        self.tr("Add face groups"),
                        self.tr("Add cell groups"),
                        self.tr("Path")]

        self.tooltip = [self.tr("Preprocessor option: --mesh"),
                        self.tr("Preprocessor sub-option: --format"),
                        self.tr("Preprocessor sub-option: --num"),
                        self.tr("Preprocessor sub-option: --reorient"),
                        self.tr("Preprocessor sub-option: --grp-fac"),
                        self.tr("Preprocessor sub-option: --grp-cel"),
                        self.tr("Preprocessor option: --mesh")]

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
                    self.setData(index, to_qvariant(self.dataMeshes[row][1]))


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
            return to_qvariant()

        col = index.column()

        if role == Qt.ToolTipRole:
            return to_qvariant(self.tooltip[col])

        elif role == Qt.DisplayRole:
            d = self.dataMeshes[index.row()][col]
            if d:
                if col == 1:
                    return to_qvariant(self.formatDict[d])
                elif col == 3:
                    return to_qvariant()
                else:
                    return to_qvariant(d)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            if col == 6:
                return to_qvariant(Qt.AlignLeft)
            else:
                return to_qvariant(Qt.AlignCenter)

        elif role == Qt.CheckStateRole:
            if col == 3:
                d = self.dataMeshes[index.row()][3]
                if d == True:
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)
            else:
                return to_qvariant()

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        self.__disableData(index.row())

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        if index.column() in [0, 6]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() in [1, 2, 4, 5]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 3:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()


        mesh = (self.dataMeshes[row][0], self.dataMeshes[row][6])

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

        elif col == 4:
            if value:
                v = from_qvariant(value, to_text_string)
            else:
                v = ''
            self.dataMeshes[row][col] = v
            if v:
                self.mdl.setMeshGroupFaces(mesh, v)

        elif col == 5:
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
            item = [mesh[0], format, "1", "", "", "", mesh[1]]
        elif format == 'cgns':
            item = [mesh[0], format, "1", "", "off", "off", mesh[1]]
        else:
            item = [mesh[0], format, "", "", "", "", mesh[1]]
        item.append(self.__isMeshPathValid(mesh))

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

        mesh = (self.dataMeshes[row][0], self.dataMeshes[row][6])
        self.dataMeshes[row][7] = self.__isMeshPathValid(mesh)

        if not self.dataMeshes[row][1] in ["cgns", "med", "ensight"]:
            if (row, 2) not in self.disabledItem:
                self.disabledItem.append((row, 2))
        else:
            if (row, 2) in self.disabledItem:
                self.disabledItem.remove((row, 2))

        if self.dataMeshes[row][1] != "cgns":
            if (row, 4) not in self.disabledItem:
                self.disabledItem.append((row, 4))
            if (row, 5) not in self.disabledItem:
                self.disabledItem.append((row, 5))
        else:
            if (row, 4) in self.disabledItem:
                self.disabledItem.remove((row, 4))
            if (row, 5) in self.disabledItem:
                self.disabledItem.remove((row, 5))


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
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(data)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 :
            return Qt.ItemIsSelectable
        else:
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
            return to_qvariant()

        if role == Qt.ToolTipRole:
            return to_qvariant()

        elif role == Qt.DisplayRole:
            data = self._data[index.row()][index.column()]
            if data:
                return to_qvariant(data)
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 :
            return Qt.ItemIsSelectable
        else:
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

        self.name_filter = str(self.tr("Imported or preprocessed meshes (mesh_input mesh_output)"))
        self.setNameFilter(self.name_filter)

        self.select_label = str(self.tr("Choose"))
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
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_SolutionDomainForm.__init__(self)
        self.setupUi(self)

        self.stbar = stbar
        self.case = case
        self.case.undoStopGlobal()
        self.mdl = SolutionDomainModel(self.case)

        # 0) Mesh Input

        self.mesh_input = self.mdl.getMeshInput()

        if self.mesh_input:
            self.radioButtonImport.setChecked(False)
            self.radioButtonExists.setChecked(True)
            self.frameMeshImport.hide()
            self.lineEditMeshInput.setText(self.mesh_input)
        else:
            self.radioButtonImport.setChecked(True)
            self.radioButtonExists.setChecked(False)
            self.frameMeshInput.hide()

        self.radioButtonImport.clicked.connect(self.slotSetImportMesh)
        self.radioButtonExists.clicked.connect(self.slotSetInputMesh)
        self.toolButtonMeshInput.pressed.connect(self.selectInputMesh)

        # 1) Meshes directory

        self.mesh_dirs = [None]

        d = self.mdl.getMeshDir()
        study_path = os.path.split(self.case['case_path'])[0]

        if d == None:
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
        self.tableViewMeshes.setItemDelegateForColumn(4, delegateGroupFaces)

        delegateGroupCells = GroupDelegate(self.tableViewMeshes)
        self.tableViewMeshes.setItemDelegateForColumn(5, delegateGroupCells)

        self.groupBoxMeshes.resizeEvent = self.MeshesResizeEvent

        # 2.2) Connections

        self.pushButtonAddMesh.clicked.connect(self.slotSearchMesh)
        self.pushButtonDeleteMesh.clicked.connect(self.slotDeleteMesh)
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

        # 5) Initialize meshes list

        # 5.1) Meshes default directory

        self.toolButtonMeshDir.pressed.connect(self.searchDir)
        self.toolButtonMeshDirClear.pressed.connect(self.clearDir)

        # 5.1) Meshes list

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


        # 5.2) Warped faces cutting

        if self.mdl.getCutStatus() == 'on':
            self.groupBoxWarp.setChecked(True)
            self.slotFacesCutting(True)
        else:
            self.groupBoxWarp.setChecked(False)
            self.slotFacesCutting(False)

        v = self.mdl.getCutAngle()
        self.warp = v
        self.lineEditWarp.setText(str(self.warp))


        # 5.3) Mesh Smoothing

        if self.mdl.getSmoothingStatus() == 'on':
            self.groupBoxMeshSmooth.setChecked(True)
            self.slotMeshSmooth(True)
        else:
            self.groupBoxMeshSmooth.setChecked(False)
            self.slotMeshSmooth(False)

        v = self.mdl.getSmoothAngle()
        self.smooth = v
        self.lineEditMeshSmooth.setText(str(self.smooth))

        # 5.4) tab Widget
        self.tabWidget.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


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

        dialog = QFileDialog(self, title, default)
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
                self.modelMeshes.dataMeshes[i][6] = mesh_list[i][1]


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
                self.modelMeshes.dataMeshes[i][6] = mesh_list[i][1]


    def selectMeshFiles(self):
        """
        Open a File Dialog in order to select mesh files.
        """
        mesh_files = []

        title    = self.tr("Select input mesh file(s)")

        default = self.mesh_dirs[0]
        if default == None:
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

        dialog = QFileDialog(self, title, default, filetypes)
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
            if row[6] != None:
                cmp_width = fm.size(Qt.TextSingleLine, str(row[6])).width()
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
            self.tableViewMeshes.setColumnHidden(4, True)
            self.tableViewMeshes.setColumnHidden(5, True)
        else:
            self.tableViewMeshes.setColumnHidden(4, False)
            self.tableViewMeshes.setColumnHidden(5, False)
            if QT_API == "PYQT4":
                self.tableViewMeshes.horizontalHeader().setResizeMode(4, QHeaderView.ResizeToContents)
                self.tableViewMeshes.horizontalHeader().setResizeMode(5, QHeaderView.ResizeToContents)
            elif QT_API == "PYQT5":
                self.tableViewMeshes.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
                self.tableViewMeshes.horizontalHeader().setSectionResizeMode(5, QHeaderView.ResizeToContents)

        # We have a bug under KDE (but not Gnome) if the line below is not commented.
        # self.tableViewMeshes.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeToContents)

        self.tableViewMeshes.horizontalHeader().setStretchLastSection(True)

        if last_width > self.tableViewMeshes.columnWidth(6):
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

            self.lineEditMeshInput.setText(mi)
            self.mdl.setMeshInput(mi)
            self.mesh_input = mi


    @pyqtSlot()
    def slotSetInputMesh(self):

        self.radioButtonImport.setChecked(False)
        self.radioButtonExists.setChecked(True)

        self.frameMeshInput.show()
        self.frameMeshImport.hide()

        if self.mesh_input == None:
            self.selectInputMesh()

        # If Dialog was canceled and no previous mesh_input was selected
        if self.mesh_input == None:
            self.slotSetImportMesh()


    @pyqtSlot()
    def slotSetImportMesh(self):

        self.radioButtonImport.setChecked(True)
        self.radioButtonExists.setChecked(False)

        self.frameMeshInput.hide()
        self.frameMeshImport.show()

        if self.mesh_input:
            self.lineEditMeshInput.setText("")
            self.mesh_input = None
            self.mdl.setMeshInput(self.mesh_input)


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
                    self.modelMeshes.dataMeshes[row][6])

            # 2) Delete mesh from view and from model

            self.modelMeshes.deleteRow(row)
            self.mdl.delMesh(mesh)

        self._tableViewLayout()


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
