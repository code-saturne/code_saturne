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
- ManageCasesView
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
from code_saturne.gui.base.QtPage import RegExpValidator, from_qvariant, to_text_string
from code_saturne.gui.studymanager_gui.ManageCasesForm import Ui_ManageCasesForm
from code_saturne.gui.studymanager_gui.ManageCasesModel import ManageCasesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ManageCasesView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------

class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, idx, name, compute, post, status, run_id, tags):
        self.index   = idx
        self.name    = name
        self.compute = compute
        self.post    = post
        self.status  = status
        self.run_id  = run_id
        self.tags    = tags

    def __repr__(self):
        return "case : %s // compute : %s // post %s // status %s // run_id %s // tags %s"\
            % (self.name, self.compute, self.post, self.status, self.run_id,
               self.tags)

#-------------------------------------------------------------------------------
# Treeitem class
#-------------------------------------------------------------------------------

class TreeItem(object):
    '''
    a python object used to return row/column data, and keep note of
    it's parents and/or children
    '''
    def __init__(self, item, header, parentItem):
        self.item = item
        self.parentItem = parentItem
        self.header = header
        self.childItems = []


    def appendChild(self, item):
        self.childItems.append(item)


    def child(self, row):
        return self.childItems[row]


    def childCount(self):
        return len(self.childItems)


    def columnCount(self):
        return 6


    def data(self, column, role):
        if self.item is None:
            if column == 0:
                return self.header
            else:
                return None
        else:
            if column == 0 and role == Qt.DisplayRole:
                return self.item.name
            elif column == 1 and role == Qt.CheckStateRole:
                value = self.item.status
                if value == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked
            elif column == 2 and role == Qt.CheckStateRole:
                value = self.item.compute
                if value == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked
            elif column == 3 and role == Qt.CheckStateRole:
                value = self.item.post
                if value == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked
            elif column == 4 and role == Qt.DisplayRole:
                return self.item.run_id
            elif column == 5 and role == Qt.DisplayRole:
                return self.item.tags
        return None


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None):
        super(LabelDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v =  RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())
            model.setData(index, p_value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class TagsDelegate(LabelDelegate):
    """
    """
    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        p_value = str(editor.text())
        model.setData(index, p_value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class TextDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None):
        super(TextDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole),
                          to_text_string)
        self.p_value = str(v)
        editor.setText(v)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        p_value = str(editor.text())
        model.setData(index, p_value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# CaseStandarItemModel class
#-------------------------------------------------------------------------------

class CaseStandardItemModel(QAbstractItemModel):

    def __init__(self, parent, case, mdl):
        """
        """
        QAbstractItemModel.__init__(self)

        self.parent = parent
        self.case   = case
        self.mdl    = mdl

        self.noderoot = {}
        self.prtlist = []
        for name in self.mdl.StudyList:
            self.prtlist.append(name)

        self.rootItem = TreeItem(None, "ALL", None)
        self.parents = {0 : self.rootItem}

        self.disabledItem = []
        self.populateModel()


    def columnCount(self, parent = None):
        if parent and parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return 6


    def data(self, index, role):
        if not index.isValid():
            return None

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return None

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return self.tr("Case name")
            elif index.column() == 1:
                return self.tr("Status")
            elif index.column() == 2:
                return self.tr("Compute")
            elif index.column() == 3:
                return self.tr("Post-processing")
            elif index.column() == 4:
                return self.tr("Run_id")
            elif index.column() == 5:
                return self.tr("Tags")

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)
        elif role == Qt.CheckStateRole:
            return item.data(index.column(), role)

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        itm = index.internalPointer()
        if index.column() == 1 or index.column() == 2 or index.column() == 3:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        elif index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return self.tr("Case name")
            elif section == 1:
                return self.tr("Status")
            elif section == 2:
                return self.tr("Compute")
            elif section == 3:
                return self.tr("Post-\nprocessing")
            elif section == 4:
                return self.tr("run_id")
            elif section == 5:
                return self.tr("Tags")
        return None


    def index(self, row, column, parent = QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        try:
            childItem = parentItem.child(row)
        except:
            childItem = None

        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()


    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        if not childItem:
            return QModelIndex()

        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)


    def rowCount(self, parent=QModelIndex()):
        if parent.column() > 0:
            return 0
        if not parent.isValid():
            p_Item = self.rootItem
        else:
            p_Item = parent.internalPointer()
        return p_Item.childCount()


    def populateModel(self):
        for name in self.prtlist:
            row = self.rowCount()
            status  = self.mdl.getStudyStatus(name)
            tags    = self.mdl.getStudyTags(name)
            item = item_class(-1, name, "off", "off", status, "", tags)
            newparent = TreeItem(item, name, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[name] = newparent

        for name in self.prtlist:
            for idx in self.mdl.getCaseList(name):
                parentItem = self.noderoot[name]
                cname = self.mdl.getCaseName(name, idx)
                compute = self.mdl.getComputeStatus(name, idx)
                post    = self.mdl.getPostStatus(name, idx)
                status  = self.mdl.getStatus(name, idx)
                run_id  = self.mdl.getRunId(name, idx)
                tags    = self.mdl.getTags(name, idx)
                item = item_class(idx, cname, compute, post, status,
                                  run_id, tags)
                new_item = TreeItem(item, "", parentItem)
                parentItem.appendChild(new_item)


    def setData(self, index, value, role=None):
        item = index.internalPointer()

        if index.column() == 1:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.status = "on"
            else:
                item.item.status = "off"
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setStatus(itm.item.name, item.item.index,
                                   item.item.status)
            else:
                self.mdl.setStudyStatus(item.item.name, item.item.status)

        elif index.column() == 2:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.compute = "on"
            else:
                item.item.compute = "off"
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setComputeStatus(itm.item.name, item.item.index,
                                          item.item.compute)

        elif index.column() == 3:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.post = "on"
            else:
                item.item.post = "off"
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setPostStatus(itm.item.name, item.item.index,
                                       item.item.post)

        elif index.column() == 4:
            run_id = str(from_qvariant(value, to_text_string))
            item.item.run_id = run_id
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setRunId(itm.item.name, item.item.index,
                                  item.item.run_id)

        elif index.column() == 5:
            tags = str(from_qvariant(value, to_text_string))
            atags = str(tags).replace(",", " ").split()
            if len(atags) > 1:
                tags = ",".join(atags)
            elif len(atags) > 1:
                tags = atags[0]
            else:
                tags = ""
            item.item.tags = tags
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setTags(itm.item.name, item.item.index, item.item.tags)
            else:
                self.mdl.setStudyTags(item.item.name, item.item.tags)

        self.dataChanged.emit(QModelIndex(), QModelIndex())

        return True


#-------------------------------------------------------------------------------
# PostScriptItemModelOutput class
#-------------------------------------------------------------------------------

class PostScriptItemModel(QStandardItemModel):

    def __init__(self, mdl, study_name, case_idx):
        """
        """
        QStandardItemModel.__init__(self)

        self.mdl = mdl
        self.study_name = ''
        self.case_idx = -1

        self.scripts = []

        self.populateModel(study_name, case_idx)

        self.headers = [self.tr("Script"),
                        self.tr("Arguments"),
                        self.tr("Status")]

        self.setColumnCount(len(self.headers))


    def data(self, index, role):
        if not index.isValid():
            return None

        # ToolTips
        if role == Qt.ToolTipRole:
            col = index.column()
            if index.column() == 0:
                return self.tr("Script file, in study POST directory")
            elif index.column() == 1:
                return self.tr("Arguments passed to script")
            elif index.column() == 2:
                return self.tr("Call this script with given arguments "
                               "when postprocessing")
            return None

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return self.tr("Script")
            elif index.column() == 1:
                return self.tr("Arguments")
            elif index.column() == 2:
                return self.tr("Status")

        # Display
        if role == Qt.DisplayRole:
            if index.column() in (0, 1):
                return self.scripts[index.row()][index.column()]
        elif role == Qt.CheckStateRole:
            if index.column() == 2:
                if self.scripts[index.row()][2] == 'off':
                    return Qt.Unchecked
                else:
                    return Qt.Checked
            else:
                return None


        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        col = index.column()
        if col == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif col == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def populateModel(self, study_name, case_idx):
        self.scripts = self.mdl.getPostScripts(study_name,
                                               case_idx)
        self.study_name = study_name
        self.case_idx = case_idx
        self.setRowCount(len(self.scripts))


    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()

        if col == 1:
            if value:
                v = from_qvariant(value, to_text_string)
            else:
                v = ''
            self.scripts[row][col] = v
            self.mdl.setPostScriptArgs(self.study_name, self.case_idx, row,
                                       self.scripts[row][col])

        elif col == 2 and role == Qt.CheckStateRole:
            state = from_qvariant(value, int)
            if state == Qt.Unchecked:
                self.scripts[row][col] = 'off'
            else:
                self.scripts[row][col] = 'on'
            self.mdl.setPostScriptStatus(self.study_name, self.case_idx, row,
                                         self.scripts[row][col])

        return True


    def addRow(self, name):
        """
        Add a row in the table.
        """
        item = [name, '', 'on']

        self.scripts.append(item)

        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        self.setRowCount(self.rowCount() - 1)
        del self.scripts[row]


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ManageCasesView(QWidget, Ui_ManageCasesForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ManageCasesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.mdl = ManageCasesModel(self.case)

        if not self.mdl.list_study:
            self.add_study()

        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.setAlternatingRowColors(True)
        self.treeViewCases.setSelectionMode(QAbstractItemView.SingleSelection)
        self.treeViewCases.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewCases.expandAll()
        self.treeViewCases.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewCases.setDragEnabled(False)

        runidDelegate = LabelDelegate(self.treeViewCases)
        self.treeViewCases.setItemDelegateForColumn(4, runidDelegate)

        tagsDelegate = TagsDelegate(self.treeViewCases)
        self.treeViewCases.setItemDelegateForColumn(5, tagsDelegate)

        self.treeViewCases.resizeColumnToContents(0)

        self.pushButtonAdd.clicked.connect(self.slotAddCase)
        self.pushButtonDelete.clicked.connect(self.slotDeleteCase)
        self.pushButtonAddStudy.clicked.connect(self.slotAddStudy)
        self.pushButtonDeleteStudy.clicked.connect(self.slotDeleteStudy)
        self.toolButtonDuplicate.clicked.connect(self.slotDuplicateCase)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.checkBoxCompare.clicked.connect(self.slotCompareStatus)
        self.pushButtonAddPostScript.clicked.connect(self.slotAddPostScriptFile)
        self.pushButtonRemovePostScript.clicked.connect(self.slotRemovePostScriptFile)
        self.pushButtonAddInput.clicked.connect(self.slotAddInputFile)
        self.pushButtonRemoveInput.clicked.connect(self.slotRemoveInputFile)

        self.listInput.currentItemChanged.connect(self.slotSelectInputRow)
        self.listInput.itemClicked.connect(self.slotClickInput)
        self.lineEditNotebookArgs.textChanged[str].connect(self.slotNotebookArgs)
        self.lineEditParametricArgs.textChanged[str].connect(self.slotParametricArgs)
        self.lineEditKwArgs.textChanged[str].connect(self.slotKwArgs)
        self.lineEditCompareArgs.textChanged[str].connect(self.slotCompareArgs)

        self.modelPostScripts = PostScriptItemModel(self.mdl, '', -1)
        self.tablePostScript.setModel(self.modelPostScripts)
        self.tablePostScript.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tablePostScript.setGridStyle(Qt.NoPen)
        self.tablePostScript.setAlternatingRowColors(True)

        hh = self.tablePostScript.horizontalHeader()
        hh.setSectionResizeMode(1, QHeaderView.Stretch)

        argsDelegate = TextDelegate(self.tablePostScript)
        self.tablePostScript.setItemDelegateForColumn(1, argsDelegate)

        self.listInput.setSelectionMode(QAbstractItemView.SingleSelection)
        self.listInput.setAlternatingRowColors(True)

        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.groupBoxCompare.hide()

        self.pushButtonDelete.setEnabled(False)
        self.pushButtonDeleteStudy.setEnabled(False)
        self.toolButtonDuplicate.setEnabled(False)
        self.pushButtonAdd.setEnabled(False)
        self.pushButtonRemoveInput.setEnabled(False)

        self.listInputClickOn = 0

    #TODO add self.disabledItem.append((row, 3)) for study nodes


    def __get_study_and_case_idx__(self):
        """
        Get study name and case index based on current selection
        """
        current = self.treeViewCases.currentIndex()

        if current == self.treeViewCases.rootIndex():
            self.pushButtonAdd.setEnabled(False)
            case_idx = -1
            study = None
        elif current.parent() == self.treeViewCases.rootIndex():
            case_idx = -1
            study = current.internalPointer().item.name
        else:
            case_idx = current.row()
            study = current.parent().internalPointer().item.name

        return study, case_idx


    def add_case(self, study):
        """
        public slot
        """
        title = self.tr("Add existing case")

        cur_path = os.getcwd()
        path = os.path.abspath(os.path.join(self.mdl.repo, study))
        try:
            os.chdir(path)
        except Exception:
            pass

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(path)
        dialog.setFileMode(QFileDialog.DirectoryOnly)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()
            dir_path = str(s[0])
            study_path, case_name = os.path.split(dir_path)
            if study != os.path.split(study_path)[1]:
                title = self.tr("WARNING")
                msg   = self.tr("This selected case is not in the directory of the study")
                QMessageBox.information(self, title, msg)
            else:
                self.mdl.addCase(study, case_name)

            log.debug("add_case -> %s" % dir_path)
        os.chdir(cur_path)


    def add_study(self):
        """
        public slot
        """
        title = self.tr("Add existing study")

        cur_path = os.getcwd()
        path = os.path.abspath(self.mdl.repo)
        os.chdir(path)

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(path)
        dialog.setFileMode(QFileDialog.DirectoryOnly)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()
            dir_path = str(s[0])
            dir_path = os.path.relpath(dir_path)
            if dir_path not in os.listdir(path):
                title = self.tr("WARNING")
                msg   = self.tr("This selected study is not in the repository directory")
                QMessageBox.information(self, title, msg)
            else:
                self.mdl.addStudy(dir_path)

                title = self.tr("Loading study/case")
                msg   = self.tr("Do you want to use automatic load ?")

                reply = QMessageBox.question(self,
                                             title,
                                             msg,
                                             QMessageBox.Yes|
                                             QMessageBox.No)
                if reply == QMessageBox.Yes:
                    self.mdl.loadCases(dir_path)

            log.debug("add_study -> %s" % dir_path)
        os.chdir(cur_path)


    def slotAddStudy(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        self.add_study()
        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.treeViewCases.expandAll()
        self.changeSelection()


    def slotDeleteStudy(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.internalPointer().item.name
        self.mdl.deleteStudy(study)
        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.treeViewCases.expandAll()
        self.changeSelection()


    def slotAddCase(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = ""
        if current.parent() == self.treeViewCases.rootIndex():
            study = current.internalPointer().item.name
        else:
            study = current.parent().internalPointer().item.name
        self.add_case(study)
        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.treeViewCases.expandAll()
        self.changeSelection()


    def slotDeleteCase(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        self.mdl.deleteCase(study, idx)
        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.treeViewCases.expandAll()
        self.changeSelection()


    def slotDuplicateCase(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        self.mdl.duplicateCase(study, idx)
        self.modelCases = CaseStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewCases.setModel(self.modelCases)
        self.treeViewCases.selectionModel().selectionChanged.connect(self.slotChangeSelection)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.treeViewCases.expandAll()
        self.changeSelection()


    def slotChangeSelection(self, new, old):
        """
        slot for change of selection
        """
        self.changeSelection()


    def changeSelection(self):
        """
        """
        self.pushButtonDelete.setEnabled(False)
        self.pushButtonDeleteStudy.setEnabled(False)
        self.toolButtonDuplicate.setEnabled(False)
        self.pushButtonAdd.setEnabled(True)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxInput.hide()
        self.groupBoxCompare.hide()

        study = None
        current = self.treeViewCases.currentIndex()
        idx = current.row()

        if current == self.treeViewCases.rootIndex():
            self.pushButtonAdd.setEnabled(False)
        elif current.parent() == self.treeViewCases.rootIndex():
            idx = -1
            study = current.internalPointer().item.name
            # study
            self.pushButtonDeleteStudy.setEnabled(True)
            idx = -1
        else:
            # case
            self.pushButtonDelete.setEnabled(True)
            self.toolButtonDuplicate.setEnabled(True)
            self.groupBoxPrepro.show()
            self.groupBoxPost.show()
            self.groupBoxInput.show()
            self.groupBoxCompare.show()
            study = current.parent().internalPointer().item.name

            # prepro
            notebook_args = self.mdl.getNotebookArgs(study, idx)
            self.lineEditNotebookArgs.setText(str(notebook_args))
            parametric_args = self.mdl.getParametricArgs(study, idx)
            self.lineEditParametricArgs.setText(str(parametric_args))
            kw_args = self.mdl.getKwArgs(study, idx)
            self.lineEditKwArgs.setText(str(kw_args))

            # compare
            status = self.mdl.getCompareStatus(study, idx)
            if status == "on":
                self.checkBoxCompare.setChecked(True)
            else:
                self.checkBoxCompare.setChecked(False)
            compare_args = self.mdl.getCompareArgs(study, idx)
            self.lineEditCompareArgs.setText(str(compare_args))

        if study is not None: # Study or case, not root
            self.groupBoxPost.show()
            self.groupBoxInput.show()

            self.modelPostScripts.populateModel(study, idx)
            self.tablePostScript.setModel(self.modelPostScripts)

            input_names = self.mdl.getPostInput(study, idx)
            self.listInput.clear()
            if input_names:
                for name in input_names:
                    self.listInput.addItem(str(name))


    @pyqtSlot()
    def slotCompareStatus(self):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        if self.checkBoxCompare.isChecked():
            self.mdl.setCompareStatus(study, idx, "on")
        else:
            self.mdl.setCompareStatus(study, idx, "off")


    @pyqtSlot(str)
    def slotNotebookArgs(self, text):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        args = str(text)
        self.mdl.setNotebookArgs(study, idx, args)


    @pyqtSlot(str)
    def slotParametricArgs(self, text):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        args = str(text)
        self.mdl.setParametricArgs(study, idx, args)


    @pyqtSlot(str)
    def slotKwArgs(self, text):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        args = str(text)
        self.mdl.setKwArgs(study, idx, args)


    @pyqtSlot(str)
    def slotCompareArgs(self, text):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        args = str(text)
        self.mdl.setCompareArgs(study, idx, args)


    @pyqtSlot(str)
    def slotPostArgs(self, text):
        """
        """
        study, idx = self.__get_study_and_case_idx__()

        args = str(text)
        self.mdl.setPostScriptArgs(study, idx, args)


    def slotPostFile(self):
        """
        public slot
        """
        study, idx = self.__get_study_and_case_idx__()

        cur_path = os.getcwd()
        base_dir = os.path.abspath(os.path.join(self.mdl.repo, study))
        rep = os.path.abspath(os.path.join(base_dir, "POST"))
        if not os.path.isdir(rep):
            rep = os.path.abspath(os.path.join(os.path.split(base_dir)[0], "POST"))
        if not os.path.isdir(rep):
            rep = base_dir
        title = self.tr("postprocess script")
        filetypes = self.tr("All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        file = str(file)

        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(rep):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the POST directory of the study")
            QMessageBox.information(self, title, msg)
        else:
            self.mdl.setPostScriptName(study, idx, file)
        return


    def slotAddPostScriptFile(self):
        """
        public slot
        """
        study, idx = self.__get_study_and_case_idx__()

        cur_path = os.getcwd()
        base_dir = os.path.abspath(os.path.join(self.mdl.repo, study))
        rep = os.path.abspath(os.path.join(base_dir, "POST"))
        if not os.path.isdir(rep):
            rep = os.path.abspath(os.path.join(os.path.split(base_dir)[0], "POST"))
        if not os.path.isdir(rep):
            rep = base_dir
        title = self.tr("input file for postprocess script")
        filetypes = self.tr("All Files (*)")
        fname = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        fname = str(fname)

        if not fname:
            return

        if rep != os.path.commonprefix([fname, rep]):
            title = self.tr("WARNING")
            msg   = self.tr("The selected file is not in the POST directory of the study")
            QMessageBox.information(self, title, msg)
        else:
            fname = os.path.relpath(fname, rep)
            self.mdl.addPostScript(study, idx, fname)
            self.modelPostScripts.populateModel(study, idx)
            self.tablePostScript.setModel(self.modelPostScripts)


    def slotRemovePostScriptFile(self):
        """
        public slot
        """
        study, idx = self.__get_study_and_case_idx__()
        name = None

        selectionModel = self.tablePostScript.selectionModel()

        # Allow handling of multiple selections. Remove last selected
        # elements first to avoid changing index of next elements to
        # remove.

        l = []
        for index in selectionModel.selectedRows():
            l.append(index.row())
        l.reverse()
        for s_idx in l:
            self.mdl.removePostScript(study, idx, s_idx)

        selectionModel.clearSelection()

        self.modelPostScripts.populateModel(study, idx)
        self.tablePostScript.setModel(self.modelPostScripts)


    def slotClickInput(self):
        """
        public slot
        """
        idx = self.listInput.currentRow();
        self.listInputClickOn -= 1
        if self.listInputClickOn < 0:
            self.listInputClickOn = 1
        if idx > -1:
            if self.listInputClickOn == 1:
                self.listInput.setCurrentRow(-1)


    def slotSelectInputRow(self, current, previous):
        """
        public slot
        """
        if self.listInput.currentRow() > -1:
            self.pushButtonRemoveInput.setEnabled(True)
            self.listInputClickOn = 3
        else:
            self.pushButtonRemoveInput.setEnabled(False)
            self.listInputClickOn = 0


    def slotAddInputFile(self):
        """
        public slot
        """
        study, idx = self.__get_study_and_case_idx__()

        cur_path = os.getcwd()
        base_dir = os.path.abspath(os.path.join(self.mdl.repo, study))
        rep = os.path.abspath(os.path.join(base_dir, "POST"))
        if not os.path.isdir(rep):
            rep = os.path.abspath(os.path.join(os.path.split(base_dir)[0], "POST"))
        if not os.path.isdir(rep):
            rep = base_dir
        title = self.tr("input file for postprocess script")
        filetypes = self.tr("All Files (*)")
        fname = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        fname = str(fname)

        if not fname:
            return

        if rep != os.path.commonprefix([fname, rep]):
            title = self.tr("WARNING")
            msg   = self.tr("The selected file is not in the POST directory of the study")
            QMessageBox.information(self, title, msg)
        else:
            fname = os.path.relpath(fname, rep)
            if self.mdl.addPostInput(study, idx, fname):
                self.listInput.addItem(str(fname))


    def slotRemoveInputFile(self):
        """
        public slot
        """
        study, idx = self.__get_study_and_case_idx__()
        name = None

        try:
            l_idx = self.listInput.currentRow()
            item = self.listInput.takeItem(l_idx)
            name = str(item.text())
        except Exception:
            return

        self.mdl.removePostInput(study, idx, name)
        self.listInput.setCurrentRow(-1)


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
