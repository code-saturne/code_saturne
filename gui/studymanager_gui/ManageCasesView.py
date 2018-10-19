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
- ManageCasesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging, os

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import RegExpValidator, to_qvariant, from_qvariant, to_text_string
from code_saturne.studymanager_gui.ManageCasesForm import Ui_ManageCasesForm
from code_saturne.studymanager_gui.ManageCasesModel import ManageCasesModel

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
    def __init__(self, idx, name, compute, post, status, run_id):
        self.index   = idx
        self.name    = name
        self.compute = compute
        self.post    = post
        self.status  = status
        self.run_id  = run_id

    def __repr__(self):
        return "case : %s // compute : %s // post %s // status %s // run_id %s"\
               % (self.name, self.compute, self.post, self.status, self.run_id)

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
        return 5


    def data(self, column, role):
        if self.item == None:
            if column == 0:
                return to_qvariant(self.header)
            else:
                return to_qvariant()
        else:
            if column == 0 and role == Qt.DisplayRole:
                return to_qvariant(self.item.name)
            elif column == 1 and role == Qt.CheckStateRole:
                value = self.item.status
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)
            elif column == 2 and role == Qt.CheckStateRole:
                value = self.item.compute
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)
            elif column == 3 and role == Qt.CheckStateRole:
                value = self.item.post
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                else:
                    return to_qvariant(Qt.Unchecked)
            elif column == 4 and role == Qt.DisplayRole:
                return to_qvariant(self.item.run_id)
        return to_qvariant()


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
    def __init__(self, parent=None, xml_model=None):
        super(LabelDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


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
            model.setData(index, to_qvariant(p_value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModelOutput class
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
            return 5


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return to_qvariant()

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return to_qvariant(self.tr("Case name"))
            elif index.column() == 1:
                return to_qvariant(self.tr("Status"))
            elif index.column() == 2:
                return to_qvariant(self.tr("Compute"))
            elif index.column() == 3:
                return to_qvariant(self.tr("Post-processing"))
            elif index.column() == 4:
                return to_qvariant(self.tr("Run_id"))

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)
        elif role == Qt.CheckStateRole:
            return item.data(index.column(), role)

        return to_qvariant()


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
                return to_qvariant(self.tr("Case name"))
            elif section == 1:
                return to_qvariant(self.tr("Status"))
            elif section == 2:
                return to_qvariant(self.tr("Compute"))
            elif section == 3:
                return to_qvariant(self.tr("Post-\nprocessing"))
            elif section == 4:
                return to_qvariant(self.tr("run_id"))
        return to_qvariant()


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
            item = item_class(-1, name, "off", "off", status, "")
            newparent = TreeItem(item, name, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[name] = newparent
#            self.disabledItem.append((row, 2))
#            self.disabledItem.append((row, 3))
#            self.disabledItem.append((row, 4))

        for name in self.prtlist:
            for idx in self.mdl.getCaseList(name):
                parentItem = self.noderoot[name]
                cname = self.mdl.getCaseName(name, idx)
                compute = self.mdl.getComputeStatus(name, idx)
                post    = self.mdl.getPostStatus(name, idx)
                status  = self.mdl.getStatus(name, idx)
                run_id  = self.mdl.getRunId(name, idx)
                item = item_class(idx, cname, compute, post, status, run_id)
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
                self.mdl.setStatus(itm.item.name, item.item.index, item.item.status)
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
                self.mdl.setComputeStatus(itm.item.name, item.item.index, item.item.compute)

        elif index.column() == 3:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.post = "on"
            else:
                item.item.post = "off"
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setPostStatus(itm.item.name, item.item.index, item.item.post)

        elif index.column() == 4:
            run_id = str(from_qvariant(value, to_text_string))
            item.item.run_id = run_id
            if item not in self.noderoot.values():
                itm = item.parentItem
                self.mdl.setRunId(itm.item.name, item.item.index, item.item.run_id)

        self.dataChanged.emit(QModelIndex(), QModelIndex())

        return True


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
        self.treeViewCases.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewCases.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewCases.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewCases.expandAll()
        self.treeViewCases.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewCases.setDragEnabled(False)

        runidDelegate = LabelDelegate(self.treeViewCases, self.mdl)
        self.treeViewCases.setItemDelegateForColumn(4, runidDelegate)

        self.treeViewCases.resizeColumnToContents(0)

        self.pushButtonAdd.clicked.connect(self.slotAddCase)
        self.pushButtonDelete.clicked.connect(self.slotDeleteCase)
        self.pushButtonAddStudy.clicked.connect(self.slotAddStudy)
        self.pushButtonDeleteStudy.clicked.connect(self.slotDeleteStudy)
        self.toolButtonDuplicate.clicked.connect(self.slotDuplicateCase)
        self.treeViewCases.clicked.connect(self.slotChangeSelection)

        self.checkBoxPrepro.clicked.connect(self.slotPreproStatus)
        self.checkBoxPost.clicked.connect(self.slotPostStatus)
        self.checkBoxCompare.clicked.connect(self.slotCompareStatus)
        self.pushButtonPrepro.clicked.connect(self.slotPreproFile)
        self.pushButtonPost.clicked.connect(self.slotPostFile)
        self.pushButtonInput.clicked.connect(self.slotInputFile)
        self.lineEditPreproArgs.textChanged[str].connect(self.slotPreproArgs)
        self.lineEditPostArgs.textChanged[str].connect(self.slotPostArgs)
        self.lineEditCompareArgs.textChanged[str].connect(self.slotCompareArgs)

        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxCompare.hide()

        self.lineEditPrepro.setEnabled(False)
        self.lineEditInput.setEnabled(False)
        self.lineEditPost.setEnabled(False)
        self.pushButtonDelete.setEnabled(False)
        self.pushButtonDeleteStudy.setEnabled(False)
        self.toolButtonDuplicate.setEnabled(False)
        self.pushButtonAdd.setEnabled(False)
    #TODO ajouter self.disabledItem.append((row, 3)) pour les noeuds study


    def add_case(self, study):
        """
        public slot
        """
        title = self.tr("Add existing case")

        cur_path = os.getcwd()
        path = os.path.abspath(os.path.join(self.mdl.repo, study))
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
                msg   = self.tr("This selected case is not in the directory of the study")
                QMessageBox.information(self, title, msg)
            else:
                self.mdl.addCase(study, dir_path)

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
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.treeViewCases.expandAll()
        self.slotChangeSelection()


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
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.treeViewCases.expandAll()
        self.slotChangeSelection()


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
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.treeViewCases.expandAll()
        self.slotChangeSelection()


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
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.treeViewCases.expandAll()
        self.slotChangeSelection()


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
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.treeViewCases.expandAll()
        self.slotChangeSelection()


    def slotChangeSelection(self):
        """
        """
        self.pushButtonDelete.setEnabled(False)
        self.pushButtonDeleteStudy.setEnabled(False)
        self.toolButtonDuplicate.setEnabled(False)
        self.pushButtonAdd.setEnabled(True)
        self.groupBoxPrepro.hide()
        self.groupBoxPost.hide()
        self.groupBoxCompare.hide()

        current = self.treeViewCases.currentIndex()
        idx = current.row()
        if current == self.treeViewCases.rootIndex():
            self.pushButtonAdd.setEnabled(False)
        elif current.parent() == self.treeViewCases.rootIndex():
            # study
            self.pushButtonDeleteStudy.setEnabled(True)
            self.groupBoxPost.show()
            self.pushButtonInput.hide()
            self.lineEditInput.hide()
            self.labelinput.hide()
            study = current.internalPointer().item.name

            status = self.mdl.getStudyPostScriptStatus(study)
            if status == "on":
                self.checkBoxPost.setChecked(True)
            else:
                self.checkBoxPost.setChecked(False)

            script_name = self.mdl.getStudyPostScriptName(study)
            self.lineEditPost.setText(str(script_name))
            if status == "on":
                if script_name != "":
                    self.pushButtonPost.setStyleSheet("background-color: green")
                else:
                    self.pushButtonPost.setStyleSheet("background-color: red")
            else:
                self.pushButtonPost.setStyleSheet("background-color: None")

            script_args = self.mdl.getStudyPostScriptArgs(study)
            self.lineEditPostArgs.setText(str(script_args))
        else:
            # case
            self.pushButtonDelete.setEnabled(True)
            self.toolButtonDuplicate.setEnabled(True)
            self.groupBoxPrepro.show()
            self.groupBoxPost.show()
            self.groupBoxCompare.show()
            self.pushButtonInput.show()
            self.lineEditInput.show()
            self.labelinput.show()
            study = current.parent().internalPointer().item.name

            # prepro
            status = self.mdl.getPreproScriptStatus(study, idx)
            if status == "on":
                self.checkBoxPrepro.setChecked(True)
            else:
                self.checkBoxPrepro.setChecked(False)

            script_name = self.mdl.getPreproScriptName(study, idx)
            self.lineEditPrepro.setText(str(script_name))
            if status == "on":
                if script_name != "":
                    self.pushButtonPrepro.setStyleSheet("background-color: green")
                else:
                    self.pushButtonPrepro.setStyleSheet("background-color: red")
            else:
                self.pushButtonPrepro.setStyleSheet("background-color: None")

            script_args = self.mdl.getPreproScriptArgs(study, idx)
            self.lineEditPreproArgs.setText(str(script_args))

            # post
            status = self.mdl.getPostScriptStatus(study, idx)
            if status == "on":
                self.checkBoxPost.setChecked(True)
            else:
                self.checkBoxPost.setChecked(False)

            script_name = self.mdl.getPostScriptName(study, idx)
            self.lineEditPost.setText(str(script_name))
            if status == "on":
                if script_name != "":
                    self.pushButtonPost.setStyleSheet("background-color: green")
                else:
                    self.pushButtonPost.setStyleSheet("background-color: red")
            else:
                self.pushButtonPost.setStyleSheet("background-color: None")

            script_args = self.mdl.getPostScriptArgs(study, idx)
            self.lineEditPostArgs.setText(str(script_args))

            input_name = self.mdl.getPostScriptInput(study, idx)
            self.lineEditInput.setText(str(input_name))

            # compare
            status = self.mdl.getCompareStatus(study, idx)
            if status == "on":
                self.checkBoxCompare.setChecked(True)
            else:
                self.checkBoxCompare.setChecked(False)
            compare_args = self.mdl.getCompareArgs(study, idx)
            self.lineEditCompareArgs.setText(str(compare_args))


    @pyqtSlot()
    def slotPreproStatus(self):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        if self.checkBoxPrepro.isChecked():
            self.mdl.setPreproScriptStatus(study, idx, "on")
            if self.mdl.getPreproScriptName(study, idx) != "":
                self.pushButtonPrepro.setStyleSheet("background-color: green")
            else:
                self.pushButtonPrepro.setStyleSheet("background-color: red")
        else:
            self.mdl.setPreproScriptStatus(study, idx, "off")
            self.pushButtonPrepro.setStyleSheet("background-color: None")


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


    @pyqtSlot()
    def slotPostStatus(self):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        if current.parent() == self.treeViewCases.rootIndex():
            study = current.internalPointer().item.name
            if self.checkBoxPost.isChecked():
                self.mdl.setStudyPostScriptStatus(study, "on")
                if self.mdl.getStudyPostScriptName(study) != "":
                    self.pushButtonPost.setStyleSheet("background-color: green")
                else:
                    self.pushButtonPost.setStyleSheet("background-color: red")
            else:
                self.mdl.setStudyPostScriptStatus(study, "off")
                self.pushButtonPost.setStyleSheet("background-color: None")
        else:
            study = current.parent().internalPointer().item.name
            if self.checkBoxPost.isChecked():
                self.mdl.setPostScriptStatus(study, idx, "on")
                if self.mdl.getPostScriptName(study, idx) != "":
                    self.pushButtonPost.setStyleSheet("background-color: green")
                else:
                    self.pushButtonPost.setStyleSheet("background-color: red")
            else:
                self.mdl.setPostScriptStatus(study, idx, "off")
                self.pushButtonPost.setStyleSheet("background-color: None")


    @pyqtSlot(str)
    def slotPreproArgs(self, text):
        """
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name
        args = str(text)
        self.mdl.setPreproScriptArgs(study, idx, args)


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
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        args = str(text)
        if current.parent() == self.treeViewCases.rootIndex():
            study = current.internalPointer().item.name
            self.mdl.setStudyPostScriptArgs(study, args)
        else:
            study = current.parent().internalPointer().item.name
            self.mdl.setPostScriptArgs(study, idx, args)


    def slotPreproFile(self):
        """
        Select a prepro script
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name

        cur_path = os.getcwd()
        rep = os.path.abspath(os.path.join(self.mdl.repo, study, "MESH"))
        os.chdir(rep)
        title = self.tr("preprocess script")
        filetypes = self.tr("(*py*);;All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        file = str(file)
        os.chdir(cur_path)

        if not file:
            return
        file = os.path.basename(file)

        if file not in os.listdir(rep):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the MESH directory of te study")
            QMessageBox.information(self, title, msg)
        else:
            self.lineEditPrepro.setText(str(file))
            self.mdl.setPreproScriptName(study, idx, file)


    def slotPostFile(self):
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

        cur_path = os.getcwd()
        rep = os.path.abspath(os.path.join(self.mdl.repo, study, "POST"))
        os.chdir(rep)
        title = self.tr("postprocess script")
        filetypes = self.tr("All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        file = str(file)
        os.chdir(cur_path)

        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(rep):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the POST directory of the study")
            QMessageBox.information(self, title, msg)
        else:
            self.lineEditPost.setText(str(file))
            if current.parent() == self.treeViewCases.rootIndex():
                self.mdl.setStudyPostScriptName(study, file)
            else:
                self.mdl.setPostScriptName(study, idx, file)
        return


    def slotInputFile(self):
        """
        public slot
        """
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        study = current.parent().internalPointer().item.name

        cur_path = os.getcwd()
        rep = os.path.abspath(os.path.join(self.mdl.repo, study, "POST"))
        os.chdir(rep)
        title = self.tr("input file for postprocess script")
        filetypes = self.tr("All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        file = str(file)
        os.chdir(cur_path)

        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(rep):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the POST directory of the study")
            QMessageBox.information(self, title, msg)
        else:
            self.lineEditInput.setText(str(file))
            self.mdl.setPostScriptInput(study, idx, file)


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
