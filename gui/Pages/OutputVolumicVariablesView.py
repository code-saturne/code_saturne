# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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
- OutputVolumicVariablesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string, logging

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
from code_saturne.Base.QtPage import PYQT_API_1
from code_saturne.Pages.OutputVolumicVariablesForm import Ui_OutputVolumicVariablesForm
from code_saturne.Pages.OutputControlModel import OutputControlModel
from code_saturne.Pages.OutputVolumicVariablesModel import OutputVolumicVariablesModel
from code_saturne.Pages.TimeStepModel import TimeStepModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputVolumicVariablesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbesValidator(QRegExpValidator):
    """
    Validator for real data.
    """
    def __init__(self, parent, xml_model):
        """
        Initialization for validator
        """
        regExp = QRegExp("^[0-9 ]*$")
        super(ProbesValidator, self).__init__(regExp, parent)
        self.parent = parent
        self.mdl = xml_model
        self.state = QValidator.Invalid


    def validate(self, stri, pos):
        """
        Validation method.

        QValidator.Invalid       0  The string is clearly invalid.
        QValidator.Intermediate  1  The string is a plausible intermediate value during editing.
        QValidator.Acceptable    2  The string is acceptable as a final result; i.e. it is valid.
        """
        state = QRegExpValidator.validate(self, stri, pos)[0]

        valid = True
        for probe in str(stri).split():
            if probe not in self.mdl.getVariableProbeList():
                valid = False

        if state == QValidator.Acceptable:
            if not valid:
                state = QValidator.Intermediate

        palette = self.parent.palette()

        if state == QValidator.Intermediate:
            palette.setColor(QPalette.Text, QColor("red"))
            self.parent.setPalette(palette)
        else:
            palette.setColor(QPalette.Text, QColor("black"))
            self.parent.setPalette(palette)

        self.state = state

        if PYQT_API_1:
            return (state, pos)
        else:
            return (state, stri, pos)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbesDelegate(QItemDelegate):
    """
    """
    def __init__(self, parent=None, xml_model=None):
        super(ProbesDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = ProbesValidator(editor, self.mdl)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = editor.text()
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, to_qvariant(value), Qt.DisplayRole)


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
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())

            if p_value in self.mdl.getLabelsList():
                default           = {}
                default['label']  = self.p_value
                default['list']   = self.mdl.getLabelsList()
                default['regexp'] = self.regExp
                log.debug("setModelData-> default = %s" % default)

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    p_value = result['label']
                    log.debug("setModelData-> result = %s" % result)
                else:
                    p_value = self.p_value

            model.setData(index, to_qvariant(p_value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------
class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, name, label, listing, post, probes):
        self.name    = name
        self.label   = label
        self.listing = listing
        self.post    = post
        self.probes  = probes

    def __repr__(self):
        return "variable - %s %s // listing %s // post %s // probes %s"\
               % (self.name, self.label, self.listing, self.post, self.probes)

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
                return to_qvariant(self.item.label)
            elif column == 1 and role == Qt.DisplayRole:
                return to_qvariant(self.item.name)
            elif column == 2 and role == Qt.CheckStateRole:
                value = self.item.listing
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                elif value == 'onoff':
                    return to_qvariant(Qt.PartiallyChecked)
                else:
                    return to_qvariant(Qt.Unchecked)
            elif column == 3 and role == Qt.CheckStateRole:
                value = self.item.post
                if value == 'on':
                    return to_qvariant(Qt.Checked)
                elif value == 'onoff':
                    return to_qvariant(Qt.PartiallyChecked)
                else:
                    return to_qvariant(Qt.Unchecked)
            elif column == 4 and role == Qt.DisplayRole:
                return to_qvariant(self.item.probes)
        return to_qvariant()


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#-------------------------------------------------------------------------------
# StandarItemModelOutput class
#-------------------------------------------------------------------------------

class VolumicOutputStandardItemModel(QAbstractItemModel):

    def __init__(self, parent, case, mdl):
        """
        """
        QAbstractItemModel.__init__(self)

        self.parent = parent
        self.case   = case
        self.mdl    = mdl

        self.noderoot = {}
        self.prtlist = []
        for (name, tpe) in self.mdl.list_name:
            if tpe not in self.prtlist:
                self.prtlist.append(tpe)

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
            if index.column() == 4:
                return to_qvariant(self.tr("Code_Saturne key word: IHISVR"))
            else:
                return to_qvariant()

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return to_qvariant(self.tr("Variable/Scalar name"))
            elif index.column() == 2:
                return to_qvariant(self.tr("Print in listing"))
            elif index.column() == 3:
                return to_qvariant(self.tr("Post-processing"))
            elif index.column() == 4:
                return to_qvariant(self.tr("Enter Probes number"))

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)
        elif role == Qt.CheckStateRole:
            return item.data(index.column(), role)
        #if role == Qt.TextAlignmentRole and index.column() > 1:
        #    return to_qvariant(Qt.AlignHCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        itm = index.internalPointer()
        if itm in self.noderoot.values():
            # traitement des categories
            if index.column() == 2 or index.column() == 3:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable | Qt.ItemIsTristate
            elif index.column() == 4:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            if index.column() == 2 or index.column() == 3:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
            elif index.column() == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("Output label"))
            elif section == 1:
                return to_qvariant(self.tr("Internal name"))
            elif section == 2:
                return to_qvariant(self.tr("Print in\nlisting"))
            elif section == 3:
                return to_qvariant(self.tr("Post-\nprocessing"))
            elif section == 4:
                return to_qvariant(self.tr("Probes"))
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
        for bs in self.prtlist:
            item = item_class("", bs, "off", "off", "")
            newparent = TreeItem(item, bs, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[bs] = newparent

        for (name, value) in self.mdl.list_name:
            row = self.rowCount()
            parentItem = self.noderoot[value]
            label = self.mdl.dicoLabelName[name]
            printing = self.mdl.getPrintingStatus(name)

            if OutputControlModel(self.case).getAssociatedWriterIdList("-1") == []:
                post = "off"
                self.mdl.setPostStatus(name, post)
                self.disabledItem.append((row, 3))
            else:
                 post = self.mdl.getPostStatus(name)

            if TimeStepModel(self.case).getTimePassing() in (0, 1):
                if name == 'local_time_step':
                    self.disabledItem.append((row, 3))
                    self.disabledItem.append((row, 4))

            if not self.mdl.getVariableProbeList():
                self.disabledItem.append((row, 4))

            listProbe = self.mdl.getProbesList(name)
            if listProbe:
                probes = " ".join(listProbe)
            else:
                probes = ""

            # StandardItemModel data
            item = item_class(name, label, printing, post, probes)
            newItem = TreeItem(item, "", parentItem)
            parentItem.appendChild(newItem)

        # update parent item
        for item in self.rootItem.childItems:
            size = len(item.childItems)
            listing = 0
            post    = 0
            probes  = 0

            for itm in item.childItems:
                if itm.item.listing == "on":
                    listing = listing + 1
                if itm.item.post == "on":
                    post = post + 1
                if itm.item.probes == "on":
                    probes = probes + 1

            if listing == 0:
                item.item.listing = "off"
            elif listing == size:
                item.item.listing = "on"
            else:
                item.item.listing = "onoff"

            if post == 0:
                item.item.post = "off"
            elif post == size:
                item.item.post = "on"
            else:
                item.item.post = "onoff"

#            if probes == 0:
#                item.item.probes = "off"
#            elif probes == size:
#                item.item.probes = "on"
#            else:
#                item.item.probes = "onoff"


    def setData(self, index, value, role=None):
        item = index.internalPointer()

        if index.column() == 0:
            label = str(from_qvariant(value, to_text_string))
            if label == "":
                label = item.label
            if item not in self.noderoot.values():
                self.mdl.setVariableLabel(item.item.label, label)
            item.item.label = label

        elif index.column() == 2:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.listing = "on"
            else:
                item.item.listing = "off"
            if item not in self.noderoot.values():
                self.mdl.setPrintingStatus(item.item.name, item.item.listing)
                # count for parent item
                size = len(item.parentItem.childItems)
                listing = 0
                for itm in item.parentItem.childItems:
                    if itm.item.listing == "on":
                        listing = listing + 1
                if listing == 0:
                    item.parentItem.item.listing = "off"
                elif listing == size:
                    item.parentItem.item.listing = "on"
                else:
                    item.parentItem.item.listing = "onoff"

            else:
                for itm in item.childItems:
                    self.mdl.setPrintingStatus(itm.item.name, item.item.listing)
                    itm.item.listing = item.item.listing

        elif index.column() == 3:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.post = "on"
            else:
                item.item.post = "off"

            if OutputControlModel(self.case).getAssociatedWriterIdList("-1") == []:
                item.item.post = "off"

            if item not in self.noderoot.values():
                self.mdl.setPostStatus(item.item.name, item.item.post)
                # count for parent item
                size = len(item.parentItem.childItems)
                post = 0
                for itm in item.parentItem.childItems:
                    if itm.item.post == "on":
                        post = post + 1
                if post == 0:
                    item.parentItem.item.post = "off"
                elif post == size:
                    item.parentItem.item.post = "on"
                else:
                    item.parentItem.item.post = "onoff"
            else:
                # set for each variable
                for itm in item.childItems:
                    self.mdl.setPrintingStatus(itm.item.name, item.item.post)
                    itm.item.post = item.item.post

        elif index.column() == 4:
            probes = str(from_qvariant(value, to_text_string))
            item.item.probes = probes
            if item not in self.noderoot.values():
                self.mdl.updateProbes(item.item.name, item.item.probes)

        self.dataChanged.emit(QModelIndex(), QModelIndex())

        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class OutputVolumicVariablesView(QWidget, Ui_OutputVolumicVariablesForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_OutputVolumicVariablesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.info_turb_name = []
        self.mdl = OutputVolumicVariablesModel(self.case)

        self.modelOutput = VolumicOutputStandardItemModel(parent, self.case, self.mdl)
        self.treeViewOutput.setModel(self.modelOutput)
        self.treeViewOutput.setAlternatingRowColors(True)
        self.treeViewOutput.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewOutput.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewOutput.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewOutput.expandAll()
        self.treeViewOutput.setDragEnabled(False)

        labelDelegate = LabelDelegate(self.treeViewOutput, self.mdl)
        self.treeViewOutput.setItemDelegateForColumn(0, labelDelegate)

        probesDelegate = ProbesDelegate(self.treeViewOutput, self.mdl)
        self.treeViewOutput.setItemDelegateForColumn(4, probesDelegate)
        self.treeViewOutput.resizeColumnToContents(0)
        self.treeViewOutput.resizeColumnToContents(1)

        self.case.undoStartGlobal()


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
