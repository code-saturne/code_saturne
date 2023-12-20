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
- OutputVolumicVariablesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

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
from code_saturne.gui.base.QtPage import ComboModel
from code_saturne.gui.case.OutputVolumicVariablesForm import Ui_OutputVolumicVariablesForm
from code_saturne.model.OutputControlModel import OutputControlModel
from code_saturne.model.OutputVolumicVariablesModel import OutputVolumicVariablesModel
from code_saturne.model.UserCalculatorModel import UserCalculatorModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputVolumicVariablesView")
log.setLevel(GuiParam.DEBUG)

_calculator_group = "User functions"
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

                from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    p_value = result['label']
                    log.debug("setModelData-> result = %s" % result)
                else:
                    p_value = self.p_value

            model.setData(index, p_value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------

class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, name, label, listing, post, monitor, value=None):
        self.name    = name
        self.label   = label
        self.status  = [listing, post, monitor]
        self.value = value

    def __repr__(self):
        return "variable - %s %s // log %s // post %s // monitor %s"\
               % (self.name, self.label,
                  self.status[0], self.status[1], self.status[2])


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
        if self.item is None:
            if column == 0:
                return self.header
            else:
                return None
        else:
            if column == 0 and role == Qt.DisplayRole:
                return self.item.label
            elif column == 1 and role == Qt.DisplayRole:
                return self.item.name
            elif column >= 2 and role == Qt.CheckStateRole:
                value = self.item.status[column - 2]
                if value == 'on':
                    return Qt.Checked
                elif value == 'onoff':
                    return Qt.PartiallyChecked
                else:
                    return Qt.Unchecked
        return None

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

    def __init__(self, parent, case, mdl, calculator=None):
        """
        """
        QAbstractItemModel.__init__(self)

        self.parent = parent
        self.case   = case
        self.mdl    = mdl
        self.calculator = calculator

        self.noderoot = {}
        self.prtlist = []
        for (name, tpe) in self.mdl.list_name:
            if tpe not in self.prtlist:
                self.prtlist.append(tpe)

        # Reordering categories order for NCFD multiphase solver
        if 'Properties' in self.prtlist:

            forced_position_labels = ['Properties',
                                      'Variables',
                                      'Time moments']

            self.prtlist.sort()
            rlist = []
            rlist.append(self.prtlist.index('Properties'))
            rlist.append(self.prtlist.index('Variables'))

            for i in range(len(self.prtlist)):
                if self.prtlist[i] not in forced_position_labels:
#                if self.prtlist[i] != 'Properties' and self.prtlist[i] != 'Variables':
                    rlist.append(i)

            if 'Time moments' in self.prtlist:
                rlist.append(self.prtlist.index('Time moments'))

            self.prtlist = [self.prtlist[j] for j in rlist]

        if self.calculator:
            if self.calculator.getNumberOfFunctions(location="cells") > 0:
                self.prtlist.append(_calculator_group)


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
            return None

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return None

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return self.tr("Variable/Scalar name")
            elif index.column() == 2:
                return self.tr("Print in listing")
            elif index.column() == 3:
                return self.tr("Post-processing")
            elif index.column() == 4:
                return self.tr("Monitoring")

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)
        elif role == Qt.CheckStateRole:
            return item.data(index.column(), role)
        #if role == Qt.TextAlignmentRole and index.column() > 1:
        #    return Qt.AlignHCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # disable item
        if (index.row(), index.column()) in self.disabledItem:
            return Qt.ItemIsEnabled

        itm = index.internalPointer()
        if itm in self.noderoot.values():
            # traitement des categories
            if index.column() >= 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable | Qt.ItemIsTristate
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            if index.column() >= 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
            elif index.column() == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return self.tr("Output label")
            elif section == 1:
                return self.tr("Internal name")
            elif section == 2:
                return self.tr("Print in\nlisting")
            elif section == 3:
                return self.tr("Post-\nprocessing")
            elif section == 4:
                return self.tr("Monitoring")
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
        for bs in self.prtlist:
            item = item_class("", bs, "off", "off", "off")
            newparent = TreeItem(item, bs, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[bs] = newparent

        _is_post_active = OutputControlModel(self.case).isVolumeWriterActive()

        for (name, value) in self.mdl.list_name:
            row = self.rowCount()
            parentItem = self.noderoot[value]
            label = self.mdl.dicoLabelName[name]
            printing = self.mdl.getPrintingStatus(name)

            if not _is_post_active:
                self.mdl.setPostStatus(name, "off")

            post = self.mdl.getPostStatus(name)

            if self.case.xmlRootNode().tagName == "Code_Saturne_GUI":
                from code_saturne.model.TimeStepModel import TimeStepModel
                if TimeStepModel(self.case).getTimePassing() in (0, 1):
                    if name == 'local_time_step':
                        self.disabledItem.append((row, 3))
                        self.disabledItem.append((row, 4))
                del TimeStepModel

            monitor = self.mdl.getMonitorStatus(name)

            # StandardItemModel data
            item = item_class(name, label, printing, post, monitor, value)
            newItem = TreeItem(item, "", parentItem)
            parentItem.appendChild(newItem)

        # Calculator functions
        for func in self.calculator.getFunctionsNamesList(location="cells"):
            row = self.rowCount()
            parentItem = self.noderoot[_calculator_group]
            printing = self.calculator.getPrintingStatus(func)

            if not _is_post_active:
                self.calculator.setPostStatus("off")
            post = self.calculator.getPostStatus(func)

            monitor = self.calculator.getMonitorStatus(func)

            item = item_class(func, func, printing, post, monitor, _calculator_group)
            newItem = TreeItem(item, "", parentItem)
            parentItem.appendChild(newItem)

        # update parent item
        for item in self.rootItem.childItems:
            size = len(item.childItems)
            listing = 0
            post    = 0
            monitor = 0

            for itm in item.childItems:
                if itm.item.status[0] == "on":
                    listing = listing + 1
                if itm.item.status[1] == "on":
                    post = post + 1
                if itm.item.status[2] == "on":
                    monitor = monitor + 1

            if listing == 0:
                item.item.status[0] = "off"
            elif listing == size:
                item.item.status[0] = "on"
            else:
                item.item.status[0] = "onoff"

            if post == 0:
                item.item.status[1] = "off"
            elif post == size:
                item.item.status[1] = "on"
            else:
                item.item.status[1] = "onoff"

            if monitor == 0:
                item.item.status[2] = "off"
            elif monitor == size:
                item.item.status[2] = "on"
            else:
                item.item.status[2] = "onoff"


    def setData(self, index, value, role=None):
        item = index.internalPointer()

        if index.column() == 0:
            label = str(from_qvariant(value, to_text_string))
            if label == "":
                label = item.label
            if item not in self.noderoot.values():
                if item.item.variable != _calculator_group:
                    self.mdl.setVariableLabel(item.item.label, label)
                else:
                    self.calculator.setName(label)
            item.item.label = label

        elif index.column() >= 2:
            c_id = index.column() - 2
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.status[c_id] = "on"
            else:
                item.item.status[c_id] = "off"
            if c_id == 1:
                if not OutputControlModel(self.case).isVolumeWriterActive():
                    item.item.status[1] = "off"
            if item not in self.noderoot.values():
                if item.item.value == _calculator_group:
                    if c_id == 0:
                        self.calculator.setPrintingStatus(item.item.name,
                                                          item.item.status[c_id])
                    elif c_id == 1:
                        self.calculator.setPostStatus(item.item.name,
                                                      item.item.status[c_id])
                    elif c_id == 2:
                        self.calculator.setMonitorStatus(item.item.name,
                                                         item.item.status[c_id])
                else:
                    if c_id == 0:
                        self.mdl.setPrintingStatus(item.item.name,
                                                   item.item.status[0])
                    elif c_id == 1:
                        self.mdl.setPostStatus(item.item.name,
                                               item.item.status[1])
                    elif c_id == 2:
                        self.mdl.setMonitorStatus(item.item.name,
                                                  item.item.status[2])
                # count for parent item
                size = len(item.parentItem.childItems)
                active = 0
                for itm in item.parentItem.childItems:
                    if itm.item.status[c_id] == "on":
                        active = active + 1
                if active == 0:
                    item.parentItem.item.status[c_id] = "off"
                elif active == size:
                    item.parentItem.item.status[c_id] = "on"
                else:
                    item.parentItem.item.status[c_id] = "onoff"

            else:
                for itm in item.childItems:
                    if itm.item.value == _calculator_group:
                        if c_id == 0:
                            self.calculator.setPrintingStatus(itm.item.name,
                                                              item.item.status[c_id])
                        elif c_id == 1:
                            self.calculator.setPostStatus(itm.item.name,
                                                          item.item.status[c_id])
                        elif c_id == 2:
                            self.calculator.setMonitorStatus(itm.item.name,
                                                             item.item.status[c_id])
                    else:
                        if c_id == 0:
                            self.mdl.setPrintingStatus(itm.item.name,
                                                       item.item.status[0])
                        elif c_id == 1:
                            self.mdl.setPostStatus(itm.item.name,
                                                   item.item.status[1])
                        elif c_id == 2:
                            self.mdl.setMonitorStatus(itm.item.name,
                                                      item.item.status[2])

                    itm.item.status[c_id] = item.item.status[c_id]

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
        self.parent = parent
        self.case.undoStopGlobal()
        self.info_turb_name = []

        self.mdl = OutputVolumicVariablesModel(self.case)
        self.calculator = UserCalculatorModel(self.case)

        self.modelOutput = VolumicOutputStandardItemModel(parent, self.case, self.mdl, self.calculator)
        self.treeViewOutput.setModel(self.modelOutput)
        self.treeViewOutput.setAlternatingRowColors(True)
        self.treeViewOutput.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewOutput.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewOutput.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewOutput.expandAll()
        self.treeViewOutput.setDragEnabled(False)

        labelDelegate = LabelDelegate(self.treeViewOutput, self.mdl)
        self.treeViewOutput.setItemDelegateForColumn(0, labelDelegate)

        self.treeViewOutput.resizeColumnToContents(0)
        self.treeViewOutput.resizeColumnToContents(1)


        if self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":

            self.groupBox_2.hide()

        elif self.case.xmlRootNode().tagName == "Code_Saturne_GUI":

            self.correctionEstimator = ComboModel(self.comboBoxIescor, 2, 1)
            self.correctionEstimator.addItem(self.tr("off"), '0')
            self.correctionEstimator.addItem(self.tr("on"), '2')

            self.driftEstimator = ComboModel(self.comboBoxIesder, 2, 1)
            self.driftEstimator.addItem(self.tr("off"), '0')
            self.driftEstimator.addItem(self.tr("on"), '2')

            self.predictionEstimator = ComboModel(self.comboBoxIespre, 2, 1)
            self.predictionEstimator.addItem(self.tr("off"), '0')
            self.predictionEstimator.addItem(self.tr("on"), '2')

            self.totalEstimator = ComboModel(self.comboBoxIestot, 2, 1)
            self.totalEstimator.addItem(self.tr("off"), '0')
            self.totalEstimator.addItem(self.tr("on"), '2')

            self.comboBoxIescor.activated[str].connect(self.slotCorrectionEstimator)
            self.comboBoxIesder.activated[str].connect(self.slotDriftEstimator)
            self.comboBoxIespre.activated[str].connect(self.slotPredictionEstimator)
            self.comboBoxIestot.activated[str].connect(self.slotTotalEstimator)

            modelIescor = self.mdl.getEstimatorModel("Correction")
            self.correctionEstimator.setItem(str_model=modelIescor)

            modelIesder = self.mdl.getEstimatorModel("Drift")
            self.driftEstimator.setItem(str_model=modelIesder)

            modelIespre = self.mdl.getEstimatorModel("Prediction")
            self.predictionEstimator.setItem(str_model=modelIespre)

            modelIestot = self.mdl.getEstimatorModel("Total")
            self.totalEstimator.setItem(str_model=modelIestot)


        self.case.undoStartGlobal()


    def initializeView(self):
        """
        """
        self.modelOutput = VolumicOutputStandardItemModel(self.parent, self.case, self.mdl)
        self.treeViewOutput.setModel(self.modelOutput)
        self.treeViewOutput.expandAll()


    @pyqtSlot(str)
    def slotCorrectionEstimator(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.correctionEstimator.dicoV2M[str(text)]
        self.mdl.setEstimatorModel("Correction", model)
        self.mdl.updateList()

        self.initializeView()


    @pyqtSlot(str)
    def slotDriftEstimator(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.driftEstimator.dicoV2M[str(text)]
        self.mdl.setEstimatorModel("Drift", model)
        self.mdl.updateList()

        self.initializeView()


    @pyqtSlot(str)
    def slotPredictionEstimator(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.predictionEstimator.dicoV2M[str(text)]
        self.mdl.setEstimatorModel("Prediction", model)
        self.mdl.updateList()

        self.initializeView()


    @pyqtSlot(str)
    def slotTotalEstimator(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.totalEstimator.dicoV2M[str(text)]
        self.mdl.setEstimatorModel("Total", model)
        self.mdl.updateList()

        self.initializeView()


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
