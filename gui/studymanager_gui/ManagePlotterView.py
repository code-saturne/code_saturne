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
This module defines the following classes:
- ManagePlotterView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *
import os

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, IntValidator
from code_saturne.Base.QtPage import RegExpValidator, to_qvariant
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.studymanager_gui.ManagePlotForm import Ui_ManagePlotForm
from code_saturne.studymanager_gui.ManagePlotterForm import Ui_ManagePlotterForm
from code_saturne.studymanager_gui.ManagePlotterSubplotForm import Ui_ManagePlotterSubplotForm
from code_saturne.studymanager_gui.ManagePlotterModel import ManagePlotterModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ManagePlotterView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Advanced dialog
#-------------------------------------------------------------------------------

class ManageSubplotDialogView(QDialog, Ui_ManagePlotterSubplotForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, mdl, study, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_ManagePlotterSubplotForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl  = mdl
        self.study = study

        title = self.tr("Select subplot(s)")
        self.setWindowTitle(title)
        self.default = default
        self.result  = self.default.copy()

        subplot_list = self.mdl.getSubplotList(self.study)

        model = QStandardItemModel(self.listViewSubplot)
        self.item_list = []
        for plotid in subplot_list:
            name = str(plotid) + ": " + self.mdl.getSubplotTitle(self.study, plotid)
            item = QStandardItem(name)
            item.setCheckable(True)
            self.item_list.append(item)
            model.appendRow(item)
            if plotid in self.default['idlist']:
                item.setCheckState(True)
            else:
                item.setCheckState(False)
        self.listViewSubplot.setModel(model)


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        self.result['idlist'] = ''

        for itm in range(len(self.item_list)):
            if self.item_list[itm].checkState():
                if self.result['idlist'] != '':
                    self.result['idlist'] = self.result['idlist'] + ' ' + str(itm)
                else:
                    self.result['idlist'] = str(itm)

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
# Advanced dialog for plot
#-------------------------------------------------------------------------------

class ManagePlotDialogView(QDialog, Ui_ManagePlotForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, mdl, study, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_ManagePlotForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl  = mdl
        self.study = study

        title = self.tr("Manage a plot")
        self.setWindowTitle(title)
        self.default = default
        self.result  = self.default.copy()

        # add les validator sur les lineEdit
        self.lineEditColor.setText(self.default['color'])
        self.lineEditFormat.setText(self.default['format'])
        self.lineEditLegend.setText(self.default['legend'])
        self.lineEditXcol.setText(self.default['xcol'])
        self.lineEditYcol.setText(self.default['ycol'])
        self.lineEditWidth.setText(self.default['width'])
        self.lineEditMarker.setText(self.default['marker'])
        self.lineEditXerr.setText(self.default['xerr'])
        self.lineEditXerrp.setText(self.default['xerrp'])
        self.lineEditYerr.setText(self.default['yerr'])
        self.lineEditYerrp.setText(self.default['yerrp'])


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        self.result['color']  = str(self.lineEditColor.text())
        self.result['format'] = str(self.lineEditFormat.text())
        self.result['legend'] = str(self.lineEditLegend.text())
        self.result['xcol']   = str(self.lineEditXcol.text())
        self.result['ycol']   = str(self.lineEditYcol.text())
        self.result['width']  = str(self.lineEditWidth.text())
        self.result['marker'] = str(self.lineEditMarker.text())
        self.result['xerr']   = str(self.lineEditXerr.text())
        self.result['xerrp']  = str(self.lineEditXerrp.text())
        self.result['yerr']   = str(self.lineEditYerr.text())
        self.result['yerrp']  = str(self.lineEditYerrp.text())

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
# Line edit delegate for float
#-------------------------------------------------------------------------------

class FloatDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(FloatDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = from_qvariant(editor.text(), float)
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, to_qvariant(value))

#-------------------------------------------------------------------------------
# Line edit delegate for integer
#-------------------------------------------------------------------------------

class IntDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(IntDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        value = from_qvariant(editor.text(), int)
        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, to_qvariant(value))


#-------------------------------------------------------------------------------
# Combo box delegate for the figure format
#-------------------------------------------------------------------------------

class FormatFigureDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None):
        super(FormatFigureDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem("pdf")
        editor.addItem("png")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"pdf": 0, "png": 1}
        row = index.row()
        string = index.model().dataFigure[row]['format']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value))


#-------------------------------------------------------------------------------
# Line edit delegate for label
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
        return editor


    def setEditorData(self, editor, index):
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        p_value = str(editor.text())
        model.setData(index, to_qvariant(p_value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QStandardItemModel for tableViewSubplot
#-------------------------------------------------------------------------------

class StandardItemModelSubplot(QStandardItemModel):
    def __init__(self, mdl, study):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = mdl
        self.study = study
        self.dataSubplot = []
        self.populateModel()

        self.headers = [self.tr("id"),
                        self.tr("title"),
                        self.tr("xlabel"),
                        self.tr("ylabel"),
                        self.tr("legend\nstatus"),
                        self.tr("legend position"),
                        self.tr("x axis\nrange"),
                        self.tr("y axis\nrange")]
        self.keys = ['id', 'title', 'xlabel', 'ylabel',\
                     'legstatus', 'legpos', 'xaxis', 'yaxis']
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                if column == 4:
                    role = Qt.CheckStateRole
                else:
                    role = Qt.DisplayRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)


    def populateModel(self):
        for idx in self.mdl.getSubplotList(self.study):
            self.addSubplot(idx)


    def addSubplot(self, idx):
        dico              = {}
        dico['id']        = idx
        dico['title']     = self.mdl.getSubplotTitle(self.study, idx)
        dico['xlabel']    = self.mdl.getSubplotXLabel(self.study, idx)
        dico['ylabel']    = self.mdl.getSubplotYLabel(self.study, idx)
        dico['legstatus'] = self.mdl.getSubplotLegStatus(self.study, idx)
        dico['legpos']    = self.mdl.getSubplotLegPos(self.study, idx)
        dico['xaxis']     = self.mdl.getSubplotXLim(self.study, idx)
        dico['yaxis']     = self.mdl.getSubplotYLim(self.study, idx)
        self.dataSubplot.append(dico)
        log.debug("populateModel-> dataSubplot = %s" % dico)
        row = self.rowCount()
        self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        row = index.row()
        column = index.column()
        dico = self.dataSubplot[row]
        key = self.keys[column]

        if dico[key] == None:
            return to_qvariant()

        if role == Qt.DisplayRole and column != 4:
            return to_qvariant(dico[key])


        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        elif role == Qt.CheckStateRole and column == 4:
            st = dico[key]
            if st == 'on':
                return to_qvariant(Qt.Checked)
            else:
                return to_qvariant(Qt.Unchecked)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() == 4:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        idx = self.dataSubplot[row]['id']

        # set title
        if column == 1:
            self.dataSubplot[row]['title'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotTitle(self.study, idx, self.dataSubplot[row]['title'])

        # set xlabel
        elif column == 2:
            self.dataSubplot[row]['xlabel'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotXLabel(self.study, idx, self.dataSubplot[row]['xlabel'])

        # set ylabel
        elif column == 3:
            self.dataSubplot[row]['ylabel'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotYLabel(self.study, idx, self.dataSubplot[row]['ylabel'])

        # set legstatus
        elif column == 4:
            v = from_qvariant(value, int)
            if v == Qt.Unchecked:
                self.dataSubplot[row]['legstatus'] = "off"
            else:
                self.dataSubplot[row]['legstatus'] = "on"
            self.mdl.setSubplotLegStatus(self.study, idx, self.dataSubplot[row]['legstatus'])

        # set ylabel
        elif column == 5:
            self.dataSubplot[row]['legpos'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotLegPos(self.study, idx, self.dataSubplot[row]['legpos'])

        # set ylabel
        elif column == 6:
            self.dataSubplot[row]['xaxis'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotXLim(self.study, idx, self.dataSubplot[row]['xaxis'])

        # set ylabel
        elif column == 7:
            self.dataSubplot[row]['yaxis'] = str(from_qvariant(value, to_text_string))
            self.mdl.setSubplotYLim(self.study, idx, self.dataSubplot[row]['yaxis'])

        self.dataChanged.emit(index, index)
        return True


#-------------------------------------------------------------------------------
# QStandardItemModel for tableViewFigure
#-------------------------------------------------------------------------------

class StandardItemModelFigure(QStandardItemModel):
    def __init__(self, mdl, study):
        """
        """
        QStandardItemModel.__init__(self)
        self.mdl = mdl
        self.study = study
        self.dataFigure = []
        self.populateModel()

        self.headers = [self.tr("id"),
                        self.tr("name"),
                        self.tr("subplot\nid list"),
                        self.tr("title"),
                        self.tr("row\n number"),
                        self.tr("column\n number"),
                        self.tr("format")]
        self.keys = ['id', 'name', 'id_list', 'title', 'row', 'column', 'format']
        self.setColumnCount(len(self.headers))

        # Initialize the flags
        for row in range(self.rowCount()):
            for column in range(self.columnCount()):
                role = Qt.DisplayRole
                index = self.index(row, column)
                value = self.data(index, role)
                self.setData(index, value)


    def populateModel(self):
        for idx in self.mdl.getFigureList(self.study):
            self.addFigure(idx)


    def addFigure(self, idx):
        dico            = {}
        dico['id']      = idx
        dico['name']    = self.mdl.getFigureName(self.study, idx)
        dico['id_list'] = self.mdl.getFigureIdList(self.study, idx)
        dico['title']   = self.mdl.getFigureTitle(self.study, idx)
        dico['row']     = self.mdl.getFigureRow(self.study, idx)
        dico['column']  = self.mdl.getFigureColumn(self.study, idx)
        dico['format']  = self.mdl.getFigureFormat(self.study, idx)
        self.dataFigure.append(dico)
        log.debug("populateModel-> dataFigure = %s" % dico)
        row = self.rowCount()
        self.setRowCount(row + 1)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        row = index.row()
        column = index.column()
        dico = self.dataFigure[row]
        key = self.keys[column]

        if dico[key] == None:
            return to_qvariant()

        if role == Qt.DisplayRole:
            return to_qvariant(dico[key])

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0 or index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role=None):
        row = index.row()
        column = index.column()
        idx = self.dataFigure[row]['id']

        self.keys = ['id', 'name', 'id_list', 'title', 'row', 'column', 'format']
        # set title
        if column == 1:
            self.dataFigure[row]['name'] = str(from_qvariant(value, to_text_string))
            self.mdl.setFigureName(self.study, idx, self.dataFigure[row]['name'])

        # set ylabel
        elif column == 3:
            self.dataFigure[row]['title'] = str(from_qvariant(value, to_text_string))
            self.mdl.setFigureTitle(self.study, idx, self.dataFigure[row]['title'])

        # set ylabel
        elif column == 4:
            self.dataFigure[row]['row'] = str(from_qvariant(value, to_text_string))
            self.mdl.setFigureRow(self.study, idx, self.dataFigure[row]['row'])

        # set ylabel
        elif column == 5:
            self.dataFigure[row]['column'] = str(from_qvariant(value, to_text_string))
            self.mdl.setFigureColumn(self.study, idx, self.dataFigure[row]['column'])

        # set ylabel
        elif column == 6:
            self.dataFigure[row]['format'] = str(from_qvariant(value, to_text_string))
            self.mdl.setFigureFormat(self.study, idx, self.dataFigure[row]['format'])

        self.dataChanged.emit(index, index)
        return True



#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------
class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, tpe, name, id_list):
        self.tpe  = tpe
        self.name = name
        self.idlist = id_list

    def __repr__(self):
        return "case : %s //  %s"\
               % (self.tpe, self.name)


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
        return 3


    def data(self, column, role):
        if self.item == None:
            if column == 0:
                return to_qvariant(self.header)
            else:
                return to_qvariant()
        else:
            if column == 0 and role == Qt.DisplayRole:
                return to_qvariant(self.item.tpe)
            elif column == 1 and role == Qt.DisplayRole:
                return to_qvariant(self.item.name)
            elif column == 2 and role == Qt.DisplayRole:
                return to_qvariant(self.item.idlist)
        return to_qvariant()


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#-------------------------------------------------------------------------------
# QStandardItemModel for treeViewMeasurement
#-------------------------------------------------------------------------------

class StandardItemModelMeasurement(QAbstractItemModel):
    def __init__(self, mdl, study):
        """
        """
        QAbstractItemModel.__init__(self)
        self.mdl = mdl
        self.study = study

        self.noderoot = {}
        self.prtlist = []
        for name in self.mdl.getMeasurementList(self.study):
            self.prtlist.append(name)

        self.rootItem = TreeItem(None, "ALL", None)
        self.parents = {0 : self.rootItem}

        self.populateModel()


    def columnCount(self, parent = None):
        if parent and parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return 3


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
                return to_qvariant(self.tr("type"))
            elif index.column() == 1:
                return to_qvariant(self.tr("identification"))
            elif index.column() == 2:
                return to_qvariant(self.tr("subplot id list"))

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("type"))
            elif section == 1:
                return to_qvariant(self.tr("identification"))
            elif section == 2:
                return to_qvariant(self.tr("subplot id list"))
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
            item = item_class("measurement", name, "")
            newparent = TreeItem(item, name, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[name] = newparent

        measurement_idx = 0
        for name in self.prtlist:
            for idx in self.mdl.getMeasurementPlotList(self.study, name):
                parentItem = self.noderoot[name]
                idlist = self.mdl.getMeasurementIdList(self.study, measurement_idx, idx)
                item = item_class("plot", str(idx), idlist)
                new_item = TreeItem(item, "", parentItem)
                parentItem.appendChild(new_item)
            measurement_idx = measurement_idx + 1

    def setData(self, index, value, role=None):
        self.dataChanged.emit(QModelIndex(), QModelIndex())
        return True


#-------------------------------------------------------------------------------
# QStandardItemModel for treeViewCases
#-------------------------------------------------------------------------------
class StandardItemModelCase(QAbstractItemModel):
    def __init__(self, mdl, study):
        """
        """
        QAbstractItemModel.__init__(self)
        self.mdl = mdl
        self.study = study

        self.noderoot = {}

        self.rootItem = TreeItem(None, "ALL", None)
        self.parents = {0 : self.rootItem}

        self.populateModel()


    def columnCount(self, parent = None):
        if parent and parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return 3


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
                return to_qvariant(self.tr("type"))
            elif index.column() == 1:
                return to_qvariant(self.tr("identification"))
            elif index.column() == 2:
                return to_qvariant(self.tr("subplot id list"))

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("type"))
            elif section == 1:
                return to_qvariant(self.tr("identification"))
            elif section == 2:
                return to_qvariant(self.tr("subplot id list"))
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
        for idx in self.mdl.getCaseList(self.study):
            name = self.mdl.getCaseName(self.study, idx)
            item = item_class("case", name, "")
            newparent = TreeItem(item, name, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[name] = newparent

            # data file for each case
            data_idx = 0
            for data in self.mdl.getCaseDataList(self.study, idx):
                item = item_class("data", data, "")
                parentItem = self.noderoot[name]
                new_item = TreeItem(item, "", parentItem)
                parentItem.appendChild(new_item)

                # plot for each data file
                for plot in self.mdl.getCasePlotList(self.study, idx, data_idx):
                    id_list = self.mdl.getCaseIdList(self.study, idx, data_idx, plot)
                    itm = item_class("plot", plot, id_list)
                    new_itm = TreeItem(itm, "", new_item)
                    new_item.appendChild(new_itm)

                data_idx = data_idx + 1


    def setData(self, index, value, role=None):
        self.dataChanged.emit(QModelIndex(), QModelIndex())
        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ManagePlotterView(QWidget, Ui_ManagePlotterForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ManagePlotterForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = ManagePlotterModel(self.case)

        if len(self.mdl.getStudyList()) == 0:
            return

        # model for study choice
        self.modelStudy = ComboModel(self.comboBoxStudy, 1, 1)
        for std in self.mdl.getStudyList():
            self.modelStudy.addItem(self.tr(std), std)
        self.current_study = self.mdl.getStudyList()[0]
        self.modelStudy.setItem(str_model = self.current_study)

        # Connections
        self.comboBoxStudy.activated[str].connect(self.slotStudy)
        self.pushButtonAddSubplot.clicked.connect(self.slotAddSubplot)
        self.pushButtonDeleteSubplot.clicked.connect(self.slotDeleteSubplot)
        self.pushButtonAssociatedSubplot.setEnabled(False)

        self.tableViewFigure.clicked.connect(self.slotChangeSelectionFigure)
        self.pushButtonAddFigure.clicked.connect(self.slotAddFigure)
        self.pushButtonDeleteFigure.clicked.connect(self.slotDeleteFigure)
        self.pushButtonAssociatedSubplot.clicked.connect(self.slotAssociatedFiguresSubplot)
        self.pushButtonDeleteFigure.setEnabled(False)

        self.treeViewMeasurement.clicked.connect(self.slotChangeSelectionMeasurement)
        self.pushButtonAddFile.clicked.connect(self.slotAddMeasurementFile)
        self.pushButtonDeleteFile.clicked.connect(self.slotDeleteMeasurementFile)
        self.pushButtonAssociatedMeasurementSubplot.clicked.connect(self.slotAssociatedMeasurementSubplot)
        self.pushButtonAddPlot.clicked.connect(self.slotMeasurementAddPlot)
        self.pushButtonModifyPlot.clicked.connect(self.slotMeasurementPlot)
        self.pushButtonDeletePlot.clicked.connect(self.slotMeasurementDeletePlot)
        self.pushButtonAssociatedMeasurementSubplot.setEnabled(False)
        self.pushButtonDeleteFile.setEnabled(False)
        self.pushButtonAddPlot.setEnabled(False)
        self.pushButtonModifyPlot.setEnabled(False)
        self.pushButtonDeletePlot.setEnabled(False)

        self.treeViewCases.clicked.connect(self.slotChangeSelectionCases)
        self.pushButtonAddData.clicked.connect(self.slotCaseAddData)
        self.pushButtonDeleteData.clicked.connect(self.slotCaseDeleteData)
        self.pushButtonAddCasePlot.clicked.connect(self.slotCaseAddPlot)
        self.pushButtonDeleteCasePlot.clicked.connect(self.slotCaseDeletePlot)
        self.pushButtonModifyCasePlot.clicked.connect(self.slotCasePlot)
        self.pushButtonAssociatedCaseSubplot.clicked.connect(self.slotCaseAssociatedSubplot)
        self.pushButtonAddData.setEnabled(False)
        self.pushButtonDeleteData.setEnabled(False)
        self.pushButtonAddCasePlot.setEnabled(False)
        self.pushButtonDeleteCasePlot.setEnabled(False)
        self.pushButtonModifyCasePlot.setEnabled(False)
        self.pushButtonAssociatedCaseSubplot.setEnabled(False)

        self.updateStudyView()


    def updateStudyView(self):
        """
        update view when study change
        """
        # model for tableViewSubplot
        self.modelSubplot = StandardItemModelSubplot(self.mdl, self.current_study)
        self.tableViewSubplot.setModel(self.modelSubplot)
        self.tableViewSubplot.resizeColumnToContents(0)
        self.tableViewSubplot.resizeRowsToContents()
        self.tableViewSubplot.setAlternatingRowColors(True)
        self.tableViewSubplot.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewSubplot.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewSubplot.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewSubplot.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSubplot.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        delegateFloat = FloatDelegate(self.tableViewSubplot)
        self.tableViewSubplot.setItemDelegateForColumn(6, delegateFloat)
        self.tableViewSubplot.setItemDelegateForColumn(7, delegateFloat)

        labelDelegate = LabelDelegate(self.tableViewSubplot, self.mdl)
        self.tableViewSubplot.setItemDelegateForColumn(1, labelDelegate)
        self.tableViewSubplot.setItemDelegateForColumn(2, labelDelegate)
        self.tableViewSubplot.setItemDelegateForColumn(3, labelDelegate)
        self.tableViewSubplot.setItemDelegateForColumn(5, labelDelegate)

        # model for tableViewFigure
        self.modelFigure = StandardItemModelFigure(self.mdl, self.current_study)
        self.tableViewFigure.setModel(self.modelFigure)
        self.tableViewFigure.resizeColumnToContents(0)
        self.tableViewFigure.resizeRowsToContents()
        self.tableViewFigure.setAlternatingRowColors(True)
        self.tableViewFigure.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewFigure.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewFigure.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewFigure.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFigure.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        delegateInt = IntDelegate(self.tableViewFigure)
        self.tableViewFigure.setItemDelegateForColumn(4, delegateInt)
        self.tableViewFigure.setItemDelegateForColumn(5, delegateInt)
        delegate_format = FormatFigureDelegate(self.tableViewFigure)
        self.tableViewFigure.setItemDelegateForColumn(6, delegate_format)
        labelFDelegate = LabelDelegate(self.tableViewFigure, self.mdl)
        self.tableViewFigure.setItemDelegateForColumn(1, labelFDelegate)
        self.tableViewFigure.setItemDelegateForColumn(3, labelFDelegate)

        # model for treeViewMeasurement
        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.setAlternatingRowColors(True)
        self.treeViewMeasurement.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewMeasurement.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewMeasurement.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewMeasurement.expandAll()
        self.treeViewMeasurement.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewMeasurement.setDragEnabled(False)

        # model for treeViewCase
        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.setAlternatingRowColors(True)
        self.treeViewCases.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewCases.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewCases.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewCases.expandAll()
        self.treeViewCases.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewCases.setDragEnabled(False)


    @pyqtSlot(str)
    def slotStudy(self, text):
        """
        INPUT choice of study
        """
        self.current_study = self.modelStudy.dicoV2M[str(text)]
        self.updateStudyView()


    def slotChangeSelectionFigure(self):
        """
        """
        self.pushButtonAssociatedSubplot.setEnabled(False)
        self.pushButtonDeleteFigure.setEnabled(False)
        current = self.tableViewFigure.currentIndex()
        idx = current.row()
        if current != self.tableViewFigure.rootIndex():
            if current.parent() == self.tableViewFigure.rootIndex():
                self.pushButtonAssociatedSubplot.setEnabled(True)
                self.pushButtonDeleteFigure.setEnabled(True)


    def slotChangeSelectionMeasurement(self):
        """
        """
        self.pushButtonAssociatedMeasurementSubplot.setEnabled(False)
        current = self.treeViewMeasurement.currentIndex()
        idx = current.row()
        if current == self.treeViewMeasurement.rootIndex():
            self.pushButtonDeleteFile.setEnabled(False)
            self.pushButtonAddPlot.setEnabled(False)
            self.pushButtonDeleteFile.setEnabled(False)
            self.pushButtonDeletePlot.setEnabled(False)
            self.pushButtonModifyPlot.setEnabled(False)
        elif current.parent() == self.treeViewMeasurement.rootIndex():
            # Measurement file selected
            self.pushButtonDeleteFile.setEnabled(True)
            self.pushButtonAddPlot.setEnabled(True)
            self.pushButtonDeleteFile.setEnabled(True)
            self.pushButtonDeletePlot.setEnabled(False)
            self.pushButtonModifyPlot.setEnabled(False)
            self.pushButtonAssociatedMeasurementSubplot.setEnabled(False)
        else:
            # plot selected
            self.pushButtonDeleteFile.setEnabled(False)
            self.pushButtonAddPlot.setEnabled(True)
            self.pushButtonDeletePlot.setEnabled(True)
            self.pushButtonModifyPlot.setEnabled(True)
            self.pushButtonAssociatedMeasurementSubplot.setEnabled(True)
            self.pushButtonDeleteFile.setEnabled(False)


    def slotChangeSelectionCases(self):
        """
        """
        self.pushButtonAssociatedCaseSubplot.setEnabled(False)
        current = self.treeViewCases.currentIndex()
        idx = current.row()
        if current == self.treeViewCases.rootIndex():
            self.pushButtonAddData.setEnabled(False)
            self.pushButtonDeleteData.setEnabled(False)
            self.pushButtonAddCasePlot.setEnabled(False)
            self.pushButtonDeleteCasePlot.setEnabled(False)
            self.pushButtonModifyCasePlot.setEnabled(False)
        elif current.parent() == self.treeViewCases.rootIndex():
            # case id selected
            self.pushButtonAddData.setEnabled(True)
            self.pushButtonDeleteData.setEnabled(False)
            self.pushButtonAddCasePlot.setEnabled(False)
            self.pushButtonDeleteCasePlot.setEnabled(False)
            self.pushButtonModifyCasePlot.setEnabled(False)
        elif current.parent().parent() == self.treeViewCases.rootIndex():
            # data file id selected
            self.pushButtonAddData.setEnabled(True)
            self.pushButtonDeleteData.setEnabled(True)
            self.pushButtonAddCasePlot.setEnabled(True)
            self.pushButtonDeleteCasePlot.setEnabled(False)
            self.pushButtonModifyCasePlot.setEnabled(False)
        else:
            # plot is selected
            self.pushButtonAddData.setEnabled(False)
            self.pushButtonDeleteData.setEnabled(False)
            self.pushButtonAddCasePlot.setEnabled(True)
            self.pushButtonDeleteCasePlot.setEnabled(True)
            self.pushButtonModifyCasePlot.setEnabled(True)
            self.pushButtonAssociatedCaseSubplot.setEnabled(True)


    def slotAddSubplot(self):
        """
        public slot
        """
        idx = self.mdl.addSubplot(self.current_study)
        self.modelSubplot.addSubplot(idx)


    def slotDeleteSubplot(self):
        """
        public slot
        """
        idx = self.tableViewSubplot.currentIndex().row()
        self.mdl.delSubplot(self.current_study, idx)

        self.tableViewSubplot.clearSelection()

        self.modelSubplot = StandardItemModelSubplot(self.mdl, self.current_study)
        self.tableViewSubplot.setModel(self.modelSubplot)


    def slotAddFigure(self):
        """
        public slot
        """
        idx = self.mdl.addFigure(self.current_study)
        self.modelFigure.addFigure(idx)


    def slotDeleteFigure(self):
        """
        public slot
        """
        idx = self.tableViewFigure.currentIndex().row()
        self.mdl.delFigure(self.current_study, idx)

        self.modelFigure = StandardItemModelFigure(self.mdl, self.current_study)
        self.tableViewFigure.setModel(self.modelFigure)
        self.tableViewFigure.clearSelection()
        self.slotChangeSelectionFigure()


    @pyqtSlot()
    def slotAssociatedFiguresSubplot(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        idx = self.tableViewFigure.currentIndex().row()

        default = {}
        default['idlist'] = self.mdl.getFigureIdList(self.current_study, idx)
        log.debug("slotAssociatedFiguresSubplot -> %s" % str(default))

        dialog = ManageSubplotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedFiguresSubplot -> %s" % str(result))
            self.mdl.setFigureIdList(self.current_study, idx, result['idlist'])

        self.modelFigure = StandardItemModelFigure(self.mdl, self.current_study)
        self.tableViewFigure.setModel(self.modelFigure)


    @pyqtSlot()
    def slotAssociatedMeasurementSubplot(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        idx = self.treeViewMeasurement.currentIndex().row()
        measurement_idx = self.treeViewMeasurement.currentIndex().parent().row()

        default = {}
        default['idlist'] = self.mdl.getMeasurementIdList(self.current_study, measurement_idx, idx)
        log.debug("slotAssociatedMeasurementSubplot -> %s" % str(default))

        dialog = ManageSubplotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedMeasurementSubplot -> %s" % str(result))
            self.mdl.setMeasurementIdList(self.current_study, measurement_idx, idx, result['idlist'])

        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.clearSelection()
        self.treeViewMeasurement.expandAll()
        self.slotChangeSelectionMeasurement()


    def slotAddMeasurementFile(self):
        """
        public slot
        """
        # file open
        cur_path = os.getcwd()
        rep = os.path.abspath(os.path.join(self.mdl.repo, self.current_study, "POST"))
        os.chdir(rep)
        title = self.tr("measurement file")
        filetypes = self.tr("All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, rep, filetypes)[0]
        file = str(file)
        os.chdir(cur_path)

        if not file:
            return
        file = os.path.basename(file)

        self.mdl.addMeasurementFile(self.current_study, file)
        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.clearSelection()
        self.treeViewMeasurement.expandAll()
        self.slotChangeSelectionMeasurement()


    def slotDeleteMeasurementFile(self):
        """
        public slot
        """
        idx = self.treeViewMeasurement.currentIndex().row()
        self.mdl.delMeasurementFile(self.current_study, idx)

        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.clearSelection()
        self.treeViewMeasurement.expandAll()
        self.slotChangeSelectionMeasurement()


    def slotMeasurementAddPlot(self):
        """
        public slot
        """
        current = self.treeViewMeasurement.currentIndex()
        measurement_idx = 0
        if current.parent() == self.treeViewMeasurement.rootIndex():
            measurement_idx = current.row()
        else:
            measurement_idx = current.parent().row()

        idx = self.mdl.addMeasurementPlot(self.current_study, measurement_idx)

        default = {}
        default['color'] = self.mdl.getMeasurementColor(self.current_study, measurement_idx, idx)
        default['format'] = self.mdl.getMeasurementFormat(self.current_study, measurement_idx, idx)
        default['legend'] = self.mdl.getMeasurementLegend(self.current_study, measurement_idx, idx)
        default['xcol'] = self.mdl.getMeasurementXcol(self.current_study, measurement_idx, idx)
        default['ycol'] = self.mdl.getMeasurementYcol(self.current_study, measurement_idx, idx)
        default['width'] = self.mdl.getMeasurementWidth(self.current_study, measurement_idx, idx)
        default['marker'] = self.mdl.getMeasurementMarker(self.current_study, measurement_idx, idx)
        default['xerr'] = self.mdl.getMeasurementXerr(self.current_study, measurement_idx, idx)
        default['xerrp'] = self.mdl.getMeasurementXerrp(self.current_study, measurement_idx, idx)
        default['yerr'] = self.mdl.getMeasurementYerr(self.current_study, measurement_idx, idx)
        default['yerrp'] = self.mdl.getMeasurementYerrp(self.current_study, measurement_idx, idx)
        log.debug("slotAssociatedMeasurementSubplot -> %s" % str(default))

        dialog = ManagePlotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedMeasurementSubplot -> %s" % str(result))
            if result['color'] != "":
                self.mdl.setMeasurementColor(self.current_study, measurement_idx, idx, result['color'])
            if result['format'] != "":
                self.mdl.setMeasurementFormat(self.current_study, measurement_idx, idx, result['format'])
            if result['legend'] != "":
                self.mdl.setMeasurementLegend(self.current_study, measurement_idx, idx, result['legend'])
            if result['xcol'] != "":
                self.mdl.setMeasurementXcol(self.current_study, measurement_idx, idx, result['xcol'])
            if result['ycol'] != "":
                self.mdl.setMeasurementYcol(self.current_study, measurement_idx, idx, result['ycol'])
            if result['width'] != "":
                self.mdl.setMeasurementWidth(self.current_study, measurement_idx, idx, result['width'])
            if result['marker'] != "":
                self.mdl.setMeasurementMarker(self.current_study, measurement_idx, idx, result['marker'])
            if result['xerr'] != "":
                self.mdl.setMeasurementXerr(self.current_study, measurement_idx, idx, result['xerr'])
            if result['xerrp'] != "":
                self.mdl.setMeasurementXerrp(self.current_study, measurement_idx, idx, result['xerrp'])
            if result['yerr'] != "":
                self.mdl.setMeasurementYerr(self.current_study, measurement_idx, idx, result['yerr'])
            if result['yerrp'] != "":
                self.mdl.setMeasurementYerrp(self.current_study, measurement_idx, idx, result['yerrp'])

        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.clearSelection()
        self.treeViewMeasurement.expandAll()
        self.slotChangeSelectionMeasurement()


    def slotMeasurementDeletePlot(self):
        """
        public slot
        """
        idx = self.treeViewMeasurement.currentIndex().row()
        measurement_idx = self.treeViewMeasurement.currentIndex().parent().row()
        self.mdl.deleteMeasurementPlot(self.current_study, measurement_idx, idx)

        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.clearSelection()
        self.treeViewMeasurement.expandAll()
        self.slotChangeSelectionMeasurement()


    @pyqtSlot()
    def slotMeasurementPlot(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        idx = self.treeViewMeasurement.currentIndex().row()
        measurement_idx = self.treeViewMeasurement.currentIndex().parent().row()

        default = {}
        default['color'] = self.mdl.getMeasurementColor(self.current_study, measurement_idx, idx)
        default['format'] = self.mdl.getMeasurementFormat(self.current_study, measurement_idx, idx)
        default['legend'] = self.mdl.getMeasurementLegend(self.current_study, measurement_idx, idx)
        default['xcol'] = self.mdl.getMeasurementXcol(self.current_study, measurement_idx, idx)
        default['ycol'] = self.mdl.getMeasurementYcol(self.current_study, measurement_idx, idx)
        default['width'] = self.mdl.getMeasurementWidth(self.current_study, measurement_idx, idx)
        default['marker'] = self.mdl.getMeasurementMarker(self.current_study, measurement_idx, idx)
        default['xerr'] = self.mdl.getMeasurementXerr(self.current_study, measurement_idx, idx)
        default['xerrp'] = self.mdl.getMeasurementXerrp(self.current_study, measurement_idx, idx)
        default['yerr'] = self.mdl.getMeasurementYerr(self.current_study, measurement_idx, idx)
        default['yerrp'] = self.mdl.getMeasurementYerrp(self.current_study, measurement_idx, idx)
        log.debug("slotAssociatedMeasurementSubplot -> %s" % str(default))

        dialog = ManagePlotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedMeasurementSubplot -> %s" % str(result))
            if result['color'] != "":
                self.mdl.setMeasurementColor(self.current_study, measurement_idx, idx, result['color'])
            if result['format'] != "":
                self.mdl.setMeasurementFormat(self.current_study, measurement_idx, idx, result['format'])
            if result['legend'] != "":
                self.mdl.setMeasurementLegend(self.current_study, measurement_idx, idx, result['legend'])
            if result['xcol'] != "":
                self.mdl.setMeasurementXcol(self.current_study, measurement_idx, idx, result['xcol'])
            if result['ycol'] != "":
                self.mdl.setMeasurementYcol(self.current_study, measurement_idx, idx, result['ycol'])
            if result['width'] != "":
                self.mdl.setMeasurementWidth(self.current_study, measurement_idx, idx, result['width'])
            if result['marker'] != "":
                self.mdl.setMeasurementMarker(self.current_study, measurement_idx, idx, result['marker'])
            if result['xerr'] != "":
                self.mdl.setMeasurementXerr(self.current_study, measurement_idx, idx, result['xerr'])
            if result['xerrp'] != "":
                self.mdl.setMeasurementXerrp(self.current_study, measurement_idx, idx, result['xerrp'])
            if result['yerr'] != "":
                self.mdl.setMeasurementYerr(self.current_study, measurement_idx, idx, result['yerr'])
            if result['yerrp'] != "":
                self.mdl.setMeasurementYerrp(self.current_study, measurement_idx, idx, result['yerrp'])

        self.modelMeasurement = StandardItemModelMeasurement(self.mdl, self.current_study)
        self.treeViewMeasurement.setModel(self.modelMeasurement)
        self.treeViewMeasurement.expandAll()
        self.treeViewMeasurement.clearSelection()
        self.slotChangeSelectionMeasurement()


    @pyqtSlot()
    def slotCasePlot(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        idx = self.treeViewCases.currentIndex().row()
        data_idx = self.treeViewCases.currentIndex().parent().row()
        case_idx = self.treeViewCases.currentIndex().parent().parent().row()
        plot_idx = self.mdl.getCasePlotId(self.current_study, case_idx, data_idx, idx)

        default = {}
        default['color'] = self.mdl.getAssociatedCaseColor(self.current_study, case_idx, data_idx, plot_idx)
        default['format'] = self.mdl.getAssociatedCaseFormat(self.current_study, case_idx, data_idx, plot_idx)
        default['legend'] = self.mdl.getAssociatedCaseLegend(self.current_study, case_idx, data_idx, plot_idx)
        default['xcol'] = self.mdl.getAssociatedCaseXcol(self.current_study, case_idx, data_idx, plot_idx)
        default['ycol'] = self.mdl.getAssociatedCaseYcol(self.current_study, case_idx, data_idx, plot_idx)
        default['width'] = self.mdl.getAssociatedCaseWidth(self.current_study, case_idx, data_idx, plot_idx)
        default['marker'] = self.mdl.getAssociatedCaseMarker(self.current_study, case_idx, data_idx, plot_idx)
        default['xerr'] = self.mdl.getAssociatedCaseXerr(self.current_study, case_idx, data_idx, plot_idx)
        default['xerrp'] = self.mdl.getAssociatedCaseXerrp(self.current_study, case_idx, data_idx, plot_idx)
        default['yerr'] = self.mdl.getAssociatedCaseYerr(self.current_study, case_idx, data_idx, plot_idx)
        default['yerrp'] = self.mdl.getAssociatedCaseYerrp(self.current_study, case_idx, data_idx, plot_idx)
        log.debug("slotAssociatedCaseSubplot -> %s" % str(default))

        dialog = ManagePlotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedCaseSubplot -> %s" % str(result))
            if result['color'] != "":
                self.mdl.setAssociatedCaseColor(self.current_study, case_idx, data_idx, plot_idx, result['color'])
            if result['format'] != "":
                self.mdl.setAssociatedCaseFormat(self.current_study, case_idx, data_idx, plot_idx, result['format'])
            if result['legend'] != "":
                self.mdl.setAssociatedCaseLegend(self.current_study, case_idx, data_idx, plot_idx, result['legend'])
            if result['xcol'] != "":
                self.mdl.setAssociatedCaseXcol(self.current_study, case_idx, data_idx, plot_idx, result['xcol'])
            if result['ycol'] != "":
                self.mdl.setAssociatedCaseYcol(self.current_study, case_idx, data_idx, plot_idx, result['ycol'])
            if result['width'] != "":
                self.mdl.setAssociatedCaseWidth(self.current_study, case_idx, data_idx, plot_idx, result['width'])
            if result['marker'] != "":
                self.mdl.setAssociatedCaseMarker(self.current_study, case_idx, data_idx, plot_idx, result['marker'])
            if result['xerr'] != "":
                self.mdl.setAssociatedCaseXerr(self.current_study, case_idx, data_idx, plot_idx, result['xerr'])
            if result['xerrp'] != "":
                self.mdl.setAssociatedCaseXerrp(self.current_study, case_idx, data_idx, plot_idx, result['xerrp'])
            if result['yerr'] != "":
                self.mdl.setAssociatedCaseYerr(self.current_study, case_idx, data_idx, plot_idx, result['yerr'])
            if result['yerrp'] != "":
                self.mdl.setAssociatedCaseYerrp(self.current_study, case_idx, data_idx, plot_idx, result['yerrp'])

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


    @pyqtSlot()
    def slotCaseAddData(self):
        """
        """
        current = self.treeViewCases.currentIndex()
        if current.parent() == self.treeViewCases.rootIndex():
            case_idx = current.row()
        else:
            case_idx = current.parent().row()

        title = self.tr("File name")
        label = self.tr("post processing file name")
        name = QInputDialog.getText(self, title, label, QLineEdit.Normal)[0]
        self.mdl.addCaseDataFile(self.current_study, case_idx, name)

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


    @pyqtSlot()
    def slotCaseDeleteData(self):
        """
        """
        idx = self.treeViewCases.currentIndex().row()
        case_idx = self.treeViewCases.currentIndex().parent().row()
        self.mdl.delCaseDataFile(self.current_study, case_idx, idx)

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


    @pyqtSlot()
    def slotCaseAddPlot(self):
        """
        """
        current = self.treeViewCases.currentIndex()
        if current.parent().parent() == self.treeViewCases.rootIndex():
            case_idx = current.parent().row()
            data_idx = current.row()
        else:
            case_idx = current.parent().parent().row()
            data_idx = current.parent().row()

        idx = self.mdl.addAssociatedCasePlot(self.current_study, case_idx, data_idx)

        default = {}
        default['color'] = self.mdl.getAssociatedCaseColor(self.current_study, case_idx, data_idx, idx)
        default['format'] = self.mdl.getAssociatedCaseFormat(self.current_study, case_idx, data_idx, idx)
        default['legend'] = self.mdl.getAssociatedCaseLegend(self.current_study, case_idx, data_idx, idx)
        default['xcol'] = self.mdl.getAssociatedCaseXcol(self.current_study, case_idx, data_idx, idx)
        default['ycol'] = self.mdl.getAssociatedCaseYcol(self.current_study, case_idx, data_idx, idx)
        default['width'] = self.mdl.getAssociatedCaseWidth(self.current_study, case_idx, data_idx, idx)
        default['marker'] = self.mdl.getAssociatedCaseMarker(self.current_study, case_idx, data_idx, idx)
        default['xerr'] = self.mdl.getAssociatedCaseXerr(self.current_study, case_idx, data_idx, idx)
        default['xerrp'] = self.mdl.getAssociatedCaseXerrp(self.current_study, case_idx, data_idx, idx)
        default['yerr'] = self.mdl.getAssociatedCaseYerr(self.current_study, case_idx, data_idx, idx)
        default['yerrp'] = self.mdl.getAssociatedCaseYerrp(self.current_study, case_idx, data_idx, idx)
        log.debug("slotAssociatedCaseSubplot -> %s" % str(default))

        dialog = ManagePlotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedCaseSubplot -> %s" % str(result))
            if result['color'] != "":
                self.mdl.setAssociatedCaseColor(self.current_study, case_idx, data_idx, idx, result['color'])
            if result['format'] != "":
                self.mdl.setAssociatedCaseFormat(self.current_study, case_idx, data_idx, idx, result['format'])
            if result['legend'] != "":
                self.mdl.setAssociatedCaseLegend(self.current_study, case_idx, data_idx, idx, result['legend'])
            if result['xcol'] != "":
                self.mdl.setAssociatedCaseXcol(self.current_study, case_idx, data_idx, idx, result['xcol'])
            if result['ycol'] != "":
                self.mdl.setAssociatedCaseYcol(self.current_study, case_idx, data_idx, idx, result['ycol'])
            if result['width'] != "":
                self.mdl.setAssociatedCaseWidth(self.current_study, case_idx, data_idx, idx, result['width'])
            if result['marker'] != "":
                self.mdl.setAssociatedCaseMarker(self.current_study, case_idx, data_idx, idx, result['marker'])
            if result['xerr'] != "":
                self.mdl.setAssociatedCaseXerr(self.current_study, case_idx, data_idx, idx, result['xerr'])
            if result['xerrp'] != "":
                self.mdl.setAssociatedCaseXerrp(self.current_study, case_idx, data_idx, idx, result['xerrp'])
            if result['yerr'] != "":
                self.mdl.setAssociatedCaseYerr(self.current_study, case_idx, data_idx, idx, result['yerr'])
            if result['yerrp'] != "":
                self.mdl.setAssociatedCaseYerrp(self.current_study, case_idx, data_idx, idx, result['yerrp'])

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


    @pyqtSlot()
    def slotCaseDeletePlot(self):
        """
        """
        idx = self.treeViewCases.currentIndex().row()
        data_idx = self.treeViewCases.currentIndex().parent().row()
        case_idx = self.treeViewCases.currentIndex().parent().parent().row()
        self.mdl.delAssociatedCasePlot(self.current_study, case_idx, data_idx, idx)

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


    @pyqtSlot()
    def slotCaseAssociatedSubplot(self):
        """
        """
        idx = self.treeViewCases.currentIndex().row()
        data_idx = self.treeViewCases.currentIndex().parent().row()
        case_idx = self.treeViewCases.currentIndex().parent().parent().row()
        plot_idx = self.mdl.getCasePlotId(self.current_study, case_idx, data_idx, idx)

        default = {}
        default['idlist'] = self.mdl.getCaseIdList(self.current_study, case_idx, data_idx, plot_idx)
        log.debug("slotAssociatedCaseSubplot -> %s" % str(default))

        dialog = ManageSubplotDialogView(self, self.case, self.mdl, self.current_study, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAssociatedCaseSubplot -> %s" % str(result))
            self.mdl.setCaseIdList(self.current_study, case_idx, data_idx, plot_idx, result['idlist'])

        self.modelCase = StandardItemModelCase(self.mdl, self.current_study)
        self.treeViewCases.setModel(self.modelCase)
        self.treeViewCases.expandAll()
        self.treeViewCases.clearSelection()
        self.slotChangeSelectionCases()


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
