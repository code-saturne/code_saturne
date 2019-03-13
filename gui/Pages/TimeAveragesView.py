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
This module contains the following classes:
- StandardItemModelAverage
- TimeAveragesView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator
from code_saturne.Base.QtPage import ComboModel, from_qvariant, to_text_string
from code_saturne.Pages.TimeAveragesForm import Ui_TimeAveragesForm
from code_saturne.model.StartRestartModel import StartRestartModel
from code_saturne.model.TimeAveragesModel import TimeAveragesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TimeAveragesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the name
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent, mdl):
        super(LabelDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl    = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_p_value = ""
        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_p_value = str(value)
        editor.setText(value)

    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())

            if p_value in self.mdl.getTimeAverageNames():
                title = self.tr("Warning")
                msg   = self.tr("This time moment field is already used.\n"\
                                "Please choose another one.")
                QMessageBox.information(self.parent, title, msg)
                return

            if p_value and p_value != self.old_p_value:
                self.mdl.setName(self.old_p_value, p_value)
                model.setData(index, p_value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QComboBox delegate for the start type
#-------------------------------------------------------------------------------

class StartTypeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(StartTypeDelegate, self).__init__(parent)
        self.parent = parent

        self.combo_items = ('time', 'time step')


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for k in self.combo_items:
            editor.addItem(k)
        editor.installEventFilter(self)
        editor.setMinimumWidth(90)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        str_model = index.model().getData(row, col)
        idx = self.combo_items.index(str_model)
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        txt   = str(comboBox.currentText())
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, txt, Qt.DisplayRole)


    def tr(self, text):
        return text

#-------------------------------------------------------------------------------
# Line edit delegate for the start time
#-------------------------------------------------------------------------------

class StartValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(StartValueDelegate, self).__init__(parent)
        self.parent = parent

    def createEditor(self, parent, option, index):
        row = index.row()
        stype = index.model().getData(row, 2)
        editor = QLineEdit(parent)
        if stype == 'time step':
            validator = IntValidator(editor, min=1)
        else:
            validator = DoubleValidator(editor, min=0.)
        editor.setValidator(validator)
        return editor

    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(str(value))

    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            row = index.row()
            stype = index.model().getData(row, 2)
            if stype == 'time step':
                value = from_qvariant(editor.text(), int)
            else:
                value = from_qvariant(editor.text(), float)
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QComboBox delegate for the boundary nature
#-------------------------------------------------------------------------------

class RestartTypeDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(RestartTypeDelegate, self).__init__(parent)
        self.parent   = parent
        self.dicoM2V = {'automatic': 'automatic',
                        'reset': 'reset'}
        self.dicoV2M = {}
        for k, v in list(self.dicoM2V.items()):
            self.dicoV2M[v] = k


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for k in list(self.dicoV2M.keys()):
            editor.addItem(k)
        editor.installEventFilter(self)
        editor.setMinimumWidth(80)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        str_model = index.model().getData(row, col)
        idx = list(self.dicoM2V.keys()).index(str_model)
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        txt   = str(comboBox.currentText())
        value = self.dicoV2M[txt]
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# StandarItemModel class for time averages
#-------------------------------------------------------------------------------

class StandardItemModelAverage(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Id"),
                        self.tr("Name"),
                        self.tr("Start type"),
                        self.tr("Start time"),
                        self.tr("Restart"),
                        self.tr("Variables")]
        self.setColumnCount(len(self.headers))
        self.mdl = mdl
        self.dataAverage = []


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.DisplayRole:
            return self.dataAverage[index.row()][index.column()]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            col = index.column()
            if col > 0 and col < 5:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            else:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        return None


    def setData(self, index, value, role):
        col = index.column()
        row = index.row()
        if col in (1, 2, 4):
            s_value = str(from_qvariant(value, to_text_string))
        if col == 1:
            self.dataAverage[row][col] = s_value
        elif col > 1 and col < 5:
            name = self.dataAverage[row][1]
            if col == 2:
                if s_value != self.dataAverage[row][2]:
                    if value == 'time step':
                        stime = self.mdl.defaultValues()['start']
                        self.mdl.setTimeStepStart(name, stime)
                        self.mdl.setTimeStart(name, -1.0)
                    else:
                        stime = self.mdl.defaultValues()['timestart']
                        self.mdl.setTimeStepStart(name, -1)
                        self.mdl.setTimeStart(name, stime)
                    self.dataAverage[row][3] = stime
                    self.dataAverage[row][col] = s_value
            elif col == 3:
                stype = self.dataAverage[row][2]
                if stype == 'time step':
                    stime = from_qvariant(value, int)
                    self.mdl.setTimeStepStart(name, stime)
                else:
                    stime = from_qvariant(value, float)
                    self.mdl.setTimeStart(name, stime)
                self.dataAverage[row][col] = stime
            elif col == 4:
                choice = s_value
                if choice == 'automatic':
                    restart = -2
                elif choice == 'reset':
                    restart = -1
                else: # if choice == "specified":
                    choice = 'automatic'
                    restart = -2
                self.mdl.setRestart(name, restart)
                self.dataAverage[row][col] = choice

        self.dataChanged.emit(index, index)
        return True


    def getData(self, row, column):
        return self.dataAverage[row][column]


    def addItem(self, label, ntdmom, ttdmom, imoold, lst):
        """
        Add a row in the table.
        """
        row = self.rowCount()
        imom = row

        item = [imom, label, None, None, imoold, lst]
        self.dataAverage.append(item)
        self.replaceItem(row, label, ntdmom, ttdmom, imoold, lst)
        self.setRowCount(row+1)


    def replaceItem(self, row, label, ntdmom, ttdmom, imoold, lst):
        """
        Replace a row in the table.
        """
        imom = row

        if ttdmom < 0:
            stype = 'time step'
            sval = ntdmom
        else:
            stype = 'time'
            sval = ttdmom
        if imoold == -2:
            rtype = 'automatic'
        elif imoold == -1:
            rtype = 'reset'
        else:
            rtype = 'automatic'
        self.dataAverage[row] = [imom, label, stype, sval, rtype, lst]


    def deleteRow(self, row):
        """
        Delete the row in the model.
        """
        del self.dataAverage[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def getItem(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataAverage[row]


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self.dataAverage = []
        self.setRowCount(0)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class TimeAveragesView(QWidget, Ui_TimeAveragesForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_TimeAveragesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = TimeAveragesModel(self.case)
        self.entriesNumber = 0

        self.label_select = None

        # Create the Page layout.

        # Models
        self.modelAverage = StandardItemModelAverage(self.mdl)
        self.treeViewAverage.setModel(self.modelAverage)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.listViewDrag.setModel(self.modelDrag)
        self.listViewDrop.setModel(self.modelDrop)
        self.listViewDrag.setAlternatingRowColors(True)
        self.listViewDrag.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.listViewDrop.setAlternatingRowColors(True)
        self.listViewDrop.setEditTriggers(QAbstractItemView.NoEditTriggers)

        # Delegates
        delegateLabel = LabelDelegate(self.treeViewAverage, self.mdl)
        delegateStartType = StartTypeDelegate(self.treeViewAverage)
        delegateStartValue = StartValueDelegate(self.treeViewAverage)
        delegateRestart = RestartTypeDelegate(self.treeViewAverage)
        self.treeViewAverage.setItemDelegateForColumn(1, delegateLabel)
        self.treeViewAverage.setItemDelegateForColumn(2, delegateStartType)
        self.treeViewAverage.setItemDelegateForColumn(3, delegateStartValue)
        self.treeViewAverage.setItemDelegateForColumn(4, delegateRestart)

        if QT_API == "PYQT4":
            self.treeViewAverage.header().setResizeMode(QHeaderView.ResizeToContents)
            self.treeViewAverage.header().setResizeMode(5, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.treeViewAverage.header().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.treeViewAverage.header().setSectionResizeMode(5, QHeaderView.Stretch)

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddAverage)
        self.pushButtonDelete.clicked.connect(self.slotdeleteTimeAverage)
        self.pushButtonAddVar.clicked.connect(self.slotAddVarAverage)
        self.pushButtonSuppressVar.clicked.connect(self.slotDeleteVarAverage)
        self.treeViewAverage.pressed[QModelIndex].connect(self.slotSelectAverage)

        # Validators

        # Initialize

        # Update list of variables, properties, scalars ...

        lst_label = [str(label) for label in list(self.mdl.dicoLabel2Name.keys())]

        self.modelDrag.setStringList(sorted(lst_label, key=str.lower))

        # Update list of averages for view from xml file

        for nb in range(self.mdl.getNumberOfTimeAverage()):
            self.entriesNumber = self.entriesNumber + 1
            label, start, timestart, restart, lst = self.mdl.getTimeAverageData(nb+1)
            self.insertAverage(label, start, timestart, restart, lst)

        self.groupBoxTimeAverage.hide()

        self.case.undoStartGlobal()


    def averageInfo(self):
        """
        Return info from the argument entry.
        """
        row = self.treeViewAverage.currentIndex().row()
        return self.modelAverage.getItem(row)


    def insertAverage(self, label, ntdmom, ttdmom, imoold, lst):
        """
        Insert values in Hlist.
        """
        if len(lst) > 0:
            idfmom = "*".join(lst)
            idfmom_view = "<" + idfmom +">"
        else:
            idfmom_view = ""

        self.modelAverage.addItem(label, ntdmom, ttdmom, imoold, idfmom_view)


    @pyqtSlot()
    def slotAddAverage(self):
        """
        Set in view
        """
        var_prop = []
        label, ntdmom, ttdmom, imoold = self.mdl.addTimeAverage()
        self.insertAverage(label, ntdmom, ttdmom, imoold, var_prop)
        self.__eraseEntries()


    @pyqtSlot()
    def slotdeleteTimeAverage(self):
        """
        Delete the selected average from the list (one by one).
        """
        row = self.treeViewAverage.currentIndex().row()
        log.debug("slotdeleteTimeAverage -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing time average")
            QMessageBox.information(self, title, msg)
        else:
            [imom, label, ntdmom, ttdmom, imoold, idfmom] = self.averageInfo()
            self.mdl.deleteTimeAverage(label)
            self.modelAverage.deleteAllData()
            for n in range(self.mdl.getNumberOfTimeAverage()):
                label, ntdmom, ttdmom, imoold, var_prop = self.mdl.getTimeAverageData(n+1)
                self.insertAverage(label, ntdmom, ttdmom, imoold, var_prop)

        self.__eraseEntries()


    @pyqtSlot("QModelIndex")
    def slotSelectAverage(self, index):
        """
        Return the selected item from the Hlist.
        """
        self.groupBoxTimeAverage.show()
        row = index.row()
        log.debug("slotSelectAverage -> %s" % (row,))

        [imom, label, ntdmom, ttdmom, imoold, idfmom] = self.averageInfo()
        self.label_select = label

        lst = [str(s) for s in idfmom.replace('>','').replace('<','').split('*')]
        if lst[0] == "":
            lst.remove(lst[0])
        self.modelDrop.setStringList(lst)


    @pyqtSlot()
    def slotAddVarAverage(self):
        """
        Add a new var from list to profile
        """
        if (self.listViewDrag.currentIndex().row() >=0) :
            lst = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.listViewDrag.currentIndex().row()]
            lst.append(var)
            self.modelDrop.setStringList(lst)
            lst = [str(s) for s in lst]
            if lst[0] == "":
                lst.remove(lst[0])
            self.label_select = self.averageInfo()[1]
            self.mdl.setVariable(self.label_select, lst)

            self.updateView()


    @pyqtSlot()
    def slotDeleteVarAverage(self):
        """
        Supress a var from profile
        """
        self.modelDrop.removeRows(self.listViewDrop.currentIndex().row(), 1)
        liste = self.modelDrop.stringList()
        liste = [str(s) for s in liste]
        self.label_select = self.averageInfo()[1]
        self.mdl.setVariable(self.label_select, liste)

        self.updateView()


    def updateView(self):
        row = self.treeViewAverage.currentIndex().row()
        if self.label_select:
            lst = self.mdl.getVariable(self.label_select)
            ntdmom = self.mdl.getTimeStepStart(self.label_select)
            ttdmom = self.mdl.getTimeStart(self.label_select)
            imoold = self.mdl.getRestart(self.label_select)
            idfmom = "*".join(lst)
            idfmom_view = "<" + idfmom +">"
            self.modelAverage.replaceItem(row, self.label_select, ntdmom, ttdmom, imoold, idfmom_view)


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.label_select = None
        self.treeViewAverage.clearSelection()
        self.groupBoxTimeAverage.hide()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
