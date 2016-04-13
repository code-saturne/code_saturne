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
This module contains the following classes:
- StandardItemModelAverage
- TimeAveragesView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string
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

from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator
from code_saturne.Base.QtPage import ComboModel, to_qvariant, from_qvariant
from code_saturne.Pages.TimeAveragesForm import Ui_TimeAveragesForm
from code_saturne.Pages.StartRestartModel import StartRestartModel
from code_saturne.Pages.TimeAveragesModel import TimeAveragesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TimeAveragesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class for time averages
#-------------------------------------------------------------------------------

class StandardItemModelAverage(QStandardItemModel):

    def __init__(self, parent):
        """
        """
        QStandardItemModel.__init__(self)
        self.parent = parent
        self.headers = [self.tr("Number"),     self.tr("Average name"), self.tr("Start"),
                        self.tr("Time start"), self.tr("Restart"),      self.tr("Variables")]
        self.setColumnCount(len(self.headers))
        self.dataAverage = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            return to_qvariant(self.dataAverage[index.row()][index.column()])
        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)
        return to_qvariant()


    def setData(self, index, value, role):
        self.dataChanged.emit(index, index)
        return True


    def addItem(self, label, ntdmom, ttdmom, imoold, lst):
        """
        Add a row in the table.
        """
        row = self.rowCount()
        imom = row + 1
        item = [imom, label, ntdmom, ttdmom, imoold, lst]
        self.dataAverage.append(item)
        if row +1 > 50:
            title = self.tr("Information")
            msg = self.tr("The maximal number of time averages cannot exceed 50. ")
            QMessageBox.information(self.parent, title, msg)
        else:
            self.setRowCount(row+1)


    def replaceItem(self, row, label, ntdmom, ttdmom, imoold, lst):
        """
        Replace a row in the table.
        """
        imom = row + 1
        self.dataAverage[row] = [imom, label, ntdmom, ttdmom, imoold, lst]


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
        self.modelAverage = StandardItemModelAverage(self)
        self.treeViewAverage.setModel(self.modelAverage)
        self.treeViewAverage.resizeColumnToContents(0)
        self.treeViewAverage.resizeColumnToContents(1)
        self.treeViewAverage.resizeColumnToContents(2)
        self.treeViewAverage.resizeColumnToContents(3)
        self.treeViewAverage.resizeColumnToContents(4)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.listViewDrag.setModel(self.modelDrag)
        self.listViewDrop.setModel(self.modelDrop)
        self.listViewDrag.setAlternatingRowColors(True)
        self.listViewDrag.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.listViewDrop.setAlternatingRowColors(True)
        self.listViewDrop.setEditTriggers(QAbstractItemView.NoEditTriggers)

        self.modelIMOOLD  = ComboModel(self.comboBoxIMOOLD, 3, 1)
        self.modelIMOOLD.addItem(self.tr('automatic'), 'automatic')
        self.modelIMOOLD.addItem(self.tr('reset'), 'reset')
        self.modelIMOOLD.addItem(self.tr('specified'), 'specified')

        self.modelStartType  = ComboModel(self.comboBoxStartType, 2, 1)
        self.modelStartType.addItem(self.tr('time'), 'time')
        self.modelStartType.addItem(self.tr('iteration'), 'iteration')

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAddAverage)
        self.pushButtonDelete.clicked.connect(self.slotdeleteTimeAverage)
        self.pushButtonAddVar.clicked.connect(self.slotAddVarAverage)
        self.pushButtonSuppressVar.clicked.connect(self.slotDeleteVarAverage)
        self.treeViewAverage.pressed[QModelIndex].connect(self.slotSelectAverage)
        self.lineEditStart.textChanged[str].connect(self.slotStart)
        self.lineEditStartTime.textChanged[str].connect(self.slotStartTime)
        self.comboBoxIMOOLD.activated[str].connect(self.slotRestartChoice)
        self.comboBoxStartType.activated[str].connect(self.slotTimeChoice)
        self.lineEditRestart.textChanged[str].connect(self.slotRestart)
        self.lineEditAverage.textChanged[str].connect(self.slotBaseName)

        # Validators
        validatorStart = IntValidator(self.lineEditStart, min=-1)
        self.lineEditStart.setValidator(validatorStart)

        validatorStartTime = DoubleValidator(self.lineEditStartTime, min=-1.)
        self.lineEditStartTime.setValidator(validatorStartTime)

        validatorRestart = IntValidator(self.lineEditRestart, min=-2, max=50)
        validatorRestart.setExclusiveValues([0])
        self.lineEditRestart.setValidator(validatorRestart)

        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorLabel =  RegExpValidator(self.lineEditAverage, QRegExp(rx))
        self.lineEditAverage.setValidator(validatorLabel)

        # Initialize

        # Update list of variables, properties, scalars ...

        lst_label = [str(label) for label in list(self.mdl.dicoLabel2Name.keys())]

        self.modelDrag.setStringList(sorted(lst_label, key=str.lower))

        # Is it a following calculation ?

        if not StartRestartModel(self.case).getRestartPath():
            self.labelRestart.setDisabled(True)
            self.comboBoxIMOOLD.setDisabled(True)
            self.lineEditRestart.setDisabled(True)
            self.treeViewAverage.hideColumn(4)
        else:
            self.slotRestartChoice("automatic")

        # Update list of averages for view from xml file

        for nb in range(self.mdl.getNumberOfTimeAverage()):
            self.entriesNumber = self.entriesNumber + 1
            label, start, timestart, restart, lst = self.mdl.getTimeAverageData(nb+1)
            self.insertAverage(label, start, timestart, restart, lst)

        self.groupBoxTimeAverage.hide()

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotStart(self, text):
        """
        Return an integer for ntdmom, value of start of calculation.
        """
        if self.lineEditStart.validator().state == QValidator.Acceptable:
            start = from_qvariant(text, int)
        else:
            start = self.mdl.defaultValues()['start']
        self.mdl.setTimeStepStart(self.label_select, start)
        self.updateView()


    @pyqtSlot(str)
    def slotStartTime(self, text):
        """
        Return an float for ttdmom, value of start of calculation.
        """
        if self.lineEditStartTime.validator().state == QValidator.Acceptable:
            start = from_qvariant(text, float)
            timestart = start
        else:
            timestart = self.mdl.defaultValues()['timestart']
        self.mdl.setTimeStart(self.label_select, timestart)
        self.updateView()


    @pyqtSlot(str)
    def slotRestartChoice(self, text):
        choice = self.modelIMOOLD.dicoV2M[str(text)]
        if choice == "automatic":
            restart = -2
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(str(restart))
        elif choice == "reset":
            restart = -1
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(str(restart))
        elif choice == "specified":
            restart = self.mdl.defaultValues()['restart']
            self.lineEditRestart.setDisabled(False)
            self.lineEditRestart.setText(str(restart))
        if self.label_select:
            self.mdl.setRestart(self.label_select, restart)
        self.updateView()


    @pyqtSlot(str)
    def slotTimeChoice(self, text):
        choice = self.modelStartType.dicoV2M[str(text)]
        self.modelStartType.setItem(str_model=choice)
        if choice == "time":
            start = -1
            self.mdl.setTimeStepStart(self.label_select, start)
            timestart = self.mdl.getTimeStart(self.label_select)
            self.lineEditStart.setText(str(start))
            self.lineEditStartTime.setText(str(timestart))
            self.lineEditStart.hide()
            self.labelStart.hide()
            self.lineEditStartTime.show()
            self.labelStartTime.show()
        elif choice == "iteration":
            timestart = -1.
            self.mdl.setTimeStart(self.label_select, timestart)
            start = self.mdl.getTimeStepStart(self.label_select)
            self.lineEditStart.setText(str(start))
            self.lineEditStartTime.setText(str(timestart))
            self.lineEditStart.show()
            self.labelStart.show()
            self.lineEditStartTime.hide()
            self.labelStartTime.hide()
        self.updateView()

    @pyqtSlot(str)
    def slotRestart(self, text):
        """
        Return an integer for imoold, value of restart of calculation.
        """
        if self.lineEditRestart.validator().state == QValidator.Acceptable:
            restart = from_qvariant(text, int)
        else:
            restart = self.mdl.defaultValues()['restart']
        if self.label_select:
            self.mdl.setRestart(self.label_select, restart)
        self.updateView()


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

        if imoold == self.mdl.defaultValues()['restart']:
            imoold = ""
        self.modelAverage.addItem(label, ntdmom, ttdmom, imoold, idfmom_view)


    @pyqtSlot()
    def slotAddAverage(self):
        """
        Set in view IMOM, NTDMOM, IMOOLD, IDFMOM
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

        self.lineEditAverage.setText(str(label))
        self.lineEditStart.setText(str(ntdmom))
        self.lineEditStartTime.setText(str(ttdmom))

        if ntdmom == -1:
            self.slotTimeChoice("time")
        else:
            self.slotTimeChoice("iteration")

        if StartRestartModel(self.case).getRestartPath():
            if not imoold or imoold == -2:
                imoold = -2
                choice = "automatic"
            elif imoold == -1:
                choice = "reset"
            else:
                choice = "specified"
            self.slotRestartChoice(str(choice))
            self.modelIMOOLD.setItem(str_model=choice)
            self.lineEditRestart.setText(str(imoold))

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
        self.mdl.setVariable(self.label_select, liste)

        self.updateView()


    @pyqtSlot(str)
    def slotBaseName(self, text):
        """
        """
        if self.lineEditAverage.validator().state == QValidator.Acceptable:
            self.mdl.setLabel(self.label_select, str(text))
            self.label_select = str(text)

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
