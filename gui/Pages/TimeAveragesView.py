# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import LABEL_LENGTH_MAX
from Base.Toolbox import GuiParam
from Base.QtPage import IntValidator, DoubleValidator, RegExpValidator
from Base.QtPage import ComboModel, to_qvariant, from_qvariant
from Pages.TimeAveragesForm import Ui_TimeAveragesForm
from Pages.StartRestartModel import StartRestartModel
from Pages.TimeAveragesModel import TimeAveragesModel

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
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
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
        self.start = 1
        self.timestart = 0.
        self.restart = -2

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
        self.connect(self.pushButtonAdd,    SIGNAL("clicked()"), self.slotAddAverage)
        self.connect(self.pushButtonEdit,   SIGNAL("clicked()"), self.slotEditAverage)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.slotdeleteTimeAverage)
        self.connect(self.pushButtonAddVar,      SIGNAL("clicked()"), self.slotAddVarAverage)
        self.connect(self.pushButtonSuppressVar, SIGNAL("clicked()"), self.slotDeleteVarAverage)
        self.connect(self.treeViewAverage,  SIGNAL("pressed(const QModelIndex &)"), self.slotSelectAverage)
        self.connect(self.lineEditStart, SIGNAL("textChanged(const QString &)"), self.slotStart)
        self.connect(self.lineEditStartTime, SIGNAL("textChanged(const QString &)"), self.slotStartTime)
        self.connect(self.comboBoxIMOOLD, SIGNAL("activated(const QString&)"), self.slotRestartChoice)
        self.connect(self.comboBoxStartType, SIGNAL("activated(const QString&)"), self.slotTimeChoice)
        self.connect(self.lineEditRestart, SIGNAL("textChanged(const QString &)"), self.slotRestart)

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

        self.modelDrag.setStringList(lst_label)

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

        self.slotTimeChoice("iteration")

        self.case.undoStartGlobal()


    def getLabel(self):
        """
        Return label of average.
        """
        label = str(self.lineEditAverage.text())
        if label in self.mdl.getTimeAverageLabels():
            default = {}
            default['label'] = label
            default['list'] = self.mdl.getTimeAverageLabels()
            rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
            default['regexp'] = QRegExp(rx)
            from Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
            dialog = VerifyExistenceLabelDialogView(self, default)
            if dialog.exec_():
                result = dialog.get_result()
                label = result['label']
                if result['label'] == default['label']:
                    label = ""
        return label


    @pyqtSignature("const QString&")
    def slotStart(self, text):
        """
        Return an integer for ntdmom, value of start of calculation.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            start = from_qvariant(text, int)
            self.start = start
        else:
            self.start = self.mdl.defaultValues()['start']


    @pyqtSignature("const QString&")
    def slotStartTime(self, text):
        """
        Return an float for ttdmom, value of start of calculation.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            start = from_qvariant(text, float)
            self.timestart = start
        else:
            self.timestart = self.mdl.defaultValues()['timestart']


    @pyqtSignature("const QString&")
    def slotRestartChoice(self, text):
        choice = self.modelIMOOLD.dicoV2M[str(text)]
        if choice == "automatic":
            self.restart = -2
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(str(self.restart))
        elif choice == "reset":
            self.restart = -1
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(str(self.restart))
        elif choice == "specified":
            self.restart = self.mdl.defaultValues()['restart']
            self.lineEditRestart.setDisabled(False)
            self.lineEditRestart.setText("")


    @pyqtSignature("const QString&")
    def slotTimeChoice(self, text):
        choice = self.modelStartType.dicoV2M[str(text)]
        self.modelStartType.setItem(str_model=choice)
        if choice == "time":
            self.start = -1
            self.lineEditStart.setText(str(self.start))
            self.lineEditStartTime.setText(str(self.timestart))
            self.lineEditStart.hide()
            self.labelStart.hide()
            self.lineEditStartTime.show()
            self.labelStartTime.show()
        elif choice == "iteration":
            self.timestart = -1.
            self.lineEditStart.setText(str(self.start))
            self.lineEditStartTime.setText(str(self.timestart))
            self.lineEditStart.show()
            self.labelStart.show()
            self.lineEditStartTime.hide()
            self.labelStartTime.hide()

    @pyqtSignature("const QString&")
    def slotRestart(self, text):
        """
        Return an integer for imoold, value of restart of calculation.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            restart = from_qvariant(text, int)
            self.restart = restart
        else:
            self.restart = self.mdl.defaultValues()['restart']


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
        idfmom = "*".join(lst)

        idfmom_view = "<" + idfmom +">"

        if imoold == self.mdl.defaultValues()['restart']:
            imoold = ""
        self.modelAverage.addItem(label, ntdmom, ttdmom, imoold, idfmom_view)


    def replaceTimeAverage(self, row, label, ntdmom, ttdmom, imoold, lst):
        """
        Insert values in Hlist.
        """
        idfmom = "*".join(lst)
        idfmom_view = "<" + idfmom + ">"

        if imoold == None:
            imoold = -1
        self.modelAverage.replaceItem(row, label, ntdmom, ttdmom, imoold, idfmom_view)


    @pyqtSignature("")
    def slotAddAverage(self):
        """
        Set in view IMOM, NTDMOM, IMOOLD, IDFMOM
        """
        var_prop = [str(s) for s in self.modelDrop.stringList()]
        log.debug("slotAddAverage -> %s" % (var_prop,))
        idfmom = "*".join(var_prop)

        if idfmom == '':
            title = self.tr("Warning")
            msg   = self.tr("You must select at least one variable or property from list")
            QMessageBox.information(self, title, msg)

        else:
            label = self.getLabel()

            self.entriesNumber = self.entriesNumber + 1
            if label == '':
                label = 'TimeAverage' + repr(self.entriesNumber)

            if label in self.mdl.getTimeAverageLabels():
                title = self.tr("Warning")
                msg = self.tr("This label already exists")
                QMessageBox.information(self, title, msg)
                return

            log.debug("slotAddAverage -> %s" % (label,))

            ntdmom = self.start
            ttdmom = self.timestart
            if StartRestartModel(self.case).getRestartPath():
                imoold = self.restart
            else:
                imoold = self.mdl.defaultValues()['restart']

            self.insertAverage(label, ntdmom, ttdmom, imoold, var_prop)
            average = idfmom.split('*')
            self.mdl.setTimeAverage(label, ntdmom, ttdmom, imoold, average)
        self.__eraseEntries()


    @pyqtSignature("")
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


    @pyqtSignature("")
    def slotEditAverage(self):
        """
        Edit average to modify the selected average from the list.
        """
        row = self.treeViewAverage.currentIndex().row()
        log.debug("slotEditAverage -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing time average")
            QMessageBox.information(self, title, msg)
        else:
            [imom, old_label, old_start, old_time_start, old_restart, old_average] = self.averageInfo()

            var_prof = [str(s) for s in self.modelDrop.stringList()]
            log.debug("slotEditAverage -> %s" % (var_prof,))
            if not var_prof:
                title = self.tr("Warning")
                msg   = self.tr("You must select at least one variable or property from list")
                QMessageBox.information(self, title, msg)

            else:
                new_label = str(self.lineEditAverage.text())

                if new_label == '':
                    new_label = old_label

                if new_label != old_label and new_label in self.mdl.getTimeAverageLabels():
                    title = self.tr("Warning")
                    msg = self.tr("This label already exists: the old label is kept")
                    QMessageBox.information(self, title, msg)
                    new_label = old_label

                log.debug("slotEditAverage -> %s" % (new_label,))

                new_start = self.start
                new_time_start = self.timestart
                if StartRestartModel(self.case).getRestartPath():
                    new_restart = self.restart
                else:
                    new_restart = self.mdl.defaultValues()['restart']

                self.replaceTimeAverage(row, new_label, new_start, new_time_start, new_restart, var_prof)
                self.mdl.replaceTimeAverage(old_label, new_label, new_start, new_time_start, new_restart, var_prof)
                idx = self.treeViewAverage.currentIndex()
                self.treeViewAverage.dataChanged(idx, idx)

                self.__eraseEntries()


    @pyqtSignature("const QModelIndex &")
    def slotSelectAverage(self, index):
        """
        Return the selected item from the Hlist.
        """
        row = index.row()
        log.debug("slotSelectAverage -> %s" % (row,))

        [imom, label, ntdmom, ttdmom, imoold, idfmom] = self.averageInfo()

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
        self.modelDrop.setStringList(lst)


    @pyqtSignature("")
    def slotAddVarAverage(self):
        """
        Add a new var from list to profile
        """
        if (self.listViewDrag.currentIndex().row() >=0) :
            lst = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.listViewDrag.currentIndex().row()]
            lst.append(var)
            self.modelDrop.setStringList(lst)
        self.start = 1
        self.timestart = 0.
        self.slotTimeChoice("iteration")


    @pyqtSignature("")
    def slotDeleteVarAverage(self):
        """
        Supress a var from profile
        """
        self.modelDrop.removeRows(self.listViewDrop.currentIndex().row(), 1)
        self.slotTimeChoice("iteration")
        self.start = 1
        self.timestart = 0.


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.lineEditAverage.setText(str(""))
        self.lineEditStart.setText(str(""))
        if not StartRestartModel(self.case).getRestartPath():
            self.lineEditRestart.setText(str(""))
        self.modelDrop.setStringList([])
        self.treeViewAverage.clearSelection()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
