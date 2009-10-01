# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
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
from TimeAveragesForm import Ui_TimeAveragesForm
from Base.QtPage import IntValidator, RegExpValidator, ComboModel
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
        self.headers = [self.tr("Number"), self.tr("Average name"), self.tr("Start"),
                        self.tr("Restart"), self.tr("Variables")]
        self.setColumnCount(len(self.headers))
        self.dataAverage = []


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataAverage[index.row()][index.column()])
        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, label, ntdmom, imoold, list):
        """
        Add a row in the table.
        """
        row = self.rowCount()
        imom = row + 1
        item = [imom, label, ntdmom, imoold, list]
        self.dataAverage.append(item)
        if row +1 > 50:
            title = self.tr("Information")
            msg = self.tr("The maximal number of time averages cannot exceed 50. ")
            QMessageBox.information(self.parent, title, msg)
        else:
            self.setRowCount(row+1)


    def replaceItem(self, row, label, ntdmom, imoold, list):
        """
        Replace a row in the table.
        """
        imom = row + 1
        self.dataAverage[row] = [imom, label, ntdmom, imoold, list]


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
        self.mdl = TimeAveragesModel(self.case)
        self.entriesNumber = 0
        self.start = 1
        self.restart = -2

        # Create the Page layout.

        # Models
        self.modelAverage = StandardItemModelAverage(self)
        self.treeViewAverage.setModel(self.modelAverage)
        self.treeViewAverage.resizeColumnToContents(0)
        self.treeViewAverage.resizeColumnToContents(1)
        self.treeViewAverage.resizeColumnToContents(2)
        self.treeViewAverage.resizeColumnToContents(3)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.listViewDrag.setModel(self.modelDrag)
        self.listViewDrop.setModel(self.modelDrop)

        # Drag ...
        self.listViewDrag.setDragDropMode(QAbstractItemView.DragOnly)
        self.listViewDrop.setDragDropOverwriteMode(False)
        self.listViewDrag.setAlternatingRowColors(True)
        self.listViewDrag.setDragEnabled(True)
        #self.listViewDrag.setAcceptDrops(True)
        # ... and Drop
        self.listViewDrop.setDragDropMode(QAbstractItemView.DragDrop)
        self.listViewDrop.setAlternatingRowColors(True)
        self.listViewDrop.setAcceptDrops(True)
        self.listViewDrop.setDragEnabled(True)
        self.listViewDrop.setDragDropOverwriteMode(False)

        self.modelIMOOLD  = ComboModel(self.comboBoxIMOOLD, 3, 1)
        self.modelIMOOLD.addItem(self.tr('automatic'), 'automatic')
        self.modelIMOOLD.addItem(self.tr('reset'), 'reset')
        self.modelIMOOLD.addItem(self.tr('specified'), 'specified')

        # Connections
        self.connect(self.pushButtonAdd,    SIGNAL("clicked()"), self.slotAddAverage)
        self.connect(self.pushButtonEdit,   SIGNAL("clicked()"), self.slotEditAverage)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.slotdeleteTimeAverage)
        self.connect(self.treeViewAverage,  SIGNAL("pressed(const QModelIndex &)"), self.slotSelectAverage)
        self.connect(self.lineEditStart, SIGNAL("textChanged(const QString &)"), self.slotStart)
        self.connect(self.comboBoxIMOOLD, SIGNAL("activated(const QString&)"), self.slotRestartChoice)
        self.connect(self.lineEditRestart, SIGNAL("textChanged(const QString &)"), self.slotRestart)

        # Validators
        validatorStart = IntValidator(self.lineEditStart, min=1)
        self.lineEditStart.setValidator(validatorStart)

        validatorRestart = IntValidator(self.lineEditRestart, min=-2, max=50)
        validatorRestart.setExclusiveValues([0])
        self.lineEditRestart.setValidator(validatorRestart)

        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorLabel =  RegExpValidator(self.lineEditAverage, QRegExp(rx))
        self.lineEditAverage.setValidator(validatorLabel)

        # Initialize

        # Update list of variables, properties, scalars ...

        liste_label = QStringList()
        for label in self.mdl.dicoLabel2Name.keys():
            liste_label.append(label)

        self.modelDrag.setStringList(liste_label)

        # Is it a following calculation ?

        if StartRestartModel(self.case).getRestart() == 'off':
            self.labelRestart.setDisabled(True)
            self.comboBoxIMOOLD.setDisabled(True)
            self.lineEditRestart.setDisabled(True)
            self.treeViewAverage.hideColumn(3)
        else:
            self.slotRestartChoice(QString("automatic"))

        # Update list of averages for view from xml file

        for nb in range(self.mdl.getNumberOfTimeAverage()):
            self.entriesNumber = self.entriesNumber + 1
            label, start, restart, list = self.mdl.getTimeAverageData(nb+1)
            self.insertAverage(label, start, restart, list)


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
            from VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
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
        start, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.start = start
        else:
            self.start = self.mdl.defaultValues()['start']


    @pyqtSignature("const QString&")
    def slotRestartChoice(self, text):
        choice = self.modelIMOOLD.dicoV2M[str(text)]
        if choice == "automatic":
            self.restart = -2
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(QString(str(self.restart)))
        elif choice == "reset":
            self.restart = -1
            self.lineEditRestart.setDisabled(True)
            self.lineEditRestart.setText(QString(str(self.restart)))
        elif choice == "specified":
            self.restart = self.mdl.defaultValues()['restart']
            self.lineEditRestart.setDisabled(False)
            self.lineEditRestart.setText(QString(""))


    @pyqtSignature("const QString&")
    def slotRestart(self, text):
        """
        Return an integer for imoold, value of restart of calculation.
        """
        restart, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.restart = restart
        else:
            self.restart = self.mdl.defaultValues()['restart']


    def averageInfo(self):
        """
        Return info from the argument entry.
        """
        row = self.treeViewAverage.currentIndex().row()
        return self.modelAverage.getItem(row)


    def insertAverage(self, label, ntdmom, imoold, list):
        """
        Insert values in Hlist.
        """
        idfmom = string.join(list,'*')
        idfmom_view = "<" + idfmom +">"

        if imoold == self.mdl.defaultValues()['restart']:
            imoold = ""
        self.modelAverage.addItem(label, ntdmom, imoold, idfmom_view)


    def replaceTimeAverage(self, row, label, ntdmom, imoold, list):
        """
        Insert values in Hlist.
        """
        idfmom = string.join(list,'*')
        idfmom_view = "<" + idfmom + ">"

        if imoold == None:
            imoold = -1
        self.modelAverage.replaceItem(row, label, ntdmom, imoold, idfmom_view)


    @pyqtSignature("")
    def slotAddAverage(self):
        """
        Set in view IMOM, NTDMOM, IMOOLD, IDFMOM
        """
        var_prop = [str(s) for s in self.modelDrop.stringList()]
        log.debug("slotAddAverage -> %s" % (var_prop,))
        idfmom = string.join(var_prop,'*')

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
            if StartRestartModel(self.case).getRestart() != 'off':
                imoold = self.restart
            else:
                imoold = self.mdl.defaultValues()['restart']

            self.insertAverage(label, ntdmom, imoold, var_prop)
            average = string.split(idfmom, '*')
            self.mdl.setTimeAverage(label, ntdmom, imoold, average)
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
            [imom, label, ntdmom, imoold, idfmom] = self.averageInfo()
            #self.modelAverage.deleteRow(row)
            self.mdl.deleteTimeAverage(label)
            self.modelAverage.deleteAllData()
            for n in range(self.mdl.getNumberOfTimeAverage()):
                label, ntdmom, imoold, var_prop = self.mdl.getTimeAverageData(n+1)
                self.insertAverage(label, ntdmom, imoold, var_prop)

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
            msg   = self.tr("You must select an existing tiem average")
            QMessageBox.information(self, title, msg)
        else:
            [imom, old_label, old_start, old_restart, old_average] = self.averageInfo()

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
                if StartRestartModel(self.case).getRestart() != 'off':
                    new_restart = self.restart
                else:
                    new_restart = self.mdl.defaultValues()['restart']

                self.replaceTimeAverage(row, new_label, new_start, new_restart, var_prof)
                self.mdl.replaceTimeAverage(old_label, new_label, new_start, new_restart, var_prof)
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

        [imom, label, ntdmom, imoold, idfmom] = self.averageInfo()

        self.lineEditAverage.setText(QString(label))
        self.lineEditStart.setText(QString(str(ntdmom)))

        if StartRestartModel(self.case).getRestart() != 'off':
            if not imoold or imoold == -2:
                imoold = -2
                choice = "automatic"
            elif imoold == -1:
                choice = "reset"
            else:
                choice = "specified"
            self.slotRestartChoice(QString(choice))
            self.modelIMOOLD.setItem(str_model=choice)
            self.lineEditRestart.setText(QString(str(imoold)))

        liste = [QString(s) for s in idfmom.replace('>','').replace('<','').split('*')]
        self.modelDrop.setStringList(liste)


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.lineEditAverage.setText(str(""))
        self.lineEditStart.setText(str(""))
        if StartRestartModel(self.case).getRestart() == 'on':
            self.lineEditRestart.setText(str(""))
        self.modelDrop.setStringList(QStringList())
        self.treeViewAverage.clearSelection()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
