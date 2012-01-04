# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
This module defines the 'Meshes and Enveloppe' page.

This module contains the following classes:
- ProfilesView
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
from Pages.ProfilesForm import Ui_ProfilesForm
from Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from Pages.ProfilesModel import ProfilesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ProfilesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelProfile(QStandardItemModel):

    def __init__(self):
        """
        """
        QStandardItemModel.__init__(self)
        self.setColumnCount(2)
        self.dataProfile = []


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataProfile[index.row()][index.column()])
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
            if section == 0:
                return QVariant(self.tr("Filename"))
            elif section == 1:
                return QVariant(self.tr("Variables"))
        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addItem(self, label, prof):
        """
        Add a row in the table.
        """
        self.dataProfile.append([label, prof])
        row = self.rowCount()
        self.setRowCount(row+1)


    def replaceItem(self, row, label, prof):
        """
        Replace a row in the table.
        """
        self.dataProfile[row] = [label, prof]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataProfile[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def getItem(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataProfile[row]


    def getLabel(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataProfile[row][0]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ProfilesView(QWidget, Ui_ProfilesForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ProfilesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = ProfilesModel(self.case)

        #  Initialize variables concerning the display of the Hlist

        self.entriesNumber = 0

        # Models
        self.modelProfile = StandardItemModelProfile()
        self.treeViewProfile.setModel(self.modelProfile)
        self.treeViewProfile.resizeColumnToContents(0)

        # QListView layout
        self.gridlayout1 = QGridLayout(self.widgetDrag)
        self.gridlayout1.setMargin(0)
        self.DragList = QListView(self.widgetDrag)
        self.gridlayout1.addWidget(self.DragList,0,0,1,1)

        self.gridlayout2 = QGridLayout(self.widgetDrop)
        self.gridlayout2.setMargin(0)
        self.DropList = QListView(self.widgetDrop)
        self.gridlayout2.addWidget(self.DropList,0,0,1,1)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.DragList.setModel(self.modelDrag)
        self.DropList.setModel(self.modelDrop)
        self.DragList.setAlternatingRowColors(True)
        self.DragList.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.DropList.setAlternatingRowColors(True)
        self.DropList.setEditTriggers(QAbstractItemView.NoEditTriggers)

        # Combo items
        self.modelFreq = ComboModel(self.comboBoxFreq, 2, 1)
        self.modelFreq.addItem(self.tr("at the end of the calculation"), "end")
        self.modelFreq.addItem(self.tr("at each 'n' time steps"), "frequency")

        self.modelFormat = ComboModel(self.comboBoxFormat, 2, 1)
        self.modelFormat.addItem(self.tr(".dat"), "DAT")
        self.modelFormat.addItem(self.tr(".csv"), "CSV")

        # Connections
        self.connect(self.treeViewProfile,       SIGNAL("pressed(const QModelIndex &)"), self.slotSelectProfile)
        self.connect(self.pushButtonAdd,         SIGNAL("clicked()"), self.slotAddProfile)
        self.connect(self.pushButtonEdit,        SIGNAL("clicked()"), self.slotEditProfile)
        self.connect(self.pushButtonDelete,      SIGNAL("clicked()"), self.slotDeleteProfile)
        self.connect(self.pushButtonAddVar,      SIGNAL("clicked()"), self.slotAddVarProfile)
        self.connect(self.pushButtonSuppressVar, SIGNAL("clicked()"), self.slotDeleteVarProfile)
        self.connect(self.comboBoxFreq,          SIGNAL("activated(const QString&)"), self.slotFrequencyType)
        self.connect(self.comboBoxFormat,        SIGNAL("activated(const QString&)"), self.slotFormatType)

        # Validators
        validatorFreq = IntValidator(self.lineEditFreq, min=0)
        validatorFreq.setExclusiveMin(True)
        self.lineEditFreq.setValidator(validatorFreq)

        validatorFloatX1 = DoubleValidator(self.lineEditX1)
        validatorFloatY1 = DoubleValidator(self.lineEditY1)
        validatorFloatZ1 = DoubleValidator(self.lineEditZ1)
        validatorFloatX2 = DoubleValidator(self.lineEditX2)
        validatorFloatY2 = DoubleValidator(self.lineEditY2)
        validatorFloatZ2 = DoubleValidator(self.lineEditZ2)

        self.lineEditX1.setValidator(validatorFloatX1)
        self.lineEditY1.setValidator(validatorFloatY1)
        self.lineEditZ1.setValidator(validatorFloatZ1)
        self.lineEditX2.setValidator(validatorFloatX2)
        self.lineEditY2.setValidator(validatorFloatY2)
        self.lineEditZ2.setValidator(validatorFloatZ2)

        rx = "[\- _A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorTitle =  RegExpValidator(self.lineEditTitle, QRegExp(rx))
        self.lineEditTitle.setValidator(validatorTitle)

        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorBaseName =  RegExpValidator(self.lineEditBaseName, QRegExp(rx))
        self.lineEditBaseName.setValidator(validatorBaseName)

        #update list of variables, properties, scalars ...
        liste_label = QStringList()
        for label in self.mdl.getVariablesAndVolumeProperties():
            liste_label.append(label)
        self.modelDrag.setStringList(liste_label)

        #update list of profiles for view from xml file
        for lab in self.mdl.getProfilesLabelsList():
            self.entriesNumber = self.entriesNumber + 1
            label, title, format, list, freq, x1, y1, z1, x2, y2, z2 = self.mdl.getProfileData(lab)
            self.__insertProfile(label, list)

        self.__eraseEntries()


    def __verifLabel(self):
        """
        Verif label.
        """
        label = str(self.lineEditBaseName.text())
        if label in self.mdl.getProfilesLabelsList():
            default = {}
            default['label'] = label
            default['list'] = self.mdl.getProfilesLabelsList()
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
    def slotFrequencyType(self, text):
        """
        Input choice for frequency for profile.
        """
        choice = self.modelFreq.dicoV2M[str(text)]

        if choice == "end":
            nfreq = -1
            self.lineEditFreq.setText(QString(str(nfreq)))
            self.lineEditFreq.setDisabled(True)
            self.labelBaseName.setText(QString("Filename"))

        if choice == "frequency":
            nfreq, ok = self.lineEditFreq.text().toInt()
            if nfreq == -1: nfreq = 1
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(QString(str(nfreq)))
            self.labelBaseName.setText(QString("Filename"))

        return choice, nfreq


    @pyqtSignature("const QString&")
    def slotFormatType(self, text):
        """
        Input choice for frequency for profile.
        """
        return self.modelFormat.dicoV2M[str(text)]


    def __infoProfile(self, row):
        """
        Return info from the argument entry.
        """
        label = self.modelProfile.getLabel(row)
        lab, title, format, list, freq, x1, y1, z1, x2, y2, z2 = self.mdl.getProfileData(label)
        return label, title, format, list, freq, x1, y1, z1, x2, y2, z2


    def __insertProfile(self, label, list):
        """
        Insert values in table view.
        """
        self.modelProfile.addItem(label, string.join(list,' ; '))



    def __replaceProfile(self, num, label, list):
        """
        Insert values in table view.
        """
        self.modelProfile.replaceItem(num-1, label, string.join(list,' ; '))


    @pyqtSignature("")
    def slotAddProfile(self):
        """
        Set in view label and variables to see on profile
        """
        var_prof = [str(s) for s in self.modelDrop.stringList()]
        log.debug("slotAddProfile -> %s" % (var_prof,))
        if not var_prof:
            title = self.tr("Warning")
            msg   = self.tr("You must select at least one variable or property from list")
            QMessageBox.information(self, title, msg)

        else:
            label = self.__verifLabel()

            self.entriesNumber = self.entriesNumber + 1
            if label == '':
                label = 'profile' + repr(self.entriesNumber)

            if label in self.mdl.getProfilesLabelsList():
                title = self.tr("Warning")
                msg = self.tr("This label already exists")
                QMessageBox.information(self, title, msg)
                return

            log.debug("slotAddProfile -> %s" % (label,))
            self.__insertProfile(label, var_prof)

            choice, freq = self.slotFrequencyType(self.comboBoxFreq.currentText())
            format = self.slotFormatType(self.comboBoxFormat.currentText())
            title = str(self.lineEditTitle.text())
            if not title: title = label
            X1, ok = self.lineEditX1.text().toDouble()
            Y1, ok = self.lineEditY1.text().toDouble()
            Z1, ok = self.lineEditZ1.text().toDouble()
            X2, ok = self.lineEditX2.text().toDouble()
            Y2, ok = self.lineEditY2.text().toDouble()
            Z2, ok = self.lineEditZ2.text().toDouble()

            self.mdl.setProfile(label, title, format, var_prof, freq, X1, Y1, Z1, X2, Y2, Z2)
            self.__eraseEntries()


    @pyqtSignature("")
    def slotDeleteProfile(self):
        """
        Delete the profile from the list (one by one).
        """
        row = self.treeViewProfile.currentIndex().row()
        log.debug("slotDeleteProfile -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing profile")
            QMessageBox.information(self, title, msg)
        else:
            label, title, format, list, freq, x1, y1, z1, x2, y2, z2 = self.__infoProfile(row)
            self.modelProfile.deleteRow(row)
            self.mdl.deleteProfile(label)
            self.__eraseEntries()


    @pyqtSignature("")
    def slotEditProfile(self):
        """
        Edit profile to modify its characteristics.
        """
        row = self.treeViewProfile.currentIndex().row()
        log.debug("slotEditProfile -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing profile")
            QMessageBox.information(self, title, msg)
        else:
            old_label, title, format, vlist, freq, x1, y1, z1, x2, y2, z2 = self.__infoProfile(row)

            var_prof = [str(s) for s in self.modelDrop.stringList()]
            log.debug("slotEditProfile -> %s" % (var_prof,))
            if not var_prof:
                title = self.tr("Warning")
                msg   = self.tr("You must select at least one variable or property from list")
                QMessageBox.information(self, title, msg)

            else:
                new_label = str(self.lineEditBaseName.text())

                if new_label == '':
                    new_label = old_label

                if new_label != old_label and new_label in self.mdl.getProfilesLabelsList():
                    title = self.tr("Warning")
                    msg = self.tr("This label already exists: the old label is kept")
                    QMessageBox.information(self, title, msg)
                    new_label = old_label

                log.debug("slotEditedProfile -> %s" % (new_label,))
                self.__replaceProfile(row+1, new_label, var_prof)
                idx = self.treeViewProfile.currentIndex()
                self.treeViewProfile.dataChanged(idx, idx)

                choice, freq = self.slotFrequencyType(self.comboBoxFreq.currentText())
                format = self.slotFormatType(self.comboBoxFormat.currentText())
                title = str(self.lineEditTitle.text())
                if not title: title = new_label
                X1, ok = self.lineEditX1.text().toDouble()
                Y1, ok = self.lineEditY1.text().toDouble()
                Z1, ok = self.lineEditZ1.text().toDouble()
                X2, ok = self.lineEditX2.text().toDouble()
                Y2, ok = self.lineEditY2.text().toDouble()
                Z2, ok = self.lineEditZ2.text().toDouble()

                self.mdl.replaceProfile(old_label, new_label, title, format, var_prof, freq, X1, Y1, Z1, X2, Y2, Z2)
                self.__eraseEntries()


    @pyqtSignature("const QModelIndex &")
    def slotSelectProfile(self, index):
        """
        Return the selected item from the list.
        """
        row = index.row()
        log.debug("slotSelectProfile -> %s" % (row,))

        label, title, format, liste, freq, x1, y1, z1, x2, y2, z2 = self.__infoProfile(row)

        self.lineEditTitle.setText(QString(str(title)))
        self.lineEditBaseName.setText(QString(str(label)))
        self.modelFormat.setItem(str_model=format)

        if freq >= 1:
            self.modelFreq.setItem(str_model='frequency')
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(QString(str(freq)))
            self.labelBaseName.setText(QString("Filename"))
        else:
            self.modelFreq.setItem(str_model='end')
            self.lineEditFreq.setText(QString(str("-1")))
            self.lineEditFreq.setDisabled(True)
            self.labelBaseName.setText(QString("Filename"))

        self.lineEditX1.setText(QString(str(x1)))
        self.lineEditY1.setText(QString(str(y1)))
        self.lineEditZ1.setText(QString(str(z1)))
        self.lineEditX2.setText(QString(str(x2)))
        self.lineEditY2.setText(QString(str(y2)))
        self.lineEditZ2.setText(QString(str(z2)))

        self.modelDrop.setStringList(QStringList())
        liste = [QString(s) for s in liste]
        self.modelDrop.setStringList(liste)


    @pyqtSignature("")
    def slotAddVarProfile(self):
        """
        Add a new var from list to profile
        """
        if (self.DragList.currentIndex().row() >=0) :
            liste = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.DragList.currentIndex().row()]
            if var not in liste :
                liste.append(var)
            self.modelDrop.setStringList(liste)


    @pyqtSignature("")
    def slotDeleteVarProfile(self):
        """
        Supress a var from profile
        """
        self.modelDrop.removeRows(self.DropList.currentIndex().row(), 1)


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.lineEditTitle.setText(QString(str("")))
        self.lineEditBaseName.setText(QString(str("")))
        self.lineEditX1.setText(QString(str("")))
        self.lineEditY1.setText(QString(str("")))
        self.lineEditZ1.setText(QString(str("")))
        self.lineEditX2.setText(QString(str("")))
        self.lineEditY2.setText(QString(str("")))
        self.lineEditZ2.setText(QString(str("")))

        self.modelFreq.setItem(str_model='end')
        self.lineEditFreq.setText(QString(str("-1")))
        self.lineEditFreq.setDisabled(True)
        self.labelBaseName.setText(QString("Filename"))
        self.modelDrop.setStringList(QStringList())
        self.treeViewProfile.clearSelection()


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
