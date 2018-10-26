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
This module contains the following classes:
- ProfilesView
"""

#-------------------------------------------------------------------------------
# Standard modules
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
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant
from code_saturne.Pages.ProfilesForm import Ui_ProfilesForm
from code_saturne.Pages.ProfilesView import StandardItemModelProfile
from code_saturne.Pages.ProfilesModelNeptune import ProfilesModel
from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ProfilesView")
log.setLevel(GuiParam.DEBUG)

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
        self.case.undoStopGlobal()
        self.mdl = ProfilesModel(self.case)

        #  Initialize variables concerning the display of the Hlist

        self.entriesNumber = 0

        # Models
        self.modelProfile = StandardItemModelProfile()
        self.treeViewProfile.setModel(self.modelProfile)
        self.treeViewProfile.resizeColumnToContents(0)

        # QListView layout
        self.gridlayout1 = QGridLayout(self.widgetDrag)
        self.gridlayout1.setContentsMargins(0, 0, 0, 0)
        self.DragList = QListView(self.widgetDrag)
        self.gridlayout1.addWidget(self.DragList,0,0,1,1)

        self.gridlayout2 = QGridLayout(self.widgetDrop)
        self.gridlayout2.setContentsMargins(0, 0, 0, 0)
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
        self.modelFreq = ComboModel(self.comboBoxFreq, 3, 1)
        self.modelFreq.addItem(self.tr("at the end of the calculation"), "end")
        self.modelFreq.addItem(self.tr("at each 'n' time steps"), "frequency")
        self.modelFreq.addItem(self.tr("Output every 'x' seconds"), 'time_value')

        self.modelFormat = ComboModel(self.comboBoxFormat, 2, 1)
        self.modelFormat.addItem(self.tr(".dat"), "DAT")
        self.modelFormat.addItem(self.tr(".csv"), "CSV")

        # Connections
        self.treeViewProfile.pressed[QModelIndex].connect(self.slotSelectProfile)
        self.pushButtonAdd.clicked.connect(self.slotAddProfile)
        self.pushButtonDelete.clicked.connect(self.slotDeleteProfile)
        self.pushButtonAddVar.clicked.connect(self.slotAddVarProfile)
        self.pushButtonSuppressVar.clicked.connect(self.slotDeleteVarProfile)
        self.comboBoxFreq.activated[str].connect(self.slotFrequencyType)
        self.comboBoxFormat.activated[str].connect(self.slotFormatType)
        self.pushButtonFormula.clicked.connect(self.slotFormula)
        self.lineEditBaseName.textChanged[str].connect(self.slotBaseName)
        self.lineEditFreq.textChanged[str].connect(self.slotFrequence)
        self.lineEditFreqTime.textChanged[str].connect(self.slotFrequenceTime)
        self.lineEditNbPoint.textChanged[str].connect(self.slotNbPoint)

        # Validators
        validatorFreq = IntValidator(self.lineEditFreq, min=0)
        validatorFreq.setExclusiveMin(True)
        self.lineEditFreq.setValidator(validatorFreq)

        validatorFreqT = DoubleValidator(self.lineEditFreqTime, min=0.)
        validatorFreqT.setExclusiveMin(True)
        self.lineEditFreqTime.setValidator(validatorFreqT)

        validatorNbPoint = IntValidator(self.lineEditNbPoint, min=0)
        self.lineEditNbPoint.setValidator(validatorNbPoint)

        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorBaseName =  RegExpValidator(self.lineEditBaseName, QRegExp(rx))
        self.lineEditBaseName.setValidator(validatorBaseName)

        #update list of variables, properties, scalars ...
        liste_label = []
        for label in self.mdl.getVariablesAndVolumeProperties():
            liste_label.append(label)
        self.modelDrag.setStringList(sorted(liste_label, key=str.lower))

        #update list of profiles for view from xml file
        for lab in self.mdl.getProfilesLabelsList():
            self.entriesNumber = self.entriesNumber + 1
            label, fmt, lst, choice, freq, formula, nb_point = self.mdl.getProfileData(lab)
            self.__insertProfile(label, lst)

        self.__eraseEntries()

        self.case.undoStartGlobal()


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
            from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
            dialog = VerifyExistenceLabelDialogView(self, default)
            if dialog.exec_():
                result = dialog.get_result()
                label = result['label']
                if result['label'] == default['label']:
                    label = ""
        return label


    @pyqtSlot(str)
    def slotFrequencyType(self, text):
        """
        Input choice for frequency for profile.
        """
        choice = self.modelFreq.dicoV2M[str(text)]

        if choice == "end":
            nfreq = -1
            self.lineEditFreq.setText(str(nfreq))
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setDisabled(True)

        elif choice == "frequency":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            nfreq = self.mdl.getOutputFrequency(self.label_select)
            if nfreq == -1:
                nfreq = 1
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(str(nfreq))

        elif choice == "time_value":
            self.lineEditFreq.hide()
            self.lineEditFreqTime.show()
            nfreq = self.mdl.getOutputFrequency(self.label_select)
            if nfreq == -1:
                nfreq = 1.0
            self.lineEditFreqTime.setText(str(nfreq))

        self.mdl.setOutputType(self.label_select, choice)
        self.mdl.setOutputFrequency(self.label_select, nfreq)


    @pyqtSlot(str)
    def slotFormatType(self, text):
        """
        Input choice for frequency for profile.
        """
        fmt = self.modelFormat.dicoV2M[str(text)]
        self.mdl.setFormat(self.label_select, fmt)


    def __infoProfile(self, row):
        """
        Return info from the argument entry.
        """
        label = self.modelProfile.getLabel(row)
        lab, fmt, lst, choice, freq, formula, nb_point = self.mdl.getProfileData(label)
        return label, fmt, lst, choice, freq, formula, nb_point


    def __insertProfile(self, label, lst):
        """
        Insert values in table view.
        """
        self.modelProfile.addItem(label, " ; ".join(lst))


    @pyqtSlot()
    def slotAddProfile(self):
        """
        Set in view label and variables to see on profile
        """
        var_prof = []
        label = self.mdl.addProfile()
        self.__insertProfile(label, var_prof)
        self.__eraseEntries()


    @pyqtSlot()
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
            label, fmt, lst, choice, freq, formula, nb_point = self.__infoProfile(row)
            self.modelProfile.deleteRow(row)
            self.mdl.deleteProfile(label)
            self.__eraseEntries()


    @pyqtSlot("QModelIndex")
    def slotSelectProfile(self, index):
        """
        Return the selected item from the list.
        """
        self.groupBoxProfile.show()

        row = index.row()
        log.debug("slotSelectProfile -> %s" % (row,))

        label, fmt, liste, choice, freq, formula, nb_point = self.__infoProfile(row)
        self.label_select = label

        self.lineEditBaseName.setText(str(label))
        self.modelFormat.setItem(str_model=fmt)

        self.modelFreq.setItem(str_model=choice)
        if choice == "end":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setText(str("-1"))
            self.lineEditFreq.setDisabled(True)

        elif choice == "frequency":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(str(freq))

        elif choice == "time_value":
            self.lineEditFreq.hide()
            self.lineEditFreqTime.show()
            self.lineEditFreqTime.setText(str(freq))

        self.lineEditNbPoint.setText(str(nb_point))

        self.modelDrop.setStringList([])
        liste = [str(s) for s in liste]

        self.modelDrop.setStringList(liste)

        exp = self.mdl.getFormula(self.label_select)
        if exp:
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(exp)
        else:
            self.pushButtonFormula.setStyleSheet("background-color: red")


    @pyqtSlot()
    def slotAddVarProfile(self):
        """
        Add a new var from list to profile
        """
        if (self.DragList.currentIndex().row() >=0) :
            liste = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.DragList.currentIndex().row()]
            if var not in liste :
                liste.append(var)
            liste = [str(s) for s in liste]
            self.modelDrop.setStringList(liste)
            self.mdl.setVariable(self.label_select, liste)

            row = self.treeViewProfile.currentIndex().row()
            liste = self.mdl.getVariable(self.label_select)
            self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    @pyqtSlot()
    def slotDeleteVarProfile(self):
        """
        Supress a var from profile
        """
        self.modelDrop.removeRows(self.DropList.currentIndex().row(), 1)
        liste = self.modelDrop.stringList()
        liste = [str(s) for s in liste]
        self.mdl.setVariable(self.label_select, liste)

        row = self.treeViewProfile.currentIndex().row()
        liste = self.mdl.getVariable(self.label_select)
        self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.groupBoxProfile.hide()
        self.label_select = None
        self.treeViewProfile.clearSelection()


    @pyqtSlot()
    def slotFormula(self):
        """
        """
        exp = self.mdl.getFormula(self.label_select)
        exa = """#example: a line segment
#(s, the parameter is always between 0 and 1)
x = 2*s + 3.2;
y = 2;
z = -0.5*s+5;"""
        req = [('x', "x formula"),
               ('y', "y formula"),
               ('z', "z formula")]
        sym = [('s', 'parameter')]

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotLineFormula -> %s" % str(result))
            self.mdl.setFormula(self.label_select, result)
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(result)


    @pyqtSlot(str)
    def slotBaseName(self, text):
        """
        """
        lst = self.mdl.getProfilesLabelsList()
        if text not in lst:
            if self.lineEditBaseName.validator().state == QValidator.Acceptable:
                self.mdl.setLabel(self.label_select, str(text))
                self.label_select = str(text)

                row = self.treeViewProfile.currentIndex().row()
                liste = self.mdl.getVariable(self.label_select)
                self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    @pyqtSlot(str)
    def slotFrequence(self, text):
        """
        """
        if self.lineEditFreq.validator().state == QValidator.Acceptable:
            self.mdl.setOutputFrequency(self.label_select, int(text))


    @pyqtSlot(str)
    def slotFrequenceTime(self, text):
        """
        """
        if self.lineEditFreqTime.validator().state == QValidator.Acceptable:
            self.mdl.setOutputFrequency(self.label_select, float(text))


    @pyqtSlot(str)
    def slotNbPoint(self, text):
        """
        """
        if self.lineEditNbPoint.validator().state == QValidator.Acceptable:
            self.mdl.setNbPoint(self.label_select, int(text))


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
