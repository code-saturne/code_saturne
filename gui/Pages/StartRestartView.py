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
This module defines the 'Start/Restart' page.

This module contains the following classes:
- StartRestartAdvancedDialogView
- StartRestartView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, IntValidator, setGreenColor
from StartRestartForm import Ui_StartRestartForm
from StartRestartAdvancedDialogForm import Ui_StartRestartAdvancedDialogForm
from Pages.StartRestartModel import StartRestartModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("StartRestartView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Popup window class
#-------------------------------------------------------------------------------

class StartRestartAdvancedDialogView(QDialog, Ui_StartRestartAdvancedDialogForm):
    """
    Building of popup window for advanced options.
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_StartRestartAdvancedDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result = self.default.copy()

        # Combo models and items
        self.modelFreq   = ComboModel(self.comboBoxFreq, 3, 1)

        self.modelFreq.addItem(self.tr("Only at the end of the calculation"), 'At the end')
        self.modelFreq.addItem(self.tr("4 restart checkpoints"), '4 output')
        self.modelFreq.addItem(self.tr("Checkpoints frequency :"), 'Frequency')

        # Connections

        self.connect(self.comboBoxFreq, SIGNAL("activated(const QString&)"), self.slotFreq)
        self.connect(self.lineEditNSUIT, SIGNAL("textChanged(const QString&)"), self.slotNsuit)

        # Validator

        validatorNSUIT = IntValidator(self.lineEditNSUIT, min=0)
        self.lineEditNSUIT.setValidator(validatorNSUIT)

        # Read of auxiliary file if calculation restart is asked

        if self.default['restart'] == "on":
            self.groupBoxRestart.show()

            if self.default['restart_with_auxiliary'] == 'on':
                self.checkBoxReadAuxFile.setChecked(True)
            else:
                self.checkBoxReadAuxFile.setChecked(False)
        else:
            self.groupBoxRestart.hide()

        # Frequency of rescue of restart file

        if self.default['restart_rescue'] == -1:
            self.nsuit = -1
            self.lineEditNSUIT.setDisabled(True)
            self.freq = 'At the end'
        elif self.default['restart_rescue'] == 0:
            self.nsuit = 0
            self.lineEditNSUIT.setDisabled(True)
            self.freq = '4 output'
        else:
            self.nsuit = self.default['restart_rescue']
            self.lineEditNSUIT.setEnabled(True)
            self.freq = 'Frequency'
        self.modelFreq.setItem(str_model=self.freq)
        self.lineEditNSUIT.setText(str(self.nsuit))


    @pyqtSignature("const QString &")
    def slotFreq(self, text):
        """
        Creation of popup window's widgets
        """
        self.freq = self.modelFreq.dicoV2M[str(text)]
        log.debug("getFreq-> %s" % self.freq)

        if self.freq == "At the end":
            self.nsuit = -1
            self.lineEditNSUIT.setText(str(self.nsuit))
            self.lineEditNSUIT.setDisabled(True)

        elif self.freq == "4 output":
            self.nsuit = 0
            self.lineEditNSUIT.setText(str(self.nsuit))
            self.lineEditNSUIT.setDisabled(True)

        elif self.freq == "Frequency":
            if self.nsuit <= 0: self.nsuit = 1
            self.lineEditNSUIT.setText(str(self.nsuit))
            self.lineEditNSUIT.setEnabled(True)


    @pyqtSignature("const QString &")
    def slotNsuit(self, text):
        n, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.nsuit = n
            log.debug("getNsuit-> nsuit = %s" % n)


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """
        if self.default['restart'] == 'on':
            if self.checkBoxReadAuxFile.isChecked():
                self.result['restart_with_auxiliary'] = 'on'
            else:
                self.result['restart_with_auxiliary'] = 'off'

        self.result['restart_rescue'] = self.nsuit
        self.result['period_rescue']  = self.freq

        QDialog.accept(self)


    def reject(self):
        """
        Method called when 'Cancel' button is clicked
        """
        QDialog.reject(self)


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class StartRestartView(QWidget, Ui_StartRestartForm):
    """
    This page is devoted to the start/restart control.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_StartRestartForm.__init__(self)
        self.setupUi(self)
        self.case = case

        self.connect(self.radioButtonYes, SIGNAL("clicked()"), self.slotStartRestart)
        self.connect(self.radioButtonNo, SIGNAL("clicked()"), self.slotStartRestart)
        self.connect(self.toolButton, SIGNAL("pressed()"), self.slotSearchRestartDirectory)
        self.connect(self.checkBox, SIGNAL("clicked()"), self.slotFrozenField)
        self.connect(self.toolButtonAdvanced, SIGNAL("pressed()"), self.slotAdvancedOptions)

        self.model = StartRestartModel(self.case)

        # Widget initialization

        if 'RESTART' in os.listdir(self.case['data_path']):
            if not os.path.isdir(os.readlink(self.case['data_path'] + '/RESTART')):
                title = self.tr("WARNING")
                msg   = self.tr("Invalid link in %s!" % self.case['data_path'])
                QMessageBox.warning(self, title, msg)

        if self.model.getRestart() == "on":
            self.radioButtonYes.setChecked(True)
            self.radioButtonNo.setChecked(False)
        else:
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)
        self.slotStartRestart()

        if self.model.getFrozenField() == 'on':
            self.checkBox.setChecked(True)
        else:
            self.checkBox.setChecked(False)


    @pyqtSignature("")
    def slotSearchRestartDirectory(self):
        """
        Search restart file (directory) in list of directories
        """
        title    = self.tr("Selection of the restart directory RESTART")
        default  = self.case['resu_path']
        dir_path = QFileDialog.getExistingDirectory(self, title, default, QFileDialog.ShowDirsOnly)

        if dir_path:
            dir_path = os.path.abspath(str(dir_path))
            dir_name = os.path.basename(dir_path)
            log.debug("slotSearchRestartDirectory-> %s" % dir_name)
            link = self.case['data_path'] + '/RESTART'

            # a link already exists
            if os.path.islink(link):
                exist_name = os.path.abspath(os.readlink(link))
                if dir_path != exist_name:
                    title = self.tr("WARNING")
                    msg   = self.tr("This symbolic link already exists in DATA:\n\n" \
                                    + exist_name + "\n\nReplace with the new directory?")
                    ans = QMessageBox.question(self, title, msg, QMessageBox.Yes, QMessageBox.No)
                    if ans == QMessageBox.Yes:
                        os.unlink(link)
                        os.symlink(dir_path, link)
                        self.model.setRestartDirectory(dir_name)
                    else:
                        self.model.setRestartDirectory(os.path.basename(exist_name))
                else:
                    self.model.setRestartDirectory(os.path.basename(exist_name))

            # a RESTART directory exists
            elif os.path.isdir(link):
                title = self.tr("WARNING")
                msg   = self.tr("The directory RESTART is already in DATA.\n\n" \
                                "Replace with the new directory?")
                ans = QMessageBox.question(self, title, msg, QMessageBox.Yes, QMessageBox.No)
                if ans == QMessageBox.Yes:
                    os.rmdir(link)
                    os.symlink(dir_path, link)
                    self.model.setRestartDirectory(dir_name)
                else:
                    self.model.setRestartDirectory('RESTART')

            # create a new link
            else:
                os.symlink(os.path.abspath(dir_path), link)
                self.model.setRestartDirectory(dir_name)

            setGreenColor(self.toolButton, False)
            self.slotStartRestart()


    @pyqtSignature("")
    def slotStartRestart(self):
        """
        Input IRESTART Code_Saturne keyword.
        """
        if self.radioButtonYes.isChecked():
            self.model.setRestart("on")
            self.toolButton.setEnabled(True)
            self.radioButtonYes.setChecked(True)
            self.radioButtonNo.setChecked(False)

            name = self.model.getRestartDirectory()
            if name:
                self.labelDir1.show()
                self.labelDir2.show()
                self.labelDir2.setText(name)
            else:
                setGreenColor(self.toolButton)

        else:
            self.model.setRestart("off")
            self.toolButton.setDisabled(True)
            self.model.setFrozenField("off")
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)

            self.labelDir1.hide()
            self.labelDir2.hide()
            setGreenColor(self.toolButton, False)


    @pyqtSignature("")
    def slotFrozenField(self):
        """
        Input if calculation on frozen velocity and pressure fields or not
        """
        if self.checkBox.isChecked():
            self.model.setFrozenField('on')
        else:
            self.model.setFrozenField('off')


    @pyqtSignature("")
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        freq, period = self.model.getRestartRescue()

        default                           = {}
        default['restart']                = self.model.getRestart()
        default['restart_with_auxiliary'] = self.model.getRestartWithAuxiliaryStatus()
        default['restart_rescue']         = freq
        default['period_rescue']          = period
##        default['main_restart']           = self.model.getMainRestartFormat()
##        default['auxiliary_restart']      = self.model.getAuxiliaryRestartFormat()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = StartRestartAdvancedDialogView(self, default)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.model.setRestart(result['restart'])
            self.model.setRestartWithAuxiliaryStatus(result['restart_with_auxiliary'])
            self.model.setRestartRescue(result['restart_rescue'])
##            self.model.setMainRestartFormat(result['main_restart'])
##            self.model.setAuxiliaryRestartFormat(result['auxiliary_restart'])


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
