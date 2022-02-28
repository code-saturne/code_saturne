# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module defines the 'Start/Restart' page.

This module contains the following classes:
- StartRestartAdvancedDialogView
- StartRestartView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types
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

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, IntValidator, from_qvariant
from code_saturne.model.SolutionDomainModel import RelOrAbsPath
from code_saturne.gui.case.StartRestartForm import Ui_StartRestartForm
from code_saturne.gui.case.StartRestartAdvancedDialogForm import Ui_StartRestartAdvancedDialogForm
from code_saturne.model.StartRestartModel import StartRestartModel, getRestartInfo

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
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_StartRestartAdvancedDialogForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result = self.default.copy()

        # Combo models and items
        self.modelFreq   = ComboModel(self.comboBoxFreq, 4, 1)

        self.modelFreq.addItem(self.tr("Never"), 'Never')
        self.modelFreq.addItem(self.tr("Only at the end of the calculation"),
                               'At the end')
        self.modelFreq.addItem(self.tr("4 restart checkpoints"), '4 output')
        self.modelFreq.addItem(self.tr("Checkpoints frequency :"), 'Frequency')

        # Connections

        self.comboBoxFreq.activated[str].connect(self.slotFreq)
        self.lineEditNSUIT.textChanged[str].connect(self.slotNsuit)

        # Validator

        validatorNSUIT = IntValidator(self.lineEditNSUIT, min=0)
        self.lineEditNSUIT.setValidator(validatorNSUIT)

        # Read of auxiliary file if calculation restart is asked

        if self.default['restart']:
            self.groupBoxRestart.show()

            if self.default['restart_with_auxiliary'] == 'on':
                self.checkBoxReadAuxFile.setChecked(True)
            else:
                self.checkBoxReadAuxFile.setChecked(False)
        else:
            self.groupBoxRestart.hide()

        # Frequency of rescue of restart file

        if self.default['restart_rescue'] == -2:
            self.nsuit = -2
            self.lineEditNSUIT.setDisabled(True)
            self.freq = 'Never'
        elif self.default['restart_rescue'] == -1:
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

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotFreq(self, text):
        """
        Creation of popup window's widgets
        """
        self.freq = self.modelFreq.dicoV2M[str(text)]
        log.debug("getFreq-> %s" % self.freq)

        if self.freq == "Never":
            self.nsuit = -2
            self.lineEditNSUIT.setText(str(self.nsuit))
            self.lineEditNSUIT.setDisabled(True)

        elif self.freq == "At the end":
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


    @pyqtSlot(str)
    def slotNsuit(self, text):
        if self.lineEditNSUIT.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, int)
            self.nsuit = n
            log.debug("getNsuit-> nsuit = %s" % n)


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """
        if self.default['restart']:
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
        self.case.undoStopGlobal()

        self.radioButtonYes.clicked.connect(self.slotStartRestart)
        self.radioButtonNo.clicked.connect(self.slotStartRestart)
        self.radioButtonAuto.clicked.connect(self.slotStartRestart)
        self.toolButton.pressed.connect(self.slotSearchRestartDirectory)
        self.toolButtonRestartMesh.pressed.connect(self.slotSearchRestartMesh)
        self.checkBox.clicked.connect(self.slotFrozenField)
        self.toolButtonAdvanced.pressed.connect(self.slotAdvancedOptions)
        self.checkBoxRestartMesh.stateChanged.connect(self.slotRestartMesh)

        self.model = StartRestartModel(self.case)

        # Widget initialization

        self.restart_path = self.model.getRestartPath()
        self.restart_mesh_path = self.model.getRestartMeshPath()

        if self.restart_path:
            if self.restart_path == '*':
                self.radioButtonNo.setChecked(False)
                self.radioButtonYes.setChecked(False)
                self.radioButtonAuto.setChecked(True)
            else:
                if not os.path.isdir(os.path.join(self.case['case_path'],
                                                  self.restart_path)):
                    title = self.tr("WARNING")
                    msg   = self.tr("Invalid path in %s!" % self.restart_path)
                    QMessageBox.warning(self, title, msg)

                self.radioButtonNo.setChecked(False)
                self.radioButtonYes.setChecked(True)
                self.radioButtonAuto.setChecked(False)

        else:
            self.radioButtonNo.setChecked(True)
            self.radioButtonYes.setChecked(False)
            self.radioButtonAuto.setChecked(False)

        self.slotStartRestart()

        if self.model.getFrozenField() == 'on':
            self.checkBox.setChecked(True)
        else:
            self.checkBox.setChecked(False)

        self.updateRestartTimes()
        self.updateRestartMeshView()

        self.case.undoStartGlobal()


    def updateRestartTimes(self):
        """
        Update information on restart times
        """

        # FIXME: ensure correct path is used in coupling situations;
        # currently, if in doubt, leave it empty (we prefer to have
        # no information than false information)
        restart_dir = None
        restart_path = None

        if self.restart_path == '*':
            d = os.path.join(os.path.split(self.case['case_path'])[0],
                             'RESU_COUPLING')
            if not os.path.isdir(d):
                restart_dir = os.path.join(self.case['case_path'], 'RESU')

            # isdir returns an error if restart_dir is None, hence the test
            if restart_dir:
                if os.path.isdir(restart_dir):
                    restart_path = '*'
                else:
                    restart_dir = None
        elif self.restart_path:
            if os.path.isabs(self.restart_path):
                restart_path = self.restart_path
            else:
                restart_path = os.path.join(self.case['case_path'],
                                            self.restart_path)

        rinfo = getRestartInfo(self.case['package'],
                               restart_dir,
                               restart_path)

        self.lineEdit.setEnabled(self.restart_path != '*')
        self.lineEdit.setFrame(self.restart_path != '*')

        if rinfo:
            if self.restart_path == '*':
                self.lineEdit.setText(rinfo[0])
            else:
                self.lineEdit.setText(self.restart_path)
            self.labelIteration.show()
            self.labelTime.show()
            self.lineEditIteration.show()
            self.lineEditTime.show()
            self.lineEditIteration.setText(str(rinfo[1]))
            self.lineEditTime.setText(str(rinfo[2]))
        else:
            self.labelIteration.hide()
            self.labelTime.hide()
            self.lineEditIteration.hide()
            self.lineEditTime.hide()


    def updateRestartMeshView(self):
        """
        Upate restart mesh path view
        """
        if self.restart_mesh_path:
            self.checkBoxRestartMesh.setChecked(True)
            self.lineEditRestartMesh.setText(self.restart_mesh_path)
            self.lineEditRestartMesh.show()
            self.toolButtonRestartMesh.show()
        else:
            self.checkBoxRestartMesh.setChecked(False)
            self.lineEditRestartMesh.setText("")
            self.lineEditRestartMesh.hide()
            self.toolButtonRestartMesh.hide()


    @pyqtSlot()
    def slotSearchRestartDirectory(self):
        """
        Search restart file (directory) in list of directories
        """

        default = None
        l_restart_dirs = []
        for d in [os.path.join(os.path.split(self.case['case_path'])[0],
                               'RESU_COUPLING'),
                  os.path.join(self.case['case_path'], 'RESU')]:
            if os.path.isdir(d):
                l_restart_dirs.append(QUrl.fromLocalFile(d))
                if not default:
                    default = d

        if not default:
            default = self.case['case_path']

        title = self.tr("Select checkpoint/restart directory")
        options = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly | QFileDialog.ShowDirsOnly

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        dialog.setOptions(options)
        dialog.setSidebarUrls(l_restart_dirs)
        dialog.setFileMode(QFileDialog.Directory)

        name_filter = str(self.tr("Checkpoint directory (checkpoint*)"))
        dialog.setNameFilter(name_filter)

        dialog.setLabelText(QFileDialog.Accept, str(self.tr("Select")))

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()

            dir_path = str(s[0])
            dir_path = os.path.abspath(dir_path)

            self.restart_path = RelOrAbsPath(dir_path, self.case['case_path'])
            self.model.setRestartPath(self.restart_path)
            self.lineEdit.setText(self.restart_path)
            self.updateRestartTimes()

            log.debug("slotSearchRestartDirectory-> %s" % self.restart_path)


    @pyqtSlot()
    def slotSearchRestartMesh(self):
        """
        Search restart mesh (file) in list of directories
        """
        title    = self.tr("Select checkpoint/restart mesh_input/output")

        default = None
        l_restart_dirs = []
        for d in [os.path.join(os.path.split(self.case['case_path'])[0],
                               'RESU_COUPLING'),
                  os.path.join(self.case['case_path'], 'RESU')]:
            if os.path.isdir(d):
                l_restart_dirs.append(QUrl.fromLocalFile(d))
                if not default:
                    default = d

        if not default:
            default = self.case['case_path']

        options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        dialog.setOptions(options)
        dialog.setSidebarUrls(l_restart_dirs)
        dialog.setFileMode(QFileDialog.ExistingFile)

        name_filter = str(self.tr("Imported or preprocessed mesh (mesh_input mesh_output *.csm)"))
        dialog.setNameFilter(name_filter)

        dialog.setLabelText(QFileDialog.Accept, str(self.tr("Select")))

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()

            path = str(s[0])
            path = os.path.abspath(path)

            self.restart_mesh_path = RelOrAbsPath(path, self.case['case_path'])
            self.model.setRestartMeshPath(self.restart_mesh_path)
            self.lineEditRestartMesh.setText(self.restart_mesh_path)

            log.debug("slotSearchRestartDirectory-> %s" % self.restart_mesh_path)


    @pyqtSlot()
    def slotStartRestart(self):
        """
        Handle restart.
        """
        if self.radioButtonYes.isChecked():
            if not self.restart_path or self.restart_path == '*':
                self.slotSearchRestartDirectory()
        elif self.radioButtonAuto.isChecked():
            self.restart_path = '*'
        else:
            self.restart_path = None

        if self.restart_path:
            self.model.setRestartPath(self.restart_path)
            if self.restart_path == '*':
                self.radioButtonYes.setChecked(False)
                self.radioButtonAuto.setChecked(True)
                self.labelRestartDir.setEnabled(False)
                self.toolButton.hide()
                self.updateRestartTimes()
            else:
                self.radioButtonYes.setChecked(True)
                self.radioButtonAuto.setChecked(False)
                self.labelRestartDir.setEnabled(True)
                self.toolButton.show()
            self.radioButtonNo.setChecked(False)
            self.frameRestart.show()
        else:
            self.model.setRestartPath(None)
            self.model.setRestartMeshPath(None)
            self.model.setFrozenField("off")
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)
            self.checkBoxRestartMesh.setChecked(False)
            self.frameRestart.hide()
            self.lineEdit.setText("")
            self.updateRestartTimes()

        self.updateRestartMeshView()


    @pyqtSlot()
    def slotRestartMesh(self):
        """
        Input different restart mesh.
        """
        if self.checkBoxRestartMesh.isChecked():
            if not self.restart_mesh_path:
                self.slotSearchRestartMesh()

        else:
            self.restart_mesh_path = None

        self.model.setRestartMeshPath(self.restart_mesh_path)

        self.updateRestartMeshView()


    @pyqtSlot()
    def slotFrozenField(self):
        """
        Input if calculation on frozen velocity and pressure fields or not
        """
        if self.checkBox.isChecked():
            self.model.setFrozenField('on')
        else:
            self.model.setFrozenField('off')


    @pyqtSlot()
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        freq, period = self.model.getRestartRescue()

        default                           = {}
        default['restart']                = self.model.getRestartPath()
        default['restart_with_auxiliary'] = self.model.getRestartWithAuxiliaryStatus()
        default['restart_rescue']         = freq
        default['period_rescue']          = period
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = StartRestartAdvancedDialogView(self, self.case, default)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.model.setRestartWithAuxiliaryStatus(result['restart_with_auxiliary'])
            self.model.setRestartRescue(result['restart_rescue'])


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
