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
This module defines the 'PerformanceTuning' page.

This module contains the following classes:
- PerformanceTuningView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, types, shutil
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

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, RegExpValidator
from code_saturne.Pages.SolutionDomainModel import RelOrAbsPath
from code_saturne.Pages.PerformanceTuningForm import Ui_PerformanceTuningForm
from code_saturne.Pages.PerformanceTuningModel import PerformanceTuningModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PerformanceTuningView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class PerformanceTuningView(QWidget, Ui_PerformanceTuningForm):
    """
    This page is devoted to the performance tuning control.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_PerformanceTuningForm.__init__(self)
        self.setupUi(self)
        self.case = case
        self.case.undoStopGlobal()

        self.mdl = PerformanceTuningModel(self.case)

        # Combo models and items

        self.modelPartType = ComboModel(self.comboBox_PartType, 4, 1)
        self.modelPartOut = ComboModel(self.comboBox_PartOutput, 3, 1)

        self.modelPartType.addItem(self.tr("Default"), 'default')
        self.modelPartType.addItem(self.tr("PT-SCOTCH / SCOTCH"), 'scotch')
        self.modelPartType.addItem(self.tr("ParMETIS / METIS"), 'metis')
        self.modelPartType.addItem(self.tr("Morton curve (bounding box)"), 'morton sfc')
        self.modelPartType.addItem(self.tr("Morton curve (bounding cube)"), 'morton sfc cube')
        self.modelPartType.addItem(self.tr("Hilbert curve (bounding box)"), 'hilbert sfc')
        self.modelPartType.addItem(self.tr("Hilbert curve (bounding cube)"), 'hilbert sfc cube')
        self.modelPartType.addItem(self.tr("Block (unoptimized)"), 'block')

        import cs_config
        cfg = cs_config.config()
        if cfg.libs['scotch'].have == "no":
            self.comboBox_PartType.setItemData(1, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['metis'].have == "no":
            self.comboBox_PartType.setItemData(2, QColor(Qt.red), Qt.TextColorRole);


        self.modelPartOut.addItem(self.tr("No"), 'no')
        self.modelPartOut.addItem(self.tr("For graph-based partitioning"), 'default')
        self.modelPartOut.addItem(self.tr("Yes"), 'yes')

        self.modelBlockIORead = ComboModel(self.comboBox_IORead, 6, 1)
        self.modelBlockIOWrite = ComboModel(self.comboBox_IOWrite, 4, 1)

        self.modelBlockIORead.addItem(self.tr("Default"), 'default')
        self.modelBlockIORead.addItem(self.tr("Standard I/O, serial"), 'stdio serial')
        self.modelBlockIORead.addItem(self.tr("Standard I/O, parallel"), 'stdio parallel')
        self.modelBlockIORead.addItem(self.tr("MPI I/O, independent"), 'mpi independent')
        self.modelBlockIORead.addItem(self.tr("MPI I/O, non-collective"), 'mpi noncollective')
        self.modelBlockIORead.addItem(self.tr("MPI I/O, collective"), 'mpi collective')

        self.modelBlockIOWrite.addItem(self.tr("Default"), 'default')
        self.modelBlockIOWrite.addItem(self.tr("Standard I/O, serial"), 'stdio serial')
        self.modelBlockIOWrite.addItem(self.tr("MPI I/O, non-collective"), 'mpi noncollective')
        self.modelBlockIOWrite.addItem(self.tr("MPI I/O, collective"), 'mpi collective')

        # Validators

        partListVd = RegExpValidator(self.lineEdit_PartList, QRegExp("[0-9- ]*"))
        self.lineEdit_PartList.setValidator(partListVd)

        # Connections

        self.radioButtonYes.clicked.connect(self.slotPartition)
        self.radioButtonNo.clicked.connect(self.slotPartition)
        self.toolButton_PartInputDir.pressed.connect(self.slotSearchPartInputDirectory)
        self.comboBox_PartOutput.activated[str].connect(self.slotPartOut)

        self.comboBox_PartType.activated[str].connect(self.slotPartType)
        self.lineEdit_PartList.textChanged[str].connect(self.slotPartitionList)
        self.spinBoxRankStep.valueChanged[int].connect(self.slotRankStep)

        self.checkBox_IgnorePerio.clicked[bool].connect(self.slotIgnorePerio)

        self.comboBox_IORead.activated[str].connect(self.slotBlockIOReadMethod)
        self.comboBox_IOWrite.activated[str].connect(self.slotBlockIOWriteMethod)

        self.spinBoxIORankStep.valueChanged[int].connect(self.slotBlockIORankStep)
        self.spinBoxIOMinBlockSize.valueChanged[int].connect(self.slotBlockIOMinSize)

        self.tabWidget.currentChanged[int].connect(self.slotchanged)

        # Widget initialization

        self.partinput_path = self.mdl.getPartitionInputPath()

        if self.partinput_path:
            if not os.path.isdir(os.path.join(self.case['case_path'],
                                              self.partinput_path)):
                title = self.tr("WARNING")
                msg   = self.tr("Invalid path in %s!" % self.partinput_path)
                QMessageBox.warning(self, title, msg)

            self.radioButtonYes.setChecked(True)
            self.radioButtonNo.setChecked(False)

        else:
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)

        self.partition_alg = str(self.mdl.getPartitionType())
        self.modelPartType.setItem(str_model=self.partition_alg)

        self.partition_out = str(self.mdl.getPartitionOut())
        self.modelPartOut.setItem(str_model=self.partition_out)

        self.partition_list = str(self.mdl.getPartitionList())
        self.lineEdit_PartList.setText(self.partition_list)

        self.rank_step = self.mdl.getPartitionRankStep()
        self.spinBoxRankStep.setValue(int(self.rank_step))

        self.slotPartition()

        if self.mdl.getIgnorePerio() == 'on':
            self.checkBox_IgnorePerio.setChecked(True)
            self.slotIgnorePerio(True)
        else:
            self.checkBox_IgnorePerio.setChecked(False)
            self.slotIgnorePerio(False)

        self.blockio_read_method = str(self.mdl.getBlockIOReadMethod())
        self.modelBlockIORead.setItem(str_model=self.blockio_read_method)

        self.blockio_write_method = str(self.mdl.getBlockIOWriteMethod())
        self.modelBlockIOWrite.setItem(str_model=self.blockio_write_method)

        self.blockio_rank_step = self.mdl.getBlockIORankStep()
        self.spinBoxIORankStep.setValue(int(self.blockio_rank_step))

        self.blockio_min_size = self.mdl.getBlockIOMinSize()
        self.spinBoxIOMinBlockSize.setValue(int(self.blockio_min_size))

        self.tabWidget.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotSearchPartInputDirectory(self):
        """
        Search for the partition input directory in list of directories
        """
        title    = self.tr("Select partition input directory")

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

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        dialog = QFileDialog(self, title, default)
        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_restart_dirs)
        dialog.setFileMode(QFileDialog.Directory)

        if dialog.exec_() == 1:

            s = dialog.selectedFiles()

            dir_path = str(s[0])
            dir_path = os.path.abspath(dir_path)

            self.partinput_path = RelOrAbsPath(dir_path, self.case['case_path'])
            self.mdl.setPartitionInputPath(self.partinput_path)
            self.lineEditPartInputDir.setText(self.partinput_path)

            log.debug("slotSearchPartInputDirectory-> %s" % self.partinput_path)


    @pyqtSlot()
    def slotPartition(self):
        """
        Determine if existing partitioning is used.
        """
        if self.radioButtonYes.isChecked():
            if not self.partinput_path:
                self.slotSearchPartInputDirectory()

        else:
            self.partinput_path = None

        if self.partinput_path:
            self.mdl.setPartitionInputPath(self.partinput_path)
            self.radioButtonYes.setChecked(True)
            self.radioButtonNo.setChecked(False)
            self.framePartInputDir.show()
            self.lineEditPartInputDir.setText(self.partinput_path)
            if not os.path.isdir(os.path.join(self.case['resu_path'],
                                              self.partinput_path)):
                self.toolButton_PartInputDir.setStyleSheet("background-color: red")
            else:
                self.toolButton_PartInputDir.setStyleSheet("background-color: green")

        else:
            self.mdl.setPartitionInputPath(None)
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)

            self.framePartInputDir.hide()
            self.lineEditPartInputDir.setText("")


    @pyqtSlot(str)
    def slotPartitionList(self, text):
        """
        Input for Partitioner.
        """
        self.partition_list = str(text)
        self.mdl.setPartitionList(self.partition_list.strip())


    @pyqtSlot(str)
    def slotPartOut(self, text):
        """
        Partitioner execution mode option.
        """
        self.partition_out = self.modelPartOut.dicoV2M[str(text)]
        self.mdl.setPartitionOut(self.partition_out)


    @pyqtSlot(str)
    def slotPartType(self, text):
        """
        Partitioner execution mode option.
        """
        self.partition_alg = self.modelPartType.dicoV2M[str(text)]
        self.mdl.setPartitionType(self.partition_alg)


    @pyqtSlot(int)
    def slotRankStep(self, text):
        """
        Input for Partitioner.
        """
        self.rank_step = self.spinBoxRankStep.value()
        self.mdl.setPartitionRankStep(self.rank_step)


    @pyqtSlot(bool)
    def slotIgnorePerio(self, checked):
        """
        Ignore periodicity.
        """
        if checked:
            self.mdl.setIgnorePerio("on")
        else:
            self.mdl.setIgnorePerio("off")


    @pyqtSlot(str)
    def slotBlockIOReadMethod(self, text):
        """
        Partitioner execution mode option.
        """
        self.blockio_read_method = self.modelBlockIORead.dicoV2M[str(text)]
        self.mdl.setBlockIOReadMethod(self.blockio_read_method)


    @pyqtSlot(str)
    def slotBlockIOWriteMethod(self, text):
        """
        Partitioner execution mode option.
        """
        self.blockio_write_method = self.modelBlockIOWrite.dicoV2M[str(text)]
        self.mdl.setBlockIOWriteMethod(self.blockio_write_method)


    @pyqtSlot(int)
    def slotBlockIORankStep(self, text):
        """
        Input for Partitioner.
        """
        self.blockio_rank_step = self.spinBoxIORankStep.value()
        self.mdl.setBlockIORankStep(self.blockio_rank_step)


    @pyqtSlot(int)
    def slotBlockIOMinSize(self, text):
        """
        Input for Partitioner.
        """
        self.blockio_min_size = self.spinBoxIOMinBlockSize.value()
        self.mdl.setBlockIOMinSize(self.blockio_min_size)


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


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
