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
This module defines the 'PerformanceTuning' page.

This module contains the following classes:
- PerformanceTuningView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types, shutil
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
from Base.QtPage import ComboModel, IntValidator, RegExpValidator, setGreenColor
from Pages.SolutionDomainModel import RelOrAbsPath
from Pages.PerformanceTuningForm import Ui_PerformanceTuningForm
from Pages.PerformanceTuningModel import PerformanceTuningModel

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
    This page is devoted to the start/restart control.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_PerformanceTuningForm.__init__(self)
        self.setupUi(self)
        self.case = case

        self.mdl = PerformanceTuningModel(self.case)

        # Combo models and items

        self.modelPartType = ComboModel(self.comboBox_PartType, 4, 1)
        self.modelPartOut = ComboModel(self.comboBox_PartOutput, 4, 1)

        self.modelPartType.addItem(self.tr("Default"), 'default')
        self.modelPartType.addItem(self.tr("PT-SCOTCH / SCOTCH"), 'scotch')
        self.modelPartType.addItem(self.tr("ParMETIS / METIS"), 'metis')
        self.modelPartType.addItem(self.tr("Morton curve (bounding box)"), 'morton sfc')
        self.modelPartType.addItem(self.tr("Morton curve (bounding cube)"), 'morton sfc cube')
        self.modelPartType.addItem(self.tr("Hilbert curve (bounding box)"), 'hilbert sfc')
        self.modelPartType.addItem(self.tr("Hilbert curve (bounding cube)"), 'hilbert sfc cube')
        self.modelPartType.addItem(self.tr("Block (unoptimized)"), 'block')

        self.modelPartOut.addItem(self.tr("No"), 'no')
        self.modelPartOut.addItem(self.tr("For graph-based partitioning"), 'default')
        self.modelPartOut.addItem(self.tr("Yes"), 'yes')

        # Validators

        rankStepVd = IntValidator(self.lineEdit_RankStep, min=1)
        self.lineEdit_RankStep.setValidator(rankStepVd)
        partListVd = RegExpValidator(self.lineEdit_PartList, QRegExp("[0-9- ]*"))
        self.lineEdit_PartList.setValidator(partListVd)

        # Connections

        self.connect(self.radioButtonYes, SIGNAL("clicked()"), self.slotPartition)
        self.connect(self.radioButtonNo, SIGNAL("clicked()"), self.slotPartition)
        self.connect(self.toolButton_PartInputDir, SIGNAL("pressed()"), self.slotSearchPartInputDirectory)
        self.connect(self.comboBox_PartOutput, SIGNAL("activated(const QString&)"), self.slotPartOut)

        self.connect(self.comboBox_PartType, SIGNAL("activated(const QString&)"), self.slotPartType)
        self.connect(self.lineEdit_PartList, SIGNAL("textChanged(const QString &)"), self.slotPartitionList)
        self.connect(self.lineEdit_RankStep, SIGNAL("textChanged(const QString &)"), self.slotRankStep)

        self.partition_alg = str(self.mdl.getPartitionType())
        self.modelPartType.setItem(str_model=self.partition_alg)

        self.partition_out = str(self.mdl.getPartitionOut())
        self.modelPartOut.setItem(str_model=self.partition_out)

        self.partition_list = str(self.mdl.getPartitionList())
        self.lineEdit_PartList.setText(QString(self.partition_list))

        self.rank_step = self.mdl.getPartitionRankStep()
        self.lineEdit_RankStep.setText(QString(str(self.rank_step)))

        self.connect(self.checkBox_IgnorePerio, SIGNAL("clicked(bool)"), self.slotIgnorePerio)

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

        self.slotPartition()

        if self.mdl.getIgnorePerio() == 'on':
            self.checkBox_IgnorePerio.setChecked(True)
            self.slotIgnorePerio(True)
        else:
            self.checkBox_IgnorePerio.setChecked(False)
            self.slotIgnorePerio(False)

    @pyqtSignature("")
    def slotSearchPartInputDirectory(self):
        """
        Search restart file (directory) in list of directories
        """
        title    = self.tr("Select checkpoint/restart directory")

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

            dir_path = str(s.first())
            dir_path = os.path.abspath(dir_path)

            self.partinput_path = RelOrAbsPath(dir_path, self.case['case_path'])
            self.mdl.setPartitionInputPath(self.partinput_path)
            self.lineEditPartInputDir.setText(self.partinput_path)

            log.debug("slotSearchPartInputDirectory-> %s" % self.partinput_path)


    @pyqtSignature("")
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
                setGreenColor(self.toolButton_PartInputDir)
            else:
                setGreenColor(self.toolButton_PartInputDir, False)

        else:
            self.mdl.setPartitionInputPath(None)
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)

            self.framePartInputDir.hide()
            self.lineEditPartInputDir.setText("")


    @pyqtSignature("const QString &")
    def slotPartitionList(self, text):
        """
        Input for Partitioner.
        """
        self.partition_list = str(text)
        self.mdl.setPartitionList(self.partition_list.strip())


    @pyqtSignature("const QString &")
    def slotPartOut(self, text):
        """
        Partitioner execution mode option.
        """
        self.partition_out = self.modelPartOut.dicoV2M[str(text)]
        self.mdl.setPartitionOut(self.partition_out)


    @pyqtSignature("const QString &")
    def slotPartType(self, text):
        """
        Partitioner execution mode option.
        """
        self.partition_alg = self.modelPartType.dicoV2M[str(text)]
        self.mdl.setPartitionType(self.partition_alg)


    @pyqtSignature("const QString &")
    def slotRankStep(self, text):
        """
        Input for Partitioner.
        """
        self.rank_step, ok = self.lineEdit_RankStep.text().toInt()
        self.mdl.setPartitionRankStep(self.rank_step)


    @pyqtSignature("bool")
    def slotIgnorePerio(self, checked):
        """
        Ignore periodicity.
        """
        if checked:
            self.mdl.setIgnorePerio("on")
        else:
            self.mdl.setIgnorePerio("off")


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
