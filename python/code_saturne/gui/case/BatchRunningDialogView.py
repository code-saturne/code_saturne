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
This module contains the following classes and function:
- BatchRunningStopByIterationDialogView
- BatchRunningListingLinesDisplayedDialogView
- ListingDialogView
- BatchRunningDialogView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys
import types
import re
import logging
import subprocess

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

from code_saturne.base import cs_case
from code_saturne.base import cs_exec_environment
from code_saturne.base import cs_batch, cs_submit

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BatchRunningDialogForm import Ui_BatchRunningDialogForm
from code_saturne.gui.case.BatchRunningDebugOptionsHelpDialogForm import Ui_BatchRunningDebugOptionsHelpDialogForm
from code_saturne.gui.case.BatchRunningStopByIterationDialogForm import Ui_BatchRunningStopByIterationDialogForm

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, IntValidator, RegExpValidator
from code_saturne.gui.base.QtPage import from_qvariant
from code_saturne.gui.base.CommandMgrDialogView import CommandMgrDialogView
from code_saturne.model.BatchRunningModel import BatchRunningModel
from code_saturne.model.ScriptRunningModel import ScriptRunningModel
from code_saturne.model.LocalizationModel import LocalizationModel, Zone

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BatchRunningDialogView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Popup help for debug tools
#-------------------------------------------------------------------------------

class BatchRunningDebugOptionsHelpDialogView(QDialog,
                                             Ui_BatchRunningDebugOptionsHelpDialogForm):
    """
    Help dialog
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningDebugOptionsHelpDialogForm.__init__(self)
        self.setupUi(self)
        self.retranslateUi(self)

        self.setWindowTitle(self.tr("Debugger syntax help"))
        self.parent = parent

        # Connections
        self.pushButtonClose.clicked.connect(self.slotDebugclose)


    @pyqtSlot()
    def slotDebugclose(self):
        """
        Close debugging page
        """

        QDialog.accept(self)
        return


#-------------------------------------------------------------------------------
# Popup window class: stop the computation at a iteration
#-------------------------------------------------------------------------------

class BatchRunningStopByIterationDialogView(QDialog,
                                            Ui_BatchRunningStopByIterationDialogForm):
    """
    Advanced dialog to stop the computation at a given iteration
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningStopByIterationDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Stop"))

        self.default = default
        self.result  = self.default.copy()

        v = IntValidator(self.lineEditStopIter, min=1)
        v.setExclusiveMin(True)
        self.lineEditStopIter.setValidator(v)

        # Previous values
        self.iter = self.default['iter']
        self.lineEditStopIter.setText(str(self.iter))

        self.lineEditStopIter.textChanged[str].connect(self.__slotStopIter)


    @pyqtSlot(str)
    def __slotStopIter(self, text):
        """
        Private slot to set an iteration number to stop the code.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            iter = from_qvariant(text, int)
            self.iter = iter


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        self.result['iter'] = self.iter
        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ListingDialogView(CommandMgrDialogView):
    def __init__(self, parent, case, title, cmd):
        self.case = case

        CommandMgrDialogView.__init__(self,
                                      parent,
                                      title,
                                      cmd,
                                      self.case['data_path'],
                                      self.case['salome'])

        self.pushButtonStop.clicked.connect(self.__slotStop)
        self.pushButtonStopAt.clicked.connect(self.__slotStopAt)
        self.pushButtonCvgTool.clicked.connect(self.__slotCvgTool)

        self.scratch_dir = ""
        self.result_dir = ""
        self.suffix   = ""
        self.listing  = "run_solver.log"
        self.n_lines = 0

        # When running under "code_saturne salome" session, we may need to
        # work around an issue due to extra (cumulative) entries in the
        # PYTHONPATH, so we use  value saved at session launch.
        pythonpath_save = None
        if self.objBr:
            pythonpath_top = os.getenv('CS_SALOME_TOP_PYTHONPATH')
            if pythonpath_top:
                pythonpath_save = os.getenv('PYTHONPATH')
                os.environ['PYTHONPATH'] = pythonpath_top

        # Start process
        self.proc.start(self.cmd)

        # Now process is launched, restore environment
        if pythonpath_save:
            os.environ['PYTHONPATH'] = pythonpath_save

        self.proc.waitForFinished(100)


    @pyqtSlot()
    def slotReadFromStdout(self):
        """
        Public slot to handle the readyReadStandardOutput signal of the process.
        """
        if self.proc is None:
            return
        self.proc.setReadChannel(QProcess.StandardOutput)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            s = (ba.data()).decode("utf-8")[:-1]
            self.logText.append(s)
            self.n_lines += 1

            # Work and result directories printed in first lines of log.
            if self.n_lines < 16:
                self.__execDir(s)


    def __execDir(self, s):
        """
        Private method. Find the directory of the code execution.
        """
        # Read directly the run directory from the sdtout of the code.

        if "Working directory" in s:
            self.scratch_dir = "Working directory"
            return
        elif self.scratch_dir == "Working directory":
            self.scratch_dir = " ".join(str(s).split())
            title = os.path.basename(self.scratch_dir)
            self.setWindowTitle(title)
            self.suffix = title
            return

        if not self.result_dir:
            if "Result directory" in s:
                self.result_dir = "Result directory"
        elif self.result_dir == "Result directory":
            self.result_dir = " ".join(str(s).split())
            title = os.path.basename(self.result_dir)
            self.setWindowTitle(title)
            self.suffix = title


    def __stopExec(self, iter, msg):
        """
        Private method. Stops the code.
        """
        line = "\n" + str(iter) + "\n\n"
        if self.scratch_dir:
            exec_dir = self.scratch_dir
        elif self.result_dir:
            exec_dir = self.result_dir
        else:
            return
        fstp = os.path.join(exec_dir, "control_file")
        f = open(fstp, 'w')
        f.write(line)
        f.close()
        QMessageBox.warning(self, self.tr("Warning"), msg)


    @pyqtSlot()
    def __slotStop(self):
        """
        Private slot. Stops the code at the end of the current iteration.
        """
        iter = 1
        msg = self.tr("Stop at the end of the current iteration.")
        self.__stopExec(iter, msg)


    @pyqtSlot()
    def __slotStopAt(self):
        """
        Private slot. Stops the code at the end of the given iteration.
        """
        default = {}
        default['iter'] = 100
        dlg = BatchRunningStopByIterationDialogView(self, default)
        if dlg.exec_():
            result = dlg.get_result()
            msg = self.tr("Stop at iteration number: %i" % result['iter'])
            self.__stopExec(result['iter'], msg)


    @pyqtSlot()
    def __slotCvgTool(self):
        """
        Private slot. Running convergence tracking tool
        """
        cmd = [os.path.join(self.case['package'].dirs['bindir'],
                            self.case['package'].name)]
        cmd.append("trackcvg")
        if self.scratch_dir:
            cmd.append("-r")
            cmd.append(self.scratch_dir)
        elif self.result_dir:
            cmd.append("-r")
            cmd.append(self.result_dir)
        self.runProcess = subprocess.Popen(cmd)


    @pyqtSlot("QProcess::ProcessState")
    def slotStateChanged(self, state):
        """
        Public slot. Handle the current status of the process.
        """
        bool = not(state == QProcess.NotRunning)
        self.pushButtonKill.setEnabled(bool)
        self.pushButtonStop.setEnabled(bool)
        self.pushButtonStopAt.setEnabled(bool)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BatchRunningDialogView(QDialog, Ui_BatchRunningDialogForm):
    """
    This class is devoted to the queue selection.
    If the batch script file name is known, informations are displayed
    in the appropiate widget.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BatchRunningDialogForm.__init__(self)
        self.setupUi(self)


        self.case = case
        self.parent = parent

        case_is_saved = not self.case.isModified()

        title = self.tr("Run computation")
        self.setWindowTitle(title)

        self.case.undoStopGlobal()

        self.mdl = ScriptRunningModel(self.case)
        self.jmdl = self.case['job_model']
        self.batch_help = None

        if self.jmdl is None:
            if self.case['data_path']:
                run_conf_path = os.path.join(self.case['data_path'], 'run.cfg')
            else:
                run_conf_path = None
            self.jmdl = BatchRunningModel(path=run_conf_path,
                                          pkg=self.case['package'])

        self.jmdl.load()

        # Remark: if the view was closed without settings being applied,
        # the job model batch dictionnary may have beed updated for display
        # purposes but not have been applied to the job header lines.
        # The batch object (and its dictionnary) could be moved from the
        # job model to the current object.

        if self.jmdl.job_header_lines != None:
            self.jmdl.batch.parse_lines(self.jmdl.job_header_lines)
            self.job_header_lines = self.jmdl.job_header_lines.copy()
        else:
            self.job_header_lines = None

        self.run_dict = {}
        for k in self.jmdl.run_dict:
            self.run_dict[k] = self.jmdl.run_dict[k]
        self.job_dict = {}
        for k in self.jmdl.job_dict:
            self.job_dict[k] = self.jmdl.job_dict[k]

        self.have_mpi = self.jmdl.have_mpi
        self.have_openmp = self.jmdl.have_openmp


        # Batch info

        self.hideBatchInfo()

        self.labelNProcs.hide()
        self.spinBoxNProcs.hide()

        self.labelNThreads.hide()
        self.spinBoxNThreads.hide()

        compute_build_id = -1
        if self.jmdl.compute_builds:
            compute_build = self.run_dict['compute_build']
            if compute_build in self.jmdl.compute_builds:
                compute_build_id = self.jmdl.compute_builds.index(compute_build)
            else:
                compute_build_id = 0
                self.run_dict['compute_build'] = self.jmdl.compute_builds[0]
        if compute_build_id < 0:
            self.labelBuildType.hide()
            self.comboBoxBuildType.hide()
        else:
            self.comboBoxBuildType.currentIndexChanged[int].connect(self.slotBuildType)
            for b in self.jmdl.compute_builds:
                build_type_label = os.path.basename(b)
                if not build_type_label:
                    build_type_label = "[default]"
                self.comboBoxBuildType.addItem(self.tr(build_type_label),
                                               build_type_label)
            self.comboBoxBuildType.setCurrentIndex(compute_build_id)

        self.class_list = None

        self.__updateRunButton__(case_is_saved)

        if self.jmdl.batch.rm_type != None:
            validatorSimpleName = RegExpValidator(self.lineEditJobName,
                                                  QRegExp("[_A-Za-z0-9]*"))
            validatorAccountName = RegExpValidator(self.lineEditJobAccount,
                                                   QRegExp("\\S+"))
            self.lineEditJobName.setValidator(validatorSimpleName)
            self.lineEditJobAccount.setValidator(validatorAccountName)
            self.lineEditJobWCKey.setValidator(validatorAccountName)

        validatorRunId = RegExpValidator(self.lineEditRunId,
                                         QRegExp("[_A-Za-z0-9]*"))
        self.lineEditRunId.setValidator(validatorRunId)

        # Connections

        if self.jmdl.batch.rm_type != None:
            self.lineEditJobName.textChanged[str].connect(self.slotJobName)
            self.spinBoxNodes.valueChanged[int].connect(self.slotJobNodes)
            self.spinBoxPpn.valueChanged[int].connect(self.slotJobPpn)
            self.spinBoxProcs.valueChanged[int].connect(self.slotJobProcs)
            self.spinBoxThreads.valueChanged[int].connect(self.slotJobThreads)
            self.spinBoxDays.valueChanged[int].connect(self.slotJobWallTime)
            self.spinBoxHours.valueChanged[int].connect(self.slotJobWallTime)
            self.spinBoxMinutes.valueChanged[int].connect(self.slotJobWallTime)
            self.spinBoxSeconds.valueChanged[int].connect(self.slotJobWallTime)
            self.comboBoxClass.activated[str].connect(self.slotClass)
            self.lineEditJobAccount.textChanged[str].connect(self.slotJobAccount)
            self.lineEditJobWCKey.textChanged[str].connect(self.slotJobWCKey)

        else:
            self.spinBoxNProcs.valueChanged[int].connect(self.slotNProcs)
            self.spinBoxNThreads.valueChanged[int].connect(self.slotNThreads)

        self.pushButtonCancel.clicked.connect(self.slotCancel)
        self.pushButtonApply.clicked.connect(self.slotApply)
        self.pushButtonRunSubmit.clicked.connect(self.slotBatchRunning)

        self.lineEditRunId.textChanged[str].connect(self.slotJobRunId)
        self.lineEdit_tool.textChanged[str].connect(self.slotDebug)
        self.toolButton_2.clicked.connect(self.slotDebugHelp)

        self.checkBoxInitOnly.stateChanged.connect(self.slotInitOnly)

        self.checkBoxTrace.stateChanged.connect(self.slotTrace)
        self.checkBoxLogParallel.stateChanged.connect(self.slotLogParallel)

        self.tabWidgetJob.setHidden(True)
        self.haveAdvparam = False

        self.pushButtonAdvparam.setCheckable(True)
        self.pushButtonAdvparam.pressed.connect(self.switchBox)

        if not self.case['data_path']:
            self.pushButtonRunSubmit.setEnabled(False)

        # initialize Widgets

        if self.jmdl.batch.rm_type != None:
            self.displayBatchInfo()

        run_id = str(self.run_dict['id'])

        if run_id != 'None':
            self.lineEditRunId.setText(run_id)
        else:
            self.lineEditRunId.setText("")

        # Advanced options

        self.debug = self.job_dict['debug_args']
        if self.debug is not None:
            self.lineEdit_tool.setText(str(self.debug))

        self.trace_iter = self.mdl.getTrace()
        self.log_parallel = self.mdl.getLogParallel()

        if self.trace_iter:
            self.checkBoxTrace.setChecked(True)
        if self.log_parallel:
            self.checkBoxLogParallel.setChecked(True)

        self.checkBoxInitOnly.setChecked(self.run_dict['initialize'] == True)

        # Script info is based on the XML model

        self.displayScriptInfo()

        # self.resize(self.minimumSizeHint())

        self.case.undoStartGlobal()


    def __updateRunButton__(self, case_is_saved):
        """
        Update push button for run
        """
        if case_is_saved:
            if self.jmdl.batch.rm_type != None:
                self.pushButtonRunSubmit.setText("Submit job")
            else:
                self.pushButtonRunSubmit.setText("Run calculation")

        else:
            if self.jmdl.batch.rm_type != None:
                self.pushButtonRunSubmit.setText("Save case and submit job")
            else:
                self.pushButtonRunSubmit.setText("Save case and run calculation")


    @pyqtSlot(str)
    def slotJobName(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobName.validator().state == QValidator.Acceptable:
            self.jmdl.batch.params['job_name'] = str(v)


    @pyqtSlot(int)
    def slotJobNodes(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        n = int(self.spinBoxNodes.text())
        self.jmdl.batch.params['job_nodes'] = str(self.spinBoxNodes.text())
        if self.jmdl.batch.params['job_ppn']:
            ppn = int(self.jmdl.batch.params['job_ppn'])
            tot_ranks = n*ppn
            self.lineEditTotMPI.setText(str(tot_ranks))


    @pyqtSlot(int)
    def slotJobPpn(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        ppn = int(self.spinBoxPpn.text())
        self.jmdl.batch.params['job_ppn']  = str(ppn)
        if self.jmdl.batch.params['job_nodes']:
            n = int(self.jmdl.batch.params['job_nodes'])
            tot_ranks = n*ppn
            self.lineEditTotMPI.setText(str(tot_ranks))
        if self.jmdl.batch.params['job_threads']:
            n_threads = int(self.jmdl.batch.params['job_threads'])
            node_threads = n_threads*ppn
            self.lineEditThreadsPerNode.setText(str(node_threads))


    @pyqtSlot(int)
    def slotJobProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.job_procs  = str(self.spinBoxProcs.text())


    @pyqtSlot(int)
    def slotJobThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        n_threads = int(self.spinBoxThreads.text())
        self.jmdl.batch.params['job_threads'] = str(n_threads)
        if self.jmdl.batch.params['job_ppn']:
            ppn = int(self.jmdl.batch.params['job_ppn'])
            node_threads = n_threads*ppn
            self.lineEditThreadsPerNode.setText(str(node_threads))


    @pyqtSlot()
    def slotJobWallTime(self):

        h_cput = self.spinBoxDays.value()*24 + self.spinBoxHours.value()
        m_cput = self.spinBoxMinutes.value()
        s_cput = self.spinBoxSeconds.value()
        self.job_walltime = h_cput*3600 + m_cput*60 + s_cput


    @pyqtSlot()
    def slotClass(self):

        self.jmdl.batch.params['job_class'] \
            = str(self.comboBoxClass.currentText())


    @pyqtSlot(str)
    def slotJobAccount(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobAccount.validator().state == QValidator.Acceptable:
            self.job_account = str(v)


    @pyqtSlot(str)
    def slotJobWCKey(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobWCKey.validator().state == QValidator.Acceptable:
            self.job_wckey = str(v)


    @pyqtSlot(int)
    def slotNProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.job_dict['n_procs'] = str(v)
        else:
            self.job_dict['n_procs'] = None


    @pyqtSlot(int)
    def slotNThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.job_dict['n_threads'] = str(v)
        else:
            self.job_dict['n_threads'] = None


    @pyqtSlot(int)
    def slotBuildType(self, v):

        if v == 0:
            self.run_dict['compute_build'] = None
        else:
            self.run_dict['compute_build'] = str(self.jmdl.compute_builds[v])
        compute_build_id = v

        pkg_compute = self.jmdl.pkg.get_alternate_version(self.jmdl.compute_builds[compute_build_id])
        config_features = pkg_compute.config.features
        if config_features['mpi'] == 'yes':
            self.have_mpi = True
        else:
            self.have_mpi = False
            self.spinBoxNProcs.setValue(1)
            self.jmdl.job_dict['n_procs'] = None

        if config_features['openmp'] == 'yes':
            self.have_openmp = True
        else:
            self.have_openmp = False
            self.spinBoxNThreads.setValue(1)
            self.jmdl.job_dict['n_threads'] = None

        self.displayScriptInfo()


    @pyqtSlot(str)
    def slotJobRunId(self, v):
        """
        """
        if self.lineEditRunId.validator().state == QValidator.Acceptable:
            self.run_dict['id'] = str(v)


    @pyqtSlot()
    def slotCancel(self):
        """
        Close dialog with no modifications
        """

        QDialog.accept(self)


    @pyqtSlot()
    def switchBox(self):
        if not self.pushButtonAdvparam.isChecked():

            if self.batch_help == None:
                self.batch_help = self.jmdl.batch.get_help_text(self.case['package'])
                self.textBrowserBatchHelp.setOpenExternalLinks(True)
                # Set fixed font for plain text (ignored by HTML)
                self.textBrowserBatchHelp.setFontFamily("monospace")
                # Set text (autodetects text/html; Markdown would need to call
                # self.textBrowserBatchHelp.setMarkdown().
                self.textBrowserBatchHelp.setText(self.batch_help)

            self.textEditJob.clear()

            self.jmdl.batch.update_lines(self.job_header_lines)

            for k in self.job_header_lines:
                self.textEditJob.append(k)

            self.tabWidgetJob.setVisible(True)
            self.haveAdvparam = True
            self.frameJob.setHidden(True)

        else:
            self.job_header_lines = self.textEditJob.toPlainText().split(os.linesep)
            self.jmdl.batch.parse_lines(self.job_header_lines)

            self.displayBatchInfo()

            self.tabWidgetJob.setHidden(True)
            self.haveAdvparam = False
            self.frameJob.setVisible(True)


    def Apply(self):
        """
        Apply changes
        """

        self.mdl.setTrace(self.trace_iter)
        self.mdl.setLogParallel(self.log_parallel)

        # Apply state
        if self.jmdl.job_header_lines != None:
            if self.pushButtonAdvparam.isChecked():
                self.job_header_lines = self.textEditJob.toPlainText().split(os.linesep)
            else:
                self.jmdl.batch.update_lines(self.job_header_lines)

            self.jmdl.job_header_lines = self.job_header_lines.copy()

        for k in self.jmdl.run_dict:
            self.jmdl.run_dict[k] = self.run_dict[k]
        for k in self.jmdl.job_dict:
            self.jmdl.job_dict[k] = self.job_dict[k]

        self.jmdl.updateComputeBuildInfo(compute_build=self.jmdl.run_dict['compute_build'])


    @pyqtSlot()
    def slotApply(self):
        """
        Apply changes and close dialog without running computation
        """

        self.Apply()
        QDialog.accept(self)


    @pyqtSlot()
    def slotBatchRunning(self):
        """
        Launch code_saturne batch running.
        """

        QDialog.accept(self)

        self.Apply()

        if self.case.isModified():
            self.parent.fileSave()
        else:
            self.jmdl.save(param=self.case['xmlfile'], force=False)

        # Ensure code is run from a case subdirectory

        prv_dir = os.getcwd()
        os.chdir(self.case['data_path'])

        # Build command line

        rm_type = self.jmdl.batch.rm_type

        cmd = sys.argv[0]

        if not cmd:
           pkg = self.case['package']
           bindir = pkg.get_dir('bindir')
           cmd = os.path.join(bindir, pkg.name)

        cmd = cs_exec_environment.enquote_arg(cmd)

        run_title = self.case.module_name()

        if rm_type is None:
            run_title += ' - Job Run'
            cmd += ' run'
        else:
            run_title += ' - Job Submission'
            cmd += ' submit'

        dlg = ListingDialogView(self.parent, self.case, run_title, cmd)
        dlg.show()

        os.chdir(prv_dir)


    @pyqtSlot(str)
    def slotDebug(self, text):
        """
        Input for Debug.
        """
        self.debug = str(text)
        self.job_dict['debug_args'] = self.debug.strip()


    @pyqtSlot()
    def slotDebugHelp(self):
        """
        Show help page for debugging
        """

        dialog = BatchRunningDebugOptionsHelpDialogView(self)
        dialog.show()


    @pyqtSlot()
    def slotInitOnly(self):
        """
        Set initialization only option
        """
        self.run_dict['initialize'] = self.checkBoxInitOnly.isChecked()


    @pyqtSlot()
    def slotTrace(self):
        """
        Update log type for trace
        """
        if self.checkBoxTrace.isChecked():
            self.trace_iter = True
        else:
            self.trace_iter = False


    @pyqtSlot()
    def slotLogParallel(self):
        """
        Update log type for trace
        """
        if self.checkBoxLogParallel.isChecked():
            self.log_parallel = True
        else:
            self.log_parallel = False


    def hideBatchInfo(self):
        """
        hide all batch info before displaying a selected subset
        """

        self.groupBoxJob.hide()

        self.labelJobName.hide()
        self.lineEditJobName.hide()
        self.labelNodes.hide()
        self.spinBoxNodes.hide()
        self.labelPpn.hide()
        self.spinBoxPpn.hide()
        self.labelProcs.hide()
        self.spinBoxProcs.hide()
        self.labelThreads.hide()
        self.spinBoxThreads.hide()
        self.labelTotMPI.hide()
        self.lineEditTotMPI.hide()
        self.labelThreadsPerNode.hide()
        self.lineEditThreadsPerNode.hide()
        self.labelClass.hide()
        self.labelWTime.hide()
        self.spinBoxDays.hide()
        self.labelDays.hide()
        self.spinBoxHours.hide()
        self.labelHours.hide()
        self.spinBoxMinutes.hide()
        self.labelMinutes.hide()
        self.spinBoxSeconds.hide()
        self.labelSeconds.hide()
        self.comboBoxClass.hide()
        self.labelJobAccount.hide()
        self.lineEditJobAccount.hide()
        self.labelJobWCKey.hide()
        self.lineEditJobWCKey.hide()


    def displayBatchInfo(self):
        """
        Layout of the second part of this page.
        """

        job_name = self.jmdl.batch.params['job_name']
        job_nodes = self.jmdl.batch.params['job_nodes']
        job_ppn = self.jmdl.batch.params['job_ppn']
        job_procs = self.jmdl.batch.params['job_procs']
        job_threads = self.jmdl.batch.params['job_threads']
        job_walltime = self.jmdl.batch.params['job_walltime']
        job_class  = self.jmdl.batch.params['job_class']
        job_account  = self.jmdl.batch.params['job_account']
        job_wckey  = self.jmdl.batch.params['job_wckey']

        if job_name != None:
            self.labelJobName.show()
            self.lineEditJobName.setText(str(job_name))
            self.lineEditJobName.show()

        if job_nodes != None:
            self.labelNodes.show()
            self.spinBoxNodes.setValue(int(job_nodes))
            self.spinBoxNodes.show()

        if job_ppn != None:
            self.labelPpn.show()
            self.spinBoxPpn.setValue(int(job_ppn))
            self.spinBoxPpn.show()

        if job_procs != None:
            self.labelProcs.show()
            self.spinBoxProcs.setValue(int(job_procs))
            self.spinBoxProcs.show()
        elif job_nodes != None and job_ppn != None:
            self.labelTotMPI.show()
            self.lineEditTotMPI.show()
            tot_ranks = int(job_nodes)*int(job_ppn)
            self.lineEditTotMPI.setText(str(tot_ranks))

        if job_threads != None:
            self.labelThreads.show()
            self.spinBoxThreads.setValue(int(job_threads))
            self.spinBoxThreads.show()
            if self.jmdl.batch.params['job_ppn'] != None:
                self.labelThreadsPerNode.show()
                self.lineEditThreadsPerNode.show()
                th_per_node = int(job_threads)*int(job_ppn)
                self.lineEditThreadsPerNode.setText(str(th_per_node))

        if job_walltime != None:
            seconds = job_walltime
            minutes = seconds / 60
            hours = minutes / 60
            days = hours / 24
            seconds = seconds % 60
            minutes = minutes % 60
            hours = hours % 24
            self.spinBoxDays.setValue(days)
            self.spinBoxHours.setValue(hours)
            self.spinBoxMinutes.setValue(minutes)
            self.spinBoxSeconds.setValue(seconds)
            self.labelWTime.show()
            self.spinBoxDays.show()
            self.labelDays.show()
            self.spinBoxHours.show()
            self.labelHours.show()
            self.spinBoxMinutes.show()
            self.labelMinutes.show()
            self.spinBoxSeconds.show()
            self.labelSeconds.show()

        if job_class != None:

            # Only one pass here
            if self.class_list is None:
                self.class_list = self.jmdl.batch.get_class_list()
                if len(self.class_list) > 0:
                    for c in self.class_list:
                        self.comboBoxClass.addItem(self.tr(c), c)
                else:
                    c = job_class
                    self.comboBoxClass.addItem(self.tr(c), c)

            # All passes
            try:
                index = self.class_list.index(job_class)
                self.comboBoxClass.setCurrentIndex(index)
            except Exception:
                if len(self.class_list) > 0:
                    job_class = self.class_list[0]
            self.labelClass.show()
            self.comboBoxClass.show()

            # update job info
            job_class = str(self.comboBoxClass.currentText())
            self.jmdl.batch.params['job_class'] = job_class

        if job_account != None:
            self.labelJobAccount.show()
            self.lineEditJobAccount.setText(str(job_account))
            self.lineEditJobAccount.show()

        if job_wckey != None:
            self.labelJobWCKey.show()
            self.lineEditJobWCKey.setText(str(job_wckey))
            self.lineEditJobWCKey.show()

        # Show Job management box

        if self.jmdl.batch.rm_type == 'LSF':
            self.groupBoxJob.setTitle("LSF job parameters")
        elif self.jmdl.batch.rm_type == 'PBS':
            self.groupBoxJob.setTitle("PBS job parameters")
        elif self.jmdl.batch.rm_type == 'SGE':
            self.groupBoxJob.setTitle("Sun Grid Engine job parameters")
        if self.jmdl.batch.rm_type == 'SLURM':
            self.groupBoxJob.setTitle("SLURM job parameters")
        else:
            self.groupBoxJob.setTitle("Batch job parameters")

        self.groupBoxJob.show()


    def displayScriptInfo(self):
        """
        Layout of the second part of this page.
        """
        if self.jmdl.batch.rm_type is None:
            if self.have_mpi:
                self.labelNProcs.show()
                self.spinBoxNProcs.show()
            else:
                self.labelNProcs.hide()
                self.spinBoxNProcs.hide()
            if self.have_openmp:
                self.labelNThreads.show()
                self.spinBoxNThreads.show()
            else:
                self.labelNThreads.hide()
                self.spinBoxNThreads.hide()
        else:
            self.labelNProcs.hide()
            self.spinBoxNProcs.hide()
            self.labelNThreads.hide()
            self.spinBoxNThreads.hide()

        if self.jmdl.batch.rm_type is None:
            if self.have_mpi:
                n_procs_s = self.jmdl.job_dict['n_procs']
                if n_procs_s:
                    n_procs = int(n_procs_s)
                else:
                    n_procs = 1
                self.spinBoxNProcs.setValue(n_procs)
            if self.have_openmp:
                n_threads_s = self.jmdl.job_dict['n_threads']
                if n_threads_s:
                    n_threads = int(n_threads_s)
                else:
                    n_threads = 1
                self.spinBoxNThreads.setValue(n_threads)
        else:
            pass


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
