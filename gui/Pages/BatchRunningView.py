# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
- BatchRunningView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys
import types
import re
import logging
import subprocess

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

import cs_case
import cs_exec_environment
import cs_runcase
import cs_submit

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BatchRunningForm import Ui_BatchRunningForm
from code_saturne.Pages.BatchRunningDebugOptionsHelpDialogForm import Ui_BatchRunningDebugOptionsHelpDialogForm
from code_saturne.Pages.BatchRunningStopByIterationDialogForm import Ui_BatchRunningStopByIterationDialogForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, RegExpValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant
from code_saturne.Base.CommandMgrDialogView import CommandMgrDialogView
from code_saturne.Pages.BatchRunningModel import BatchRunningModel
from code_saturne.Pages.ScriptRunningModel import ScriptRunningModel
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BatchRunningView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Popup help for debug tools
#-------------------------------------------------------------------------------

class BatchRunningDebugOptionsHelpDialogView(QDialog, Ui_BatchRunningDebugOptionsHelpDialogForm):
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

class BatchRunningStopByIterationDialogView(QDialog, Ui_BatchRunningStopByIterationDialogForm):
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

        CommandMgrDialogView.__init__(self, parent, title, cmd, self.case['scripts_path'], self.case['salome'])

        self.pushButtonStop.clicked.connect(self.__slotStop)
        self.pushButtonStopAt.clicked.connect(self.__slotStopAt)
        self.pushButtonCvgTool.clicked.connect(self.__slotCvgTool)

        self.scratch_dir = ""
        self.result_dir = ""
        self.suffix   = ""
        self.listing  = "listing"
        self.n_lines = 0

        self.proc.start(self.cmd)
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
        import subprocess

        cmd = [self.case['package'].config.python]
        cmd.append(os.path.join(self.case['package'].dirs['bindir'][1],
                                self.case['package'].name))
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

class BatchRunningView(QWidget, Ui_BatchRunningForm):
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

        Ui_BatchRunningForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent

        self.case.undoStopGlobal()

        self.mdl = ScriptRunningModel(self.case)

        # Check if the script file name is already defined

        if self.case['scripts_path']:
            if not self.case['runcase']:
                if self.case['package'].runcase in os.listdir(self.case['scripts_path']):
                    runcase_path = os.path.join(self.case['scripts_path'],
                                                self.case['package'].runcase)
                    self.case['runcase'] = cs_runcase.runcase(runcase_path,
                                                              package=self.case['package'])

        self.jmdl = BatchRunningModel(parent, self.case)

        # Get MPI and OpenMP features

        self.have_mpi = False
        self.have_openmp = False
        config_features = self.case['package'].config.features

        config = configparser.ConfigParser()
        config.read(self.case['package'].get_configfiles())

        compute_build_id = -1
        self.compute_versions = None
        if config.has_option('install', 'compute_versions'):
            self.compute_versions = config.get('install', 'compute_versions').split(':')
            compute_build_id = 0
            if len(self.compute_versions) > 1:
                run_build = self.jmdl.dictValues['run_build']
                if self.compute_versions.count(run_build) > 0:
                    compute_build_id = self.compute_versions.index(run_build)
                elif run_build:
                    compute_build_id = -1
                    self.jmdl.dictValues['run_build'] = None
                    self.jmdl.updateBatchFile('run_build')

        if compute_build_id >= 0:
            pkg_compute = self.case['package'].get_alternate_version(self.compute_versions[compute_build_id])
            config_features = pkg_compute.config.features
        if config_features['mpi'] == 'yes':
            self.have_mpi = True
        if config_features['openmp'] == 'yes':
            self.have_openmp = True

        # Batch info

        self.hideBatchInfo()

        self.labelNProcs.hide()
        self.spinBoxNProcs.hide()

        self.labelNThreads.hide()
        self.spinBoxNThreads.hide()

        if compute_build_id < 0:
            self.labelBuildType.hide()
            self.comboBoxBuildType.hide()
        else:
            self.comboBoxBuildType.currentIndexChanged[int].connect(self.slotBuildType)
            for b in self.compute_versions:
                build_type_label = os.path.basename(b)
                if not build_type_label:
                    build_type_label = "[default]"
                self.comboBoxBuildType.addItem(self.tr(build_type_label),
                                               to_qvariant(build_type_label))
            self.comboBoxBuildType.setCurrentIndex(compute_build_id)

        self.class_list = None

        self.__updateRunButton__()

        if self.jmdl.batch.rm_type != None:

            validatorSimpleName = RegExpValidator(self.lineEditJobName,
                                                  QRegExp("[_A-Za-z0-9]*"))
            validatorAccountName = RegExpValidator(self.lineEditJobAccount,
                                                   QRegExp("\\S+"))
            self.lineEditJobName.setValidator(validatorSimpleName)
            self.lineEditJobAccount.setValidator(validatorAccountName)
            self.lineEditJobWCKey.setValidator(validatorAccountName)

        else:

            if self.jmdl.dictValues['run_nprocs'] == None:
                try:
                    # For backwards compatibility
                    # (this is a specific case, as we move information from
                    # the XML model to the batch script)
                    self.jmdl.dictValues['run_nprocs'] = self.mdl.getString('n_procs')
                    self.mdl.setString('n_procs', None)
                except Exception:
                    pass

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

        self.comboBoxRunType.activated[str].connect(self.slotArgRunType)
        self.pushButtonRunSubmit.clicked.connect(self.slotBatchRunning)
        self.lineEditRunId.textChanged[str].connect(self.slotJobRunId)

        self.lineEdit_tool.textChanged[str].connect(self.slotDebug)
        self.toolButton_2.clicked.connect(self.slotDebugHelp)

        self.checkBoxTrace.stateChanged.connect(self.slotTrace)
        self.checkBoxLogParallel.stateChanged.connect(self.slotLogParallel)

        # Combomodels

        self.modelArg_cs_verif = ComboModel(self.comboBoxRunType, 2, 1)

        self.modelArg_cs_verif.addItem(self.tr("Import mesh only"), 'none')
        self.modelArg_cs_verif.addItem(self.tr("Mesh preprocessing"), 'mesh preprocess')
        self.modelArg_cs_verif.addItem(self.tr("Mesh quality criteria"), 'mesh quality')
        if self.case['prepro'] == True:
            self.modelArg_cs_verif.setItem(str_model=self.mdl.getRunType(self.case['prepro']))
            self.labelRunType.show()
            self.comboBoxRunType.show()
        else:
            self.labelRunType.hide()
            self.comboBoxRunType.hide()

        # initialize Widgets

        if self.jmdl.batch.rm_type != None and self.case['runcase']:
            self.displayBatchInfo()

        #rm_type = self.jmdl.batch.rm_type
        run_id = str(self.jmdl.dictValues['run_id'])

        #if run_id == 'None' and self.case['scripts_path']:
        #    run_id, run_title = self.__suggest_run_id()
        #    self.__updateRuncase(run_id)

        if run_id != 'None':
            self.lineEditRunId.setText(run_id)
        else:
            self.lineEditRunId.setText("")

        # Advanced options

        self.debug = self.mdl.getString('debug')
        if self.debug is not None:
            self.lineEdit_tool.setText(str(self.debug))

        if self.mdl.getTrace():
            self.checkBoxTrace.setChecked(True)
        if self.mdl.getLogParallel():
            self.checkBoxLogParallel.setChecked(True)

        # Script info is based on the XML model

        self.displayScriptInfo()

        self.case.undoStartGlobal()


    def __caseIsSaved__(self):
        """
        Launch Code_Saturne batch running.
        """
        # Is the file saved?

        xml_current = os.path.basename(self.case['xmlfile'])
        xml_param = None
        if self.case['runcase']:
            xml_param = self.case['runcase'].get_parameters()
        if not xml_current or xml_current != xml_param:
            self.case['saved'] = "no"

        is_saved = True
        if self.case['saved'] == "no":
            is_saved = False
        if len(self.case['undo']) > 0 or len(self.case['redo']) > 0:
            is_saved = False

        return is_saved


    def __updateRunButton__(self):
        """
        Update push button for run
        """
        if self.__caseIsSaved__():
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
            self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_name')


    @pyqtSlot(int)
    def slotJobNodes(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        n = int(self.spinBoxNodes.text())
        self.jmdl.batch.params['job_nodes'] = str(self.spinBoxNodes.text())
        self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_nodes')
        if self.job_ppn:
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
        self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_ppn')
        if self.job_nodes:
            n = int(self.jmdl.batch.params['job_nodes'])
            tot_ranks = n*ppn
            self.lineEditTotMPI.setText(str(tot_ranks))
        if self.job_threads:
            n_threads = int(self.jmdl.batch.params['job_threads'])
            node_threads = n_threads*ppn
            self.lineEditThreadsPerNode.setText(str(node_threads))


    @pyqtSlot(int)
    def slotJobProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.jmdl.batch.params['job_procs']  = str(self.spinBoxProcs.text())
        self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_procs')


    @pyqtSlot(int)
    def slotJobThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        n_threads = int(self.spinBoxThreads.text())
        self.jmdl.batch.params['job_threads']  = str(n_threads)
        self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_threads')
        if self.job_ppn:
            ppn = int(self.jmdl.batch.params['job_ppn'])
            node_threads = n_threads*ppn
            self.lineEditThreadsPerNode.setText(str(node_threads))


    @pyqtSlot()
    def slotJobWallTime(self):

        h_cput = self.spinBoxDays.value()*24 + self.spinBoxHours.value()
        m_cput = self.spinBoxMinutes.value()
        s_cput = self.spinBoxSeconds.value()
        self.jmdl.batch.params['job_walltime'] = h_cput*3600 + m_cput*60 + s_cput
        self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_walltime')


    @pyqtSlot()
    def slotClass(self):

        self.jmdl.batch.params['job_class'] = str(self.comboBoxClass.currentText())
        if len(self.jmdl.batch.params['job_class']) > 0:
            self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_class')


    @pyqtSlot(str)
    def slotJobAccount(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobAccount.validator().state == QValidator.Acceptable:
            self.jmdl.batch.params['job_account'] = str(v)
            self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_account')


    @pyqtSlot(str)
    def slotJobWCKey(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobWCKey.validator().state == QValidator.Acceptable:
            self.jmdl.batch.params['job_wckey'] = str(v)
            self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_wckey')


    @pyqtSlot(str)
    def slotArgRunType(self, text):
        """
        Input run type option.
        """
        run_type = self.modelArg_cs_verif.dicoV2M[str(text)]
        self.mdl.setRunType(run_type)


    @pyqtSlot(int)
    def slotNProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.jmdl.dictValues['run_nprocs'] = str(v)
        else:
            self.jmdl.dictValues['run_nprocs'] = None
        self.jmdl.updateBatchFile('run_nprocs')


    @pyqtSlot(int)
    def slotNThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.jmdl.dictValues['run_nthreads'] = str(v)
        else:
            self.jmdl.dictValues['run_nthreads'] = None
        self.jmdl.updateBatchFile('run_nthreads')


    @pyqtSlot(int)
    def slotBuildType(self, v):

        if v == 0:
            self.jmdl.dictValues['run_build'] = None
        else:
            self.jmdl.dictValues['run_build'] = str(self.compute_versions[v])
        self.jmdl.updateBatchFile('run_build')
        compute_build_id = v

        pkg_compute = self.case['package'].get_alternate_version(self.compute_versions[compute_build_id])
        config_features = pkg_compute.config.features
        if config_features['mpi'] == 'yes':
            self.have_mpi = True
        else:
            self.have_mpi = False
            self.spinBoxNProcs.setValue(1)
            self.jmdl.dictValues['run_nprocs'] = None
            self.jmdl.updateBatchFile('run_nprocs')

        if config_features['openmp'] == 'yes':
            self.have_openmp = True
        else:
            self.have_openmp = False
            self.spinBoxNThreads.setValue(1)
            self.jmdl.dictValues['run_nthreads'] = None
            self.jmdl.updateBatchFile('run_nthreads')

        self.displayScriptInfo()


    @pyqtSlot(str)
    def slotJobRunId(self, v):
        """
        """
        if self.lineEditRunId.validator().state == QValidator.Acceptable:
            self.jmdl.dictValues['run_id'] = str(v)
            self.jmdl.updateBatchFile('run_id')
            if v != 'None':
                self.__updateRuncase(str(v))


    @pyqtSlot()
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        log.debug("slotAdvancedOptions")

        dialog = BatchRunningAdvancedOptionsDialogView(self)

        if dialog.exec_():
            log.debug("slotAdvancedOptions validated")


    @pyqtSlot()
    def slotBatchRunning(self):
        """
        Launch Code_Saturne batch running.
        """
        # Is the file saved?

        if not self.__caseIsSaved__():
            self.parent.fileSave(renew_page=False)
            self.__updateRunButton__

        if not self.__caseIsSaved__():
            title = self.tr("Warning")
            msg   = self.tr("The current case must be saved before "\
                            "running the computation script.")
            QMessageBox.information(self, title, msg)
            return

        # Ensure code is run from a case subdirectory

        prv_dir = os.getcwd()
        os.chdir(self.case['scripts_path'])

        # Do we have a mesh ?

        have_mesh = False
        node_ecs = self.case.xmlGetNode('solution_domain')
        if node_ecs.xmlGetNode('meshes_list'):
            if node_ecs.xmlGetNode('meshes_list').xmlGetNodeList('mesh'):
                have_mesh = True
        if node_ecs.xmlGetNode('mesh_input', 'path'):
            have_mesh = True
        if not have_mesh:
            title = self.tr("Warning")
            msg   = self.tr("You have to select a mesh.\n\n")
            QMessageBox.information(self, title, msg)
            return

        # Verify if boundary condition definitions exist
        if self.case['prepro'] == False:
            bd = LocalizationModel('BoundaryZone', self.case)
            if not bd.getZones():
                if self.case['no_boundary_conditions'] == False:
                    title = self.tr("Warning")
                    msg   = self.tr("No boundary definition declared.\n\n")
                    QMessageBox.warning(self, title, msg)
                    self.case['no_boundary_conditions'] = True

        # Build command line

        rm_type = self.jmdl.batch.rm_type

        batch = self.case['runcase'].path

        cmd = None
        run_title = None

        run_id = self.jmdl.dictValues['run_id']

        if rm_type == None:
            if not run_id:
                tmp_run_id, run_title = self.__suggest_run_id()
                self.__updateRuncase(tmp_run_id)
            else:
                run_title = self.case['package'].code_name + run_id
            cmd = cs_exec_environment.enquote_arg(batch)
        else:
            run_title = self.case['package'].code_name + ' - Job Submission'
            cmd = cs_exec_environment.enquote_arg(sys.argv[0]) \
                  + ' submit ' +  cs_exec_environment.enquote_arg(batch)

        dlg = ListingDialogView(self.parent, self.case, run_title, cmd)
        dlg.show()

        if rm_type == None and not run_id:
            self.__updateRuncase('')  # remove --id <id> from runcase

        os.chdir(prv_dir)


    @pyqtSlot(str)
    def slotDebug(self, text):
        """
        Input for Debug.
        """
        self.debug = str(text)
        self.mdl.setString('debug', self.debug.strip())


    @pyqtSlot()
    def slotDebugHelp(self):
        """
        Show help page for debugging
        """

        dialog = BatchRunningDebugOptionsHelpDialogView(self)
        dialog.show()

        return


    @pyqtSlot()
    def slotTrace(self):
        """
        Update log type for trace
        """
        if self.checkBoxTrace.isChecked():
            trace = True
        else:
            trace = False
        self.mdl.setTrace(trace)


    @pyqtSlot()
    def slotLogParallel(self):
        """
        Update log type for trace
        """
        if self.checkBoxLogParallel.isChecked():
            logp = True
        else:
            logp = False
        self.mdl.setLogParallel(logp)


    def __suggest_run_id(self):
        """
        Return an id.
        """
        cmd = os.path.join(self.case['package'].get_dir('bindir'),
                           self.case['package'].name)
        cmd = cs_exec_environment.enquote_arg(cmd) + " run --suggest-id"
        xmlfile = self.case['xmlfile']
        if xmlfile:
            cmd += " --param " + cs_exec_environment.enquote_arg(xmlfile)

        r_title = subprocess.Popen(cmd,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   universal_newlines=True).stdout.read()[:-1]
        r_id = os.path.join(self.case['resu_path'], r_title)

        run_id = r_id
        run_title = r_title

        return os.path.basename(run_id), run_title


    def __updateRuncase(self, run_id):
        """
        Update the command line in the launcher C{runcase}.
        """
        runcase = self.case['runcase']

        runcase.set_run_id(run_id=run_id)
        runcase.save()


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

        self.job_name  = self.jmdl.batch.params['job_name']
        self.job_nodes = self.jmdl.batch.params['job_nodes']
        self.job_ppn  = self.jmdl.batch.params['job_ppn']
        self.job_procs = self.jmdl.batch.params['job_procs']
        self.job_threads = self.jmdl.batch.params['job_threads']
        self.job_walltime = self.jmdl.batch.params['job_walltime']
        self.job_class  = self.jmdl.batch.params['job_class']
        self.job_account  = self.jmdl.batch.params['job_account']
        self.job_wckey  = self.jmdl.batch.params['job_wckey']

        if self.job_name != None:
            self.labelJobName.show()
            self.lineEditJobName.setText(str(self.job_name))
            self.lineEditJobName.show()

        if self.job_nodes != None:
            self.labelNodes.show()
            self.spinBoxNodes.setValue(int(self.job_nodes))
            self.spinBoxNodes.show()

        if self.job_ppn != None:
            self.labelPpn.show()
            self.spinBoxPpn.setValue(int(self.job_ppn))
            self.spinBoxPpn.show()

        if self.job_procs != None:
            self.labelProcs.show()
            self.spinBoxProcs.setValue(int(self.job_procs))
            self.spinBoxProcs.show()
        elif self.job_nodes != None and self.job_ppn != None:
            self.labelTotMPI.show()
            self.lineEditTotMPI.show()
            tot_ranks = int(self.job_nodes)*int(self.job_ppn)
            self.lineEditTotMPI.setText(str(tot_ranks))

        if self.job_threads != None:
            self.labelThreads.show()
            self.spinBoxThreads.setValue(int(self.job_threads))
            self.spinBoxThreads.show()
            if self.job_ppn != None:
                self.labelThreadsPerNode.show()
                self.lineEditThreadsPerNode.show()
                th_per_node = int(self.job_threads)*int(self.job_ppn)
                self.lineEditThreadsPerNode.setText(str(th_per_node))

        if self.job_walltime != None:
            seconds = self.job_walltime
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

        if self.job_class != None:

            # Only one pass here
            if self.class_list == None:
                self.class_list = self.jmdl.batch.get_class_list()
                if len(self.class_list) > 0:
                    for c in self.class_list:
                        self.comboBoxClass.addItem(self.tr(c), to_qvariant(c))
                else:
                    c = self.job_class
                    self.comboBoxClass.addItem(self.tr(c), to_qvariant(c))

            # All passes
            try:
                index = self.class_list.index(self.job_class)
                self.comboBoxClass.setCurrentIndex(index)
            except Exception:
                if len(self.class_list) > 0:
                    self.job_class = self.class_list[0]
            self.labelClass.show()
            self.comboBoxClass.show()

            # update runcase (compute class specific to ivanoe)
            self.jmdl.batch.params['job_class'] = str(self.comboBoxClass.currentText())
            if len(self.jmdl.batch.params['job_class']) > 0:
                self.jmdl.batch.update_lines(self.case['runcase'].lines, 'job_class')

        if self.job_account != None:
            self.labelJobAccount.show()
            self.lineEditJobAccount.setText(str(self.job_account))
            self.lineEditJobAccount.show()

        if self.job_wckey != None:
            self.labelJobWCKey.show()
            self.lineEditJobWCKey.setText(str(self.job_wckey))
            self.lineEditJobWCKey.show()

        # Show Job management box

        if self.jmdl.batch.rm_type == 'LOADL':
            self.groupBoxJob.setTitle("Load Leveler job parameters")
            self.labelJobAccount.setText(str("Group"))
            self.lineEditJobAccount.setToolTip("To obtain a list of defined groups, run <b><tt>xloadl</tt></b>, then select <i>File -> Build a Job</i>, and check the group names in the <i>Group</i> field")
        elif self.jmdl.batch.rm_type == 'LSF':
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

        # Update file

        self.jmdl.batch.update_lines(self.case['runcase'].lines)


    def displayScriptInfo(self):
        """
        Layout of the second part of this page.
        """
        if self.jmdl.batch.rm_type == None:
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

        if self.jmdl.batch.rm_type == None:
            if self.have_mpi:
                n_procs_s = self.jmdl.dictValues['run_nprocs']
                if n_procs_s:
                    n_procs = int(n_procs_s)
                else:
                    n_procs = 1
                self.spinBoxNProcs.setValue(n_procs)
            if self.have_openmp:
                n_threads_s = self.jmdl.dictValues['run_nthreads']
                if n_threads_s:
                    n_threads = int(n_threads_s)
                else:
                    n_threads = 1
                self.spinBoxNThreads.setValue(n_threads)
        else:
            pass


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
