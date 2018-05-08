# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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
- BatchRunningAdvancedOptionsDialogView
- BatchRunningStopByIterationDialogView
- BatchRunningListingLinesDisplayedDialogView
- ListingDialogView
- BatchRunningView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys
import string, types
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

import cs_case
import cs_exec_environment
import cs_runcase

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BatchRunningForm import Ui_BatchRunningForm
from code_saturne.Pages.BatchRunningAdvancedOptionsDialogForm import Ui_BatchRunningAdvancedOptionsDialogForm
from code_saturne.Pages.BatchRunningStopByIterationDialogForm import Ui_BatchRunningStopByIterationDialogForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, RegExpValidator, setGreenColor
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
# Popup advanced options
#-------------------------------------------------------------------------------


class BatchRunningAdvancedOptionsDialogView(QDialog, Ui_BatchRunningAdvancedOptionsDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_BatchRunningAdvancedOptionsDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Advanced options"))
        self.parent = parent

        # Combo models
        self.modelCSOUT1       = ComboModel(self.comboBox_6, 2, 1)
        self.modelCSOUT2       = ComboModel(self.comboBox_7, 3, 1)

        # Combo items
        self.modelCSOUT1.addItem(self.tr("to standard output"), 'stdout')
        self.modelCSOUT1.addItem(self.tr("to listing"), 'listing')

        self.modelCSOUT2.addItem(self.tr("no output"), 'null')
        self.modelCSOUT2.addItem(self.tr("to standard output"), 'stdout')
        self.modelCSOUT2.addItem(self.tr("to listing_n<p>"), 'listing')

        # Connections
        self.connect(self.toolButton_2, SIGNAL("clicked()"), self.slotSearchFile)
        self.connect(self.lineEdit_3, SIGNAL("textChanged(const QString &)"), self.slotValgrind)
        self.connect(self.comboBox_6, SIGNAL("activated(const QString&)"), self.slotLogType)
        self.connect(self.comboBox_7, SIGNAL("activated(const QString&)"), self.slotLogType)

        # Previous values
        self.valgrind = self.parent.mdl.getString('valgrind')
        if self.valgrind is not None:
            self.lineEdit_3.setText(str(self.valgrind))

        self.setLogType()


    @pyqtSignature("const QString &")
    def slotValgrind(self, text):
        """
        Input for Valgrind.
        """
        self.valgrind = str(text)


    def setLogType(self):
        """
        Set logging arguments.
        """
        self.log_type = self.parent.mdl.getLogType()
        self.modelCSOUT1.setItem(str_model=self.log_type[0])
        self.modelCSOUT2.setItem(str_model=self.log_type[1])


    @pyqtSignature("const QString &")
    def slotLogType(self, text):
        """
        Input logging options.
        """
        self.log_type = [self.modelCSOUT1.dicoV2M[str(self.comboBox_6.currentText())],
                         self.modelCSOUT2.dicoV2M[str(self.comboBox_7.currentText())]]


    @pyqtSignature("")
    def slotSearchFile(self):
        """
        Choice temporary directory for batch
        """
        file_name = ""

        title = self.tr("Select file for use VALGRIND option")
        path  = os.getcwd()
        filetypes = self.tr("All Files (*)")
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes)
        file_name = str(file_name)

        # TO CHECK ...
        if file_name:
            self.valgrind = str(self.lineEdit_3.text())
            if not self.valgrind:
                new = file_name + " --tool=memcheck"
            else:
                new = ""
                for i in self.valgrind.split():
                    if i == self.valgrind.split()[0]:
                        i = file_name
                        new = new + i + ' '
            self.valgrind = new
            self.lineEdit_3.setText(str(self.valgrind))


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """

        self.parent.mdl.setString('valgrind', self.valgrind.strip())

        self.parent.mdl.setLogType(self.log_type)

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text

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

        self.connect(self.lineEditStopIter,
                     SIGNAL("textChanged(const QString &)"),
                     self.__slotStopIter)


    @pyqtSignature("const QString &")
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

        self.connect(self.pushButtonStop,   SIGNAL('clicked()'), self.__slotStop)
        self.connect(self.pushButtonStopAt, SIGNAL('clicked()'), self.__slotStopAt)

        self.scratch_dir = ""
        self.result_dir = ""
        self.suffix   = ""
        self.listing  = "listing"
        self.n_lines = 0

        self.proc.start(self.cmd)


    @pyqtSignature("")
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


    @pyqtSignature("")
    def __slotStop(self):
        """
        Private slot. Stops the code at the end of the current iteration.
        """
        iter = 1
        msg = self.tr("Stop at the end of the current iteration.")
        self.__stopExec(iter, msg)


    @pyqtSignature("")
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


    @pyqtSignature("QProcess::ProcessState")
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
    This class is devoted to the Computer selection.
    When a new computer is selected, The old data frame is deleted and
    a new apropriate frame is open.
    If the batch script file name is known, informations are display
    in the apropiate widget.
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

        # Get batch type

        config = configparser.ConfigParser()
        config.read(self.case['package'].get_configfiles())

        cs_batch_type = None
        if config.has_option('install', 'batch'):
            cs_batch_type = config.get('install', 'batch')
            if os.path.isabs(cs_batch_type):
                i = cs_batch_type.rfind(".")
                if i > -1:
                    cs_batch_type = cs_batch_type[i+1:]

        self.case['batch_type'] = cs_batch_type

        # Get MPI and OpenMP features

        self.have_mpi = False
        self.have_openmp = False
        config_features = self.case['package'].config.features
        if config.has_option('install', 'compute_versions'):
            compute_versions = config.get('install', 'compute_versions').split(':')
            if compute_versions[0]:
                pkg_compute = self.case['package'].get_alternate_version(compute_versions[0])
                config_features = pkg_compute.config.features
        if config_features['mpi'] == 'yes':
            self.have_mpi = True
        if config_features['openmp'] == 'yes':
            self.have_openmp = True

        self.jmdl = BatchRunningModel(parent, self.case)

        # Batch info

        self.hideBatchInfo()

        self.labelNProcs.hide()
        self.spinBoxNProcs.hide()

        self.labelNThreads.hide()
        self.spinBoxNThreads.hide()

        self.class_list = None

        if self.case['batch_type'] != None:

            self.groupBoxArchi.setTitle("Job and script files")
            self.labelBatch.show()
            self.toolButtonSearchBatch.show()

            validatorSimpleName = RegExpValidator(self.lineEditJobName,
                                                  QRegExp("[_A-Za-z0-9]*"))
            validatorAccountName = RegExpValidator(self.lineEditJobAccount,
                                                   QRegExp("\\S+"))
            self.lineEditJobName.setValidator(validatorSimpleName)
            self.lineEditJobAccount.setValidator(validatorAccountName)
            self.lineEditJobWCKey.setValidator(validatorAccountName)
            self.pushButtonRunSubmit.setText("Submit job")

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


        # Connections

        if self.case['batch_type'] != None:
            self.connect(self.lineEditJobName, SIGNAL("textChanged(const QString &)"),
                         self.slotJobName)
            self.connect(self.spinBoxNodes, SIGNAL("valueChanged(int)"),
                         self.slotJobNodes)
            self.connect(self.spinBoxPpn, SIGNAL("valueChanged(int)"),
                         self.slotJobPpn)
            self.connect(self.spinBoxProcs, SIGNAL("valueChanged(int)"),
                         self.slotJobProcs)
            self.connect(self.spinBoxThreads, SIGNAL("valueChanged(int)"),
                         self.slotJobThreads)
            self.connect(self.spinBoxDays, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxHours, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxMinutes, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.spinBoxSeconds, SIGNAL("valueChanged(int)"),
                         self.slotJobWallTime)
            self.connect(self.comboBoxClass, SIGNAL("activated(const QString&)"),
                         self.slotClass)
            self.connect(self.lineEditJobAccount, SIGNAL("textChanged(const QString &)"),
                         self.slotJobAccount)
            self.connect(self.lineEditJobWCKey, SIGNAL("textChanged(const QString &)"),
                         self.slotJobWCKey)

        else:
            self.connect(self.spinBoxNProcs, SIGNAL("valueChanged(int)"), self.slotNProcs)
            self.connect(self.spinBoxNThreads, SIGNAL("valueChanged(int)"), self.slotNThreads)

        self.connect(self.toolButtonSearchBatch, SIGNAL("clicked()"), self.slotSearchBatchFile)
        self.connect(self.comboBoxRunType, SIGNAL("activated(const QString&)"), self.slotArgRunType)
        self.connect(self.toolButtonAdvanced, SIGNAL("clicked()"), self.slotAdvancedOptions)
        self.connect(self.pushButtonRunSubmit, SIGNAL("clicked()"), self.slotBatchRunning)

        # Combomodels

        self.modelArg_cs_verif = ComboModel(self.comboBoxRunType, 2, 1)

        self.modelArg_cs_verif.addItem(self.tr("Import mesh only"), 'none')
        self.modelArg_cs_verif.addItem(self.tr("Mesh preprocessing"), 'mesh preprocess')
        self.modelArg_cs_verif.addItem(self.tr("Mesh quality criteria"), 'mesh quality')
        self.modelArg_cs_verif.addItem(self.tr("Standard"), 'standard')
        self.modelArg_cs_verif.setItem(str_model=self.mdl.getRunType())

        # initialize Widgets

        # Check if the script file name is already defined

        if self.case['runcase']:
            name = os.path.basename(self.case['runcase'].path)
            self.labelBatchName.setText(str(name))
            setGreenColor(self.toolButtonSearchBatch, False)
        else:
            setGreenColor(self.toolButtonSearchBatch, True)

        if self.case['batch_type'] != None and self.case['runcase']:
            self.displayBatchInfo()

        # Script info is based on the XML model

        self.displayScriptInfo()

        self.case.undoStartGlobal()


    @pyqtSignature("const QString &")
    def slotJobName(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobName.validator().state == QValidator.Acceptable:
            self.jmdl.dictValues['job_name'] = str(v)
            self.jmdl.updateBatchFile('job_name')


    @pyqtSignature("int")
    def slotJobNodes(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.jmdl.dictValues['job_nodes'] = str(self.spinBoxNodes.text())
        self.jmdl.updateBatchFile('job_nodes')


    @pyqtSignature("int")
    def slotJobPpn(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.jmdl.dictValues['job_ppn']  = str(self.spinBoxPpn.text())
        self.jmdl.updateBatchFile('job_ppn')


    @pyqtSignature("int")
    def slotJobProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.jmdl.dictValues['job_procs']  = str(self.spinBoxProcs.text())
        self.jmdl.updateBatchFile('job_procs')


    @pyqtSignature("int")
    def slotJobThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        self.jmdl.dictValues['job_threads']  = str(self.spinBoxThreads.text())
        self.jmdl.updateBatchFile('job_threads')


    @pyqtSignature("")
    def slotJobWallTime(self):

        h_cput = self.spinBoxDays.value()*24 + self.spinBoxHours.value()
        m_cput = self.spinBoxMinutes.value()
        s_cput = self.spinBoxSeconds.value()
        self.jmdl.dictValues['job_walltime'] = h_cput*3600 + m_cput*60 + s_cput
        self.jmdl.updateBatchFile('job_walltime')


    @pyqtSignature("")
    def slotClass(self):

        self.jmdl.dictValues['job_class'] = str(self.comboBoxClass.currentText())
        if len(self.jmdl.dictValues['job_class']) > 0:
            self.jmdl.updateBatchFile('job_class')


    @pyqtSignature("const QString &")
    def slotJobAccount(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobAccount.validator().state == QValidator.Acceptable:
            self.jmdl.dictValues['job_account'] = str(v)
            self.jmdl.updateBatchFile('job_account')


    @pyqtSignature("const QString &")
    def slotJobWCKey(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if self.lineEditJobWCKey.validator().state == QValidator.Acceptable:
            self.jmdl.dictValues['job_wckey'] = str(v)
            self.jmdl.updateBatchFile('job_wckey')


    @pyqtSignature("const QString &")
    def slotArgRunType(self, text):
        """
        Input run type option.
        """
        self.run_type = self.modelArg_cs_verif.dicoV2M[str(text)]
        self.mdl.setRunType(self.run_type)


    @pyqtSignature("int")
    def slotNProcs(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.jmdl.dictValues['run_nprocs'] = str(v)
        else:
            self.jmdl.dictValues['run_nprocs'] = None
        self.jmdl.updateBatchFile('run_nprocs')


    @pyqtSignature("int")
    def slotNThreads(self, v):
        """
        Increment, decrement and colorize the input argument entry
        """
        if v > 1:
            self.jmdl.dictValues['run_nthreads'] = str(v)
        else:
            self.jmdl.dictValues['run_nthreads'] = None
        self.jmdl.updateBatchFile('run_nthreads')


    @pyqtSignature("")
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        log.debug("slotAdvancedOptions")

        dialog = BatchRunningAdvancedOptionsDialogView(self)

        if dialog.exec_():
            log.debug("slotAdvancedOptions validated")


    @pyqtSignature("")
    def slotBatchRunning(self):
        """
        Launch Code_Saturne batch running.
        """
        # Is the file saved?

        if self.case['new'] == "yes" or len(self.case['undo']) > 0 or len(self.case['redo']) > 0:

            title = self.tr("Warning")
            msg   = self.tr("The current case must be saved before "\
                            "running the ") + self.tr(self.case['package']).code_name + self.tr(" script.")
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

        bd = LocalizationModel('BoundaryZone', self.case)
        if not bd.getZones():
            if self.case['no_boundary_conditions'] == False:
                title = self.tr("Warning")
                msg   = self.tr("No boundary definition declared.\n\n")
                QMessageBox.warning(self, title, msg)
                self.case['no_boundary_conditions'] = True

        # Build command line

        key = self.case['batch_type']

        batch = cs_exec_environment.enquote_arg(self.case['runcase'].path)

        if key == None:
            run_id, run_title = self.__suggest_run_id()
            self.__updateRuncase(run_id)
            cmd = batch
            key = 'localhost'
        elif key[0:3] == 'CCC':
            cmd = 'msub ' + batch
        elif key[0:5] == 'LOADL':
            cmd = 'llsubmit ' + batch
        elif key[0:3] == 'LSF':
            cmd = 'bsub < ' + batch
        elif key[0:3] == 'PBS' or key[0:3] == 'SGE':
            cmd = 'qsub ' + batch
        elif key[0:5] == 'SLURM':
            cmd = 'sbatch ' + batch
        else:
            pass

        if self.case['salome'] or key == 'localhost':
            dlg = ListingDialogView(self.parent, self.case, run_title, cmd)
            dlg.show()
        else:
            cs_exec_environment.run_command(cmd)

        if self.case['salome'] or key == 'localhost':
            self.__updateRuncase('')  # remove --id <id> from runcase

        os.chdir(prv_dir)


    def __suggest_run_id(self):
        """
        Return an id.
        """
        cmd = os.path.join(self.case['package'].get_dir('bindir'),
                           self.case['package'].name)
        cmd = cs_exec_environment.enquote_arg(cmd) + " run --suggest-id"
        r_title = subprocess.Popen(cmd,
                                   shell=True,
                                   stdout=subprocess.PIPE).stdout.read()[:-1]
        r_id = os.path.join(self.case['resu_path'], r_title)

        run_id = r_id
        run_title = r_title
        j = 1
        while os.path.isdir(run_id):
            j += 1
            run_id = r_id + "_" + str(j)
            run_title = r_title + "(" + str(j) + ")"

        return os.path.basename(run_id), run_title


    def __updateRuncase(self, run_id):
        """
        Update the command line in the launcher C{runcase}.
        """
        runcase = self.case['runcase']

        runcase.set_run_id(run_id=run_id)
        runcase.save()


    def getCommandOutput(self, cmd):
        """
        Run a command and return it's standard output.
        """
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        lines = []
        while True:
            l = p.stdout.readline()
            lines.append(l.strip())
            if len(l) == 0 and p.poll() != None:
                break
        output = p.communicate()

        if p.returncode == 0:
            return lines


    def getClassList(self):
        """
        Layout of the second part of this page.
        """

        self.class_list = []

        try:

            if self.case['batch_type'][0:3] == 'CCC':
                output = self.getCommandOutput('class')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:5] == 'LOADL':
                output = self.getCommandOutput('llclass')
                ignore = True
                for l in output:
                    if l[0:3] == '---':
                        ignore = not ignore
                    elif ignore == False:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'LSF':
                output = self.getCommandOutput('bqueues')
                ignore = True
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'PBS':
                output = self.getCommandOutput('qstat -q')
                ignore = True
                for l in output:
                    if l[0:3] == '---':
                        ignore = not ignore
                    elif ignore == False:
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:3] == 'SGE':
                output = self.getCommandOutput('qconf -sc')
                for l in output:
                    if l[0:1] != '#':
                        self.class_list.append(l.split(' ')[0])

            elif self.case['batch_type'][0:5] == 'SLURM':
                output = self.getCommandOutput('sinfo -s')
                for l in output[1:]:
                    if len(l) == 0:
                        break
                    else:
                        name = l.split(' ')[0]
                        if name[-1:] == '*':
                            name = name[:-1]
                        self.class_list.append(name)

        except Exception:
            pass


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

        self.job_name  = self.jmdl.dictValues['job_name']
        self.job_nodes = self.jmdl.dictValues['job_nodes']
        self.job_ppn  = self.jmdl.dictValues['job_ppn']
        self.job_procs = self.jmdl.dictValues['job_procs']
        self.job_threads = self.jmdl.dictValues['job_threads']
        self.job_walltime = self.jmdl.dictValues['job_walltime']
        self.job_class  = self.jmdl.dictValues['job_class']
        self.job_account  = self.jmdl.dictValues['job_account']
        self.job_wckey  = self.jmdl.dictValues['job_wckey']

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

        if self.job_threads != None:
            self.labelThreads.show()
            self.spinBoxThreads.setValue(int(self.job_threads))
            self.spinBoxThreads.show()

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
                self.getClassList()
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
            self.jmdl.dictValues['job_class'] = str(self.comboBoxClass.currentText())
            if len(self.jmdl.dictValues['job_class']) > 0:
                self.jmdl.updateBatchFile('job_class')

        if self.job_account != None:
            self.labelJobAccount.show()
            self.lineEditJobAccount.setText(str(self.job_account))
            self.lineEditJobAccount.show()

        if self.job_wckey != None:
            self.labelJobWCKey.show()
            self.lineEditJobWCKey.setText(str(self.job_wckey))
            self.lineEditJobWCKey.show()

        # Show Job management box

        if self.case['batch_type'][0:5] == 'LOADL':
            self.groupBoxJob.setTitle("Load Leveler job parameters")
            self.labelJobAccount.setText(str("Group"))
            self.lineEditJobAccount.setToolTip("To obtain a list of defined groups, run <b><tt>xloadl</tt></b>, then select <i>File -> Build a Job</i>, and check the group names in the <i>Group</i> field")
        elif self.case['batch_type'][0:3] == 'LSF':
            self.groupBoxJob.setTitle("LSF job parameters")
        elif self.case['batch_type'][0:3] == 'PBS':
            self.groupBoxJob.setTitle("PBS job parameters")
        elif self.case['batch_type'][0:3] == 'SGE':
            self.groupBoxJob.setTitle("Sun Grid Engine job parameters")
        if self.case['batch_type'][0:5] == 'SLURM':
            self.groupBoxJob.setTitle("SLURM job parameters")
        else:
            self.groupBoxJob.setTitle("Batch job parameters")

        self.groupBoxJob.show()

        # Update file

        self.jmdl.updateBatchFile()


    def displayScriptInfo(self):
        """
        Layout of the second part of this page.
        """

        if self.case['batch_type'] == None:
            if self.have_mpi:
                self.labelNProcs.show()
                self.spinBoxNProcs.show()
            if self.have_openmp:
                self.labelNThreads.show()
                self.spinBoxNThreads.show()
        else:
            self.labelNProcs.hide()
            self.spinBoxNProcs.hide()
            self.labelNThreads.hide()
            self.spinBoxNThreads.hide()

        if self.case['batch_type'] == None:
            if self.have_mpi:
                n_procs_s = self.jmdl.dictValues['run_nprocs']
                if n_procs_s:
                    n_procs = int(n_procs_s)
                else:
                    n_procs = 1
                self.spinBoxNProcs.setValue(n_procs)
            if self.have_openmp == 'yes':
                n_threads_s = self.jmdl.dictValues['run_nthreads']
                if n_threads_s:
                    n_threads = int(n_threads_s)
                else:
                    n_threads = 1
                self.spinBoxNThreads.setValue(n_threads)
        else:
            pass


    @pyqtSignature("")
    def slotSearchBatchFile(self):
        """
        Open a FileDialog in order to search the batch command file
        in the system file.
        """
        file_name = ""
        if self.case['scripts_path'] and os.path.isdir(self.case['scripts_path']):
            path = self.case['scripts_path']
        else:
            path = os.getcwd()
        title = self.tr("Select the batch script")
        filetypes = self.tr("All Files (*)")
        file_name = QFileDialog.getOpenFileName(self, title, path, filetypes)
        file_name = str(file_name)

        if file_name:

            launcher = os.path.basename(file_name)
            setGreenColor(self.toolButtonSearchBatch, False)

            if self.case['scripts_path'] == os.path.dirname(file_name):
                self.case['runcase'] = cs_runcase.runcase(os.path.join(self.case['scripts_path'],
                                                                       launcher),
                                                          package=self.case['package'])
                self.labelBatchName.setText(str(launcher))
                self.hideBatchInfo()
                if self.case['batch_type'] != None:
                    self.displayBatchInfo()
            else:
                title = self.tr("Warning")
                msg   = self.tr("The new batch file is not in scripts "\
                                "directory given in the 'Identity and paths' "\
                                "section.\n\n" + \
                                "Verify the existence and location of these files, "\
                                "and the 'Identity and Pathes' section")
                QMessageBox.warning(self, title, msg)


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
