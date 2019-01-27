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
- OpenTurnsView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

import os
import subprocess
#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator, from_qvariant
from code_saturne.Pages.OpenTurnsForm import Ui_OpenTurnsForm
from code_saturne.model.OpenTurnsModel import OpenTurnsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OpenTurnsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class OpenTurnsView(QWidget, Ui_OpenTurnsForm):
    """
    OpenTurns Page viewer class
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_OpenTurnsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.mdl = OpenTurnsModel(case)
        if not self.mdl.getHostName():
            self.mdl.setDefaultOptions()

        # Combo model
        config = configparser.ConfigParser()
        config.read(self.case['package'].get_configfiles())
        self.nmodes = 1

        dist_hosts = None
        if config.has_section('distant_hosts'):
            if len(config.options('distant_hosts')) != 0:
                dist_hosts = config.options('distant_hosts')
                self.nmodes += len(dist_hosts)

        self.modelOtStudyHosts = ComboModel(self.comboBoxStudyMode,self.nmodes,1)
        self.modelOtStudyHosts.addItem(self.tr("Localhost"), 'localhost')

        self.hosts_bmgr = {}
        if config.has_option('install', 'batch'):
            self.hosts_bmgr['localhost'] = (config.get('install', 'batch')).lower()
        else:
            self.hosts_bmgr['localhost'] = 'none'

        # Check for distant builds:
        # Hosts are stored in the form <batch_rm>_<host_name> hence the split
        # used hereafter to determine the "real" host name
        self.hosts_binpath =  {}

        self.distant_host_builds = None
        if dist_hosts != None:
            self.distant_host_builds = {}
            for key in dist_hosts:
                host_name = key.split('_')[1]
                self.hosts_bmgr[host_name] = key.split('_')[0]

                self.hosts_binpath[host_name] = config.get('distant_hosts',
                                                           key)

                self.addDistantBuilds(host_name)

                dh_not_found = False
                if self.distant_host_builds[host_name] == None:
                    dh_not_found = True

                host_tag = 'distant : ' + host_name
                self.modelOtStudyHosts.addItem(self.tr(host_tag),
                                               host_name,
                                               warn=dh_not_found)
                if dh_not_found:
                    self.modelOtStudyHosts.disableItem(str_model=host_name)

        # ---------------------------------------
        # Connections:
        self.comboBoxStudyMode.activated[str].connect(self.slotOtStudyMode)
        self.comboBoxDistantBuilds.activated[str].connect(self.slotBuildChoice)

        self.lineEditDistWorkDir.textChanged[str].connect(self.slotOtDistWdir)

        self.spinBoxNumberNodes.valueChanged[int].connect(self.slotUpdateNodesNumber)
        self.spinBoxNumberTasks.valueChanged[int].connect(self.slotUpdateTasksNumber)
        self.spinBoxNumberThreads.valueChanged[int].connect(self.slotUpdateThreadsNumber)

        self.spinBoxNumberDays.valueChanged[int].connect(self.slotUpdateWCDays)
        self.spinBoxNumberHours.valueChanged[int].connect(self.slotUpdateWCHours)
        self.spinBoxNumberMinutes.valueChanged[int].connect(self.slotUpdateWCMinutes)
        self.spinBoxNumberSeconds.valueChanged[int].connect(self.slotUpdateWCSeconds)

        self.pushButtonLaunchOT.clicked.connect(self.slotLaunchCsOt)

        # ---------------------------------------
        # Hide/Show initial elements
        if self.nmodes == 1:
            self.groupBoxLocalLaunch.show()
            self.groupBoxDistantLaunch.hide()
            self.setAvailableBuildsList('localhost')
        else:
            self.setAvailableBuildsList(self.mdl.getHostName())
            if self.mdl.getHostName() == "localhost":
                self.groupBoxLocalLaunch.show()
                self.groupBoxDistantLaunch.hide()
            else:
                self.groupBoxLocalLaunch.hide()
                self.groupBoxDistantLaunch.show()

        # ---------------------------------------
        # Initial values
        self.lineEditOutputFile.setText(self.mdl.resfile_name)
        if dist_hosts != None:
            self.modelOtStudyHosts.setItem(str_model=self.mdl.host_name)
        else:
            self.modelOtStudyHosts.setItem(str_model='localhost')

        self.spinBoxLocalProcs.setValue(int(self.mdl.nprocs))
        self.spinBoxLocalThreads.setValue(1)

        self.spinBoxNumberNodes.setValue(int(self.mdl.nnodes))
        self.spinBoxNumberTasks.setValue(int(self.mdl.ntasks))
        self.spinBoxNumberThreads.setValue(int(self.mdl.nthreads))

        wct = self.mdl.getWallClockTime()
        self.spinBoxNumberDays.setValue(int(wct[0]))
        self.spinBoxNumberHours.setValue(int(wct[1]))
        self.spinBoxNumberMinutes.setValue(int(wct[2]))
        self.spinBoxNumberSeconds.setValue(int(wct[3]))

        self.lineEditWCKEY.setText(self.mdl.wckey)


    @pyqtSlot(str)
    def slotOtStudyMode(self, text):
        """
        Host type: localhost or a distant one (defined in code_saturne.cfg).
        """
        host_name = self.modelOtStudyHosts.dicoV2M[str(text)]

        self.mdl.setHostName(host_name)
        self.mdl.setBatchManager(self.hosts_bmgr[host_name])

        if host_name != 'localhost':
            self.groupBoxDistantLaunch.show()
        else:
            self.groupBoxDistantLaunch.hide()

        self.setAvailableBuildsList(host_name)

        self.mdl.arch_path = self.hosts_binpath[host_name]


    @pyqtSlot(str)
    def slotBuildChoice(self, text):
        """
        Sets the hostname
        """
        self.mdl.setBuildName(text)


    @pyqtSlot(str)
    def slotOtDistWdir(self, text):
        """
        Set the distant workdir path
        """
        self.mdl.setDistWorkdir(text)


    @pyqtSlot(int)
    def slotUpdateNodesNumber(self, v):
        """
        Update the number of required computation nodes
        """

        n = int(self.spinBoxNumberNodes.text())
        self.mdl.setClusterParams(nnodes=n)

    @pyqtSlot(int)
    def slotUpdateNprocs(self, v):
        """
        Update the number of required processes
        """

        n = int(self.spinBoxLocalProcs.text())
        self.mdl.setNprocs(n)

    @pyqtSlot(int)
    def slotUpdateTasksNumber(self, v):
        """
        Update the number of required mpi tasks per node
        """

        n = int(self.spinBoxNumberTasks.text())
        self.mdl.setClusterParams(ntasks=n)


    @pyqtSlot(int)
    def slotUpdateThreadsNumber(self, v):
        """
        Update the number of required threads per processor
        """

        n = int(self.spinBoxNumberThreads.text())
        self.mdl.setClusterParams(nthreads=n)


    @pyqtSlot(int)
    def slotUpdateWCDays(self, v):
        """
        Update the wall clock days value
        """

        d, h, m, s = self.mdl.getWallClockTime()
        d = str(int(self.spinBoxNumberDays.text()))

        self.mdl.setWallClockTime(d, h, m, s)


    @pyqtSlot(int)
    def slotUpdateWCHours(self, v):
        """
        Update the wall clock hours value
        """

        d, h, m, s = self.mdl.getWallClockTime()
        h = str(int(self.spinBoxNumberHours.text()))

        self.mdl.setWallClockTime(d, h, m, s)


    @pyqtSlot(int)
    def slotUpdateWCMinutes(self, v):
        """
        Update the wall clock minutes value
        """

        d, h, m, s = self.mdl.getWallClockTime()
        m = str(int(self.spinBoxNumberMinutes.text()))

        self.mdl.setWallClockTime(d, h, m, s)


    @pyqtSlot(int)
    def slotUpdateWCSeconds(self, v):
        """
        Update the wall clock seconds value
        """

        d, h, m, s = self.mdl.getWallClockTime()
        s = str(int(self.spinBoxNumberSeconds.text()))

        self.mdl.setWallClockTime(d, h, m, s)


    @pyqtSlot()
    def slotLaunchCsOt(self):
        """
        Translate the Code_Sature reference case and study into an OpenTurs
        physical model and study
        """

        if self.case['salome']:
            import salome_ot

            # Update the study cfg file
            self.mdl.update_cfg_file()

            # Generating OpenTurns _exec function
            self.create_cs_exec_function()

            # Load the Code_Saturne cas as Physical model and launch
            # OpenTurns GUI
            cs_exec_script_name = os.path.join(self.mdl.otstudy_path, 'cs_execute_job.py')
            salome_ot.loadYacsPyStudy(cs_exec_script_name)

        else:
            print("This option is only available within the SALOME_CFD platform")


    def tr(self, text):
        """
        Translation
        """
        return text

    def addDistantBuilds(self, host_name):
        """
        Search for the distant builds of Code_Saturne for the given distant
        host.
        """

        host_path = self.hosts_binpath[host_name]

        builds_list = __getListOfDistantBuilds__(host_name,
                                                 host_path)

        self.distant_host_builds[host_name] = builds_list


    def setAvailableBuildsList(self, host_name):
        """
        Set the list of available builds per host
        """

        if host_name == 'localhost':

            self.groupBoxDistantLaunch.hide()
            self.groupBoxLocalLaunch.show()

            self.comboBoxDistantBuilds.hide()
            self.labelDistantBuilds.hide()
            return

        else:
            self.groupBoxDistantLaunch.show()
            self.groupBoxLocalLaunch.hide()

            self.comboBoxDistantBuilds.show()
            self.labelDistantBuilds.show()

            dist_builds_list = self.distant_host_builds[host_name]
            self.modelOtDistantBuilds = ComboModel(self.comboBoxDistantBuilds,
                                                   len(dist_builds_list),
                                                   1)
            for db in dist_builds_list:
                self.modelOtDistantBuilds.addItem(self.tr(db), db)





    def create_cs_exec_function(self):
        """
        This function generates the _exec function needed by OpenTurns for a study
        using distant launching on clusters.
        Takes as input:
            - Code_Saturne study path
            - OT_params.cfg name
            - list of OTurns input variables
            - list of OTurns output variables
            - the requested cluster name
            - the result file name which contains the output values
        """

        cluster = self.mdl.host_name

        exec_file_name = os.path.join(self.mdl.otstudy_path, 'cs_execute_job.py')

        f = open(exec_file_name, 'wa')

        script_cmd = "\n"
        script_cmd += "# ============================================================================== \n"
        script_cmd += "# OPENTURNS EXEC FUNCTION WHICH LAUNCHES CODE_SATURNE ON A DISTANT CLUSTER \n"
        script_cmd += "# ============================================================================== \n"

        script_cmd += "\n\n"


        nvars  = len(self.mdl.input_variables)

        script_cmd = 'def _exec('
        vars_dict = '{'

        toffset    = '    '
        cmd1 = 'cfd_eval = cfd_openturns_study('
        loffset1   = toffset
        for i in range(len(cmd1)):
            loffset1 += ' '

        loffset2   = '           '

        iv = -1
        for i in range(nvars):
            if i == 0:
                script_cmd += self.mdl.input_variables[i]
                vars_dict += '"' + self.mdl.input_variables[i] + '":'
                vars_dict += self.mdl.input_variables[i]
            else:
                script_cmd += ", " + self.mdl.input_variables[i]
                vars_dict += ', \n'
                vars_dict += loffset1 + loffset2 + '   '
                vars_dict += '"' + self.mdl.input_variables[i] + '":'
                vars_dict += self.mdl.input_variables[i]

        script_cmd += '):\n\n'
        vars_dict += '}'

        script_cmd += toffset + "import sys\n"

        salome_pydir = os.path.join(self.case['package'].dirs['pythondir'][1],
                                    'salome')
        script_cmd += toffset
        script_cmd += "sys.path.insert(-1, '%s')\n\n" % salome_pydir
        script_cmd += toffset + "from CFDSTUDYOTURNS_StudyInterface import cfd_openturns_study"
        script_cmd += "\n\n"

        script_cmd += toffset + cmd1 + 'study_path = "' + self.mdl.otstudy_path + '",\n'
        script_cmd += loffset1 + 'study_cfg  = "openturns_study.cfg",\n'
        script_cmd += loffset1 + 'vars_dico  = ' + vars_dict + ')\n\n'

        script_cmd += toffset + 'cfd_eval.study2code() \n\n'
        script_cmd += toffset + 'cfd_eval.run() \n\n'

        vals_list = ''
        for i in range(len(self.mdl.output_variables)):
            if i > 0:
                vals_list += ', '
            vals_list += self.mdl.output_variables[i]

        script_cmd += toffset + vals_list + ' = cfd_eval.code2study()\n\n'

        script_cmd += toffset + 'return '
        script_cmd += vals_list + '\n'

        f.write(script_cmd)
        f.close()

#-------------------------------------------------------------------------------
# Utility functions
#-------------------------------------------------------------------------------

def __getListOfDistantBuilds__(host_name, search_path):
    """
    This functions retrieve the list of Code_Saturne builds in a given
    directory on a distant cluster.
    Returns None if no builds are found.
    """

    # Constructing the search command
    search_cmd = ''
    search_cmd += "cd "+ search_path + "\n"
    search_cmd += "build_list=\n"
    search_cmd += "find . -mindepth 1 -maxdepth 1 -type d |\n"
    search_cmd += "while read -r d; do\n"
    search_cmd += "  echo ${d:2}\n"
    search_cmd += "done\n"
    search_cmd += "'"

    ssh_cmd = subprocess.Popen(["ssh", host_name, search_cmd],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    vr = ssh_cmd.stdout.readlines()
    ve = ssh_cmd.stderr.readlines()

    if vr == [] or vr == None or 'Could not resolve hostname' in ve:
        dist_versions = None
    else:
        vr.sort()
        dist_versions = []
        for iv in range(len(vr)):
            dist_versions.append(str(vr[iv][:-1]))


    return dist_versions

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
