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
Command Manager
===============
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_Commons, CFDSTUDYGUI_SolverGUI
from CFDSTUDYGUI_Commons import sgPyQt, LoggingMgr
from ui_CFDSTUDYGUI_QProcessDialog import Ui_CFDSTUDYGUI_QProcessDialog
import CFDSTUDYGUI_DataModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_CommandMgr")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Classes definitions
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_QProcessDialog(QDialog, Ui_CFDSTUDYGUI_QProcessDialog):
    """
    Advanced dialog.
    """
    def __init__(self, parent, title, cmd_list, obj_directory="", start_directory=""):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_CFDSTUDYGUI_QProcessDialog.__init__(self)
        self.setupUi(self)
        self.setWindowTitle(title)
        self.pushButton.setEnabled(False)

        self.objBr = None
        if obj_directory != None and obj_directory != "":
            self.objBr = obj_directory
        self.start_directory = start_directory
        self.proc = QProcess()

        self.proc.readyReadStandardOutput.connect(self.__readFromStdout)
        self.proc.readyReadStandardError.connect(self.__readFromStderr)

        self.procErrorFlag = False

        self.cmd_list = cmd_list
        self.cmd = self.cmd_list.pop(0)
        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)
        self.__process()


    def __process(self):
        if self.proc.exitStatus() == QProcess.NormalExit and not self.procErrorFlag:
            if self.start_directory != "":
                self.proc.setWorkingDirectory(self.start_directory)
            self.proc.start(self.cmd)
            if self.cmd_list:
                self.cmd = self.cmd_list.pop(0)
                self.proc.finished['int', 'QProcess::ExitStatus'].connect(self.__process)
            else:
                self.proc.finished['int', 'QProcess::ExitStatus'].connect(self.__finished)


    def __readFromStdout(self):
        """
        Private slot to handle the readyReadStandardOutput signal of the process.
        """
        if self.proc is None:
            return
        self.proc.setReadChannel(QProcess.StandardOutput)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            s = (ba.data()).decode("utf-8")[:-1]
            self.logText.append(s)


    def __readFromStderr(self):
        """
        Private slot to handle the readyReadStandardError signal of the process.
        """
        if self.proc is None:
            return
#MP        self.proc.setReadChannel(QProcess.Exception)
        self.proc.setReadChannel(QProcess.StandardError)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            s = (ba.data()).decode("utf-8")[:-1]
            self.logText.append('<font color="red">' + s + '</font>')
            self.procErrorFlag = True


    def __finished(self):
        if self.objBr:
            CFDSTUDYGUI_DataModel.UpdateSubTree(self.objBr)
        QApplication.restoreOverrideCursor()
        self.pushButton.setEnabled(True)


    def event(self, e):
        if e.type() == 9999:
            return QDialog.event(self, QCloseEvent())

        if e.type() == 9998:
            return QDialog.event(self, QCloseEvent())
        if e.type() == QEvent.Close:
            e.ignore()
            return True

        return QDialog.event(self, e)


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def runCommand(cmd, start_directory, prefix, *args):
    """
    Run command cmd and asynchronize put it to LogWindow
    Each string logged with prefix
    """
    import subprocess
    import os
    try:
        pipe = subprocess.Popen(cmd, bufsize = 0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:
            text = pipe.stdout.readline().decode()
            if not text:
                break

            sgPyQt.message( prefix + text, False )

    except OSError as e:
        sgPyQt.message( prefix + "Exception had occured during script execution " + e.__str__(), True )
