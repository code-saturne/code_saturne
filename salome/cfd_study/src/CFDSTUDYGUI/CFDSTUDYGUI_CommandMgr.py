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

from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QApplication, QCursor, QDialog, QLabel, QGridLayout, \
                        QCloseEvent, QTextEdit, QTextOption, QDockWidget, QWidget, QFont
from PyQt4.QtCore import Qt, QObject, QVariant, SIGNAL, QEvent, QProcess, QString

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_Commons, CFDSTUDYGUI_SolverGUI
from CFDSTUDYGUI_Commons import sgPyQt, LoggingMgr
#from CFDSTUDYGUI_SolverGUI import _d_DockWindowsRuncase
from CFDSTUDYGUI_Management import _d_DockWindowsRuncase
from ui_CFDSTUDYGUI_QProcessDialog import Ui_CFDSTUDYGUI_QProcessDialog

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_CommandMgr")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------

StopEventType = -1000

#-------------------------------------------------------------------------------
# Classes definitions
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_CommandMgr(QObject):
    def __init__(self):
        QObject.__init__(self)
        self._LoggingMgr = LoggingMgr()
        self.textEdit = CFDSTUDYGUI_MyWidgetListing()


    def runFunctionDlg(self, function, Message, RedirectOut, **kwargs):
        """
        Calls custom function with list of arguments in background mode.
        @type function: C{function}
        @param function: popup menu from the Object Browser.
        @type Message: C{String}
        @param Message:
        @type RedirectOut: C{True} or C{False}
        @param RedirectOut:
        @type kwargs: C{Dictionary}
        @param kwargs:
        """
        self._dlg = CFDSTUDYGUI_CommandDlg(Message)
        self.kw = kwargs

        import thread
        thread.start_new_thread( self._runFunction,\
                                ( self._dlg, function, RedirectOut ) )

        self._dlg.show()


    def _runFunction( self, dlg, function, RedirectOut ):

        if RedirectOut:
            import sys
            self._LoggingMgr.start( sys )

        function(**self.kw)
        if RedirectOut:
            self._LoggingMgr.finish(sys)
        #end of operation -> close the dialog
            QApplication.postEvent(dlg, QCloseEvent())

        QApplication.postEvent(dlg, QEvent(9999))

    def runCommandDlg(self, sObjRep, Message, cmd, start_directory = "", prefix = ""):
        """
        Executing of custom shell command in background mode.
        All output information catched by LogWindow.
        """
        #self._dlg = CFDSTUDYGUI_CommandDlg(Message)
        self.sObjR = sObjRep
        if "str" in str(type(cmd)):
            # run runcase (we are in a runcase case)
            self.runTextEdit("", cmd)
        import thread
        #thread.start_new_thread( self._runCommand, ( self._dlg, cmd, start_directory, prefix) )
        thread.start_new_thread( self._runCommand, ( cmd, start_directory, prefix) )
        #self._dlg.show()


    def _runCommand(self, cmd, start_directory = '', prefix = '', log_file = ''):
        """
        Run command cmd and asynchronize put it to LogWindow and to log file.
        Each string is logged with prefix.
        """
        import subprocess
        import os
        import re

        if start_directory != None and start_directory != "":
            os.chdir(start_directory)

        aLogFile = None
        if log_file != '':
            aLogFile = open( log_file, 'w' )

        #try:
        pipe = subprocess.Popen(cmd, bufsize = 0, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, close_fds=True)
        try:
            while True:
                text = pipe.stdout.readline()
                if not text:
                    break
                if aLogFile:
                    aLogFile.write( text )
                if re.search( "\n$", text ):
                    text = text[:len(text)-1]
                ev=QtCore.QEvent(QtCore.QEvent.User)
                ev.text= prefix + text
                QtGui.qApp.postEvent(self.textEdit,ev)
                #sgPyQt.message( prefix + text, False )

        except OSError, e:
            sgPyQt.message( prefix + "Exception had occured during script execution " + e.__str__(), True )
        if aLogFile:
            aLogFile.close()
        if pipe:
            pipe.stdout.close()
            #QApplication.postEvent(dlg, QEvent(9998))
            import CFDSTUDYGUI_DataModel
            CFDSTUDYGUI_DataModel.UpdateSubTree(self.sObjR)


    def runTextEdit(self,texte,aTitleCase):
        """
        """
        import string
        dsk = sgPyQt.getDesktop()
        runDockExist = False
        ldockWindows = dsk.findChildren(QDockWidget)

        for dock in ldockWindows:
            dockTitle = dock.windowTitle()

            if string.strip(str(dock.windowTitle())) == string.strip(aTitleCase):

                runDockExist = True
                [self.textEdit] = dock.findChildren(CFDSTUDYGUI_MyWidgetListing)

                self.textEdit.textEditListing.clear()

                runDock = dock
        if runDockExist == False:
            dock = QDockWidget(aTitleCase)
            dock.setWidget(self.textEdit)

            dsk.addDockWidget(Qt.RightDockWidgetArea,dock)
            runDock = dock
            studyId = sgPyQt.getStudyId()
            if studyId not in _d_DockWindowsRuncase.keys():
                _d_DockWindowsRuncase[studyId] = []
            _d_DockWindowsRuncase[studyId].append(runDock)

        for dockw in ldockWindows:
            titleDock = dockw.windowTitle()
            if "Browser" not in titleDock and "xml" in titleDock:
                dsk.tabifyDockWidget(dockw,runDock)

        runDock.setVisible(True)
        runDock.show()
        runDock.raise_()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_MyWidgetListing(QWidget):
    def __init__(self, parent = None ):
        QWidget.__init__(self, parent)
        self.textEditListing = QTextEdit()
        self.textEditListing.setWordWrapMode(QTextOption.NoWrap)
        self.textEditListing.setCurrentFont(QFont(QString("Courier"), 10))
        layout = QGridLayout()
        layout.addWidget(self.textEditListing,0,0)
        self.setLayout(layout)


    def event(self, ev):
        if ev.type() == QEvent.User:
            self.textEditListing.append(ev.text)
        return QWidget.event(self, ev)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_CommandDlg(QDialog):
    def __init__(self, text,parent = None ):
        QDialog.__init__(self, parent)
        self.setWindowTitle(" ")
        layout = QGridLayout()
        label = QLabel(text, self)
        layout.addWidget(label)
        self.setLayout(layout)
        self.setMaximumSize(self.minimumSize())
        self.setSizeGripEnabled(False)

    def event(self, e):
        if e.type() == 9999:
            return QDialog.event(self, QCloseEvent())

        if e.type() == 9998:
            return QDialog.event(self, QCloseEvent())
        if e.type() == QEvent.Close:
            e.ignore()
            return True

        return QDialog.event(self, e)


    def show(self):
        QDialog.show(self)


    def accept( self ):
        QDialog.accept(self)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_QProcessDialog(QDialog, Ui_CFDSTUDYGUI_QProcessDialog):
    """
    Advanced dialog.
    """
    def __init__(self, parent, title, cmd_list, start_directory=""):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_CFDSTUDYGUI_QProcessDialog.__init__(self)
        self.setupUi(self)
        self.setWindowTitle(title)
        self.pushButton.setEnabled(False)

        if start_directory != None and start_directory != "":
            os.chdir(start_directory)

        self.proc = QProcess()
        self.connect(self.proc, SIGNAL('readyReadStandardOutput()'), self.__readFromStdout)
        self.connect(self.proc, SIGNAL('readyReadStandardError()'),  self.__readFromStderr)
        self.procErrorFlag = False

        self.cmd_list = cmd_list
        self.cmd = self.cmd_list.pop(0)
        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)
        self.__process()


    def __process(self):
        if self.proc.exitStatus() == QProcess.NormalExit and not self.procErrorFlag:
            self.proc.start(self.cmd)
            if self.cmd_list:
                self.cmd = self.cmd_list.pop(0)
                self.connect(self.proc,
                             SIGNAL('finished(int, QProcess::ExitStatus)'),
                             self.__process)
            else:
                self.connect(self.proc,
                             SIGNAL('finished(int, QProcess::ExitStatus)'),
                             self.__finished)


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
            str = QString()
            s = QString(str.fromUtf8(ba.data()))[:-1]
            self.logText.append(s)


    def __readFromStderr(self):
        """
        Private slot to handle the readyReadStandardError signal of the process.
        """
        if self.proc is None:
            return
        self.proc.setReadChannel(QProcess.StandardError)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            str = QString()
            s = QString(str.fromUtf8(ba.data()))[:-1]
            self.logText.append(s.prepend('<font color="red">').append('</font>'))
            self.procErrorFlag = True


    def __finished(self):
        QApplication.restoreOverrideCursor()
        self.pushButton.setEnabled(True)

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

    if start_directory != None and start_directory != "":
        os.chdir(start_directory)

    try:
        pipe = subprocess.Popen(cmd, bufsize = 0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        while True:
            text = pipe.stdout.readline()
            if not text:
                break

            sgPyQt.message( prefix + text, False )

    except OSError, e:
        sgPyQt.message( prefix + "Exception had occured during script execution " + e.__str__(), True )
