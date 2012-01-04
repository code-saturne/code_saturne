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
Process Manager
===============

A process represent a run of the CFD code. The Process Manager is able to manage
several runs at the same time. All out put tho the console is redirected to the
SALOME Log Window. When the listing of the CFD code is redirected to the
standard output, it is display in a new specific View in thz SALME workspace.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import thread

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtGui import QApplication, QMessageBox, QTextEdit, QFont
from PyQt4.QtCore import QObject, QEvent, QString

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_Commons
import CFDSTUDYGUI_DataModel
import CFDSTUDYGUI_SolverGUI
from CFDSTUDYGUI_Commons import Trace, sgPyQt, LoggingMgr
from CFDSTUDYGUI_Commons import CaseInProcessStart, CaseInProcessEnd, UpdateScriptFolder

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

#Global variables
CurrentProcID = 1

#Process events
PE_Progress   = -1
PE_Finished   = -2
PE_Prepare    = -3

#Process statuses
PS_INIT                = 0
PS_RUN                 = 1
PS_STOP                = 2
PS_EXCEPTION           = 3

#Process results
PR_NotFinished         = 0
PR_Successed           = 1
PR_SuccessedWithErrors = 2
PR_InterruptedByUser   = 3
PR_NotSuccessed        = 4

#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_ProcessMgr(QObject):
    """
    Class for management of solver process
    _ProcessMap structure:
        [proc_id] -> [process temp directory, process status, process result, case path]
    _ProcessLogMap structure:
        [proc_id] -> [log lines count, log widget, log_file, log_message]
    """
    def __init__(self, dskAgent):
        """
        Constructor.
        """
        QObject.__init__(self)
        self._LoggingMgr = LoggingMgr()
        self._ProcessMap = {}
        self._ProcessLogMap = {}
        self._Lock = thread.allocate_lock()
        self._DskAgent = dskAgent


    def addProcess(self, command, case_path, xml_file, tmp_prefix, arg_output):
        """
        Adds command for execution of new process.
        This method is called by the CFD solver GUI, when the user starts his calculus.
        """
        from time import gmtime, strftime
        import os
        import os.path
        import re
        global CurrentProcID

        #tmp dir suffix
        tmp_suffix = strftime("%m%d%H%M")

        proc_dir = tmp_prefix + '.' + tmp_suffix
        case_name = os.path.basename(case_path)

        #init a status map
        self._setProcessInfo(CurrentProcID, proc_dir, case_path)

        # start the CFD code
        thread.start_new_thread(self._runProcess, (command, CurrentProcID))

        # start the monitoring of the running calculus
        thread.start_new_thread(self._monitorProcess, (case_path, proc_dir, CurrentProcID, tmp_suffix))

        info_tab = None
        #check ARG_CS_OUTPUT
        if re.search('--log 0', arg_output): #standard output
            info_tab = CFDSTUDYGUI_ProcessTab(self._DskAgent.workspace())

        #init log status
        self._setLogInfo(CurrentProcID, info_tab)

        CurrentProcID += 1 #increment the future process id

        event = QEvent(CaseInProcessStart)
        event.setData(case_path)
        QApplication.postEvent(self, event)


    def _runProcess(self, command, proc_id):
        """
        Starts the CFD code launcher.
        """
        import subprocess
        import os
        import re
        import time

        self._setProcessStatus(proc_id, PS_RUN)

        try:
            pipe = subprocess.Popen(command, bufsize = 0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            logCount = 1
            while True:
                text = pipe.stdout.readline()
                if not text:
                    break
                log_widget = self._getLogWidget(proc_id)
                if log_widget:
                    event = QEvent(PE_Progress)
                    event.setData(text)
                    QApplication.postEvent(log_widget, event)
                else:
                    self._appendLogMsg(proc_id, text)

                #dublicate msg to the Message Window
                #remove '\n' symbol from end of the line
                #if re.search("\n$", text):
                #    text = text[:len(text)-1]
                #if re.search("Execution\s+$", text):
                #    sgPyQt.message(text, False)
                #    sgPyQt.message("  ********************************************", False)
                #    logCount = 0
                #if re.search("Fin normale du calcul", text) or re.search("ERREUR", text):
                #    sgPyQt.message("  ********************************************", False)
                #    logCount = 1
                #if logCount:
                #    sgPyQt.message(text, False)

            #check for compilation errors
            comp_log = os.path.join(self._getProcessInfo(proc_id)[0], 'compil.log')
            if os.path.exists(comp_log):
                fd = open(comp_log)
                content = fd.readlines()
                fd.close()

                error = False
                for i in content:
                    if re.match(".* Error .*", i):
                        error = True
                        break

                if error:
                    if Trace(): print "compilation error!!!!"
                    self._setProcessStatus(proc_id, PS_EXCEPTION)
                else:
                    self._setProcessStatus(proc_id, PS_STOP)
            else:
                self._setProcessStatus(proc_id, PS_STOP)
        except OSError, e:
            sgPyQt.message("Exception had occured during script execution " + e.__str__(), True)
            self._setProcessStatus(proc_id, PS_EXCEPTION)

        #waiting that monitor of process allow remove proces from manager
        while True:
            if not self._getProcessResult(proc_id) == PR_NotFinished:
                break
            time.sleep(1)

        #update obj browser
        anEvent = QEvent(UpdateScriptFolder)
        anEvent.setData(self._getLogInfo(proc_id)[2])
        QApplication.postEvent(self, anEvent)

        #remove process information from the manager
        self._Lock.acquire()
        del self._ProcessMap[proc_id]
        del self._ProcessLogMap[proc_id]
        self._Lock.release()


    def _monitorProcess(self, case_path, tmp_path, proc_id, tmp_suffix):
        """
        Starts the monitoring of the running calculus.
        """
        import os
        import os.path
        import re
        import time

        script_path = os.path.join(case_path, 'SCRIPTS')
        case_name = os.path.basename(case_path)

        log_file = None

        names = ["runningext." + tmp_suffix, "runningstd." + tmp_suffix]

        # searching output file
        while log_file == None:
            lst = os.listdir(script_path)
            for f in lst:
                if f in names:
                    #file is found
                    if Trace(): print "file is found", f
                    log_file = os.path.join(script_path, f)
                    break
            time.sleep(0.5)

        self._setLogFile(proc_id, log_file)
        log_widget = self._getLogWidget(proc_id)
        if log_widget:
            log_widget.setWindowTitle(os.path.basename(case_path) + "::" + os.path.basename(log_file))
            log_widget.show()

        #update obj browser
        anEvent = QEvent(UpdateScriptFolder)
        anEvent.setData(self._getLogInfo(proc_id)[2])
        QApplication.postEvent(self, anEvent)

        #scanning log file
        while True:
            try:
                fd = open(log_file)
                content = fd.readlines()
                fd.close()
                #back order for optimize detection
                lines_count = len(content)
                log_lines_count = self._getLogCount(proc_id)

                new_log = ""

                log_widget = self._getLogWidget(proc_id)
                if log_widget and lines_count > log_lines_count:
                    for i in range (log_lines_count, lines_count):
                        new_log += content[i]

                    # log needs to update
                    if log_widget:
                        event = QEvent(PE_Progress)
                        event.setData(new_log)
                        QApplication.postEvent(log_widget, event)
                    else:
                        self._appendLogMsg(proc_id, new_log)

                    #dublicate msg to the Message Window
                    #remove '\n' symbol from end of the line
                    #if re.search("\n$", new_log):
                    #    new_log = new_log[:len(new_log)-1]
                    #sgPyQt.message(new_log, False)

                    self._setLogCount(proc_id, lines_count)

            except IOError:
                if Trace(): print "no such file"

            #remove process from the active process
            status = self._getProcessStatus(proc_id)
            if status == PS_STOP or status == PS_EXCEPTION:
                if Trace(): print "process finished!!!"
                break

            #sleepping few seconds
            time.sleep(2)

        #checking process status
        event = QEvent(PE_Finished)

        #if status successiful
        if status == PS_STOP:
            #searching "stop file"
            if os.path.exists(os.path.join(tmp_path, 'ficstp.mod')):
                event.setData([proc_id, PR_InterruptedByUser])
            else:
                #searching error file
                if os.path.exists(os.path.join(tmp_path, 'erreur')) or \
                       os.path.exists(os.path.join(tmp_path, 'error')):
                    event.setData([proc_id, PR_SuccessedWithErrors])
                else:
                    #successed without errors
                    event.setData([proc_id, PR_Successed])
        elif status == PS_EXCEPTION:
            #system exception during execution
            event.setData([proc_id, PR_NotSuccessed])

        #post finish event to process widget
        if log_widget:
            QApplication.postEvent(log_widget, event)

        #post finish event to itself for ActionHandler
        anEvent = QEvent(CaseInProcessEnd)
        anEvent.setData(self._getProcessInfo(proc_id)[3])
        QApplication.postEvent(self, anEvent)

        #needs to update process result -> process can be removed
        self._setProcessResult(proc_id, event.data()[1])


    def _setProcessStatus(self, proc_id, status):
        """
        Locks process map and update status of process
        """
        self._Lock.acquire()

        if proc_id in self._ProcessMap:
            self._ProcessMap[proc_id][1] = status
        self._Lock.release()


    def _getProcessStatus(self, proc_id):
        """
        Locks process map and get status of process
        """
        status = None
        self._Lock.acquire()
        status = self._ProcessMap[proc_id][1]
        self._Lock.release()

        return status


    def _setProcessResult(self, proc_id, result):
        """
        Locks process map and update result of process
        """
        self._Lock.acquire()

        if proc_id in self._ProcessMap:
            self._ProcessMap[proc_id][2] = result
        self._Lock.release()


    def _getProcessResult(self, proc_id):
        """
        Locks process map and get result of process
        """
        self._Lock.acquire()
        result = self._ProcessMap[proc_id][2]
        self._Lock.release()

        return result


    def _setProcessInfo(self, proc_id, proc_dir, case_path):
        """
        Locks process map and set information about process
        """
        self._Lock.acquire()

        if not proc_id in self._ProcessMap:
            self._ProcessMap[proc_id] = [proc_dir, PS_INIT, PR_NotFinished, case_path]
        self._Lock.release()


    def _getProcessInfo(self, proc_id):
        """
        Locks process map and get information about process
        """
        self._Lock.acquire()

        res = ""

        if proc_id in self._ProcessMap:
            res = self._ProcessMap[proc_id]
        self._Lock.release()

        return res


    def _setLogInfo(self, proc_id, log_widget = None):
        """
        Locks process log map and get information about log of process
        """
        self._Lock.acquire()
        self._ProcessLogMap[proc_id] = [0, log_widget, "", ""]
        self._Lock.release()


    def _getLogInfo(self, proc_id):
        """
        Locks process log map and get information about log of process
        """
        res = None

        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            res = self._ProcessLogMap[proc_id]
        self._Lock.release()

        return res


    def _getLogWidget(self, proc_id):
        """
        Locks process log map and gets widget for log text
        """
        res = None

        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            res = self._ProcessLogMap[proc_id][1]
        self._Lock.release()

        return res


    def _setLogWidget(self, proc_id, log_win):
        """
        Locks process log map and sets widget for log text
        """
        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            self._ProcessLogMap[proc_id][1] = log_win
        self._Lock.release()


    def _getLogCount(self, proc_id):
        """
        Locks process log map and get line count in log text
        """
        res = None

        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            res = self._ProcessLogMap[proc_id][0]
        self._Lock.release()

        return res


    def _setLogCount(self, proc_id, count):
        """
        Locks process log map and get information about log of process
        """
        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            self._ProcessLogMap[proc_id][0] = count
        self._Lock.release()


    def _setLogFile(self, proc_id, log_file):
        """
        Sets log file path
        """
        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            self._ProcessLogMap[proc_id][2] = log_file
        self._Lock.release()


    def _getLogMsg(self, proc_id):
        """
        Locks process log map and get log message before activation of view tab
        """
        res = None

        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            res = self._ProcessLogMap[proc_id][3]
        self._Lock.release()

        return res


    def _appendLogMsg(self, proc_id, text):
        """
        Locks process log map and append text to the log message before activation of view tab
        """
        self._Lock.acquire()
        if proc_id in self._ProcessLogMap:
            self._ProcessLogMap[proc_id][3] += text
        self._Lock.release()


    def _findProcId(self, logObject):
        """
        Returns process id for log file object, if solver running
        """
        if not logObject: return False

        log_path = CFDSTUDYGUI_DataModel._GetPath(logObject)

        proc_id = None
        #find process id
        self._Lock.acquire()
        for index in self._ProcessLogMap:
            if self._ProcessLogMap[index][2] == log_path:
                proc_id = index
                break
        self._Lock.release()

        return proc_id


    def isActiveProcess(self, logObject):
        """
        Returns true if solver for log object running, else - false
        """
        proc_id = self._findProcId(logObject)
        return proc_id != None


    def stopCurrentProcess(self, logObject):
        """
        Stops solver by log file object
        """
        proc_id = self._findProcId(logObject)

        if proc_id:
            self._stopProcess(proc_id)
        else:
            if Trace(): print "Can't find process identificator"


    def showCurrentProcess(self, logObject):
        """
        Creates log window solver by log file object
        """
        import os.path

        proc_id = self._findProcId(logObject)
        if proc_id:
            info_tab = self._getLogWidget(proc_id)

            if info_tab == None:
                #create log widget
                info_tab = CFDSTUDYGUI_ProcessTab(self._DskAgent.workspace())
                case_name = os.path.basename(self._getProcessInfo(proc_id)[3])
                info_tab.setWindowTitle(case_name + "::" + logObject.GetName())
                info_tab.show()
                info_tab.append(self._getLogMsg(proc_id))
                self._setLogWidget(proc_id, info_tab)
            else:
                QApplication.postEvent(info_tab, QEvent(QEvent.FocusIn))
        else:
            if Trace(): print "Can't find process identificator"


    def _stopProcess(self, proc_id):
        """
        Stops the process with given process id
        """
        import os.path
        proc_dir = self._getProcessInfo(proc_id)[0]
        fd = open(os.path.join(proc_dir, 'ficstp'), 'w')
        fd.write("\n1")
        fd.close()


    def event(self, e):
        """
        Describes new events I{e} added to the C{QApplication}. These eents allox to
        refresh the Object Browser, during the run of the CFD code.
        """
        if e.type() == CaseInProcessStart:
            CFDSTUDYGUI_DataModel.setCaseInProcess(e.data(), True)
            CFDSTUDYGUI_SolverGUI.updateObjectBrowser()
            return True

        elif e.type() == CaseInProcessEnd:
            # update the icon of th case in the Object Browser (case path in: e.data())
            CFDSTUDYGUI_DataModel.setCaseInProcess(e.data(), False)

            # update the view in the Object Browser with the new folder results
            case = CFDSTUDYGUI_DataModel.findMaxDeepObject(e.data())
            #case = CFDSTUDYGUI_DataModel.GetCase(obj)
            if case:
                lst = CFDSTUDYGUI_DataModel.ScanChildren(case, "RESU")
                if len(lst) == 1:
                    # update RESU folder
                    CFDSTUDYGUI_DataModel._RebuildTreeRecursively(lst[0])
                else:
                    raise ValueError, "Invalid case folders"

            CFDSTUDYGUI_SolverGUI.updateObjectBrowser()
            return True

        elif e.type() == UpdateScriptFolder:
            # get SCRIPTS folder (e.data() contens the path of the file runningstd.MMDDHHMM or runningext.MMDDHHMM)
            obj = CFDSTUDYGUI_DataModel.findMaxDeepObject(e.data())
            if obj:
                # update SCRIPTS folder
                CFDSTUDYGUI_DataModel._RebuildTreeRecursively(obj)
                CFDSTUDYGUI_SolverGUI.updateObjectBrowser()
            else:
                raise ValueError,  "Invalid object path"

            return True

        return QObject.event(self, e)


class CFDSTUDYGUI_ProcessTab(QTextEdit):
    """
    Defines a new Tab window C{QTextEdit} in order to display the listing of the
    CFD code, when the listing is redirected to the standard output.
    """
    def __init__(self, parent):
        """
        Constructor.
        """
        QTextEdit.__init__(self, parent)
        self.setCurrentFont(QFont(QString("Courier"), 10))


    def event(self, e):
        """
        Describes new events I{e} added to the C{QApplication}. It is allow to
        display the right C{QMessageBox} after the end of the run of the CFD code.
        """
        if e.type() == PE_Progress:
            #update log
            new_log = e.data()
            self.append(new_log)
            return True

        elif e.type() == PE_Prepare:
            caption = e.data()
            self.setWindowTitle(caption)
            self.show()
            return True

        elif e.type() == PE_Finished:
            #close process widget
            proc_id, result = e.data()
            mess = "Unknown result of process!"

            if result == PR_Successed:
                mess = self.tr("PROCESS_DLG_PROCESS_SUCCESS")
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
            elif result == PR_SuccessedWithErrors:
                mess = self.tr("PROCESS_DLG_PROCESS_ERROR")
                QMessageBox.critical(None, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)
            elif result == PR_InterruptedByUser:
                mess = self.tr("PROCESS_DLG_PROCESS_STOP")
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
            elif result == PR_NotSuccessed:
                mess = self.tr("PROCESS_DLG_PROCESS_NOT_SUCCESS")
                QMessageBox.critical(None, "Error", mess, QMessageBox.Ok, QMessageBox.NoButton)

            return True

        return QTextEdit.event(self, e)
