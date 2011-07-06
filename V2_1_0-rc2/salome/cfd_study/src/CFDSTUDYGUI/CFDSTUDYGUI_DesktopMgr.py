# -*- coding: utf-8 -*-
#============================================================================
#
#     This file is part of CFDSTUDY the plug-in for Salome
#     of Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     CFDSTUDY is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     CFDSTUDY is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#============================================================================

"""
Desktop Manager
===============

Desktop represents a main frame (QMainWindow) of a SALOME application.
It contains a menu bar, tool bars, and central area for GUI controls of
components: Object Browser, Python console, 3D/2D viewers, etc.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import QObject, SIGNAL

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from CFDSTUDYGUI_Commons import CFD_Code, Trace, CFD_Saturne, CFD_Neptune, sgPyQt, sg
from CFDSTUDYGUI_ActionsHandler import CFDSTUDYGUI_ActionsHandler
from CFDSTUDYGUI_ProcessMgr import CFDSTUDYGUI_ProcessMgr
from CFDSTUDYGUI_Agents import *

#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class CFDSTUDYGUI_DesktopMgr(QObject):
    """
    Auxilliary class for GUI management amoung opened SALOME studies.
    It helps to destroy objects with corresponding Desktop.
    """
    def __init__(self):
        """
        Constructor.
        """
        QObject.__init__(self, None)
        self._ActionHandlerMap = {}


    def slotDeleteDsk(self):
        """
        Destroys objects with corresponding Desktop.
        """
        dsk = self.sender()
        if Trace(): print "CFDSTUDYGUI_DesktopMgr::slotDeleteDsk() ", dsk
        if dsk in self._ActionHandlerMap:
            del self._ActionHandlerMap[dsk]


    def getActionHandler(self, dsk):
        """
        Returns existing or creates new ActionHandler associated to a dekstop.

        @type dsk: C{QMainWindow}
        @param dsk: main window of a SALOME application
        @return: ActionHandler associated to the current SALOME study.
        @rtype: C{CFDSTUDYGUI_ActionsHandler}
        """
        if not dsk in self._ActionHandlerMap:
            ah = CFDSTUDYGUI_ActionsHandler()
            ah.createActions()
            self._ActionHandlerMap[dsk] = ah
            self.connect(dsk, SIGNAL("destroyed(QObject*)"), self.slotDeleteDsk)

        return self._ActionHandlerMap[dsk]


    def getProcessMgr(self, dsk):
        """
        Returns existing or creates new Process Manager. Usefull for the CFD code ruuning.

        @type dsk: C{QMainWindow}
        @param dsk: main window of a SALOME application
        @return: Process Manager.
        @rtype: C{CFDSTUDYGUI_ProcessMgr}
        """
        ah = self.getActionHandler(dsk)
        return ah.processMgr()


    def setWorkspace(self, dsk, ws):
        """
        Stores a workspace I{ws} to an associated desktop I{dsk}.

        @type dsk: C{QMainWindow}
        @param dsk: main window of a SALOME application.
        @type ws: C{QWidget}
        @param ws: workspace.
        """
        ah = self.getActionHandler(dsk)
        ah.dskAgent().setWorkspace(ws)
        #updates in ProcessMgr automatically
