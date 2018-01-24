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
Desktop Manager
===============

Desktop represents a main frame (QMainWindow) of a SALOME application.
It contains a menu bar, tool bars, and central area for GUI controls of
components: Object Browser, Python console, 3D/2D viewers, etc.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

from __future__ import print_function

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from CFDSTUDYGUI_Commons import Trace
from CFDSTUDYGUI_ActionsHandler import CFDSTUDYGUI_ActionsHandler
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
        if Trace(): print("CFDSTUDYGUI_DesktopMgr::slotDeleteDsk() ", dsk)
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
            dsk.destroyed.connect(self.slotDeleteDsk)

        return self._ActionHandlerMap[dsk]


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

#-------------------------------------------------------------------------------
