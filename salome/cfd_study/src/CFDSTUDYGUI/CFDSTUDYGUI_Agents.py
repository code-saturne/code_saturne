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
Desktop Agent
=============

Container for the SALOME workspace. The SALOME workspace is a C{QWidget} that
is the main window's central widget.
"""


#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------
import logging
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_Agents")
log.setLevel(logging.NOTSET)

class Desktop_Agent:
    """
    Container for the SALOME workspace.
    """
    def __init__(self):
        """
        Constructor.
        """
        log.debug("CFDSTUDYGUI_Agents.Desktop_Agent.__init__: ")
        self._WORKSPACE = None


    def setWorkspace(self, ws):
        """
        Stores the SALOME Workspace I{ws} into the Desktop Manager.

        @type ws: C{QWidget}
        @param ws: main window's central widget.
        """
        log.debug("CFDSTUDYGUI_Agents.Desktop_Agent.setWorkspace: ")
        self._WORKSPACE = ws


    def workspace(self):
        """
        Returns the SALOME Workspace I{ws} into the Desktop Manager.

        @return: main window's central widget.
        @rtype: C{QWidget}
        """
        log.debug("CFDSTUDYGUI_Agents.Desktop_Agent.workspace: ")
        return self._WORKSPACE

#-------------------------------------------------------------------------------
