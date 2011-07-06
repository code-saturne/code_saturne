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

class Desktop_Agent:
    """
    Container for the SALOME workspace.
    """
    def __init__(self):
        """
        Constructor.
        """
        self._WORKSPACE = None


    def setWorkspace(self, ws):
        """
        Stores the SALOME Workspace I{ws} into the Desktop Manager.

        @type ws: C{QWidget}
        @param ws: main window's central widget.
        """
        self._WORKSPACE = ws


    def workspace(self):
        """
        Returns the SALOME Workspace I{ws} into the Desktop Manager.

        @return: main window's central widget.
        @rtype: C{QWidget}
        """
        return self._WORKSPACE


#class SolverProcess_Agent:
#    """
#    Usefull for a CFD code run management.
#    """
#    def __init__(self):
#        """
#        Constructor.
#        """
#        self._Solver = None
#
#
#    def setProcessMgr(self, pm):
#        self._Solver = solver
#
#
#    def isEmpty(self):
#        return self._Solver == None
#
#
#    def stopSolver(self):
#        if self._Solver:
#            self._Solver.Stop()
