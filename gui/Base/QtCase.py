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
This module extends a case with undo/redo functionnality.

This module defines the following classes:
- QtCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, unittest, logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.XMLengine import *

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *

#-------------------------------------------------------------------------------
# XML utility functions
#-------------------------------------------------------------------------------

class QtCase(Case, QObject):
    undo_signal = pyqtSignal()

    def __init__(self, package=None, file_name="", studymanager=False):
        """
        Instantiate a new dico and a new xml doc
        """
        Case.__init__(self, package, file_name, studymanager)
        QObject.__init__(self)


    def undoStop(self):
        self.record_local = True


    def undoStart(self):
        self.record_local = False


    def undoStopGlobal(self):
        self.record_global = False


    def undoStartGlobal(self):
        self.record_global = True


    def undoGlobal(self, f, c):
        if self['current_page'] != '' and self.record_local == False and self.record_global == True:
            if sys.version[0] == '2':
                self['dump_python'].append([f.__module__, f.func_name, c])
            else:
                self['dump_python'].append([f.__module__, f.__name__, c])
            if self.xml_prev != self.toString() or self.xml_prev == "":
                # control if function have same arguments
                # last argument is value
                same = True
                if self.record_argument_prev == None:
                    same = False
                elif (len(c) == len(self.record_argument_prev) and len(c) >= 2):
                    for i in range(0, len(c)-1):
                        if c[i] != self.record_argument_prev[i]:
                            same = False

                if same:
                    pass
                else:
                    self['undo'].append([self['current_page'], self.toString(), self['current_index'], self['current_tab']])
                    self.xml_prev = self.toString()
                    self.record_func_prev = None
                    self.record_argument_prev = c
                    self.undo_signal.emit()


    def undo(self, f, c):
        if self['current_page'] != '' and self.record_local == False and self.record_global == True:
            if sys.version[0] == '2':
                self['dump_python'].append([f.__module__, f.func_name, c])
            else:
                self['dump_python'].append([f.__module__, f.__name__, c])
            if self.xml_prev != self.toString():
                # control if function have same arguments
                # last argument is value
                same = True
                if self.record_argument_prev == None:
                    same = False
                elif (len(c) == len(self.record_argument_prev) and len(c) >= 2):
                    for i in range(0, len(c)-1):
                        if c[i] != self.record_argument_prev[i]:
                            same = False

                if self.record_func_prev == f and same:
                    pass
                else:
                    self.record_func_prev = f
                    self.record_argument_prev = c
                    self['undo'].append([self['current_page'], self.toString(), self['current_index'], self['current_tab']])
                    self.xml_prev = self.toString()
                    self.undo_signal.emit()


#-------------------------------------------------------------------------------
# End of QtCase
#-------------------------------------------------------------------------------
