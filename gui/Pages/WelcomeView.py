# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
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
#-------------------------------------------------------------------------------

"""
This module defines the welcome page.
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from PyQt4 import QtGui

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from WelcomeForm import Ui_WelcomeForm

#-------------------------------------------------------------------------------
# This class defines the welcome page
#-------------------------------------------------------------------------------


class WelcomeView(QtGui.QWidget, Ui_WelcomeForm):
    """
    Class for the welcome page
    """
    
    def __init__(self):
        """
        Constructor
        """
        QtGui.QWidget.__init__(self)
        Ui_WelcomeForm.__init__(self)
        self.setupUi(self)
        

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End 
#-------------------------------------------------------------------------------
