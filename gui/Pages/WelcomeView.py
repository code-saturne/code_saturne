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
This module defines the welcome page.
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

from code_saturne.Base           import QtGui
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

try:
    from code_saturne.Pages.WelcomeForm import Ui_WelcomeForm
except Exception:
    import os, sys
    sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
    from code_saturne.Pages.WelcomeForm import Ui_WelcomeForm

#-------------------------------------------------------------------------------
# This class defines the welcome page
#-------------------------------------------------------------------------------

class WelcomeView(QWidget, Ui_WelcomeForm):
    """
    Class for the welcome page
    """
    def __init__(self):
        """
        Constructor
        """
        QWidget.__init__(self)
        Ui_WelcomeForm.__init__(self)
        self.setupUi(self)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
