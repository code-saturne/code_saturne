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
This module contains the following classes and function:
- MatisseRangeDescriptionView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from MatisseRangeDescriptionForm import Ui_MatisseRangeDescriptionForm
import Base.QtPage as QtPage

import Pages.MatisseGeomModel as MatisseGeom
from Pages.MatisseRangeDescriptionModel import MatisseRangeDescriptionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseRangeDescriptionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseRangeDescriptionView(QWidget, Ui_MatisseRangeDescriptionForm):
    """
    """
    def __init__(self, parent, case, rangeType):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseRangeDescriptionForm.__init__(self)
        self.setupUi(self)

        self.rangeType = rangeType
        self.case = case

        # Create the Page layout.
        if self.rangeType == "inlet_range":
            self.widgetLine.initWidget(self.case, "inlet_range_line")
            self.widgetHeight.initWidget(self.case, "inlet_range_height")
        if self.rangeType == "outlet_range":
            self.widgetLine.initWidget(self.case, "outlet_range_line")
            self.widgetHeight.initWidget(self.case, "outlet_range_height")


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------