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
This module contains the following classes:
- BoundaryConditionsMeteoView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.BoundaryConditionsMeteoForm import Ui_BoundaryConditionsMeteoForm
from Pages.AtmosphericFlowsModel import AtmosphericFlowsModel

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsMeteoView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsMeteoView(QWidget, Ui_BoundaryConditionsMeteoForm):
    """
    Boundary condifition for the velocity part
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsMeteoForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup teh widget
        """
        self.__case = case
        self.__boundary = None
        self.__model = AtmosphericFlowsModel(self.__case)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        if self.__model.getAtmosphericFlowsModel() != "off" \
            and self.__model.getMeteoDataStatus() == "on":
            self.show()
            self.__boundary = boundary
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
