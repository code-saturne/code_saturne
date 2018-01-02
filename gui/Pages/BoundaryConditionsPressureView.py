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
This module contains the following classes:
- BoundaryConditionsPressureView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant

from code_saturne.Pages.BoundaryConditionsPressureForm import \
     Ui_BoundaryConditionsPressureForm
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsPressureView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsPressureView(QWidget, Ui_BoundaryConditionsPressureForm):
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsPressureForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        # Connections
        self.lineEditPressure.textChanged[str].connect(self.slotPressureValue)

        # Validators
        validatorP = DoubleValidator(self.lineEditPressure, min = 0.0)

        # Apply validators
        self.lineEditPressure.setValidator(validatorP)

        self.__case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        label = boundary.getLabel()
        self.nature  = boundary.getNature()
        self.__boundary = Boundary(self.nature, label, self.__case)
        self.initialize()


    def initialize(self):

        # Initialize thermodynamic value
        pressure = self.__boundary.getPressureValue()
        self.lineEditPressure.setText(str(pressure))

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(str)
    def slotPressureValue(self, text):
        """
        INPUT outlet pressure
        """
        if self.sender().validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.__boundary.setPressureValue(t)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
