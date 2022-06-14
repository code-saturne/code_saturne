# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsPressure import Ui_BoundaryConditionsPressure

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, from_qvariant, ComboModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsPressureView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsPressureView(QWidget, Ui_BoundaryConditionsPressure) :
    """
    Boundary condition for energy
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsPressure.__init__(self)
        self.setupUi(self)

        # Connections
        self.lineEditPressure.textChanged[str].connect(self.__slotPressure)

        validatorPressure = DoubleValidator(self.lineEditPressure)

        self.lineEditPressure.setValidator(validatorPressure)

        self.comboBoxPressureMode.activated[str].connect(self._slotChoicePressure)
        self._modelPressure = ComboModel(self.comboBoxPressureMode, 2, 1);


    def setup(self, case, fieldId = None):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        self.lineEditPressure.show()
        val = boundary.getReferencePressure()
        self.lineEditPressure.setText(str(val))

        self._modelPressure.addItem(self.tr("Mean outlet pressure"), "dpdndtau")
        self._modelPressure.addItem(self.tr("Homogeneous outlet pressure"), "dirichlet")
        pressureChoice = self.__boundary.getPressureChoice()
        self._modelPressure.setItem(str_model=pressureChoice)

        self.show()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotPressure(self, text):
        """
        INPUT fraction value
        """
        if self.lineEditPressure.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setReferencePressure(value)

    @pyqtSlot(str)
    def _slotChoicePressure(self, text):
        """
        Input choice of pressure boundary condition.
        """

        p_choice = self._modelPressure.dicoV2M[str(text)]

        self.__boundary.setPressureChoice(p_choice)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
