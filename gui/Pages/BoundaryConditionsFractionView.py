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
- BoundaryConditionsFractionView
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

from code_saturne.Pages.BoundaryConditionsFraction import Ui_BoundaryConditionsFraction

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, from_qvariant, ComboModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsFractionView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsFractionView(QWidget, Ui_BoundaryConditionsFraction) :
    """
    Boundary condition for energy
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsFraction.__init__(self)
        self.setupUi(self)

        self.__modelFraction = ComboModel(self.comboBoxFraction, 2, 1)
        self.__modelFraction.addItem(self.tr("Imposed value"), 'dirichlet')
        self.__modelFraction.addItem(self.tr("Automatic value"), 'automatic')

        # Connections
        self.lineEditFraction.textChanged[str].connect(self.__slotFraction)
        self.comboBoxFraction.activated[str].connect(self.__slotChoiceFraction)

        validatorFraction = DoubleValidator(self.lineEditFraction, min = 0., max = 1.)
        validatorFraction.setExclusiveMin(False)
        validatorFraction.setExclusiveMax(False)

        self.lineEditFraction.setValidator(validatorFraction)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None
        self.__currentField = fieldId


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        if self.__boundary.getNature() == "outlet":
            self.comboBoxFraction.show()
            mdl = boundary.getFractionChoice(self.__currentField)
            self.__modelFraction.setItem(str_model=mdl)
            if mdl == "dirichlet":
                self.lineEditFraction.show()
                val = boundary.getFraction(self.__currentField)
                self.lineEditFraction.setText(str(val))
            elif mdl == "automatic":
                self.lineEditFraction.hide()
        else:
            self.comboBoxFraction.hide()
            self.lineEditFraction.show()
            val = boundary.getFraction(self.__currentField)
            self.lineEditFraction.setText(str(val))
        self.show()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotFraction(self, text):
        """
        INPUT fraction value
        """
        if self.lineEditFraction.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setFraction(self.__currentField, value)


    @pyqtSlot(str)
    def __slotChoiceFraction(self, text):
        """
        INPUT choice of method of calculation of fraction for outlet
        """
        fraction_choice = self.__modelFraction.dicoV2M[str(text)]
        self.__boundary.setFractionChoice(self.__currentField, fraction_choice)

        if fraction_choice == "dirichlet":
            self.lineEditFraction.show()
            val = self.__boundary.getFraction(self.__currentField)
            self.lineEditFraction.setText(str(val))
        elif fraction_choice == "automatic":
            self.lineEditFraction.hide()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
