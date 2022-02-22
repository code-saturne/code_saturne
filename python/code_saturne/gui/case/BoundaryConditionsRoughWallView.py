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
- BoundaryConditionsRoughWallView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant

from code_saturne.gui.case.BoundaryConditionsRoughWallForm import Ui_BoundaryConditionsRoughWallForm

from code_saturne.model.TurbulenceModel import TurbulenceModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsRoughWallView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsRoughWallView(QWidget, Ui_BoundaryConditionsRoughWallForm):
    """
    Boundary condition for smooth or rough wall.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsRoughWallForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None

        self.case.undoStopGlobal()

        self.radioButtonSmooth.clicked.connect(self.__slotRoughness)
        self.radioButtonRough.clicked.connect(self.__slotRoughness)

        self.lineEditRoughCoef.textChanged[str].connect(self.__slotRoughnessHeight)

        validatorRoughCoef = DoubleValidator(self.lineEditRoughCoef)
        self.lineEditRoughCoef.setValidator(validatorRoughCoef)

        turb_mdl = TurbulenceModel(self.case)
        if turb_mdl.getWallFunction() in [0, 1, 7]:
            self.radioButtonRough.setEnabled(False)

        self.case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.show()
        self.__boundary = boundary

        if self.__boundary.getRoughnessChoice() == "on":
            self.radioButtonSmooth.setChecked(False)
            self.radioButtonRough.setChecked(True)
        else:
            self.radioButtonSmooth.setChecked(True)
            self.radioButtonRough.setChecked(False)

        self.__slotRoughness()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot()
    def __slotRoughness(self):
        """
        Private slot.

        Selects if the wall is rough or smooth.
        """
        if self.radioButtonSmooth.isChecked():
            self.frameRoughness.hide()
            self.__boundary.setRoughnessChoice('off')

        elif self.radioButtonRough.isChecked():
            self.frameRoughness.show()
            self.__boundary.setRoughnessChoice('on')
            r = self.__boundary.getRoughness()
            self.lineEditRoughCoef.setText(str(r))


    @pyqtSlot(str)
    def __slotRoughnessHeight(self, text):
        """
        Private slot.

        Input the roughness height for the selected wall.

        @type text: C{QString}
        @param text: roughness height.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            r = from_qvariant(text, float)
            self.__boundary.setRoughness(r)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
