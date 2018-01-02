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
This module contains the following classes and function:
- BodyForcesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

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
from code_saturne.Base.QtPage import DoubleValidator, from_qvariant
from code_saturne.Pages.BodyForcesForm import Ui_BodyForcesForm
from code_saturne.Pages.BodyForcesModel import BodyForcesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BodyForcesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BodyForcesView(QWidget, Ui_BodyForcesForm):
    """
    Class to open the Body Forces (gravity) Page.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BodyForcesForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.case.undoStopGlobal()

        self.mdl = BodyForcesModel(self.case)

        # Connections

        self.lineEditX.textChanged[str].connect(self.slotGravityX)
        self.lineEditY.textChanged[str].connect(self.slotGravityY)
        self.lineEditZ.textChanged[str].connect(self.slotGravityZ)

        # Validators
        validatorX = DoubleValidator(self.lineEditX)
        validatorY = DoubleValidator(self.lineEditY)
        validatorZ = DoubleValidator(self.lineEditZ)

        self.lineEditX.setValidator(validatorX)
        self.lineEditY.setValidator(validatorY)
        self.lineEditZ.setValidator(validatorZ)

        # Initialization

        gravity_x = self.mdl.getGravity(self.mdl.nodes[0])
        gravity_y = self.mdl.getGravity(self.mdl.nodes[1])
        gravity_z = self.mdl.getGravity(self.mdl.nodes[2])

        self.lineEditX.setText(str(gravity_x))
        self.lineEditY.setText(str(gravity_y))
        self.lineEditZ.setText(str(gravity_z))

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotGravityX(self, text):
        """
        Input GX
        """
        if self.lineEditX.validator().state == QValidator.Acceptable:
            gravity_x = from_qvariant(text, float)
            self.mdl.setGravity('gravity_x', gravity_x)


    @pyqtSlot(str)
    def slotGravityY(self, text):
        """
        Input GY
        """
        if self.lineEditY.validator().state == QValidator.Acceptable:
            gravity_y = from_qvariant(text, float)
            self.mdl.setGravity('gravity_y', gravity_y)


    @pyqtSlot(str)
    def slotGravityZ(self, text):
        """
        Input GZ
        """
        if self.lineEditZ.validator().state == QValidator.Acceptable:
            gravity_z = from_qvariant(text, float)
            self.mdl.setGravity('gravity_z', gravity_z)


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
