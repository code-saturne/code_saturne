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
- CoriolisSourceTermsView
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

from code_saturne.Pages.CoriolisSourceTermsForm import Ui_CoriolisSourceTermsForm
from code_saturne.Pages.CoriolisSourceTermsModel import CoriolisSourceTermsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CoriolisSourceTermsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class CoriolisSourceTermsView(QWidget, Ui_CoriolisSourceTermsForm):
    """
    Class to open the Coriolis Source Terms Page.
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_CoriolisSourceTermsForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.case.undoStopGlobal()

        self.mdl = CoriolisSourceTermsModel(self.case)

        # Connections

        self.lineEditOMEGAX.textChanged[str].connect(self.slotOmegaX)
        self.lineEditOMEGAY.textChanged[str].connect(self.slotOmegaY)
        self.lineEditOMEGAZ.textChanged[str].connect(self.slotOmegaZ)

        # Validators
        validatorOmegaX = DoubleValidator(self.lineEditOMEGAX)
        validatorOmegaY = DoubleValidator(self.lineEditOMEGAY)
        validatorOmegaZ = DoubleValidator(self.lineEditOMEGAZ)

        self.lineEditOMEGAX.setValidator(validatorOmegaX)
        self.lineEditOMEGAY.setValidator(validatorOmegaY)
        self.lineEditOMEGAZ.setValidator(validatorOmegaZ)

        # Initialization
        omega_x = self.mdl.getOmega(self.mdl.nodes[0])
        omega_y = self.mdl.getOmega(self.mdl.nodes[1])
        omega_z = self.mdl.getOmega(self.mdl.nodes[2])

        self.lineEditOMEGAX.setText(str(omega_x))
        self.lineEditOMEGAY.setText(str(omega_y))
        self.lineEditOMEGAZ.setText(str(omega_z))

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotOmegaX(self, text):
        """
        Input OMEGAX
        """
        if self.lineEditOMEGAX.validator().state == QValidator.Acceptable:
            omega_x = from_qvariant(text, float)
            self.mdl.setOmega('omega_x', omega_x)


    @pyqtSlot(str)
    def slotOmegaY(self, text):
        """
        Input OMEGAY
        """
        if self.lineEditOMEGAY.validator().state == QValidator.Acceptable:
            omega_y = from_qvariant(text, float)
            self.mdl.setOmega('omega_y', omega_y)


    @pyqtSlot(str)
    def slotOmegaZ(self, text):
        """
        Input OmegaZ
        """
        if self.lineEditOMEGAZ.validator().state == QValidator.Acceptable:
            omega_z = from_qvariant(text, float)
            self.mdl.setOmega('omega_z', omega_z)


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
