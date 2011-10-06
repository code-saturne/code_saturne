# -*- coding: iso-8859-1 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from CoriolisSourceTermsForm import Ui_CoriolisSourceTermsForm
import Base.QtPage as QtPage
from Pages.CoriolisSourceTermsModel import CoriolisSourceTermsModel

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
        self.mdl = CoriolisSourceTermsModel(self.case)

        # Connections

        self.connect(self.lineEditOMEGAX, SIGNAL("textChanged(const QString &)"), self.slotOmegaX)
        self.connect(self.lineEditOMEGAY, SIGNAL("textChanged(const QString &)"), self.slotOmegaY)
        self.connect(self.lineEditOMEGAZ, SIGNAL("textChanged(const QString &)"), self.slotOmegaZ)

        # Validators
        validatorOmegaX = QtPage.DoubleValidator(self.lineEditOMEGAX)
        validatorOmegaY = QtPage.DoubleValidator(self.lineEditOMEGAY)
        validatorOmegaZ = QtPage.DoubleValidator(self.lineEditOMEGAZ)

        self.lineEditOMEGAX.setValidator(validatorOmegaX)
        self.lineEditOMEGAY.setValidator(validatorOmegaY)
        self.lineEditOMEGAZ.setValidator(validatorOmegaZ)

        # Initialization
        omega_x = self.mdl.getOmega(self.mdl.nodes[0])
        omega_y = self.mdl.getOmega(self.mdl.nodes[1])
        omega_z = self.mdl.getOmega(self.mdl.nodes[2])

        self.lineEditOMEGAX.setText(QString(str(omega_x)))
        self.lineEditOMEGAY.setText(QString(str(omega_y)))
        self.lineEditOMEGAZ.setText(QString(str(omega_z)))


    @pyqtSignature("const QString&")
    def slotOmegaX(self, text):
        """
        Input OMEGAX
        """
        omega_x, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOmega('omega_x', omega_x)


    @pyqtSignature("const QString&")
    def slotOmegaY(self, text):
        """
        Input OMEGAY
        """
        omega_y, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setOmega('omega_y', omega_y)


    @pyqtSignature("const QString&")
    def slotOmegaZ(self, text):
        """
        Input OmegaZ
        """
        omega_z, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
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
