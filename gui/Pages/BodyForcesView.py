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
from code_saturne.Pages.CoriolisSourceTermsModel import CoriolisSourceTermsModel

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

        self.gmdl = BodyForcesModel(self.case)
        self.cmdl = None

        if case.xmlRootNode().tagName != "NEPTUNE_CFD_GUI":
            node_pm = self.case.xmlGetNode('thermophysical_models')
            if node_pm:
                node = node_pm.xmlGetNode('groundwater_model',  'model')
                if not node or node['model'] == 'off':
                    self.cmdl = CoriolisSourceTermsModel(self.case)
        if self.cmdl is None:
            self.groupBoxCoriolis.hide()

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
        gravity_x = self.gmdl.getGravity(self.gmdl.nodes[0])
        gravity_y = self.gmdl.getGravity(self.gmdl.nodes[1])
        gravity_z = self.gmdl.getGravity(self.gmdl.nodes[2])

        self.lineEditX.setText(str(gravity_x))
        self.lineEditY.setText(str(gravity_y))
        self.lineEditZ.setText(str(gravity_z))

        # Coriolis source terms view

        if self.cmdl:
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
            omega_x = self.cmdl.getOmega(self.cmdl.nodes[0])
            omega_y = self.cmdl.getOmega(self.cmdl.nodes[1])
            omega_z = self.cmdl.getOmega(self.cmdl.nodes[2])

            self.lineEditOMEGAX.setText(str(omega_x))
            self.lineEditOMEGAY.setText(str(omega_y))
            self.lineEditOMEGAZ.setText(str(omega_z))

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotGravityX(self, text):
        """
        Input GX
        """
        if self.lineEditX.validator().state == QValidator.Acceptable:
            gravity_x = from_qvariant(text, float)
            self.gmdl.setGravity('gravity_x', gravity_x)


    @pyqtSlot(str)
    def slotGravityY(self, text):
        """
        Input GY
        """
        if self.lineEditY.validator().state == QValidator.Acceptable:
            gravity_y = from_qvariant(text, float)
            self.gmdl.setGravity('gravity_y', gravity_y)


    @pyqtSlot(str)
    def slotGravityZ(self, text):
        """
        Input GZ
        """
        if self.lineEditZ.validator().state == QValidator.Acceptable:
            gravity_z = from_qvariant(text, float)
            self.gmdl.setGravity('gravity_z', gravity_z)


    @pyqtSlot(str)
    def slotOmegaX(self, text):
        """
        Input OMEGAX
        """
        if self.lineEditOMEGAX.validator().state == QValidator.Acceptable:
            omega_x = from_qvariant(text, float)
            self.cmdl.setOmega('omega_x', omega_x)


    @pyqtSlot(str)
    def slotOmegaY(self, text):
        """
        Input OMEGAY
        """
        if self.lineEditOMEGAY.validator().state == QValidator.Acceptable:
            omega_y = from_qvariant(text, float)
            self.cmdl.setOmega('omega_y', omega_y)


    @pyqtSlot(str)
    def slotOmegaZ(self, text):
        """
        Input OmegaZ
        """
        if self.lineEditOMEGAZ.validator().state == QValidator.Acceptable:
            omega_z = from_qvariant(text, float)
            self.cmdl.setOmega('omega_z', omega_z)


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
