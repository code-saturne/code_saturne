# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from BodyForcesForm import Ui_BodyForcesForm
import Base.QtPage as QtPage
from Pages.BodyForcesModel import BodyForcesModel

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
        self.mdl = BodyForcesModel(self.case)

        # Connections

        self.connect(self.lineEditX, SIGNAL("textChanged(const QString &)"), self.slotGravityX)
        self.connect(self.lineEditY, SIGNAL("textChanged(const QString &)"), self.slotGravityY)
        self.connect(self.lineEditZ, SIGNAL("textChanged(const QString &)"), self.slotGravityZ)
        self.connect(self.radioButtonYes, SIGNAL("clicked()"), self.slotHydrostaticPressure)
        self.connect(self.radioButtonNo, SIGNAL("clicked()"), self.slotHydrostaticPressure)

        # Validators
        validatorX = QtPage.DoubleValidator(self.lineEditX)
        validatorY = QtPage.DoubleValidator(self.lineEditY)
        validatorZ = QtPage.DoubleValidator(self.lineEditZ)

        self.lineEditX.setValidator(validatorX)
        self.lineEditY.setValidator(validatorY)
        self.lineEditZ.setValidator(validatorZ)

        # Initialization

        if self.mdl.getHydrostaticPressure() == "on":
            self.radioButtonYes.setChecked(True)
            self.radioButtonNo.setChecked(False)
        else:
            self.radioButtonYes.setChecked(False)
            self.radioButtonNo.setChecked(True)

        gravity_x = self.mdl.getGravity(self.mdl.nodes[0])
        gravity_y = self.mdl.getGravity(self.mdl.nodes[1])
        gravity_z = self.mdl.getGravity(self.mdl.nodes[2])

        self.lineEditX.setText(QString(str(gravity_x)))
        self.lineEditY.setText(QString(str(gravity_y)))
        self.lineEditZ.setText(QString(str(gravity_z)))


    @pyqtSignature("const QString&")
    def slotGravityX(self, text):
        """
        Input GX
        """
        gravity_x, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setGravity('gravity_x', gravity_x)


    @pyqtSignature("const QString&")
    def slotGravityY(self, text):
        """
        Input GY
        """
        gravity_y, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setGravity('gravity_y', gravity_y)


    @pyqtSignature("const QString&")
    def slotGravityZ(self, text):
        """
        Input GZ
        """
        gravity_z, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setGravity('gravity_z', gravity_z)


    @pyqtSignature("")
    def slotHydrostaticPressure(self):
        """
        Input IHYDPR.
        """
        if self.radioButtonYes.isChecked():
            hpr = 'on'
        else:
            hpr = 'off'
        self.mdl.setHydrostaticPressure(hpr)


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
