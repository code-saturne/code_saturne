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
This module defines the values of reference.

This module contains the following classes and function:
- ReferenceValuesView
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
from Pages.ReferenceValuesForm import Ui_ReferenceValuesForm
import Base.QtPage as QtPage
from Pages.ReferenceValuesModel import ReferenceValuesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ReferenceValuesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ReferenceValuesView(QWidget, Ui_ReferenceValuesForm):
    """
    Class to open Reference Pressure Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ReferenceValuesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = ReferenceValuesModel(self.case)

        # Connections

        self.connect(self.lineEditP0, SIGNAL("textChanged(const QString &)"), self.slotPressure)
        self.connect(self.lineEditT0, SIGNAL("textChanged(const QString &)"), self.slotTemperature)
        self.connect(self.lineEditMassMolar, SIGNAL("textChanged(const QString &)"), self.slotMassemol)

        # Validators

        validatorP0 = QtPage.DoubleValidator(self.lineEditP0, min=0.0)
        self.lineEditP0.setValidator(validatorP0)

        validatorT0 = QtPage.DoubleValidator(self.lineEditT0,  min=0.0)
        validatorT0.setExclusiveMin(True)
        self.lineEditT0.setValidator(validatorT0)

        validatorMM = QtPage.DoubleValidator(self.lineEditMassMolar, min=0.0)
        validatorMM.setExclusiveMin(True)
        self.lineEditMassMolar.setValidator(validatorMM)

        # Display

        model, node = self.mdl.getParticularPhysical()

        if model == "atmo":
            self.groupBoxTemperature.show()
            self.labelInfoT0.hide()
            self.groupBoxMassMolar.hide()
        elif model != "off":
            self.groupBoxTemperature.show()
            self.groupBoxMassMolar.show()
        else:
            self.groupBoxTemperature.hide()
            self.groupBoxMassMolar.hide()

        # Initialization

        p = self.mdl.getPressure()
        self.lineEditP0.setText(QString(str(p)))

        model, node = self.mdl.getParticularPhysical()
        if model == "atmo":
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(QString(str(t)))
        elif model != "off":
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(QString(str(t)))
            m = self.mdl.getMassemol()
            self.lineEditMassMolar.setText(QString(str(m)))


    @pyqtSignature("const QString&")
    def slotPressure(self,  text):
        """
        Input PRESS.
        """
        p, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setPressure(p)


    @pyqtSignature("const QString&")
    def slotTemperature(self,  text):
        """
        Input TEMPERATURE.
        """
        t, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setTemperature(t)


    @pyqtSignature("const QString&")
    def slotMassemol(self,  text):
        """
        Input Mass molar.
        """
        m, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setMassemol(m)


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
