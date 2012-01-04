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
This module defines the values of reference.

This module contains the following classes and function:
- SteadyManagementView
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
from Pages.SteadyManagementForm import Ui_SteadyManagementForm
import Base.QtPage as QtPage
from Base.XMLvariables import Variables
from Pages.SteadyManagementModel import SteadyManagementModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SteadyManagementView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class SteadyManagementView(QWidget, Ui_SteadyManagementForm):
    """
    Class to open steady management Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_SteadyManagementForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = SteadyManagementModel(self.case)

        # Connections

        self.connect(self.lineEditRELXST, SIGNAL("textChanged(const QString &)"), self.slotRelaxCoef)
        self.connect(self.lineEditNTMABS, SIGNAL("textChanged(const QString &)"), self.slotNbIter)
        self.connect(self.checkBoxINPDT0, SIGNAL("clicked()"), self.slotZeroIteration)

        # Validators

        validatorRELXST = QtPage.DoubleValidator(self.lineEditRELXST, min=0.0, max=1.0)
        validatorRELXST.setExclusiveMin(True)
        self.lineEditRELXST.setValidator(validatorRELXST)

        validatorNTMABS = QtPage.IntValidator(self.lineEditNTMABS, min=0)
        self.lineEditNTMABS.setValidator(validatorNTMABS)

        # Initialization

        relax_coef = self.mdl.getRelaxCoefficient()
        self.lineEditRELXST.setText(QString(str(relax_coef)))

        nb_iter = self.mdl.getNbIter()
        self.lineEditNTMABS.setText(QString(str(nb_iter)))

        if self.mdl.getZeroIteration() == 'on':
            self.checkBoxINPDT0.setChecked(True)
        else:
            self.checkBoxINPDT0.setChecked(False)


    @pyqtSignature("const QString&")
    def slotRelaxCoef(self, text):
        """
        Input relaxation coefficient.
        """
        relax_coef, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setRelaxCoefficient(relax_coef)


    @pyqtSignature("const QString&")
    def slotNbIter(self, text):
        """
        Input itarations number.
        """
        nb_iter, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setNbIter(nb_iter)


    @pyqtSignature("")
    def slotZeroIteration(self):
        """
        Input zero iteration number.
        """
        if self.checkBoxINPDT0.isChecked():
            zero_iter = 'on'
        else:
            zero_iter = 'off'
        self.mdl.setZeroIteration(zero_iter)


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
