# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Pages.SteadyManagementForm import Ui_SteadyManagementForm
from code_saturne.Base.QtPage import DoubleValidator, IntValidator, from_qvariant
from code_saturne.Base.XMLvariables import Variables
from code_saturne.Pages.SteadyManagementModel import SteadyManagementModel

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
        self.case.undoStopGlobal()
        self.mdl = SteadyManagementModel(self.case)

        # Connections

        self.lineEditRELXST.textChanged[str].connect(self.slotRelaxCoef)
        self.lineEditNTMABS.textChanged[str].connect(self.slotNbIter)
        self.checkBoxINPDT0.clicked.connect(self.slotZeroIteration)

        # Validators

        validatorRELXST = DoubleValidator(self.lineEditRELXST, min=0.0, max=1.0)
        validatorRELXST.setExclusiveMin(True)
        self.lineEditRELXST.setValidator(validatorRELXST)

        validatorNTMABS = IntValidator(self.lineEditNTMABS, min=0)
        self.lineEditNTMABS.setValidator(validatorNTMABS)

        # Initialization

        relax_coef = self.mdl.getRelaxCoefficient()
        self.lineEditRELXST.setText(str(relax_coef))

        nb_iter = self.mdl.getNbIter()
        self.lineEditNTMABS.setText(str(nb_iter))

        if self.mdl.getZeroIteration() == 'on':
            self.checkBoxINPDT0.setChecked(True)
        else:
            self.checkBoxINPDT0.setChecked(False)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotRelaxCoef(self, text):
        """
        Input relaxation coefficient.
        """
        if self. lineEditRELXST.validator().state == QValidator.Acceptable:
            relax_coef = from_qvariant(text, float)
            self.mdl.setRelaxCoefficient(relax_coef)


    @pyqtSlot(str)
    def slotNbIter(self, text):
        """
        Input itarations number.
        """
        if self. lineEditNTMABS.validator().state == QValidator.Acceptable:
            nb_iter = from_qvariant(text, int)
            self.mdl.setNbIter(nb_iter)


    @pyqtSlot()
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
