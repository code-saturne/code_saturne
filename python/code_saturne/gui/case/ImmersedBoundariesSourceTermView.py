# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the 'Main fields initialization' page.

This module contains the following classes:
- ImmersedBoundariesSourceTermView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator
from code_saturne.gui.case.ImmersedBoundariesSourceTerm import Ui_ImmersedBoundariesSourceTerm
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesSourceTermView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# MainFieldsInitialization class
#-------------------------------------------------------------------------------

class ImmersedBoundariesSourceTermView(QWidget, Ui_ImmersedBoundariesSourceTerm):

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesSourceTerm.__init__(self)
        self.setupUi(self)

        self.case = None
        self.zone = None
        self.zone_id = None


    def setup(self, case, ibm, current_obj):
        self.case = case
        self.ibm = ibm
        self.current_obj = current_obj
        self.case.undoStopGlobal()

        #Temperature initialization
        self.pushButtonThermalSourceTerm.clicked.connect(self.slotThermalSourceTerm)

        self.update()
        self.case.undoStartGlobal()


    def update(self):

        # Thermal source term
        exp = self.ibm.getObjectThermalSourceTermFormula(self.current_obj-1)
        if exp:
            self.pushButtonThermalSourceTerm.setStyleSheet("background-color: green")
            self.pushButtonThermalSourceTerm.setToolTip(exp)
        else:
            self.pushButtonThermalSourceTerm.setStyleSheet("background-color: red")


    @pyqtSlot()
    def slotThermalSourceTerm(self):
        """
        Formula for Temperature
        """
        exp, req, sym = self.ibm.getFormulaThermalSourceTerm(self.current_obj-1)

        exa = (
            "#Explicit part (W/m^3)\n"
            "st_exp = 0.;\n\n"
            "#Implicit part (W/m^3/K)\n"
            "st_imp = 0.;"
        )

        object_name = self.ibm.getObjectName(self.current_obj)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = object_name,
                                variable_name = "porous_thermal_st",
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermalSourceTerm -> %s" % str(result))
            self.ibm.setObjectThermalSourceTermFormula(self.current_obj-1, result)
            self.pushButtonThermalSourceTerm.setStyleSheet("background-color: green")
            self.pushButtonThermalSourceTerm.setToolTip(result)
