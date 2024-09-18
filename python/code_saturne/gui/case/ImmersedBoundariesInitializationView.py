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
- ImmersedBoundariesInitializationView
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
from code_saturne.gui.case.ImmersedBoundariesInitialization import Ui_ImmersedBoundariesInitialization
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesInitializationView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# MainFieldsInitialization class
#-------------------------------------------------------------------------------

class ImmersedBoundariesInitializationView(QWidget, Ui_ImmersedBoundariesInitialization):

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesInitialization.__init__(self)
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
        self.pushButtonTemperature.clicked.connect(self.slotTemperature)

        self.update()
        self.case.undoStartGlobal()


    def update(self):

        exp = self.ibm.getObjectTemperatureFormula(self.current_obj-1)
        if exp:
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(exp)
        else:
            self.pushButtonTemperature.setStyleSheet("background-color: red")


    @pyqtSlot()
    def slotTemperature(self):
        """
        Formula for Temperature
        """
        exp, req, sym = self.ibm.getFormulaTemperature(self.current_obj-1)

        exa = """#example :
temperature = 273.15;"""

        object_name = self.ibm.getObjectName(self.current_obj)

        dialog = QMegEditorView(parent        = self,
                                #function_type = 'ibm_ini',
                                function_type = 'ibm_vol',
                                zone_name     = object_name,
                                variable_name = "porous_temperature",
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaTemperatureInit -> %s" % str(result))
            self.ibm.setObjectTemperatureFormula(self.current_obj-1, result)
            self.pushButtonTemperature.setStyleSheet("background-color: green")
            self.pushButtonTemperature.setToolTip(result)
