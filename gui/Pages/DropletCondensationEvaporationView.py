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
This module defines the 'Droplet Condensation-Evaporation' page.

This module contains the following classes:
- DropletCondensationEvaporationView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, from_qvariant
from DropletCondensationEvaporation import Ui_DropletCondensationEvaporation
from DropletCondensationEvaporationModel import DropletCondensationEvaporationModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("DropletCondensationEvaporationView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
#  class InterfacialForces
#-------------------------------------------------------------------------------

class DropletCondensationEvaporationView(QWidget, Ui_DropletCondensationEvaporation):
    """
    Droplet Condensation-Evaporation model layout.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_DropletCondensationEvaporation.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = DropletCondensationEvaporationModel(self.case)

        self.modelYPlus = ComboModel(self.comboBoxYPlus, 3, 1)
        self.modelYPlus.addItem(self.tr("Boundary cell center"), "center")
        self.modelYPlus.addItem(self.tr("Y+ = "), "Yplus_value")
        self.modelYPlus.addItem(self.tr("Droplets diameter"), "diameter")

        # Validators

        validatorYplus = DoubleValidator(self.lineEditYPlus, min = 0.0)

        validatorYplus.setExclusiveMin(True)

        self.lineEditYPlus.setValidator(validatorYplus)

        # Connect signals to slots
        self.comboBoxYPlus.activated[str].connect(self.slotYPlus)
        self.lineEditYPlus.textChanged[str].connect(self.slotYPlusValue)

        isYPlus = self.mdl.getYPlusModel()
        self.modelYPlus.setItem(str_model=isYPlus)

        if isYPlus == "Yplus_value" :
           self.lineEditYPlus.show()
           self.lineEditYPlus.setText(str(self.mdl.getYPlusValue()))
        else :
           self.lineEditYPlus.hide()

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotYPlus(self, text):
        """
        configure Y Plus model
        """
        value = self.modelYPlus.dicoV2M[text]
        log.debug("slotYPlus -> %s" % value)
        self.mdl.setYPlusModel(value)

        if value == "Yplus_value" :
           self.lineEditYPlus.show()
           self.lineEditYPlus.setText(str(self.mdl.getYPlusValue()))
        else :
           self.lineEditYPlus.hide()


    @pyqtSlot(str)
    def slotYPlusValue(self, text):
        """
        Update the Yplus value
        """
        if self.lineEditYPlus.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.mdl.setYPlusValue(value)


