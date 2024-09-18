# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Standard modules
# -------------------------------------------------------------------------------

import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.case.ImmersedBoundariesBoundaryConditionsForm import Ui_ImmersedBoundariesBoundaryConditionsForm
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel

# -------------------------------------------------------------------------------
# Widgets import
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesBoundaryConditionsView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class ImmersedBoundariesBoundaryConditionsView(QWidget, Ui_ImmersedBoundariesBoundaryConditionsForm):
    """ Display available boundary process for a given zone (boundary of the immersed object) """

    def __init__(self, parent, case, zone_name):
        QWidget.__init__(self, parent)
        Ui_ImmersedBoundariesBoundaryConditionsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.object_name = zone_name
        self.ibm = ImmersedBoundariesModel(self.case)
        self.current_obj = self.ibm.getObjectNumFromName(self.object_name)
        self.mfm = MainFieldsModel(self.case)
        self.case.undoStopGlobal()

        self.IBMneptuneWallWidget.hide()
        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            groupBoxIBMBoundary.hide()
            return

        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":

            #Wall boundary condition for immersed object
            if (self.ibm.getObjectBoundaryConditionNature(self.current_obj) == 'Wall'):
                self.IBMneptuneWallWidget.show()
                self.IBMneptuneWallWidget.setup(self.case,
                                                self.ibm,
                                                self.current_obj)



#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
