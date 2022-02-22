# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
# Library modules import
# -------------------------------------------------------------------------------

import logging

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *
from code_saturne.gui.case.CouplingParametersForm import Ui_CouplingParametersForm
from code_saturne.model.MobileMeshModel import MobileMeshModel
from code_saturne.model.ConjugateHeatTransferModel import ConjugateHeatTransferModel
from code_saturne.model.InternalCouplingModel import InternalCouplingModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("Boundary")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class CouplingParametersView(QWidget, Ui_CouplingParametersForm):

    def __init__(self, parent, case):
        QWidget.__init__(self, parent)
        Ui_CouplingParametersForm.__init__(self)
        self.setupUi(self)
        self.case = case
        self.list_of_tabs = ["internal", "cht", "fsi", "nepcat"]  # keep track of tab indices

        if self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            if "internal" in self.list_of_tabs:
                int_index = self.list_of_tabs.index("internal")
                self.couplingModelsTabWidget.removeTab(int_index)
                self.list_of_tabs.pop(int_index)
            if "fsi" in self.list_of_tabs:
                fsi_index = self.list_of_tabs.index("fsi")
                self.couplingModelsTabWidget.removeTab(fsi_index)
                self.list_of_tabs.pop(fsi_index)
            if "nepcat" not in self.list_of_tabs:
                nepcat_index = 1
                self.list_of_tabs.insert(nepcat_index, "nepcat")
                self.couplingModelsTabWidget.insertTab(nepcat_index, self.cathareCouplingTab, "Cathare coupling")
            self.cathareCouplingTab.setup(self.case)
        else:
            if "internal" not in self.list_of_tabs:
                int_index = 0
                self.list_of_tabs.insert(int_index, "internal")
                self.couplingModelsTabWidget.insertTab(int_index,
                                                       self.internalCouplingTab,
                                                       "Internal coupling")
            if "nepcat" in self.list_of_tabs:
                nepcat_index = self.list_of_tabs.index("nepcat")
                self.couplingModelsTabWidget.removeTab(nepcat_index)
                self.list_of_tabs.pop(nepcat_index)
            if "fsi" not in self.list_of_tabs:
                fsi_index = 1
                self.list_of_tabs.insert(fsi_index, "fsi")
                self.couplingModelsTabWidget.insertTab(fsi_index, self.fluidStructureInteractionTab,
                                                       "Fluid-structure interaction")

        if InternalCouplingModel(self.case).getZonesList():
            self.internalCouplingTab.setEnabled(True)
            self.internalCouplingTab.setup(self.case)
        else:
            self.internalCouplingTab.setEnabled(False)

        if MobileMeshModel(self.case).getMethod() == "off":
            self.fluidStructureInteractionTab.setEnabled(False)
        else:
            self.fluidStructureInteractionTab.setEnabled(True)
            self.fluidStructureInteractionTab.setup(self.case)

        if ConjugateHeatTransferModel(self.case).getNumberOfSyrthesCoupling() > 0:
            self.conjugateHeatTransferTab.setEnabled(True)
            self.conjugateHeatTransferTab.setup(self.case)
        else:
            self.conjugateHeatTransferTab.setEnabled(False)

        self.case.undoStopGlobal()
