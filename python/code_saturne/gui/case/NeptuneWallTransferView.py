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

"""
This module defines the 'Nucleate boiling' page.

This module contains the following classes:
- NucleateBoilingView
"""

# -------------------------------------------------------------------------------
# Library modules import
# -------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.base.QtPage import ComboModel, from_qvariant
from code_saturne.gui.case.NeptuneWallTransferForm import Ui_NeptuneWallTransferForm
from code_saturne.model.NeptuneWallTransferModel import NeptuneWallTransferModel
from code_saturne.model.MainFieldsModel import MainFieldsModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NeptuneWallTransferView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
#  class InterfacialForces
# -------------------------------------------------------------------------------

class NeptuneWallTransferView(QWidget, Ui_NeptuneWallTransferForm):

    def __init__(self, parent, case):
        QWidget.__init__(self, parent)
        Ui_NeptuneWallTransferForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.xml_model = NeptuneWallTransferModel(self.case)

        self.set_connections()
        # Connections must be set before ComboModel is declared.
        self.combomodel_wallmodel = ComboModel(self.comboBoxWallTransferType)
        self.fill_widgets()
        self.initialize_widgets()
        self.update_view()

    def initialize_widgets(self):
        predefined_flow = MainFieldsModel(self.case).getPredefinedFlow()
        if predefined_flow in ["free_surface", "boiling_flow"]:
            self.combomodel_wallmodel.setItem(str_model="nucleate_boiling")
            self.setEnabled(False)
        elif predefined_flow == "droplet_flow":
            self.combomodel_wallmodel.setItem(str_model="droplet_evaporation_condensation")
            self.setEnabled(False)
        elif predefined_flow == "multiregime":
            self.combomodel_wallmodel.disableItem(str_model="none")
            self.combomodel_wallmodel.setItem(str_model=self.xml_model.wall_transfer_type)
            self.boilingWidget.setEnabled(False)
            self.dropletWidget.setEnabled(False)
        else:
            self.combomodel_wallmodel.setItem(str_model=self.xml_model.wall_transfer_type)

    def update_view(self):
        model = self.xml_model.wall_transfer_type
        self.boilingWidget.hide()
        self.dropletWidget.hide()
        if model == "nucleate_boiling":
            self.boilingWidget.setup(self.case)
            self.boilingWidget.show()
        elif model == "droplet_evaporation_condensation":
            self.dropletWidget.setup(self.case)
            self.dropletWidget.show()
        elif model == "none":
            pass
        else:
            raise ValueError("Unknown model {0}".format(model))
        return

    def fill_widgets(self):
        self.combomodel_wallmodel.addItem(self.tr("None"), "none")
        self.combomodel_wallmodel.addItem(self.tr("Nucleate boiling"), "nucleate_boiling")
        self.combomodel_wallmodel.addItem(self.tr("Droplet evaporation/condensation"),
                                          "droplet_evaporation_condensation")
        return

    def set_connections(self):
        self.comboBoxWallTransferType.currentTextChanged[str].connect(self.slot_set_wall_model)

    @pyqtSlot(str)
    def slot_set_wall_model(self, text):
        model = self.combomodel_wallmodel.dicoV2M[text]
        self.xml_model.wall_transfer_type = model
        self.update_view()
        return
