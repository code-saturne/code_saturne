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

import sys
from code_saturne.model.XMLvariables import Model
from code_saturne.model.XMLengine import *
from code_saturne.model.XMLmodel import *


# -------------------------------------------------------------------------------
# Constructor
# -------------------------------------------------------------------------------

class NeptuneWallTransferModel(Variables, Model):
    AVAILABLE_MODELS = ["none", "nucleate_boiling", "droplet_evaporation_condensation"]

    def __init__(self, case):
        self.case = case
        self._wall_transfer_type = "none"

        xml_node_closure = self.case.xmlGetNode('closure_modeling')
        self._xml_node_mass_transfer = xml_node_closure.xmlInitNode('mass_transfer_model')
        self.load_xml()

    @property
    def wall_transfer_type(self):
        return self._wall_transfer_type

    @wall_transfer_type.setter
    def wall_transfer_type(self, model):
        if model not in self.AVAILABLE_MODELS:
            raise ValueError("Unknown model '{0}'. Model should be in {1}".format(model, self.AVAILABLE_MODELS))
        self._wall_transfer_type = model
        self.update_xml()

    def clear(self):
        self.wall_transfer_type = "none"
        self._xml_node_mass_transfer.xmlRemoveChildren()
        return

    def update_xml(self):
        self._xml_node_mass_transfer.xmlSetData('wall_transfer_type', self.wall_transfer_type)
        return

    def load_xml(self):
        value = self._xml_node_mass_transfer.xmlGetString('wall_transfer_type')
        if value == "":
            value = "none"
        self.wall_transfer_type = value
        return
