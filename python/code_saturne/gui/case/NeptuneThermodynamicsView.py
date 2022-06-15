# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel

from NeptuneThermodynamics import Ui_NeptuneThermodynamics


class NeptuneThermodynamicsView(QWidget, Ui_NeptuneThermodynamics):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        Ui_NeptuneThermodynamics.__init__(self)
        self.setupUi(self)

        self.case = None
        self.zone_name = None
        self.tabWidget.currentChanged[int].connect(self.slotRefresh)

    def setup(self, case, zone_name):
        self.case = case
        self.zone_name = zone_name
        self.tabFieldProperties.setup(self.case, self.zone_name)
        self.tabSaturationProperties.setup(self.case, self.zone_name)
        self.tabInteractionProperties.setup(self.case, self.zone_name)

    @pyqtSlot(int)
    def slotRefresh(self, index):
        if index < 2:
            self.tabWidget.widget(index).setup(self.case, self.zone_name)
        else:
            self.tabWidget.widget(index).refresh()
