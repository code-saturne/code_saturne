# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
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
# Standard modules
# -------------------------------------------------------------------------------

import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.Base.QtCore import *
from code_saturne.Base.QtGui import *
from code_saturne.Base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Pages.VolumicConditionsForm import Ui_VolumicConditionsForm

# -------------------------------------------------------------------------------
# Widgets import
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("VolumicConditionsView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class VolumicConditionsView(QWidget, Ui_VolumicConditionsForm):
    """ Display available volumic treatments for a given zone """

    def __init__(self, parent, case, zone_name):
        QWidget.__init__(self, parent)
        Ui_VolumicConditionsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.zone_name = zone_name
        self.case.undoStopGlobal()

        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            self.saturneInitializationWidget.hide()
            self.neptuneInitializationWidget.setup(self.case, self.zone_name)
        else:
            self.neptuneInitializationWidget.hide()
            self.saturneInitializationWidget.setup(self.case, self.zone_name)
        self.porosityPage.setup(self.case, self.zone_name)
        self.headLossesPage.setup(self.case, self.zone_name)
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            self.saturneSourceTermsWidget.hide()
            self.neptuneSourceTermsWidget.setup(self.case, self.zone_name)
        else:
            self.neptuneSourceTermsWidget.hide()
            self.saturneSourceTermsWidget.setup(self.case, self.zone_name)
        self.groundwaterLawPage.setup(self.case, self.zone_name)
