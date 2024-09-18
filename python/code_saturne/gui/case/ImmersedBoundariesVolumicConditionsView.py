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
from code_saturne.gui.case.ImmersedBoundariesVolumicConditionsForm import Ui_ImmersedBoundariesVolumicConditionsForm
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel

# -------------------------------------------------------------------------------
# Widgets import
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesVolumicConditionsView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class ImmersedBoundariesVolumicConditionsView(QWidget, Ui_ImmersedBoundariesVolumicConditionsForm):
    """ Display available volumic treatments for a given zone """

    def __init__(self, parent, case, zone_name):
        QWidget.__init__(self, parent)
        Ui_ImmersedBoundariesVolumicConditionsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent
        self.object_name = zone_name
        self.ibm = ImmersedBoundariesModel(self.case)
        self.current_obj = self.ibm.getObjectNumFromName(self.object_name)
        self.case.undoStopGlobal()

        #print("object_name, num, list = ", self.zone_name, self.current_obj, self.ibm.getObjectsNameList())

        self.IBMneptuneFSIWidget.hide()
        self.IBMneptunePropertiesWidget.hide()
        self.IBMneptuneInitializationWidget.hide()
        self.IBMneptuneSourceTermWidget.hide()

        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":

            #Object Physical Properties
            if (self.ibm.getObjectPhysicalProperties(self.current_obj) == 'on'):
                self.IBMneptunePropertiesWidget.show()
                self.IBMneptunePropertiesWidget.setup(self.case,
                                                      self.ibm,
                                                      self.current_obj)

            #FSI Object Physical Properties
            if (self.ibm.getObjectFSI(self.current_obj) == 'on'):
                self.IBMneptuneFSIWidget.show()
                self.IBMneptuneFSIWidget.setup(self.case,
                                               self.ibm,
                                               self.current_obj)
            #Object Initialization
            if (self.ibm.getObjectInit(self.current_obj) == 'on'):
                self.IBMneptuneInitializationWidget.show()
                self.IBMneptuneInitializationWidget.setup(self.case,
                                                          self.ibm,
                                                          self.current_obj)

            #Object source term
            if (self.ibm.getObjectCHT(self.current_obj) == 'on'):
                self.IBMneptuneSourceTermWidget.show()
                self.IBMneptuneSourceTermWidget.setup(self.case,
                                                      self.ibm,
                                                      self.current_obj)

            #Remove tab
            if (self.ibm.getObjectThermalSourceTerm(self.current_obj) == 'off'):
                for i in range(self.IBMtabWidget.count()):
                    if self.IBMtabWidget.tabText(i) == "Thermal source term":
                        self.IBMtabWidget.removeTab(i)
                        break

            if (self.ibm.getObjectInit(self.current_obj) == 'off'):
                for i in range(self.IBMtabWidget.count()):
                    if self.IBMtabWidget.tabText(i) == "Initialization":
                        self.IBMtabWidget.removeTab(i)
                        break

            if (self.ibm.getObjectPhysicalProperties(self.current_obj) == 'off'):
                for i in range(self.IBMtabWidget.count()):
                    if self.IBMtabWidget.tabText(i) == "Physical properties":
                        self.IBMtabWidget.removeTab(i)
                        break

            if (self.ibm.getObjectFSI(self.current_obj) == 'off'):
                for i in range(self.IBMtabWidget.count()):
                    if self.IBMtabWidget.tabText(i) == "Fluid structure interaction":
                        self.IBMtabWidget.removeTab(i)
                        break


        self.case.undoStartGlobal()

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
