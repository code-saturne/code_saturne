# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module contains the following classes:
- BoundaryConditionsWallView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsWall import Ui_BoundaryConditionsWall

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel

from code_saturne.model.MainFieldsModel import MainFieldsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsWallView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsWallView(QWidget, Ui_BoundaryConditionsWall) :
    """
    Boundary condition for energy
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsWall.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxWallModel.activated[str].connect(self.__slotWall)

        self.__WallModel = ComboModel(self.comboBoxWallModel, 5, 1)
        self.__WallModel.addItem(self.tr("adherence"), "adherence")
        self.__WallModel.addItem(self.tr("friction"), "friction")
        self.__WallModel.addItem(self.tr("dU2/dn = 0"), "du2_dn")
        self.__WallModel.addItem(self.tr("dVR/dn = 0"), "dvr_dn")
        self.__WallModel.addItem(self.tr("droplet friction"), "droplet_friction")


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None
        self.__currentField = fieldId


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        if len(MainFieldsModel(self.case).getSolidFieldIdList()) > 0:
            mdl = boundary.getWallModel(self.__currentField)
            self.__WallModel.setItem(str_model = mdl)
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotWall(self, text):
        """
        INPUT wall model
        """
        mdl =  self.__WallModel.dicoV2M[str(text)]
        self.__boundary.setWallModel(self.__currentField, mdl)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
