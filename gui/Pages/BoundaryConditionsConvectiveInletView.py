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
This module contains the following classes:
- BoundaryConditionsConvectiveInletView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BoundaryConditionsConvectiveInletForm import Ui_BoundaryConditionsConvectiveInletForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary
from code_saturne.Pages.CompressibleModel import CompressibleModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsConvectiveInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsConvectiveInletView(QWidget, Ui_BoundaryConditionsConvectiveInletForm):
    """
    Boundary condition for velocity in inlet, without particular physics.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsConvectiveInletForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        # Connections
        self.groupBoxConvectiveInlet.clicked[bool].connect(self.__slotConvectiveInlet)

        # Validators

        # Apply validators

        self.__case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        checked = False
        hide = False

        mdl = CompressibleModel(self.__case)
        if mdl.getCompressibleModel() != "off":
            hide = True

        if not hide and self.__boundary.getConvectiveInletStatus() == "on":
            checked = True

        if checked:
            self.groupBoxConvectiveInlet.setChecked(True)
        else:
            self.groupBoxConvectiveInlet.setChecked(False)
        self.__slotConvectiveInlet(checked)

        if hide:
            self.hide()
        else:
            self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(bool)
    def __slotConvectiveInlet(self, checked):
        """
        Private slot.

        Activates convective inlet boundary condition.

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox convective inlet parameters.
        """
        self.groupBoxConvectiveInlet.setFlat(not checked)

        if checked:
            self.__boundary.setConvectiveInletStatus("on")
        else:
            self.__boundary.setConvectiveInletStatus("off")


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
