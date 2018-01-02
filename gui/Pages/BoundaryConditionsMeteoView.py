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
- BoundaryConditionsMeteoView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

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

from code_saturne.Pages.BoundaryConditionsMeteoForm import Ui_BoundaryConditionsMeteoForm
from code_saturne.Pages.AtmosphericFlowsModel import AtmosphericFlowsModel

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsMeteoView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsMeteoView(QWidget, Ui_BoundaryConditionsMeteoForm):
    """
    Boundary condifition for the velocity part
    """
    def __init__(self, parent):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsMeteoForm.__init__(self)
        self.setupUi(self)


    def setup(self, case, velocityWidget, turbulenceWidget, scalarsWidget):
        """
        Setup the widget.
        """
        self.__case = case
        self.velocityWidget = velocityWidget
        self.turbulenceWidget = turbulenceWidget
        self.scalarsWidget = scalarsWidget
        self.__boundary = None

        sca_mo  = DefineUserScalarsModel(self.__case)
        self.species_list = sca_mo.getUserScalarNameList()

        self.__case.undoStopGlobal()

        self.__model = AtmosphericFlowsModel(self.__case)

        self.checkBoxReadData.clicked[bool].connect(self.__slotReadData)
        self.checkBoxAutoNature.clicked[bool].connect(self.__slotAutoNature)

        self.__case.undoStartGlobal()


    def showWidget(self, b):
        """
        Show the widget.
        """
        self.__b = b
        if self.__model.getAtmosphericFlowsModel() != "off" \
            and self.__model.getMeteoDataStatus() == "on":
            self.show()

            label = b.getLabel()
            nature = "meteo_" + b.getNature()
            self.__boundary = Boundary(nature, label, self.__case)

            if self.__boundary.getMeteoDataStatus() == 'on':
                self.checkBoxReadData.setChecked(True)
                self.checkBoxAutoNature.setEnabled(True)
                self.velocityWidget.hideWidget()
                self.turbulenceWidget.hideWidget()
            else:
                self.checkBoxReadData.setChecked(False)
                self.checkBoxAutoNature.setEnabled(False)
                if nature == "meteo_inlet":
                    self.velocityWidget.showWidget(b)
                    self.turbulenceWidget.showWidget(b)
                else:
                    self.velocityWidget.hideWidget()
                    self.turbulenceWidget.hideWidget()

            if self.__boundary.getAutomaticNatureStatus() == 'on':
                self.checkBoxAutoNature.setChecked(True)
            else:
                self.checkBoxAutoNature.setChecked(False)

            self.scalarsWidget.showWidget(self.__b)

        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all.
        """
        self.hide()


    def __slotReadData(self, bool):
        """
        Input if the meteo data must be read.
        """
        if bool == True:
            self.__boundary.setMeteoDataStatus('on')
            self.checkBoxAutoNature.setEnabled(True)
            self.velocityWidget.hideWidget()
            self.turbulenceWidget.hideWidget()
        else:
            self.__boundary.setMeteoDataStatus('off')
            self.checkBoxAutoNature.setChecked(False)
            self.__boundary.setAutomaticNatureStatus('off')
            self.checkBoxAutoNature.setEnabled(False)
            if self.__boundary.getNature() == "meteo_inlet":
                self.velocityWidget.showWidget(self.__b)
                self.turbulenceWidget.showWidget(self.__b)
            else:
                self.velocityWidget.hideWidget()
                self.turbulenceWidget.hideWidget()

        self.scalarsWidget.showWidget(self.__b)


    def __slotAutoNature(self, bool):
        """
        Input if the nature of the boundary must be detected automaticaly.
        """
        if bool == True:
            self.__boundary.setAutomaticNatureStatus('on')
        else:
            self.__boundary.setAutomaticNatureStatus('off')


    def tr(self, text):
        """
        Translation.
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
