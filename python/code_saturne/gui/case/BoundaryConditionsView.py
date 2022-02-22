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
- StandardItemModelBoundaries
- BoundaryConditionsView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsForm import Ui_BoundaryConditionsForm

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary
from code_saturne.model.MobileMeshModel import MobileMeshModel
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.model.LagrangianModel import LagrangianModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsView")
log.setLevel(GuiParam.DEBUG)


class BoundaryConditionsView(QWidget, Ui_BoundaryConditionsForm):
    """
    Main boundary conditions view class
    """
    def __init__(self, parent, case, zone_name):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsForm.__init__(self)
        self.setupUi(self)
        self.case = case

        self.case.undoStopGlobal()

        d = LocalizationModel('BoundaryZone', self.case)
        for zone in d.getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone

        _title = zone_name
        if self.zone:
            _title+= " ["
            _title+= self.zone.getModel2ViewDictionary()[self.zone.getNature()]
            _title+= "]"

        self.groupBoxMain.setTitle(_title)

        # Set the case for custom widgets
        self.roughWidget.setup(self.case)
        self.slidingWidget.setup(self.case)
        self.convectiveInletWidget.setup(self.case)
        self.mappedInletWidget.setup(self.case)
        self.velocityWidget.setup(self.case)
        self.turbulenceWidget.setup(self.case)
        self.compressibleOutletWidget.setup(self.case)
        self.coalWidget.setup(self.case)
        self.scalarsWidget.setup(self.case)
        self.meteoWidget.setup(self.case, self.velocityWidget,
                               self.turbulenceWidget,
                               self.scalarsWidget)
        self.mobileMeshWidget.setup(self.case)
        self.radiativeWidget.setup(self.case)
        self.electricalWidget.setup(self.case)
        self.hydraulicheadWidget.setup(self.case)
        self.pressureWidget.setup(self.case)
        self.externalHeadLossesWidget.setup(self.case)
        self.particlesWidget.setup(self.case)

        self.__hideAllWidgets()
        self.__selectBoundary()

        self.case.undoStartGlobal()

    def __selectBoundary(self):
        """
        Select a boundary in the QTreeView.
        """
        label = self.zone.getLabel()
        nature = self.zone.getNature()
        log.debug("slotSelectBoundary label %s (%s)" % (label, nature))

        self.__hideAllWidgets()
        boundary = Boundary(nature, label, self.case)

        if LagrangianModel(self.case).getLagrangianModel() != 'off':
            self.particlesWidget.showWidget(self.zone)
        if nature == 'wall':
            self.__selectWallBoundary(boundary)
        elif nature == 'inlet':
            self.__selectInletBoundary(boundary)
        elif nature == 'outlet':
            self.__selectOutletBoundary(boundary)
        elif nature == 'symmetry':
            self.__selectSymmetryBoundary(boundary)
        elif nature == 'free_inlet_outlet':
            self.__selectInletOutletBoundary(boundary)
        elif nature == 'imposed_p_outlet':
            self.__selectImposedPressureOutletBoundary(boundary)
        elif nature == 'groundwater':
            self.__selectGroundwaterBoundary(boundary)


    def __selectInletBoundary(self, boundary):
        """
        Shows widgets for inlet.
        """

        if GroundwaterModel(self.case).getGroundwaterModel() != 'off':
            # TODO is this the desired behaviour ?
            return

        self.hydraulicheadWidget.hideWidget()
        if self.coalWidget.getCoalNumber() == 0:
            # TODO how should the widgets behave when groundwater model is activated ?
            if GroundwaterModel(self.case).getGroundwaterModel() == "off":
                self.velocityWidget.showWidget(boundary)
            else:
                self.velocityWidget.hideWidget()
                self.hydraulicheadWidget.showWidget(boundary)
            self.coalWidget.hideWidget()
        else:
            self.velocityWidget.hideWidget()
            self.coalWidget.showWidget(boundary)

        self.turbulenceWidget.showWidget(boundary)
        self.meteoWidget.showWidget(boundary)
        self.scalarsWidget.showWidget(boundary)
        self.mobileMeshWidget.showWidget(boundary)
        self.electricalWidget.showWidget(boundary)
        self.externalHeadLossesWidget.hideWidget()
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.showWidget(boundary)
        self.mappedInletWidget.showWidget(boundary)


    def __selectWallBoundary(self, boundary):
        """
        Shows widgets for wall.
        """
        self.slidingWidget.showWidget(boundary)
        self.roughWidget.showWidget(boundary)

        self.scalarsWidget.showWidget(boundary)
        self.mobileMeshWidget.showWidget(boundary)
        self.radiativeWidget.showWidget(boundary)
        self.electricalWidget.showWidget(boundary)
        self.externalHeadLossesWidget.hideWidget()
        self.hydraulicheadWidget.hideWidget()
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()


    def __selectOutletBoundary(self, boundary):
        """
        Shows widgets for wall.
        """
        self.scalarsWidget.showWidget(boundary)
        self.mobileMeshWidget.showWidget(boundary)
        self.meteoWidget.showWidget(boundary)
        if self.compressibleOutletWidget.getCompressibleModel() != "off":
            self.compressibleOutletWidget.showWidget(boundary)
        else:
            self.compressibleOutletWidget.hideWidget()
        self.electricalWidget.showWidget(boundary)
        self.externalHeadLossesWidget.hideWidget()
        if GroundwaterModel(self.case).getGroundwaterModel() == "off":
            self.hydraulicheadWidget.hideWidget()
        else:
            self.hydraulicheadWidget.showWidget(boundary)
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()


    def __selectInletOutletBoundary(self, boundary):
        """
        Shows widgets for free inlet outlet.
        """
        self.coalWidget.hideWidget()
        self.velocityWidget.hideWidget()
        self.turbulenceWidget.hideWidget()
        self.meteoWidget.hideWidget()
        self.scalarsWidget.hideWidget()
        self.mobileMeshWidget.showWidget(boundary)
        self.electricalWidget.hideWidget()
        self.externalHeadLossesWidget.showWidget(boundary)
        self.hydraulicheadWidget.hideWidget()
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()


    def __selectImposedPressureOutletBoundary(self, boundary):
        """
        Shows widgets for imposed pressure outlet.
        """
        self.coalWidget.hideWidget()
        self.velocityWidget.hideWidget()
        self.turbulenceWidget.hideWidget()
        self.meteoWidget.hideWidget()
        self.scalarsWidget.hideWidget()
        self.mobileMeshWidget.showWidget(boundary)
        self.electricalWidget.hideWidget()
        self.externalHeadLossesWidget.hideWidget()
        self.hydraulicheadWidget.hideWidget()
        self.pressureWidget.showWidget(boundary)
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()


    def __selectGroundwaterBoundary(self, boundary):
        """
        Shows widgets for groundwater flow.
        """
        self.coalWidget.hideWidget()
        self.velocityWidget.hideWidget()
        self.turbulenceWidget.hideWidget()
        self.meteoWidget.hideWidget()
        self.scalarsWidget.showWidget(boundary)
        self.mobileMeshWidget.showWidget(boundary)
        self.electricalWidget.hideWidget()
        self.externalHeadLossesWidget.hideWidget()
        self.hydraulicheadWidget.showWidget(boundary)
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()


    def __selectSymmetryBoundary(self, boundary):
        """
        Shows widgets for wall.
        """
        self.mobileMeshWidget.showWidget(boundary)


    def __hideAllWidgets(self):
        """
        Hides all promoted QWidgets, for all nature.
        """
        self.slidingWidget.hideWidget()
        self.roughWidget.hideWidget()
        self.velocityWidget.hideWidget()
        self.turbulenceWidget.hideWidget()
        self.coalWidget.hideWidget()
        self.compressibleOutletWidget.hideWidget()
        self.meteoWidget.hideWidget()
        self.scalarsWidget.hideWidget()
        self.mobileMeshWidget.hideWidget()
        self.radiativeWidget.hideWidget()
        self.electricalWidget.hideWidget()
        self.externalHeadLossesWidget.hideWidget()
        self.hydraulicheadWidget.hideWidget()
        self.pressureWidget.hideWidget()
        self.convectiveInletWidget.hideWidget()
        self.mappedInletWidget.hideWidget()
        self.particlesWidget.hideWidget()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
