# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This module contains the following class:
- NumericalParamGlobalView
"""

#-------------------------------------------------------------------------------
# Library modules import
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

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator, from_qvariant
from code_saturne.gui.case.NumericalParamGlobalForm import Ui_NumericalParamGlobalForm
from code_saturne.model.NumericalParamGlobalModel import NumericalParamGlobalModel
from code_saturne.model.GroundwaterModel import GroundwaterModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("NumericalParamGlobalView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class NumericalParamGlobalView(QWidget, Ui_NumericalParamGlobalForm):
    """
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_NumericalParamGlobalForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = NumericalParamGlobalModel(self.case)
        self.browser = tree

        self.labelSRROM.hide()
        self.lineEditSRROM.hide()

        # Combo models
        self.modelGradientType = ComboModel(self.comboBoxGradientType, 5, 1)
        self.modelGradientType.addItem(self.tr("Automatic"), 'default')
        self.modelGradientType.addItem(self.tr("Green-Gauss with iterative face values"),
                                       'green_iter')
        self.modelGradientType.addItem(self.tr("Least squares"), 'lsq')
        self.modelGradientType.addItem(self.tr("Green-Gauss with least squares gradient face values"),
                                       'green_lsq')
        self.modelGradientType.addItem(self.tr("Green-Gauss with vertex interpolated face values"),
                                       'green_vtx')

        self.modelExtNeighbors = ComboModel(self.comboBoxExtNeighbors, 8, 1)
        self.modelExtNeighbors.addItem(self.tr("Automatic"), 'default')
        self.modelExtNeighbors.addItem(self.tr("None (face adjacent only)"), 'none')
        self.modelExtNeighbors.addItem(self.tr("Boundary only"), 'boundary')
        self.modelExtNeighbors.addItem(self.tr("Optimized (heuristics)"), 'optimized')
        self.modelExtNeighbors.addItem(self.tr("Optimized (heuristics), complete on boundary"),
                                       'optimized_with_boundary')
        self.modelExtNeighbors.addItem(self.tr("Opposite adjacent cell centers"),
                                       'cell_center_opposite')
        self.modelExtNeighbors.addItem(self.tr("Opposite adjacent cell centers, complete on boundary"),
                                       'cell_center_opposite_with_boundary')
        self.modelExtNeighbors.addItem(self.tr("Full (all vertex adjacent)"), 'complete')

        if str(self.model.getExtendedNeighborType()) == 'non_ortho_max':
            self.modelExtNeighbors.addItem(self.tr("Non-orthogonal faces threshold (legacy)"),
                                           'non_ortho_max')

        self.modelDensityVar = ComboModel(self.comboBoxDensityVar, 6, 1)
        self.modelDensityVar.addItem(self.tr("Automatic"), 'default')
        self.modelDensityVar.addItem(self.tr("Boussinesq approximation (rho constant except in the buoyant term)"),
                                       'boussi')
        self.modelDensityVar.addItem(self.tr("Dilatable steady algorithm"), 'dilat_std')
        self.modelDensityVar.addItem(self.tr("Dilatable unsteady algorithm"),
                                       'dilat_unstd')
        self.modelDensityVar.addItem(self.tr("Low-Mach algorithm"),
                                       'low_mach')
        self.modelDensityVar.addItem(self.tr("Algorithm for fire"),
                                       'algo_fire')

        # Connections
        self.checkBoxIVISSE.clicked.connect(self.slotIVISSE)

        self.checkBoxIPUCOU.clicked.connect(self.slotIPUCOU)
        self.checkBoxImprovedPressure.clicked.connect(self.slotImprovedPressure)
        self.checkBoxICFGRP.clicked.connect(self.slotICFGRP)
        self.lineEditRELAXP.textChanged[str].connect(self.slotRELAXP)
        self.comboBoxGradientType.activated[str].connect(self.slotGradientType)
        self.comboBoxExtNeighbors.activated[str].connect(self.slotExtNeighbors)
        self.lineEditSRROM.textChanged[str].connect(self.slotSRROM)
        self.comboBoxDensityVar.activated[str].connect(self.slotDensityVar)

        # Validators
        validatorRELAXP = DoubleValidator(self.lineEditRELAXP, min=0., max=1.)
        validatorRELAXP.setExclusiveMin(True)
        validatorSRROM = DoubleValidator(self.lineEditSRROM, min=0., max=1.)
        validatorSRROM.setExclusiveMax(True)
        self.lineEditRELAXP.setValidator(validatorRELAXP)
        self.lineEditSRROM.setValidator(validatorSRROM)

        if self.model.getTransposedGradient() == 'on':
            self.checkBoxIVISSE.setChecked(True)
        else:
            self.checkBoxIVISSE.setChecked(False)

        if self.model.getVelocityPressureCoupling() == 'on':
            self.checkBoxIPUCOU.setChecked(True)
        else:
            self.checkBoxIPUCOU.setChecked(False)

        from code_saturne.model.FluidCharacteristicsModel import FluidCharacteristicsModel
        fluid = FluidCharacteristicsModel(self.case)
        modl_atmo, modl_joul, modl_thermo, modl_gas, modl_coal, modl_comp, modl_hgn = \
            fluid.getThermoPhysicalModel()
        modl_gwf  = GroundwaterModel(self.case).getGroundwaterModel()

        if self.model.getHydrostaticPressure() == 'on':
            self.checkBoxImprovedPressure.setChecked(True)
        else:
            self.checkBoxImprovedPressure.setChecked(False)

        self.lineEditRELAXP.setText(str(self.model.getPressureRelaxation()))
        self.modelGradientType.setItem(str_model=str(self.model.getGradientReconstruction()))
        self.modelExtNeighbors.setItem(str_model=str(self.model.getExtendedNeighborType()))
        self.modelDensityVar.setItem(str_model=str(self.model.getDensityVar()))

        if self.model.getGradientReconstruction() == 'green_iter':
            self.labelExtNeighbors.hide()
            self.comboBoxExtNeighbors.hide()

        from code_saturne.model.TimeStepModel import TimeStepModel
        idtvar = TimeStepModel(self.case).getTimePassing()

        if modl_joul != 'off' or modl_gas != 'off' or modl_coal != 'off':
            self.labelSRROM.show()
            self.lineEditSRROM.show()
            self.lineEditSRROM.setText(str(self.model.getDensityRelaxation()))

        if modl_comp != 'off':
            self.checkBoxICFGRP.show()
            if self.model.getHydrostaticEquilibrium() == 'on':
                self.checkBoxICFGRP.setChecked(True)
            else:
                self.checkBoxICFGRP.setChecked(False)
            self.checkBoxIPUCOU.hide()
            self.lineEditRELAXP.hide()
            self.labelRELAXP.hide()
            self.checkBoxImprovedPressure.hide()
        else:
            self.checkBoxICFGRP.hide()
            if idtvar == -1:
                self.checkBoxIPUCOU.setEnabled(False)
            else:
                self.checkBoxIPUCOU.show()
            self.lineEditRELAXP.show()
            self.labelRELAXP.show()
            self.checkBoxImprovedPressure.show()

        self.modelDensityVar.disableItem(str_model = 'algo_fire')
        if modl_gas != 'off':
            self.modelDensityVar.enableItem(str_model = 'algo_fire')

        # For the moment, the Low Mach algorithm is disabled in the GUI
        self.modelDensityVar.disableItem(str_model = 'low_mach')

        # Hide non relevant info for ground water flows

        if modl_gwf != 'off':
            self.groupBoxTerms.hide()
            self.groupBoxVelPres.hide()
            self.groupBoxDensityVar.hide()

        # Update the Tree files and folders
        self.browser.configureTree(self.case)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotIVISSE(self):
        """
        Set value for parameter IVISSE
        """
        if self.checkBoxIVISSE.isChecked():
            self.model.setTransposedGradient("on")
        else:
            self.model.setTransposedGradient("off")


    @pyqtSlot()
    def slotIPUCOU(self):
        """
        Set value for parameter IPUCOU
        """
        if self.checkBoxIPUCOU.isChecked():
            self.model.setVelocityPressureCoupling("on")
        else:
            self.model.setVelocityPressureCoupling("off")


    @pyqtSlot()
    def slotICFGRP(self):
        """
        Set value for parameter IPUCOU
        """
        if self.checkBoxICFGRP.isChecked():
            self.model.setHydrostaticEquilibrium("on")
        else:
            self.model.setHydrostaticEquilibrium("off")


    @pyqtSlot()
    def slotImprovedPressure(self):
        """
        Input IHYDPR.
        """
        if self.checkBoxImprovedPressure.isChecked():
            self.model.setHydrostaticPressure("on")
        else:
            self.model.setHydrostaticPressure("off")


    @pyqtSlot(str)
    def slotRELAXP(self, text):
        """
        Set value for parameter RELAXP
        """
        if self.lineEditRELAXP.validator().state == QValidator.Acceptable:
            relaxp = from_qvariant(text, float)
            self.model.setPressureRelaxation(relaxp)
            log.debug("slotRELAXP-> %s" % relaxp)


    @pyqtSlot(str)
    def slotSRROM(self, text):
        """
        Set value for parameter SRROM
        """
        if self.lineEditSRROM.validator().state == QValidator.Acceptable:
            srrom = from_qvariant(text, float)
            self.model.setDensityRelaxation(srrom)
            log.debug("slotSRROM-> %s" % srrom)


    @pyqtSlot(str)
    def slotGradientType(self, text):
        """
        Set value for parameter GradientType
        """
        grd_type = self.modelGradientType.dicoV2M[str(text)]
        self.model.setGradientReconstruction(grd_type)
        log.debug("slotGradientType-> %s" % grd_type)

        if grd_type == 'green_iter':
            self.labelExtNeighbors.hide()
            self.comboBoxExtNeighbors.hide()
        else:
            self.labelExtNeighbors.show()
            self.comboBoxExtNeighbors.show()


    @pyqtSlot(str)
    def slotExtNeighbors(self, text):
        """
        Set extended neighborhood type
        """
        enh_type = self.modelExtNeighbors.dicoV2M[str(text)]
        self.model.setExtendedNeighborType(enh_type)
        log.debug("slotExtNeighbors-> %s" % enh_type)

    @pyqtSlot(str)
    def slotDensityVar(self, text):
        """
        Set algorithm for density variation in time
        """
        dv_type = self.modelDensityVar.dicoV2M[str(text)]
        self.model.setDensityVar(dv_type)
        log.debug("slotDensityVar-> %s" % dv_type)

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
