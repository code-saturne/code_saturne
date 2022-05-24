# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2015 EDF S.A.
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
This module defines the values of reference.

This module contains the following classes and function:
- GroundwaterView
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

from code_saturne.gui.base.Toolbox           import GuiParam
from code_saturne.gui.base.QtPage            import ComboModel, from_qvariant, DoubleValidator
from code_saturne.gui.case.GroundwaterForm  import Ui_GroundwaterForm
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("GroundwaterView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class GroundwaterView(QWidget, Ui_GroundwaterForm):
    """
    Class to open Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_GroundwaterForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = GroundwaterModel(self.case)

        self.list_scalars = []
        self.m_sca = DefineUserScalarsModel(self.case)
        for s in self.m_sca.getUserScalarNameList():
            self.list_scalars.append((s, self.tr("Additional scalar")))

        # ComboBox
        self.modelPermeability = ComboModel(self.comboBoxPermeability,2,1)
        self.modelFlowType = ComboModel(self.comboBoxFlowType,2,1)
        self.modelUnsaturated = ComboModel(self.comboBoxUnsaturated,2,1)
        self.modelChemistryModel = ComboModel(self.comboBoxChemistryModel,2,1)
        self.modelSpeciesName = ComboModel(self.comboBoxSpeciesName,1,1)

        self.modelPermeability.addItem(self.tr("isotropic"), 'isotropic')
        self.modelPermeability.addItem(self.tr("anisotropic"), 'anisotropic')
        self.modelFlowType.addItem(self.tr("steady"), 'steady')
        self.modelFlowType.addItem(self.tr("unsteady"), 'unsteady')
        self.modelUnsaturated.addItem(self.tr("True"), 'true')
        self.modelUnsaturated.addItem(self.tr("False"), 'false')
        self.modelChemistryModel.addItem(self.tr("Kd"), 'Kd')
        self.modelChemistryModel.addItem(self.tr("EK"), 'EK')

        self.scalar = ""
        scalar_list = self.m_sca.getUserScalarNameList()
        for s in self.m_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)

        if scalar_list != []:
            self.scalar = scalar_list[0]
            for scalar in scalar_list:
                self.modelSpeciesName.addItem(scalar)
        else:
            self.groupBoxTransport.hide()

        # Set up validators
        self.lineEditDecayRate.setValidator(DoubleValidator(self.lineEditDecayRate))

        # Connections
        self.comboBoxPermeability.activated[str].connect(self.slotPermeabilityType)
        self.comboBoxFlowType.activated[str].connect(self.slotFlowType)
        self.comboBoxUnsaturated.activated[str].connect(self.slotUnsaturated)
        self.checkBoxGravity.clicked.connect(self.slotGravity)
        self.comboBoxChemistryModel.activated[str].connect(self.slotChemistryModel)
        self.comboBoxSpeciesName.activated[str].connect(self.slotSpeciesName)
        self.lineEditDecayRate.textChanged[str].connect(self.slotDecayRate)

        self.initializeWidget(scalar_list)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def initializeWidget(self, scalar_list):
        """
        """
        value = self.mdl.getPermeabilityType()
        self.modelPermeability.setItem(str_model=value)

        value = self.mdl.getFlowType()
        self.modelFlowType.setItem(str_model=value)

        value = self.mdl.getUnsaturatedZone()
        self.modelUnsaturated.setItem(str_model=value)

        if self.mdl.getGravity() == 'on':
            self.checkBoxGravity.setChecked(True)
        else:
            self.checkBoxGravity.setChecked(False)

        if scalar_list != []:
            value = self.mdl.getDecayRate(self.scalar)
            self.lineEditDecayRate.setText(str(value))
            value = self.mdl.getChemistryModel(self.scalar)
            self.modelChemistryModel.setItem(str_model=value)


    @pyqtSlot(str)
    def slotPermeabilityType(self, text):
        """
        Input permeability type : isotrop or anisotrop.
        """
        mdl = self.modelPermeability.dicoV2M[str(text)]
        self.mdl.setPermeabilityType(mdl)


    @pyqtSlot(str)
    def slotFlowType(self, text):
        """
        Input flow type : steady or unsteady.
        """
        mdl = self.modelFlowType.dicoV2M[str(text)]
        self.mdl.setFlowType(mdl)


    @pyqtSlot(str)
    def slotUnsaturated(self, text):
        """
        Input flow type : steady or unsteady.
        """
        mdl = self.modelUnsaturated.dicoV2M[str(text)]
        self.mdl.setUnsaturatedZone(mdl)


    @pyqtSlot()
    def slotGravity(self):
        """
        Input if gravity is taken into account or not
        """
        if self.checkBoxGravity.isChecked():
            self.mdl.setGravity('on')
        else:
            self.mdl.setGravity('off')


    @pyqtSlot(str)
    def slotSpeciesName(self, text):
        """
        Method to choose the scalar which properties shall be changed
        """
        mdl = self.modelSpeciesName.dicoV2M[str(text)]
        self.scalar = str(text)
        scal = self.scalar
        value = self.mdl.getDecayRate(scal)
        self.lineEditDecayRate.setText(str(value))
        value = self.mdl.getChemistryModel(scal)
        self.modelChemistryModel.setItem(str_model=value)


    @pyqtSlot(str)
    def slotChemistryModel(self, text):
        """
        Input chemistry model for soil-water partition : Kd or EK.
        """
        choice = self.modelChemistryModel.dicoV2M[str(text)]
        scal = self.scalar
        self.mdl.setChemistryModel(scal, choice)


    @pyqtSlot(str)
    def slotDecayRate(self, text):
        """
        """
        if self.lineEditDecayRate.validator().state == QValidator.Acceptable:
            val = float(text)
            scal = self.scalar
            self.mdl.setDecayRate(scal, val)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
