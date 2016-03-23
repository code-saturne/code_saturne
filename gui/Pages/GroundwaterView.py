# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox           import GuiParam
from code_saturne.Base.QtPage            import ComboModel, from_qvariant, DoubleValidator
from code_saturne.Pages.GroundwaterForm  import Ui_GroundwaterForm
from code_saturne.Pages.GroundwaterModel import GroundwaterModel

from code_saturne.Pages.QMeiEditorView   import QMeiEditorView

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

        # ComboBox
        self.modelPermeability = ComboModel(self.comboBoxPermeability,2,1)
        self.modelDispersion = ComboModel(self.comboBoxDispersion,2,1)
        self.modelFlowType = ComboModel(self.comboBoxFlowType,2,1)
        self.modelUnsaturated = ComboModel(self.comboBoxUnsaturated,2,1)

        self.modelPermeability.addItem(self.tr("isotropic"), 'isotropic')
        self.modelPermeability.addItem(self.tr("anisotropic"), 'anisotropic')
        self.modelDispersion.addItem(self.tr("isotropic"), 'isotropic')
        self.modelDispersion.addItem(self.tr("anisotropic"), 'anisotropic')
        self.modelFlowType.addItem(self.tr("steady"), 'steady')
        self.modelFlowType.addItem(self.tr("unsteady"), 'unsteady')
        self.modelUnsaturated.addItem(self.tr("True"), 'true')
        self.modelUnsaturated.addItem(self.tr("False"), 'false')

        # Connections
        self.comboBoxPermeability.activated[str].connect(self.slotPermeabilityType)
        self.comboBoxDispersion.activated[str].connect(self.slotDispersionType)
        self.comboBoxFlowType.activated[str].connect(self.slotFlowType)
        self.comboBoxUnsaturated.activated[str].connect(self.slotUnsaturated)
        self.checkBoxGravity.clicked.connect(self.slotGravity)

        self.initializeWidget()

        self.case.undoStartGlobal()


    @pyqtSlot()
    def initializeWidget(self):
        """
        """
        value = self.mdl.getPermeabilityType()
        self.modelPermeability.setItem(str_model=value)

        value = self.mdl.getDispersionType()
        self.modelDispersion.setItem(str_model=value)

        value = self.mdl.getFlowType()
        self.modelFlowType.setItem(str_model=value)

        value = self.mdl.getUnsaturatedZone()
        self.modelUnsaturated.setItem(str_model=value)

        if self.mdl.getGravity() == 'on':
            self.checkBoxGravity.setChecked(True)
        else:
            self.checkBoxGravity.setChecked(False)


    @pyqtSlot(str)
    def slotPermeabilityType(self, text):
        """
        Input permeability type : isotrop or anisotrop.
        """
        mdl = self.modelPermeability.dicoV2M[str(text)]
        self.mdl.setPermeabilityType(mdl)


    @pyqtSlot(str)
    def slotDispersionType(self, text):
        """
        Input viscosity type : isotrop or anisotrop.
        """
        mdl = self.modelDispersion.dicoV2M[str(text)]
        self.mdl.setDispersionType(mdl)


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


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
