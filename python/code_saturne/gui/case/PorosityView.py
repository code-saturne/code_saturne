# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------

"""
This module defines the Porosity model data management.

This module contains the following classes:
- Porosity
- PorosityView
"""

# -------------------------------------------------------------------------------
# Library modules import
# -------------------------------------------------------------------------------

import sys, logging

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
from code_saturne.gui.base.QtPage import ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.PorosityForm import Ui_PorosityForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.PorosityModel import PorosityModel
from code_saturne.model.NotebookModel import NotebookModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PorosityView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main view class
# -------------------------------------------------------------------------------

class PorosityView(QWidget, Ui_PorosityForm):

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_PorosityForm.__init__(self)
        self.setupUi(self)

        self.case = None
        self.model = None
        self.notebook = None
        self.zone = None

    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()
        self.model = PorosityModel(self.case)
        self.notebook = NotebookModel(self.case)
        localization_model = LocalizationModel("VolumicZone", self.case)
        for zone in localization_model.getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone

        if self.zone.isNatureActivated("porosity"):
            self.setViewFromCase()
        else:  # TODO check if content of tab should still be visible or not
            self.displayDefaultView()
            self.setEnabled(False)

        self.case.undoStartGlobal()

    def setViewFromCase(self):
        # TODO should method be renamed ?
        if self.case.module_name() == 'code_saturne':
            self.modelPorosityType = ComboModel(self.comboBoxType, 3, 1)
            self.modelPorosityType.addItem(self.tr("isotropic"), 'isotropic')
            self.modelPorosityType.addItem(self.tr("anisotropic"), 'anisotropic')
            self.modelPorosityType.addItem(self.tr("integral model"), 'integral')
        else:
            self.modelPorosityType = ComboModel(self.comboBoxType, 1, 1)
            self.modelPorosityType.addItem(self.tr("isotropic"), 'isotropic')
            self.modelPorosityType.disableItem(index=0)
        self.comboBoxType.activated[str].connect(self.slotPorosity)
        self.pushButtonPorosity.clicked.connect(self.slotFormulaPorosity)
        self.selectPorosityZones()

    def selectPorosityZones(self):
        zone_label = self.zone.getLabel()
        zone_id = self.zone.getCodeNumber()

        if hasattr(self, "modelScalars"): del self.modelScalars
        log.debug("slotSelectPorosityZones label %s " % zone_label)
        self.groupBoxType.show()
        self.groupBoxDef.show()

        choice = self.model.getPorosityModel(zone_id)
        self.modelPorosityType.setItem(str_model=choice)

        exp = self.model.getPorosityFormula(zone_id)
        if exp:
            self.pushButtonPorosity.setToolTip(exp)
            self.pushButtonPorosity.setStyleSheet("background-color: green")
        else:
            self.pushButtonPorosity.setStyleSheet("background-color: red")

    @pyqtSlot(str)
    def slotPorosity(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        zone_id = self.zone.getCodeNumber()
        choice = self.modelPorosityType.dicoV2M[str(text)]

        self.model.setPorosityModel(zone_id, choice)

    @pyqtSlot()
    def slotFormulaPorosity(self):
        """
        User formula for density
        """
        zone_label = self.zone.getLabel()
        zone_id = self.zone.getCodeNumber()

        choice = self.model.getPorosityModel(zone_id)
        fname = 'porosity'
        if choice == 'anisotropic':
            fname += '+tensorial_porosity'

        exp, req, sca, sym = self.model.getPorosityFormulaComponents(zone_id)

        if exp is None:
            exp = self.getDefaultPorosityFormula(choice)

        exa = """#example: \n""" + self.model.getDefaultPorosityFormula(choice)

        dialog = QMegEditorView(parent=self,
                                function_type='vol',
                                zone_name=zone_label,
                                variable_name=fname,
                                expression=exp,
                                required=req,
                                symbols=sym,
                                known_fields=sca,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorosity -> %s" % str(result))
            self.model.setPorosityFormula(zone_id, str(result))
            self.pushButtonPorosity.setToolTip(result)
            self.pushButtonPorosity.setStyleSheet("background-color: green")

    def displayDefaultView(self):
        self.modelPorosityType = ComboModel(self.comboBoxType, 1, 1)
        self.modelPorosityType.addItem(self.tr("isotropic"), 'isotropic')
        self.modelPorosityType.setItem(str_model="isotropic")
        self.groupBoxType.show()
        self.groupBoxDef.show()


# -------------------------------------------------------------------------------
# Testing part
# -------------------------------------------------------------------------------


if __name__ == "__main__":
    pass

# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
