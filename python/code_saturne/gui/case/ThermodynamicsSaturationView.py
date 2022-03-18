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
This module define the 'ThermodynamicsSaturationView' page.
This module contains the following classes:
- ThermodynamicsSaturationView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel
from code_saturne.gui.base.QtPage import to_text_string
from code_saturne.gui.case.ThermodynamicsSaturation import Ui_ThermodynamicsSaturation
from code_saturne.model.ThermodynamicsModel import *
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.SpeciesModel import SpeciesModel
from code_saturne.model.OutputFieldsModel import OutputFieldsModel
from code_saturne.model.NonCondensableModel import NonCondensableModel
from code_saturne.model.LocalizationModel import LocalizationModel

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.InterfacialEnthalpyModel import InterfacialEnthalpyModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermodynamicsFieldView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Combo box delegate for the material
#-------------------------------------------------------------------------------

_ok_str  = "background-color: green"
_nok_str = "background-color: red"

#-------------------------------------------------------------------------------
# MainFieldsView class
#-------------------------------------------------------------------------------

class ThermodynamicsSaturationView(QWidget, Ui_ThermodynamicsSaturation):
    """
    Thermodynamics layout.
    """

    def __init__(self, parent = None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ThermodynamicsSaturation.__init__(self)
        self.setupUi(self)

        self.case      = None
        self.mdl       = None
        self.notebook  = None
        self.zone      = None
        self.zone_name = None
        self.zone_id   = None

    def setup(self,case, zone_name):
        self.case = case
        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_name = zone.getLabel()
                self.zone_id   = zone.getCodeNumber()
        self.case.undoStopGlobal()

        self.mdl      = ThermodynamicsModel(self.case)
        self.notebook = NotebookModel(self.case)
        self.ncond    = NonCondensableModel(self.case)
        self.interf   = InterfacialEnthalpyModel(self.case)
        self.mfm      = MainFieldsModel(self.case)

        # Set fields string:
        field_ids = self.interf.getEnthalpyCoupleFieldId()
        couple_str = ""
        if field_ids:
            couple_str = self.mfm.getLabel(field_ids[0])
            couple_str += " (liquid)/"
            couple_str += self.mfm.getLabel(field_ids[1])
            couple_str += " (gas)"

        self.lineEditCoupledFields.setText(couple_str)
        self.lineEditCoupledFields.setEnabled(False)

        # Check if user properties:
        self.user_method = False
        if field_ids:
            if self.mdl.getMethod(field_ids[0]) == "user_properties" \
                    and self.mdl.getMethod(field_ids[1]) == "user_properties":
                self.user_method = True

        is_main_zone = (zone_name == "all_cells")
        # Dico
        self.m_out = OutputFieldsModel(self.case)
        self.currentFluid = 0

        # Connections
        self.sat_ppts = {"SaturationTemperature":"Tsat",
                         "d_Tsat_d_P":"dTsatdp",
                         "LatentHeat":"Hlat",
                         "SaturationEnthalpyLiquid":"HsatL",
                         "SaturationEnthalpyGas":"HsatG",
                         "d_Hsat_d_P_Liquid":"dHsatLdp",
                         "d_Hsat_d_P_Gas":"dHsatGdp"}


        for key in self.sat_ppts.keys():
            _k = self.sat_ppts[key]
            _button  = getattr(self, "pushButton" + _k)

            # WARNING: when connecting a pushButton, a 'state' signal
            # is emitted. If you do not provide state to lambda,
            # input will be False/True instead of the variable
            _button.clicked.connect(lambda state, name=key: self.slotFormula(name))


        # load Field

        self.initializeWidget()

        self.case.undoStartGlobal()

    def initializeWidget(self):

        # hide groupBoxConstantProperties
        self.groupBoxEauvap.setVisible(self.user_method)

        if self.user_method:
            for key in self.sat_ppts.keys():
                _k = self.sat_ppts[key]
                _button = getattr(self, "pushButton" + _k)

                exp = self.mdl.getFormula("none", key, self.zone_id)
                if exp:
                    _button.setStyleSheet(_ok_str)
                    _button.setToolTip(exp)
                else:
                    _button.setStyleSheet(_nok_str)

    @pyqtSlot()
    def slotFormula(self, name):
        """
        User formula for a given saturation property.
        """

        exp, req, sca, symbols = self.mdl.getFormulaComponents('none', name, zone=self.zone_id)

        if not exp:
            exp = self.mdl.getDefaultFormula('none', name)

        exa = self.mdl.getExampleFormula('none', name)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = symbols,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            _k = self.sat_ppts[name]
            _button  = getattr(self, "pushButton" + _k)
            log.debug("slotFormula%s -> %s" % (_k, str(result)))
            self.mdl.setFormula('none', name, result, zone=self.zone_id)
            _button.setStyleSheet(_ok_str)
            _button.setToolTip(result)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

