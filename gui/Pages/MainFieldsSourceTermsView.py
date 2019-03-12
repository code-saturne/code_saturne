# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
- MainFieldsSourceTermsView
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

from code_saturne.Pages.MainFieldsSourceTerms import Ui_MainFieldsSourceTerms

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.model.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.MainFieldsSourceTermsModel import MainFieldsSourceTermsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MainFieldsSourceTermsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MainFieldsSourceTermsView(QWidget, Ui_MainFieldsSourceTerms):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MainFieldsSourceTerms.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.parent = parent

        self.mdl      = MainFieldsSourceTermsModel(self.case)
        self.mfm      = MainFieldsModel(self.case)
        self.notebook = NotebookModel(self.case)
        self.volzone  = LocalizationModel('VolumicZone', self.case)

        # 0/ Read label names from XML file

        # Velocity

        # Thermal scalar
        self.th_sca_name = 'enthalpy'

        # 1/ Combo box models
        self.modelZone     = ComboModel(self.comboBoxZone, 1, 1)

        self.zone = ""
        zones = self.volzone.getZones()
        for zone in zones:
            active = 0
            if ('thermal_source_term' in zone.getNature().keys()):
                if (zone.getNature()['thermal_source_term']  == "on"):
                    active = 1

            if (active):
                label = zone.getLabel()
                name = str(zone.getCodeNumber())
                self.modelZone.addItem(self.tr(label), name)
                if label == "all_cells":
                    self.zone = name
                if not self.zone:
                    self.zone = name

        self.modelZone.setItem(str_model = self.zone)

        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for fieldId in self.mfm.getFieldIdList() :
            label = self.mfm.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)

        self.currentId = -1
        if len(self.mfm.getFieldIdList()) > 0:
            self.currentId = self.mfm.getFieldIdList()[0]
            self.modelField.setItem(str_model = self.currentId)


        # 2/ Connections
        self.comboBoxZone.activated[str].connect(self.slotZone)
        self.comboBoxField.activated[str].connect(self.slotField)
        self.pushButtonThermal.clicked.connect(self.slotThermalFormula)

        # 3/ Initialize widget

        self.initialize(self.zone)

        self.case.undoStartGlobal()


    def initialize(self, zone_num):
        """
        Initialize widget when a new volumic zone is chosen
        """
        zone = self.case.xmlGetNode("zone", id=zone_num)

        if zone['thermal_source_term']  == "on":
            self.pushButtonThermal.show()
            self.th_sca_name = 'enthalpy'
            exp = self.mdl.getThermalFormula(self.zone,
                                             self.currentId,
                                             self.th_sca_name)
            if exp:
                self.pushButtonThermal.setToolTip(exp)
                self.pushButtonThermal.setStyleSheet("background-color: green")
            else:
                self.pushButtonThermal.setStyleSheet("background-color: red")
        else:
            self.pushButtonThermal.hide()
            self.labelThermal.hide()



    @pyqtSlot(str)
    def slotZone(self, text):
        """
        INPUT label for choice of zone
        """
        self.zone = self.modelZone.dicoV2M[str(text)]
        self.initialize(self.zone)


    @pyqtSlot()
    def slotThermalFormula(self):
        """
        Input the initial formula of thermal scalar
        """
        exa = """#example: """

        exp, req, sym = self.mdl.getThermalFormulaComponents(self.zone,
                                                             self.currentId,
                                                             self.th_sca_name)

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermal -> %s" % str(result))
            self.mdl.setThermalFormula(self.zone,
                                       self.currentId,
                                       self.th_sca_name,
                                       str(result))
            self.pushButtonThermal.setToolTip(result)
            self.pushButtonThermal.setStyleSheet("background-color: green")


    @pyqtSlot(str)
    def slotField(self, text):
        """
        INPUT label for choice of field
        """
        self.currentId = self.modelField.dicoV2M[str(text)]
#        self.initializeVariables(self.zone, self.currentid)

        exp = self.mdl.getThermalFormula(self.zone,
                                         self.currentId,
                                         self.th_sca_name)
        if exp:
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(exp)
        else:
            self.pushButtonThermal.setStyleSheet("background-color: red")


#    def initializeVariables(self, zone, fieldId):
#        """
#        Initialize variables when a new volumic zone or fieldId is choosen
#        """
#        # Thermal Initialization



    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
