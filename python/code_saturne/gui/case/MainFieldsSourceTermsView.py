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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.MainFieldsSourceTerms import Ui_MainFieldsSourceTerms

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.model.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
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

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MainFieldsSourceTerms.__init__(self)
        self.parent = parent
        self.setupUi(self)

        self.case = None
        self.zone = None
        self.zone_id = None
        self.mdl = None
        self.mfm = None
        self.notebook = None
        self.th_sca_name = 'enthalpy'

        self.defineConnections()

    def defineConnections(self):
        self.comboBoxField.activated[str].connect(self.slotField)
        self.pushButtonThermal.clicked.connect(self.slotThermalFormula)

    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()
        self.mdl = MainFieldsSourceTermsModel(self.case)
        self.mfm = MainFieldsModel(self.case)
        self.notebook = NotebookModel(self.case)
        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_id = str(zone.getCodeNumber())
        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for fieldId in self.mfm.getFieldIdList():
            label = self.mfm.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)
        self.currentId = -1
        if len(self.mfm.getFieldIdList()) > 0:
            self.currentId = self.mfm.getFieldIdList()[0]
            self.modelField.setItem(str_model=self.currentId)

        if self.zone.isNatureActivated("source_term"):
            self.setViewFromCase()
        else:
            self.setEnabled(False)
        self.case.undoStartGlobal()

    def setViewFromCase(self):
        """
        Initialize widget when a new volumic zone_info is chosen
        """
        if self.isSourceTermActived("thermal"):
            self.pushButtonThermal.show()
            self.th_sca_name = 'enthalpy'
            exp = self.mdl.getThermalFormula(self.zone_id,
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

    @pyqtSlot()
    def slotThermalFormula(self):
        """
        Input the initial formula of thermal scalar
        """
        exa = """#example: """

        exp, req, sym, knf = self.mdl.getThermalFormulaComponents(self.zone_id,
                                                                  self.currentId,
                                                                  self.th_sca_name)

        name = 'enthalpy_%s' % (str(self.currentId))
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'src',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = knf,
                                examples      = exa,
                                source_type   = 'thermal_source_term')

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaThermal -> %s" % str(result))
            self.mdl.setThermalFormula(self.zone_id,
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

        exp = self.mdl.getThermalFormula(self.zone_id,
                                         self.currentId,
                                         self.th_sca_name)
        if exp:
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(exp)
        else:
            self.pushButtonThermal.setStyleSheet("background-color: red")

    def isSourceTermActived(self, nature: str) -> bool:
        source_term = nature + "_source_term"
        zone_info = self.zone.getNature()
        if source_term in zone_info:
            return zone_info[source_term] == "on"
        else:
            return False

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
