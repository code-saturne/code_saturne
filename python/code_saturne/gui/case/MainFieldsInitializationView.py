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
This module defines the 'Main fields initialization' page.

This module contains the following classes:
- MainFieldsInitializationView
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator
from MainFieldsInitialization import Ui_MainFieldsInitialization
from code_saturne.model.MainFieldsInitializationModel import MainFieldsInitializationModel
from code_saturne.model.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from code_saturne.model.NonCondensableModel import NonCondensableModel
from code_saturne.model.SpeciesModel import SpeciesModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel

from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MainFieldsInitializationView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# MainFieldsInitialization class
#-------------------------------------------------------------------------------

class MainFieldsInitializationView(QWidget, Ui_MainFieldsInitialization):
    """
    Main fields layout.
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MainFieldsInitialization.__init__(self)
        self.setupUi(self)

        self.case = None
        self.mdl = None
        self.volzone = None
        self.NonCondensable = None
        self.SpeciesModel = None
        self.zone = None
        self.zone_id = None

        self.hideAllWidgets()
        self.defineConnections()

    def defineConnections(self):
        self.comboBoxField.activated[str].connect(self.slotField)
        self.comboBoxEnergy.activated[str].connect(self.slotEnergyModel)
        self.comboBoxNonCondensable.activated[str].connect(self.slotNonCondensableType)
        self.comboBoxScalar.activated[str].connect(self.slotScalarName)
        self.pushButtonPressure.clicked.connect(self.slotPressure)
        self.pushButtonVelocity.clicked.connect(self.slotVelocity)
        self.pushButtonFraction.clicked.connect(self.slotFraction)
        self.pushButtonEnergy.clicked.connect(self.slotEnergy)
        self.pushButtonNonCondensable.clicked.connect(self.slotNonCondensable)
        self.pushButtonScalar.clicked.connect(self.slotScalar)

    def hideAllWidgets(self):
        self.labelEnergy.hide()
        self.comboBoxEnergy.hide()
        self.pushButtonEnergy.hide()
        self.labelNonCondensable.hide()
        self.comboBoxNonCondensable.hide()
        self.pushButtonNonCondensable.hide()
        self.labelScalar.hide()
        self.comboBoxScalar.hide()
        self.pushButtonScalar.hide()

    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()

        self.mdl = MainFieldsInitializationModel(self.case)
        self.volzone = LocalizationModel('VolumicZone', self.case)
        self.NonCondensable = NonCondensableModel(self.case)
        self.SpeciesModel = SpeciesModel(self.case)

        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_id = str(zone.getCodeNumber())

        self.modelField = ComboModel(self.comboBoxField, 1, 1)
        for fieldId in self.mdl.getFieldIdList():
            label = self.mdl.getLabel(fieldId)
            name = str(fieldId)
            self.modelField.addItem(self.tr(label), name)

        need_none = (SpeciesModel(self.case).getScalarByFieldId("none")!=[])
        if need_none:
            self.modelField.addItem(self.tr('Non-convected scalars'), 'none')

        self.currentid = -1
        if len(self.mdl.getFieldIdList()) > 0:
            self.currentid = self.mdl.getFieldIdList()[0]
            self.modelField.setItem(str_model = self.currentid)
        self.modelEnergy = ComboModel(self.comboBoxEnergy, 3, 1)
        self.modelEnergy.addItem(self.tr("Enthalpy"), "enthalpy")
        self.modelEnergy.addItem(self.tr("Temperature"), "temperature")
        self.modelEnergy.addItem(self.tr("Saturation enthalpy"), "hsat_P")

        if int(self.currentid) > 0:
            if ThermodynamicsModel(self.case).getMaterials(self.currentid) == 'user_material' :
                self.modelEnergy.disableItem(1)
                self.modelEnergy.disableItem(2)

        self.modelNonCondensable = ComboModel(self.comboBoxNonCondensable, 1, 1)
        self.currentNonCond = ""
        self.currentNonCondLabel = ""

        self.modelScalar = ComboModel(self.comboBoxScalar, 1, 1)
        self.currentScalar = ""
        self.currentScalarLabel = ""

        if self.zone.isNatureActivated("initialization"):
            self.setViewFromCase()
        else:
            self.setEnabled(False)
        self.case.undoStartGlobal()

    def setViewFromCase(self):
        exp = self.mdl.getFormulaPressure(self.zone_id)
        if exp:
            self.pushButtonPressure.setStyleSheet("background-color: green")
            self.pushButtonPressure.setToolTip(exp)
        else:
            self.pushButtonPressure.setStyleSheet("background-color: red")
        if (len(self.mdl.getFieldIdList()) > 0):
            self.groupBoxDefinition.show()
            self.initializeVariables(self.zone_id, self.currentid)

            exp = self.mdl.getFormula(self.zone_id, self.currentid, 'velocity')
            if exp:
                self.pushButtonVelocity.setStyleSheet("background-color: green")
                self.pushButtonVelocity.setToolTip(exp)
            else:
                self.pushButtonVelocity.setStyleSheet("background-color: red")
            exp = self.mdl.getFormula(self.zone_id, self.currentid, 'volume_fraction')
            if exp:
                self.pushButtonFraction.setStyleSheet("background-color: green")
                self.pushButtonFraction.setToolTip(exp)
            else:
                self.pushButtonFraction.setStyleSheet("background-color: red")

            if self.mdl.getEnergyResolution(self.currentid) == "on":
                exp = self.mdl.getFormula(self.zone_id, self.currentid, 'enthalpy')
                if exp:
                    self.pushButtonEnergy.setStyleSheet("background-color: green")
                    self.pushButtonEnergy.setToolTip(exp)
                else:
                    self.pushButtonEnergy.setStyleSheet("background-color: red")

            lst = self.NonCondensable.getNonCondensableByFieldId(self.currentid)
            if len(lst) > 0:
                exp = self.mdl.getFormulaNonCondensable(self.zone_id, self.currentid, self.currentNonCond)
                if exp:
                    self.pushButtonNonCondensable.setStyleSheet("background-color: green")
                    self.pushButtonNonCondensable.setToolTip(exp)
                else:
                    self.pushButtonNonCondensable.setStyleSheet("background-color: red")

            lst = self.SpeciesModel.getScalarByFieldId(self.currentid)
            if len(lst) > 0:
                exp = self.mdl.getFormulaScalar(self.zone_id, self.currentid, self.currentScalar)
                if exp:
                    self.pushButtonScalar.setStyleSheet("background-color: green")
                    self.pushButtonScalar.setToolTip(exp)
                else:
                    self.pushButtonScalar.setStyleSheet("background-color: red")
        else:
            self.groupBoxDefinition.hide()

    @pyqtSlot(str)
    def slotField(self, text):
        """
        INPUT label for choice of field
        """
        self.currentid = self.modelField.dicoV2M[str(text)]
        self.initializeVariables(self.zone_id, self.currentid)

        if self.currentid != 'none':
            # Velocity
            exp = self.mdl.getFormula(self.zone_id,
                                      self.currentid,
                                      'velocity')

            if exp:
                self.pushButtonVelocity.setStyleSheet("background-color: green")
                self.pushButtonVelocity.setToolTip(exp)
            else:
                self.pushButtonVelocity.setStyleSheet("background-color: red")

            # Volume Fraction
            exp = self.mdl.getFormula(self.zone_id,
                                      self.currentid,
                                      'volume_fraction')
            if exp:
                self.pushButtonFraction.setStyleSheet("background-color: green")
                self.pushButtonFraction.setToolTip(exp)
            else:
                self.pushButtonFraction.setStyleSheet("background-color: red")

            # Energy
            if self.mdl.getEnergyResolution(self.currentid) == "on":
                exp = self.mdl.getFormula(self.zone_id,
                                          self.currentid,
                                          'enthalpy')
                if exp:
                    self.pushButtonEnergy.setStyleSheet("background-color: green")
                    self.pushButtonEnergy.setToolTip(exp)
                else:
                    self.pushButtonEnergy.setStyleSheet("background-color: red")

            # Non condensable gases
            lst = self.NonCondensable.getNonCondensableByFieldId(self.currentid)
            if len(lst) > 0 :
                exp = self.mdl.getFormulaNonCondensable(self.zone_id,
                                                        self.currentid,
                                                        self.currentNonCond)
                if exp:
                    self.pushButtonNonCondensable.setStyleSheet("background-color: green")
                    self.pushButtonNonCondensable.setToolTip(exp)
                else:
                    self.pushButtonNonCondensable.setStyleSheet("background-color: red")

        # Scalars (can exist for 'none')
        lst = self.SpeciesModel.getScalarByFieldId(self.currentid)
        if len(lst) > 0 :
            exp = self.mdl.getFormulaScalar(self.zone_id,
                                            self.currentid,
                                            self.currentScalar)
            if exp:
                self.pushButtonScalar.setStyleSheet("background-color: green")
                self.pushButtonScalar.setToolTip(exp)
            else:
                self.pushButtonScalar.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def slotEnergyModel(self, text):
        """
        INPUT label for choice of energy model
        """
        model = self.modelEnergy.dicoV2M[str(text)]
        self.mdl.setEnergyModel(self.zone_id, self.currentid, model)

        if model != "hsat_P":
            self.pushButtonEnergy.setEnabled(True)
            exp = self.mdl.getFormula(self.zone_id, self.currentid, 'enthalpy')
            if exp:
                if (model == 'enthalpy' and 'temperature' in exp):
                    self.pushButtonEnergy.setStyleSheet("background-color: red")
                elif (model == 'temperature' and 'enthalpy' in exp):
                    self.pushButtonEnergy.setStyleSheet("background-color: red")
                else:
                    self.pushButtonEnergy.setStyleSheet("background-color: green")
                    self.pushButtonEnergy.setToolTip(exp)
            else:
                self.pushButtonEnergy.setStyleSheet("background-color: red")
        else:
            self.pushButtonEnergy.setEnabled(False)
            self.pushButtonEnergy.setStyleSheet("background-color: None")


    @pyqtSlot(str)
    def slotNonCondensableType(self, text):
        """
        INPUT label for choice of non condensable model
        """
        self.currentNonCond = self.modelNonCondensable.dicoV2M[str(text)]
        self.currentNonCondLabel = str(text)
        exp = self.mdl.getFormulaNonCondensable(self.zone_id, self.currentid, self.currentNonCond)
        if exp:
            self.pushButtonNonCondensable.setStyleSheet("background-color: green")
            self.pushButtonNonCondensable.setToolTip(exp)
        else:
            self.pushButtonNonCondensable.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def slotScalarName(self, text):
        """
        INPUT label for choice of scalar
        """
        self.currentScalar = self.modelScalar.dicoV2M[str(text)]
        self.currentScalarLabel = str(text)
        exp = self.mdl.getFormulaScalar(self.zone_id, self.currentid, self.currentScalar)
        if exp:
            self.pushButtonScalar.setStyleSheet("background-color: green")
            self.pushButtonScalar.setToolTip(exp)
        else:
            self.pushButtonScalar.setStyleSheet("background-color: red")


    @pyqtSlot()
    def slotVelocity(self):
        """
        Formula for velocity
        """
        exp, req, sym = self.mdl.getFormulaComponents(self.zone_id,
                                                      self.currentid,
                                                      'velocity')
        if not exp:
            exp = "u = 0.;\nv = 0.;\nw = 0.;"

        exa = "u = 3.0;\nv = 1.0;\nw = 0.0;\n"

        name = 'velocity_%s' % (str(self.currentid))
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula(self.zone_id, self.currentid, 'velocity', result)
            self.pushButtonVelocity.setStyleSheet("background-color: green")
            self.pushButtonVelocity.setToolTip(result)


    @pyqtSlot()
    def slotFraction(self):
        """
        Formula for fraction
        """
        exp, req, sym = self.mdl.getFormulaComponents(self.zone_id,
                                                      self.currentid,
                                                      'volume_fraction')
        if not exp:
            if self.currentid == "1":
                exp = "vol_f = 1.;\n"
            else:
                exp = "vol_f = 0.;\n"

        exa = "vol_f = 1.0;\n"

        name = 'volume_fraction_%s' % (str(self.currentid))
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula(self.zone_id, self.currentid, 'volume_fraction', result)
            self.pushButtonFraction.setStyleSheet("background-color: green")
            self.pushButtonFraction.setToolTip(result)


    @pyqtSlot()
    def slotEnergy(self):
        """
        Formula for energy
        """
        exp, req, sym = self.mdl.getFormulaComponents(self.zone_id,
                                                      self.currentid,
                                                      'enthalpy')

        th_sca_label = self.mdl.getEnergyModel(self.zone_id, self.currentid)
        if not exp:
            if str(th_sca_label) == 'enthalpy':
                exp = th_sca_label + """ = 50000.0;\n"""

            elif str(th_sca_label) == 'temperature':
                exp = th_sca_label + """ = 293.15;\n"""

        elif ('enthalpy' in exp and str(th_sca_label) == 'temperature'):
            exp = th_sca_label + """ = 293.15;\n"""

        elif ('temperature' in exp and str(th_sca_label) == 'enthalpy'):
            exp = th_sca_label + """ = 50000.0;\n"""


        exa = """#example: """

        name = 'enthalpy_%s' % (str(self.currentid))
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula(self.zone_id, self.currentid, 'enthalpy', result)
            self.pushButtonEnergy.setStyleSheet("background-color: green")
            self.pushButtonEnergy.setToolTip(result)


    @pyqtSlot()
    def slotPressure(self):
        """
        Formula for pressure
        """
        exp, req, sym = self.mdl.getPressureFormulaComponents(self.zone_id)

        if not exp:
            exp = """pressure = 101325.;\n"""

        exa = """#example :
rho0 = 1.8;
zmax = 1.2;
Po = 1.e6;
pressure = P0 + rho0 * g * (zmax - z);"""

        name = 'pressure'
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormulaPressure(self.zone_id, result)
            self.pushButtonPressure.setStyleSheet("background-color: green")
            self.pushButtonPressure.setToolTip(result)


    @pyqtSlot()
    def slotNonCondensable(self):
        """
        Formula for non condensable
        """
        exp, req, sym = self.mdl.getNonCondensableFormulaComponents(self.zone_id,
                                                                    self.currentid,
                                                                    self.currentNonCond)
        if not exp:
            exp = self.currentNonCondLabel + """ = 0;\n"""

        exa = """#example: """

        name = self.currentNonCondLabel
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormulaNonCondensable(self.zone_id, self.currentid, self.currentNonCond, result)
            self.pushButtonNonCondensable.setStyleSheet("background-color: green")
            self.pushButtonNonCondensable.setToolTip(result)


    @pyqtSlot()
    def slotScalar(self):
        """
        Formula for species
        """
        exp, req, sym = self.mdl.getScalarFormulaComponents(self.zone_id,
                                                            self.currentid,
                                                            self.currentScalar)
        if not exp:
            exp = self.currentScalarLabel + """ = 0;\n"""

        exa = """#example: """

        name = self.currentScalarLabel
        zone_name = self.zone.getLabel()

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ini',
                                zone_name     = zone_name,
                                variable_name = name,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormulaScalar(self.zone_id, self.currentid, self.currentScalar, result)
            self.pushButtonScalar.setStyleSheet("background-color: green")
            self.pushButtonScalar.setToolTip(result)


    def initializeVariables(self, zone, fieldId):
        """
        Initialize variables when a new volumic zone or fieldId is choosen
        """
        # Energy initialization
        self.labelEnergy.hide()
        self.comboBoxEnergy.hide()
        self.pushButtonEnergy.hide()

        self.labelVelocity.setVisible(fieldId!='none')
        self.pushButtonVelocity.setVisible(fieldId!='none')
        self.labelFraction.setVisible(fieldId!='none')
        self.pushButtonFraction.setVisible(fieldId!='none')

        if fieldId != 'none' and self.mdl.getEnergyResolution(fieldId) == "on":
            self.labelEnergy.show()
            self.comboBoxEnergy.show()
            self.pushButtonEnergy.show()

            model = self.mdl.getEnergyModel(zone, fieldId)
            self.modelEnergy.setItem(str_model=model)

            if model == "enthalpy" or model == "temperature":
                self.pushButtonEnergy.setEnabled(True)
            elif model == "hsat_P":
                self.pushButtonEnergy.setEnabled(False)
                self.pushButtonEnergy.setStyleSheet("background-color: None")

            if ThermodynamicsModel(self.case).getMaterials(self.currentid) == 'user_material' :
                self.modelEnergy.disableItem(1)
                self.modelEnergy.disableItem(2)
            else :
                self.modelEnergy.enableItem(1)
                self.modelEnergy.enableItem(2)

        # Non-condensable initialization
        if fieldId != 'none':
            lst = self.NonCondensable.getNonCondensableByFieldId(fieldId)
            if len(lst) > 0 :
                self.labelNonCondensable.show()
                self.comboBoxNonCondensable.show()
                self.pushButtonNonCondensable.show()

                if len(self.modelNonCondensable.getItems()) != 0 :
                    for nb in range(len(self.modelNonCondensable.getItems())):
                        self.modelNonCondensable.delItem(0)

                for var in lst :
                    label = self.NonCondensable.getNonCondLabel(var)
                    self.modelNonCondensable.addItem(self.tr(label), var)

                self.currentNonCond = lst[0]
                self.currentNonCondLabel = self.modelNonCondensable.dicoM2V[lst[0]]
                self.modelNonCondensable.setItem(str_model = self.currentNonCond)
            else:
                self.labelNonCondensable.hide()
                self.comboBoxNonCondensable.hide()
                self.pushButtonNonCondensable.hide()

        # species initialization
        lst = self.SpeciesModel.getScalarByFieldId(fieldId)
        if len(lst) > 0 :
            self.labelScalar.show()
            self.comboBoxScalar.show()
            self.pushButtonScalar.show()

            if len(self.modelScalar.getItems()) != 0 :
                for nb in range(len(self.modelScalar.getItems())):
                    self.modelScalar.delItem(0)

            for var in lst :
                label = self.SpeciesModel.getScalarLabelByName(var)
                self.modelScalar.addItem(self.tr(label), var)

            self.currentScalar = lst[0]
            self.currentScalarLabel = self.modelScalar.dicoM2V[lst[0]]
            self.modelScalar.setItem(str_model = self.currentScalar)
        else:
            self.labelScalar.hide()
            self.comboBoxScalar.hide()
            self.pushButtonScalar.hide()


