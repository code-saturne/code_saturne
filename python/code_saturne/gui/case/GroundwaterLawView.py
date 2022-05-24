# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------

"""
This module defines the GroundwaterLaw model data management.

This module contains the following classes:
- GroundwaterLawView
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator
from code_saturne.gui.case.GroundwaterLawForm import Ui_GroundwaterLawForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.GroundwaterLawModel import GroundwaterLawModel
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("GroundwaterLawView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main view class
# -------------------------------------------------------------------------------

class GroundwaterLawView(QWidget, Ui_GroundwaterLawForm):

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_GroundwaterLawForm.__init__(self)
        self.setupUi(self)
        self.case = None
        self.mdl = None
        self.scalar = ""
        self.list_scalars = []
        self.zone = None
        self.modelGroundwaterLawType = None
        self.modelNameDiff = None

    def setup(self, case, zone_name):
        self.case = case
        self.case.undoStopGlobal()

        for zone in LocalizationModel("VolumicZone", self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone

        if self.zone.isNatureActivated("groundwater_law"):
            self.mdl = GroundwaterLawModel(self.case)

            # Combo model
            self.modelGroundwaterLawType = ComboModel(self.comboBoxType, 2, 1)
            self.modelGroundwaterLawType.addItem(self.tr("User law"), 'user')
            self.modelGroundwaterLawType.addItem(self.tr("Van Genuchten law"), 'VanGenuchten')
            self.modelNameDiff = ComboModel(self.comboBoxNameDiff, 1, 1)
            scalar_model = DefineUserScalarsModel(self.case)
            for s in scalar_model.getUserScalarNameList():
                self.list_scalars.append((s, self.tr("Additional scalar")))
            scalar_list = scalar_model.getUserScalarNameList()
            for s in scalar_model.getScalarsVarianceList():
                if s in scalar_list: scalar_list.remove(s)
            if scalar_list != []:
                self.scalar = scalar_list[0]
                for scalar in scalar_list:
                    self.modelNameDiff.addItem(scalar)
            self.setValidators()
            self.setConnections()
            self.forgetStandardWindows()
            self.selectGroundwaterLawZones()
        else:  # TODO check if content of tab should remain visible or not
            self.displayDefaultView()
            self.setEnabled(False)
        self.case.undoStartGlobal()

    def displayDefaultView(self):
        self.forgetStandardWindows()
        self.groupBoxGroundProperties.show()
        self.groupBoxSoilDensity.show()
        self.groupBoxType.show()
        self.groupBoxVanGenuchten.show()
        self.groupBoxUser.hide()
        self.groupBoxSoluteProperties.hide()

    def setConnections(self):
        self.comboBoxType.activated[str].connect(self.slotGroundwaterLaw)
        self.lineEditKs.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks"))
        self.lineEditKsXX.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xx"))
        self.lineEditKsYY.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_yy"))
        self.lineEditKsZZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_zz"))
        self.lineEditKsXY.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xy"))
        self.lineEditKsXZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xz"))
        self.lineEditKsYZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_yz"))
        self.lineEditThetas.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "thetas"))
        self.lineEditKsSaturated.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks"))
        self.lineEditKsSaturatedXX.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xx"))
        self.lineEditKsSaturatedYY.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_yy"))
        self.lineEditKsSaturatedZZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_zz"))
        self.lineEditKsSaturatedXY.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xy"))
        self.lineEditKsSaturatedXZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_xz"))
        self.lineEditKsSaturatedYZ.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "ks_yz"))
        self.lineEditThetasSaturated.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "thetas"))
        self.lineEditThetar.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "thetar"))
        self.lineEditN.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "n"))
        self.lineEditL.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "l"))
        self.lineEditAlpha.textChanged[str].connect(lambda x: self.slotSetModelValue(x, "alpha"))
        self.lineEditSoilDensity.textChanged[str].connect(self.slotSetSoilDensity)
        self.pushButtonUserLaw.clicked.connect(self.slotFormula)
        self.comboBoxNameDiff.activated[str].connect(self.slotNameDiff)
        self.lineEditDiffusivity.textChanged[str].connect(lambda x: self.slotSetScalarProperty(x, "diffusivity"))
        self.lineEditKd.textChanged[str].connect(lambda x: self.slotSetScalarProperty(x, "kd"))
        self.lineEditkplus.textChanged[str].connect(lambda x: self.slotSetScalarProperty(x, "kplus"))
        self.lineEditkminus.textChanged[str].connect(lambda x: self.slotSetScalarProperty(x, "kminus"))

    def setValidators(self):
        self.lineEditKs.setValidator(DoubleValidator(self.lineEditKs))
        self.lineEditKsXX.setValidator(DoubleValidator(self.lineEditKsXX))
        self.lineEditKsYY.setValidator(DoubleValidator(self.lineEditKsYY))
        self.lineEditKsZZ.setValidator(DoubleValidator(self.lineEditKsZZ))
        self.lineEditKsXY.setValidator(DoubleValidator(self.lineEditKsXY))
        self.lineEditKsXZ.setValidator(DoubleValidator(self.lineEditKsXZ))
        self.lineEditKsYZ.setValidator(DoubleValidator(self.lineEditKsYZ))
        self.lineEditThetas.setValidator(DoubleValidator(self.lineEditThetas))
        self.lineEditThetar.setValidator(DoubleValidator(self.lineEditThetar))
        self.lineEditN.setValidator(DoubleValidator(self.lineEditN))
        self.lineEditL.setValidator(DoubleValidator(self.lineEditL))
        self.lineEditAlpha.setValidator(DoubleValidator(self.lineEditAlpha))
        self.lineEditLongitudinal.setValidator(DoubleValidator(self.lineEditLongitudinal))
        self.lineEditTransverse.setValidator(DoubleValidator(self.lineEditTransverse))
        self.lineEditKsSaturated.setValidator(DoubleValidator(self.lineEditKsSaturated))
        self.lineEditKsSaturatedXX.setValidator(DoubleValidator(self.lineEditKsSaturatedXX))
        self.lineEditKsSaturatedYY.setValidator(DoubleValidator(self.lineEditKsSaturatedYY))
        self.lineEditKsSaturatedZZ.setValidator(DoubleValidator(self.lineEditKsSaturatedZZ))
        self.lineEditKsSaturatedXY.setValidator(DoubleValidator(self.lineEditKsSaturatedXY))
        self.lineEditKsSaturatedXZ.setValidator(DoubleValidator(self.lineEditKsSaturatedXZ))
        self.lineEditKsSaturatedYZ.setValidator(DoubleValidator(self.lineEditKsSaturatedYZ))
        self.lineEditThetasSaturated.setValidator(DoubleValidator(self.lineEditThetasSaturated))
        self.lineEditSoilDensity.setValidator(DoubleValidator(self.lineEditSoilDensity))
        self.lineEditDiffusivity.setValidator(DoubleValidator(self.lineEditDiffusivity))
        self.lineEditKd.setValidator(DoubleValidator(self.lineEditKd))
        self.lineEditkplus.setValidator(DoubleValidator(self.lineEditkplus))
        self.lineEditkminus.setValidator(DoubleValidator(self.lineEditkminus))

    @pyqtSlot("QModelIndex")
    def selectGroundwaterLawZones(self):
        label = self.zone.getLabel()
        name = self.zone.getCodeNumber()

        if hasattr(self, "modelScalars"): del self.modelScalars
        log.debug("slotSelectGroundwaterLawZones label %s " % label)

        self.groupBoxGroundProperties.show()

        # ground properties
        self.groupBoxSoilDensity.show()
        value = self.mdl.getSoilDensity(name)
        self.lineEditSoilDensity.setText(str(value))

        if GroundwaterModel(self.case).getUnsaturatedZone() == "true":
            self.groupBoxType.show()

            choice = self.mdl.getGroundwaterLawModel(name)
            self.modelGroundwaterLawType.setItem(str_model=choice)

            if choice == "user":
                self.groupBoxVanGenuchten.hide()
                self.groupBoxUser.show()
                exp = self.mdl.getGroundwaterLawFormula(name)
                if exp:
                    self.pushButtonUserLaw.setStyleSheet("background-color: green")
                    self.pushButtonUserLaw.setToolTip(exp)
                else:
                    self.pushButtonUserLaw.setStyleSheet("background-color: red")
            else:
                self.groupBoxVanGenuchten.show()
                self.groupBoxUser.hide()
                self.initializeVanGenuchten(name)
        else:
            self.groupBoxSaturated.show()
            if GroundwaterModel(self.case).getPermeabilityType() == "anisotropic":
                self.labelKsSaturated.hide()
                self.lineEditKsSaturated.hide()
                value = self.mdl.getValue(name, "ks_xx")
                self.lineEditKsSaturatedXX.setText(str(value))
                value = self.mdl.getValue(name, "ks_yy")
                self.lineEditKsSaturatedYY.setText(str(value))
                value = self.mdl.getValue(name, "ks_zz")
                self.lineEditKsSaturatedZZ.setText(str(value))
                value = self.mdl.getValue(name, "ks_xy")
                self.lineEditKsSaturatedXY.setText(str(value))
                value = self.mdl.getValue(name, "ks_xz")
                self.lineEditKsSaturatedXZ.setText(str(value))
                value = self.mdl.getValue(name, "ks_yz")
                self.lineEditKsSaturatedYZ.setText(str(value))
            else:
                self.labelKsSaturatedXX.hide()
                self.labelKsSaturatedYY.hide()
                self.labelKsSaturatedZZ.hide()
                self.labelKsSaturatedXY.hide()
                self.labelKsSaturatedXZ.hide()
                self.labelKsSaturatedYZ.hide()
                self.lineEditKsSaturatedXX.hide()
                self.lineEditKsSaturatedYY.hide()
                self.lineEditKsSaturatedZZ.hide()
                self.lineEditKsSaturatedXY.hide()
                self.lineEditKsSaturatedXZ.hide()
                self.lineEditKsSaturatedYZ.hide()
                value = self.mdl.getValue(name, "ks")
                self.lineEditKsSaturated.setText(str(value))
            value = self.mdl.getValue(name, "thetas")
            self.lineEditThetasSaturated.setText(str(value))

        # solute properties
        scal = self.scalar
        if scal == "":
            self.groupBoxSoluteProperties.hide()
        else:
            self.groupBoxSoluteProperties.show()

            self.modelNameDiff.setItem(str_model=str(scal))
            value = self.mdl.getGroundWaterScalarPropertyByZone(scal, name, 'diffusivity')
            self.lineEditDiffusivity.setText(str(value))
            value = self.mdl.getGroundWaterScalarPropertyByZone(scal, name, 'kd')
            self.lineEditKd.setText(str(value))
            # chemistry model
            if GroundwaterModel(self.case).getChemistryModel(scal) == "EK":
                self.lineEditkplus.show()
                self.label_kplus.show()
                self.lineEditkminus.show()
                self.label_kminus.show()
                value = self.mdl.getGroundWaterScalarPropertyByZone(scal, name, 'kplus')
                self.lineEditkplus.setText(str(value))
                value = self.mdl.getGroundWaterScalarPropertyByZone(scal, name, 'kminus')
                self.lineEditkminus.setText(str(value))
            else:
                self.lineEditkplus.hide()
                self.label_kplus.hide()
                self.lineEditkminus.hide()
                self.label_kminus.hide()

        value = self.mdl.getDispersionCoefficient(name, "longitudinal")
        self.lineEditLongitudinal.setText(str(value))
        value = self.mdl.getDispersionCoefficient(name, "transverse")
        self.lineEditTransverse.setText(str(value))

    def initializeVanGenuchten(self, name):
        """
        initialize variables for groupBoxVanGenuchten
        """
        value = self.mdl.getValue(name, "thetas")
        self.lineEditThetas.setText(str(value))
        value = self.mdl.getValue(name, "thetar")
        self.lineEditThetar.setText(str(value))
        value = self.mdl.getValue(name, "n")
        self.lineEditN.setText(str(value))
        value = self.mdl.getValue(name, "l")
        self.lineEditL.setText(str(value))
        value = self.mdl.getValue(name, "alpha")
        self.lineEditAlpha.setText(str(value))
        if GroundwaterModel(self.case).getPermeabilityType() == "anisotropic":
            self.labelKs.hide()
            self.lineEditKs.hide()
            value = self.mdl.getValue(name, "ks_xx")
            self.lineEditKsXX.setText(str(value))
            value = self.mdl.getValue(name, "ks_yy")
            self.lineEditKsYY.setText(str(value))
            value = self.mdl.getValue(name, "ks_zz")
            self.lineEditKsZZ.setText(str(value))
            value = self.mdl.getValue(name, "ks_xy")
            self.lineEditKsXY.setText(str(value))
            value = self.mdl.getValue(name, "ks_xz")
            self.lineEditKsXZ.setText(str(value))
            value = self.mdl.getValue(name, "ks_yz")
            self.lineEditKsYZ.setText(str(value))
        else:
            self.labelKsXX.hide()
            self.labelKsYY.hide()
            self.labelKsZZ.hide()
            self.labelKsXY.hide()
            self.labelKsXZ.hide()
            self.labelKsYZ.hide()
            self.lineEditKsXX.hide()
            self.lineEditKsYY.hide()
            self.lineEditKsZZ.hide()
            self.lineEditKsXY.hide()
            self.lineEditKsXZ.hide()
            self.lineEditKsYZ.hide()
            value = self.mdl.getValue(name, "ks")
            self.lineEditKs.setText(str(value))

    def forgetStandardWindows(self):
        """
        For forget standard windows
        """
        self.groupBoxGroundProperties.hide()
        self.groupBoxType.hide()
        self.groupBoxUser.hide()
        self.groupBoxVanGenuchten.hide()
        self.groupBoxSaturated.hide()
        self.groupBoxSoluteProperties.hide()

    @pyqtSlot(str)
    def slotGroundwaterLaw(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        name = self.zone.getCodeNumber()
        choice = self.modelGroundwaterLawType.dicoV2M[str(text)]

        self.mdl.setGroundwaterLawModel(name, choice)

        if choice == "user":
            self.groupBoxVanGenuchten.hide()
            self.groupBoxUser.show()
            exp = self.mdl.getGroundwaterLawFormula(name)
            if exp:
                self.pushButtonUserLaw.setStyleSheet("background-color: green")
                self.pushButtonUserLaw.setToolTip(exp)
            else:
                self.pushButtonUserLaw.setStyleSheet("background-color: red")
        else:
            self.groupBoxVanGenuchten.show()
            self.groupBoxUser.hide()
            self.initializeVanGenuchten(name)

    def slotSetModelValue(self, text, field_name):

        # Avoid crash when change is induced by the combobox which
        # has no validator...
        _v = self.sender().validator()
        if _v and _v.state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(self.zone.getCodeNumber(), field_name, val)

    def slotSetDispersionCoefficient(self, text, field_name):

        # Avoid crash when change is induced by the combobox which
        # has no validator...
        _v = self.sender().validator()
        if _v and _v.state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setDispersionCoefficient(self.zone.getCodeNumber(), field_name, val)

    def slotSetScalarProperty(self, text, field_name):

        # Avoid crash when change is induced by the combobox which
        # has no validator...
        _v = self.sender().validator()
        if _v and _v.state == QValidator.Acceptable:
            value = float(text)
            scalar = self.scalar
            self.mdl.setGroundWaterScalarPropertyByZone(scalar,
                                                        self.zone.getCodeNumber(),
                                                        field_name,
                                                        value)

    @pyqtSlot(str)
    def slotSetSoilDensity(self, text):
        """
        """
        if self.lineEditSoilDensity.validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setSoilDensity(self.zone.getCodeNumber(), val)

    @pyqtSlot()
    def slotFormula(self):
        """
        User formula for Groundwater functions
        """
        label = self.zone.getLabel()
        name = self.zone.getCodeNumber()

        exp, req, sym = self.mdl.getGroundwaterLawFormulaComponents(name)

        exa = """#example: \n""" + self.mdl.getDefaultGroundwaterLawFormula()

        dialog = QMegEditorView(parent=self,
                                function_type='vol',
                                zone_name=label,
                                variable_name='capacity+saturation+permeability',
                                expression=exp,
                                required=req,
                                symbols=sym,
                                examples=exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormula -> %s" % str(result))
            self.mdl.setGroundwaterLawFormula(name, str(result))
            self.pushButtonUserLaw.setStyleSheet("background-color: green")
            self.pushButtonUserLaw.setToolTip(result)

    @pyqtSlot(str)
    def slotNameDiff(self, text):
        """
        Method to choose the scalar which properties shall be changed
        """
        name = self.zone.getCodeNumber()
        log.debug("slotNameDiff -> %s" % (text))
        self.scalar = str(text)
        scalar = self.scalar
        value = self.mdl.getGroundWaterScalarPropertyByZone(scalar, name, 'diffusivity')
        self.lineEditDiffusivity.setText(str(value))
        value = self.mdl.getGroundWaterScalarPropertyByZone(scalar, name, 'kd')
        self.lineEditKd.setText(str(value))
        if GroundwaterModel(self.case).getChemistryModel(scalar) == "EK":
            self.lineEditkplus.show()
            self.label_kplus.show()
            self.lineEditkminus.show()
            self.label_kminus.show()
            value = self.mdl.getGroundWaterScalarPropertyByZone(scalar, name, 'kplus')
            self.lineEditkplus.setText(str(value))
            value = self.mdl.getGroundWaterScalarPropertyByZone(scalar, name, 'kminus')
            self.lineEditkminus.setText(str(value))
        else:
            self.lineEditkplus.hide()
            self.label_kplus.hide()
            self.lineEditkminus.hide()
            self.label_kminus.hide()


# -------------------------------------------------------------------------------
# Testing part
# -------------------------------------------------------------------------------


if __name__ == "__main__":
    pass

# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
