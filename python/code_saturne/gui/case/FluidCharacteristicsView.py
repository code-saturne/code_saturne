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
This module defines the Page for the physical properties of the fluid.
These properties can be reference value or initial value

This module contains the following classes and function:
- FluidCharacteristicsView
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
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.gui.case.FluidCharacteristicsForm import Ui_FluidCharacteristicsForm
from code_saturne.model.AtmosphericFlowsModel import AtmosphericFlowsModel
from code_saturne.model.FluidCharacteristicsModel import FluidCharacteristicsModel
from code_saturne.model.DefineUserScalarsModel import DefineUserScalarsModel
from code_saturne.model.ThermalScalarModel import ThermalScalarModel
from code_saturne.model.GroundwaterModel import GroundwaterModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.InternalCouplingModel import InternalCouplingModel
from code_saturne.model.LocalizationModel import LocalizationModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FluidCharacteristicsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class FluidCharacteristicsView(QWidget, Ui_FluidCharacteristicsForm):

    """
    Class to open Molecular Properties Page.
    """
    density = """# Density of air
density = 1.293*(273.15 / temperature);


# density for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

rho1 = 1.25051;
rho2 = 1.7832;
A = (Y1 / rho1) + (Y2 /rho2);
density = 1.0 / A;
"""

    density_h = """# Density
density = enthalpy / 1040. * 1.29;

# density for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

rho1 = 1.25051;
rho2 = 1.7832;
A = (Y1 / rho1) + (Y2 /rho2);
density = 1.0 / A;
"""

    density_wo = """density = 1.25051;

"""
    molecular_viscosity="""# Sutherland's Formula
# Gas             Cst    T0      mu0
# air             120    291.15  18.27e-6
# nitrogen        111    300.55  17.81e-6
# oxygen          127    292.25  20.18e-6
# carbon dioxide  240    293.15  14.8e-6
# carbon monoxide 118    288.15  17.2e-6
# hydrogen        72     293.85  8.76e-6
# ammonia         370    293.15  9.82e-6
# sulfur dioxide  416    293.65  12.54e-6
# helium          79.4   273     19e-6

CST = 120;
T0 = 291.15;
mu_ref = 18.27e-6;

if ( temperature > 0 && temperature < 555) {
molecular_viscosity = mu_ref * ((T0+CST) / (temperature+CST)) * (temperature/T0)^(3./2.);
} else {
molecular_viscosity = -999.0;
}
"""

    molecular_viscosity_h="""CST = 120;
T0 = 291.15;
mu_ref = 18.27e-6;
temperature = enthalpy / 1040.;

if ( enthalpy > 0) {
molecular_viscosity = mu_ref * (T0+CST / temperature+CST) * (temperature/T0)^(3./2.);
} else {
molecular_viscosity = -999.0;
}
"""

    molecular_viscosity_wo="""CST = 120;
T0 = 291.15;
mu_ref = 18.27e-6;
molecular_viscosity = mu_ref * (T0+CST);

"""

    specific_heat="""# specific heat for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

Cp1 = 520.3;
Cp2 = 1040.0;
specific_heat = Y1 * Cp1 + Y2 *Cp2;
"""

    volume_viscosity="""# volume_viscosity
"""

    thermal_conductivity="""# oxygen
thermal_conductivity = 6.2e-5 * temperature + 8.1e-3;

# nitrogen
thermal_conductivity = 6.784141e-5 * temperature + 5.564317e-3;

# hydrogen
thermal_conductivity = 4.431e-4 * temperature + 5.334e-2;
"""

    thermal_conductivity_h="""temperature = enthalpy / 1040.;
thermal_conductivity = 6.2e-5 * temperature + 8.1e-3;
"""

    def __init__(self, parent=None): #, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_FluidCharacteristicsForm.__init__(self)
        self.setupUi(self)

        self.case     = None
        self.mdl      = None
        self.notebook = None
        self.zone     = None
        self.zone_name = None
        self.zone_id   = None
        self.is_main_zone = False
        self.is_solid  = False


    def setup(self, case, zone_name):

        self.case = case

        for zone in LocalizationModel('VolumicZone', self.case).getZones():
            if zone.getLabel() == zone_name:
                self.zone = zone
                self.zone_name = zone.getLabel()
                self.zone_id   = str(zone.getCodeNumber())
                if self.zone_name == "all_cells":
                    self.is_main_zone = True
                elif not zone.isNatureActivated("physical_properties"):
                    return
                if zone.isNatureActivated("solid"):
                    self.is_solid = True

        self.case.undoStopGlobal()

        self.mdl = FluidCharacteristicsModel(self.case)
        self.notebook = NotebookModel(self.case)

        self.internal_cpl = InternalCouplingModel(self.case)
        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarName()
        tsm = self.mdl.tsm

        # Particular Widget init. taking into account chosen fluid model
        mdl_atmo, mdl_joule, mdl_thermal, mdl_gas, mdl_coal, mdl_comp, mdl_hgn=\
            self.mdl.getThermoPhysicalModel()

        # Combo models

        self.modelRho       = ComboModel(self.comboBoxRho,      2, 1)
        self.modelMu        = ComboModel(self.comboBoxMu,       2, 1)
        self.modelSigma     = ComboModel(self.comboBoxSigma,    1, 1)
        self.modelCp        = ComboModel(self.comboBoxCp,       2, 1)
        self.modelAl        = ComboModel(self.comboBoxAl,       2, 1)
        self.modelDiff      = ComboModel(self.comboBoxDiff,     2, 1)
        self.modelNameDiff  = ComboModel(self.comboBoxNameDiff, 1, 1)
        self.modelViscv0    = ComboModel(self.comboBoxViscv0,   3, 1)
        self.modelDiftl0    = ComboModel(self.comboBoxDiftl0,   3, 1)
        self.modelMaterial  = ComboModel(self.comboBoxMaterial, 1, 1)
        self.modelMethod    = ComboModel(self.comboBoxMethod,   1, 1)
        self.modelReference = ComboModel(self.comboBoxReference,  1, 1)

        if mdl_joule == 'off':
            self.modelRho.addItem(self.tr('constant'), 'constant')
            self.modelRho.addItem(self.tr('material law'), 'thermal_law')
        if mdl_atmo != 'off':
            self.modelRho.addItem(self.tr('predefined law'), 'predefined_law')
        elif mdl_hgn != "off":
            self.modelRho.addItem(self.tr('linear mixture law'), 'predefined_law')
        elif mdl_joule != 'off':
            self.modelRho.addItem(self.tr('predefined law'), 'predefined_law')
        elif mdl_comp != 'off':
            self.modelRho.addItem(self.tr('predefined law'), 'predefined_law')
        elif mdl_gas != 'off' or mdl_coal != 'off':
            self.modelRho.addItem(self.tr('predefined law'), 'predefined_law')
        self.modelRho.addItem(self.tr('user law'), 'user_law')

        if mdl_joule == 'off':
            self.modelMu.addItem(self.tr('constant'), 'constant')
            self.modelMu.addItem(self.tr('material law'), 'thermal_law')
        if mdl_joule != 'off':
            self.modelMu.addItem(self.tr('predefined law'), 'predefined_law')
        if mdl_hgn != "off":
            self.modelMu.addItem(self.tr('linear mixture law'), 'predefined_law')
        self.modelMu.addItem(self.tr('user law'), 'user_law')

        if mdl_hgn != "off":
            self.modelSigma.addItem(self.tr('constant'), 'constant')

        if mdl_joule == 'off':
            self.modelCp.addItem(self.tr('constant'), 'constant')
            self.modelCp.addItem(self.tr('material law'), 'thermal_law')
        if mdl_joule != 'off':
            self.modelCp.addItem(self.tr('predefined law'), 'predefined_law')
        self.modelCp.addItem(self.tr('user law'), 'user_law')

        if mdl_joule == 'off':
            self.modelAl.addItem(self.tr('constant'), 'constant')
            self.modelAl.addItem(self.tr('material law'), 'thermal_law')
        if mdl_joule != 'off':
            self.modelAl.addItem(self.tr('predefined law'), 'predefined_law')
        self.modelAl.addItem(self.tr('user law'), 'user_law')

        self.modelDiff.addItem(self.tr('constant'), 'constant')
        self.modelDiff.addItem(self.tr('user law'), 'user_law')

        self.modelViscv0.addItem(self.tr('constant'), 'constant')
        self.modelViscv0.addItem(self.tr('user law'), 'user_law')
        self.modelViscv0.addItem(self.tr('material law'), 'thermal_law')

        self.modelDiftl0.addItem(self.tr('constant'), 'constant')
        self.modelDiftl0.addItem(self.tr('user law'), 'user_law')
        self.modelDiftl0.addItem(self.tr('material law'), 'thermal_law')
        if mdl_gas != 'off' or mdl_coal != 'off':
            self.modelDiftl0.addItem(self.tr('predefined law'), 'predefined_law')

        self.scalar = ""
        scalar_list = self.mdl.m_sca.getUserScalarNameList()
        for s in self.mdl.m_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)

        if scalar_list != []:
            self.scalar = scalar_list[0]
            for scalar in scalar_list:
                self.modelNameDiff.addItem(scalar)

        # Validators

        validatorP0 = DoubleValidator(self.lineEditP0, min=0.0)
        self.lineEditP0.setValidator(validatorP0)

        validatorT0 = DoubleValidator(self.lineEditT0,  min=0.0)
        validatorT0.setExclusiveMin(True)
        self.lineEditT0.setValidator(validatorT0)

        validatorOxydant = DoubleValidator(self.lineEditOxydant,  min=0.0)
        validatorOxydant.setExclusiveMin(True)
        self.lineEditOxydant.setValidator(validatorOxydant)

        validatorFuel = DoubleValidator(self.lineEditFuel,  min=0.0)
        validatorFuel.setExclusiveMin(True)
        self.lineEditFuel.setValidator(validatorFuel)

        validatorMM = DoubleValidator(self.lineEditMassMolar, min=0.0)
        validatorMM.setExclusiveMin(True)
        self.lineEditMassMolar.setValidator(validatorMM)

        validatorRho    = DoubleValidator(self.lineEditRho,    min = 0.0)
        validatorRho1   = DoubleValidator(self.lineEditRho1,   min = 0.0)
        validatorRho2   = DoubleValidator(self.lineEditRho2,   min = 0.0)
        validatorMu     = DoubleValidator(self.lineEditMu,     min = 0.0)
        validatorMu1    = DoubleValidator(self.lineEditMu1,    min = 0.0)
        validatorMu2    = DoubleValidator(self.lineEditMu2,    min = 0.0)
        validatorSigma  = DoubleValidator(self.lineEditSigma,  min = 0.0)
        validatorCp     = DoubleValidator(self.lineEditCp,     min = 0.0)
        validatorAl     = DoubleValidator(self.lineEditAl,     min = 0.0)
        validatorDiff   = DoubleValidator(self.lineEditDiff,   min = 0.0)
        validatorViscv0 = DoubleValidator(self.lineEditViscv0, min = 0.0)
        validatorDiftl0 = DoubleValidator(self.lineEditDiftl0, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorRho1.setExclusiveMin(True)
        validatorRho2.setExclusiveMin(True)
        validatorMu.setExclusiveMin(True)
        validatorMu1.setExclusiveMin(True)
        validatorMu2.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)
        validatorDiff.setExclusiveMin(True)
        validatorDiftl0.setExclusiveMin(True)

        self.lineEditRho.setValidator(validatorRho)
        self.lineEditRho1.setValidator(validatorRho1)
        self.lineEditRho2.setValidator(validatorRho2)
        self.lineEditMu.setValidator(validatorMu)
        self.lineEditMu1.setValidator(validatorMu1)
        self.lineEditMu2.setValidator(validatorMu2)
        self.lineEditSigma.setValidator(validatorSigma)
        self.lineEditCp.setValidator(validatorCp)
        self.lineEditAl.setValidator(validatorAl)
        self.lineEditDiff.setValidator(validatorDiff)
        self.lineEditViscv0.setValidator(validatorViscv0)
        self.lineEditDiftl0.setValidator(validatorDiftl0)

        # Connections

        self.lineEditP0.textChanged[str].connect(self.slotPressure)
        self.lineEditT0.textChanged[str].connect(self.slotTemperature)
        self.lineEditOxydant.textChanged[str].connect(self.slotTempOxydant)
        self.lineEditFuel.textChanged[str].connect(self.slotTempFuel)
        self.lineEditMassMolar.textChanged[str].connect(self.slotMassemol)

        self.comboBoxRho.currentIndexChanged[str].connect(self.slotStateRho)
        self.comboBoxMu.currentIndexChanged[str].connect(self.slotStateMu)
        self.comboBoxSigma.currentIndexChanged[str].connect(self.slotStateSigma)
        self.comboBoxCp.currentIndexChanged[str].connect(self.slotStateCp)
        self.comboBoxAl.currentIndexChanged[str].connect(self.slotStateAl)
        self.comboBoxDiff.currentIndexChanged[str].connect(self.slotStateDiff)
        self.comboBoxNameDiff.activated[str].connect(self.slotNameDiff)
        self.comboBoxViscv0.currentIndexChanged[str].connect(self.slotStateViscv0)
        self.comboBoxMaterial.activated[str].connect(self.slotMaterial)
        self.comboBoxMethod.activated[str].connect(self.slotMethod)
        self.comboBoxReference.activated[str].connect(self.slotReference)
        self.lineEditRho.textChanged[str].connect(self.slotRho)
        self.lineEditRho1.textChanged[str].connect(self.slotRho1)
        self.lineEditRho2.textChanged[str].connect(self.slotRho2)
        self.lineEditMu.textChanged[str].connect(self.slotMu)
        self.lineEditMu1.textChanged[str].connect(self.slotMu1)
        self.lineEditMu2.textChanged[str].connect(self.slotMu2)
        self.lineEditSigma.textChanged[str].connect(self.slotSigma)
        self.lineEditCp.textChanged[str].connect(self.slotCp)
        self.lineEditAl.textChanged[str].connect(self.slotAl)
        self.lineEditDiff.textChanged[str].connect(self.slotDiff)
        self.lineEditDiftl0.textChanged[str].connect(self.slotDiftl0)
        self.lineEditViscv0.textChanged[str].connect(self.slotViscv0)
        self.pushButtonRho.clicked.connect(self.slotFormulaRho)
        self.pushButtonMu.clicked.connect(self.slotFormulaMu)
        self.pushButtonCp.clicked.connect(self.slotFormulaCp)
        self.pushButtonAl.clicked.connect(self.slotFormulaAl)
        self.pushButtonDiff.clicked.connect(self.slotFormulaDiff)
        self.pushButtonViscv0.clicked.connect(self.slotFormulaViscv0)

        self.initializeWidget()

        self.case.undoStartGlobal()


    def initializeWidget(self):
        """
        """
        mdls = self.mdl.getThermoPhysicalModel()
        mdl_atmo, mdl_joule, mdl_thermal, mdl_gas, mdl_coal, mdl_comp, mdl_hgn = mdls

        if self.is_solid:
            mdl_atmo = "off"
            mdl_joule = "off"
            mdl_gas = "off"
            mdl_coal = "off"
            mdl_comp = "off"
            mdl_hgn = "off"

        is_main_zone = self.is_main_zone

        self.groupBoxPressure.setVisible(is_main_zone)
        self.groupBoxTemperature.setVisible(is_main_zone)

        self.groupBoxMassMolar.hide()

        isLargeScaleMeteoChecked = AtmosphericFlowsModel(self.case).getLargeScaleMeteoStatus() == 'on'
        if mdl_atmo != "off" and isLargeScaleMeteoChecked:
            self.groupBoxTemperature.setEnabled(False)
            t = AtmosphericFlowsModel(self.case).getMeteoT0()
            self.lineEditT0.setText(str(t))

        elif mdl_comp != "off" or mdl_coal != "off":
            m = self.mdl.getMassemol()
            self.lineEditMassMolar.setText(str(m))
            self.groupBoxMassMolar.setVisible(is_main_zone)
            if mdl_comp != "off":
                t = self.mdl.getTemperature()
                self.lineEditT0.setText(str(t))

        elif mdl_thermal != "off" or mdl_gas == 'd3p':
            t = self.mdl.getTemperature()
            self.lineEditT0.setText(str(t))
            if mdl_thermal == "temperature_celsius":
                self.labelUnitT0.setText("\xB0 C")

            if self.mdl.getMaterials() != "user_material":
                self.labelInfoT0.hide()
        else:
            self.groupBoxTemperature.hide()

        if mdl_gas == 'd3p':
            self.groupBoxTempd3p.setVisible(is_main_zone)
            t_oxy  = self.mdl.getTempOxydant()
            t_fuel = self.mdl.getTempFuel()
            self.lineEditOxydant.setText(str(t_oxy))
            self.lineEditFuel.setText(str(t_fuel))
        else:
            self.groupBoxTempd3p.hide()

        darc = GroundwaterModel(self.case).getGroundwaterModel()
        if darc != 'off':
            self.groupBoxPressure.hide()
        elif mdl_atmo != 'off' and isLargeScaleMeteoChecked:
            self.groupBoxPressure.setEnabled(False)
            p = AtmosphericFlowsModel(self.case).getMeteoPsea()
            self.lineEditP0.setText(str(p))
            self.groupBoxPressure.setEnabled(False)

        else:
            p = self.mdl.getPressure()
            self.lineEditP0.setText(str(p))

        if self.mdl.tables == 0 or mdl_joule != 'off' or mdl_comp != 'off':
            self.groupBoxTableChoice.hide()
        else:
            self.groupBoxTableChoice.setVisible(is_main_zone)

            fluids = self.mdl.getLibPropertiesDict()
            o_keys = list(fluids.keys())
            o_keys.sort()
            for f in o_keys:
                if fluids[f] != 0:
                    self.modelMaterial.addItem(self.tr(f), f)
                else:
                    self.modelMaterial.addItem(self.tr(f), True)

            material = self.mdl.getMaterials()
            self.modelMaterial.setItem(str_model=material)
            self.updateMethod()

        # VoF
        if mdl_hgn != 'off':
            self.widgetRefRho.hide()
            self.widgetVofRho.setVisible(is_main_zone)
            self.widgetRefMu.hide()
            self.widgetVofMu.setVisible(is_main_zone)
            self.groupBoxSigma.setVisible(is_main_zone)
        elif mdl_atmo != 'off' and isLargeScaleMeteoChecked:
            self.widgetRefRho.hide()
            self.widgetVofRho.hide()
            self.widgetRefMu.setVisible(is_main_zone)
            self.widgetVofMu.hide()
            self.groupBoxSigma.hide()
        else:
            if mdl_comp != "off":
                self.widgetRefRho.hide()
            else:
                self.widgetRefRho.setVisible(is_main_zone)
            self.widgetVofRho.hide()
            self.widgetRefMu.setVisible(is_main_zone)
            self.widgetVofMu.hide()
            self.groupBoxSigma.hide()

        # compressible
        self.groupBoxViscv0.hide()

        # combustion
        self.groupBoxDiftl0.hide()

        if self.scalar == "":
            self.groupBoxDiff.hide()
        else :
            self.groupBoxDiff.show()
            self.lineEditDiff.setText(str(self.mdl.m_sca.getScalarDiffusivityInitialValue(self.scalar)))

            diff_choice =  self.mdl.m_sca.getScalarDiffusivityChoice(self.scalar)
            self.modelDiff.setItem(str_model=diff_choice)
            self.modelNameDiff.setItem(str_model=str(self.scalar))
            if diff_choice  != 'user_law':
                self.pushButtonDiff.hide()
                self.pushButtonDiff.setEnabled(False)
                self.pushButtonDiff.setStyleSheet("background-color: None")
            else:
                self.pushButtonDiff.show()
                self.pushButtonDiff.setEnabled(True)
                name = self.mdl.m_sca.getScalarDiffusivityName(self.scalar)
                exp = self.mdl.m_sca.getDiffFormula(self.scalar, self.zone_id)
                if exp:
                    self.pushButtonDiff.setStyleSheet("background-color: green")
                    self.pushButtonDiff.setToolTip(exp)
                else:
                    self.pushButtonDiff.setStyleSheet("background-color: red")

            self.lineEditDiff.setVisible(is_main_zone)
            self.labelDiff.setVisible(is_main_zone)
            self.labelUnitDiff.setVisible(is_main_zone)
            self.comboBoxDiff.setEnabled(is_main_zone)

        # Solid zone

        if self.is_solid:
            self.groupBoxTempd3p.hide()
            self.groupBoxMu.hide()
            self.groupBoxDiftl0.hide()
            self.groupBoxDiff.hide()

        # Standard Widget initialization
        for tag, symbol in self.mdl.lst:
            __model  = getattr(self, "model"      + symbol)
            __line   = getattr(self, "lineEdit"   + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label"      + symbol)
            __labelu = getattr(self, "labelUnit"  + symbol)

            try:
                __combo  = getattr(self, "comboBox" + symbol)
                __combo.setEnabled(is_main_zone)
            except:
                pass

            if tag != 'dynamic_diffusion':
                __labelv = getattr(self, "labelVar"   + symbol)
                c = self.mdl.getPropertyMode(tag)
                if c not in __model.getItems():
                    c = 'constant'
                    self.mdl.setPropertyMode(tag, c)
                if self.is_solid:
                    c = "user_law"

                __model.setItem(str_model=c)
                if c == 'user_law':
                    __button.show()
                    __button.setEnabled(True)
                    __label.setText(self.tr("Reference value"))
                else:
                    __button.hide()
                    __button.setEnabled(False)
                    __label.setText(self.tr("Reference value"))
                if c == 'thermal_law':
                    __line.hide()
                    __label.hide()
                    __labelu.hide()
                    __labelv.hide()
                else:
                    __line.setVisible(is_main_zone)
                    __label.setVisible(is_main_zone)
                    __labelu.setVisible(is_main_zone)
                    __labelv.setVisible(is_main_zone)
                try:
                    if self.mdl.getMaterials() == "user_material":
                        __model.disableItem(str_model='thermal_law')
                    else:
                        __model.enableItem(str_model='thermal_law')
                except Exception:
                    pass

            else:
                __label.setText(self.tr("Reference value"))

            __line.setText(str(self.mdl.getInitialValue(tag)))


        # If no thermal scalar, hide specific heat and thermal conductivity.
        if mdl_thermal == 'off':
            self.groupBoxCp.hide()
            self.groupBoxAl.hide()

        if mdl_gas != 'off' or mdl_coal != 'off':
            self.groupBoxDiftl0.setVisible(is_main_zone)

        for tag, symbol in self.mdl.lst:
            __model  = getattr(self, "model" + symbol)
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label" + symbol)
            __combo  = getattr(self, "comboBox" + symbol)

            # Gas or coal combustion
            if mdl_gas != 'off' or mdl_coal != 'off':
                if tag == 'density':
                    __model.setItem(str_model='predefined_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()
                    self.mdl.setPropertyMode(tag, 'predefined_law')
                    __label.setText(self.tr("Calculation by\n perfect gas law"))
                    __line.setText(str(""))
                    __line.setEnabled(False)
                elif tag == 'dynamic_diffusion':
                    __model.setItem(str_model='predefined_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()
                else:
                    choice_p = self.mdl.getPropertyMode(tag)
                    __model.setItem(str_model=choice_p)

            # Joule
            if mdl_joule == 'arc':
                __model.setItem(str_model='predefined_law')
                __combo.setEnabled(False)
                __button.setEnabled(False)
                __button.hide()
                self.mdl.setPropertyMode(tag, 'predefined_law')

            # Atmospheric Flows
            if mdl_atmo != 'off':
                if tag == 'density':
                    __model.disableItem(str_model='constant')
                    __model.disableItem(str_model='predefined_law')
                    __model.setItem(str_model='predefined_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()

            # VoF
            if mdl_hgn != 'off':
                if tag == 'molecular_viscosity' or tag == 'density':
                    __model.disableItem(str_model='constant')
                    __model.disableItem(str_model='predefined_law')
                    __model.setItem(str_model='predefined_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()
                if tag == 'surface_tension':
                    __model.disableItem(str_model='constant')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()

            # Compressible Flows
            if mdl_comp != 'off':
                self.groupBoxViscv0.setVisible(is_main_zone)
                if tag == 'density':
                    __model.setItem(str_model='predefined_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()
                    self.mdl.setPropertyMode(tag, 'predefined_law')
                    __line.setEnabled(True)
                if tag == 'specific_heat':
                    __model.setItem(str_model='constant')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __button.hide()
                    self.mdl.setPropertyMode(tag, 'constant')
                    self.groupBoxCp.setTitle('Isobaric specific heat')

                if tag == 'volume_viscosity':
                    __combo.setEnabled(True)
                    c = self.mdl.getPropertyMode(tag)
                    if c == 'user_law':
                        __button.setEnabled(True)
                    else:
                        __button.setEnabled(False)
            else:
                if tag == 'specific_heat':
                    self.groupBoxCp.setTitle('Specific heat')


        if mdl_hgn != 'off':
            self.lineEditRho1.setText(str(self.mdl.getVofValueDensity(0)))
            self.lineEditRho2.setText(str(self.mdl.getVofValueDensity(1)))

            self.lineEditMu1.setText(str(self.mdl.getVofValueViscosity(0)))
            self.lineEditMu2.setText(str(self.mdl.getVofValueViscosity(1)))

            self.lineEditSigma.setText(str(self.mdl.getInitialValueSurfaceTension()))

    def updateTypeChoice(self, old_choice):
        """
        add/suppress thermo tables for each properties
        """
        for tag, symbol in self.mdl.lst:
            __model  = getattr(self, "model" + symbol)
            if self.mdl.getMaterials() == "user_material":
                __model.disableItem(str_model='thermal_law')
            else:
                __model.enableItem(str_model='thermal_law')
                if old_choice == "user_material":
                    self.mdl.setPropertyMode(tag, 'thermal_law')
                self.__changeChoice(str("material law"), symbol, tag)
            c = self.mdl.getPropertyMode(tag)
            __model.setItem(str_model=c)


    def updateMethod(self):
        """
        update method list with material choice
        """

        for nb in range(len(self.modelMethod.getItems())):
            self.modelMethod.delItem(0)

        material = self.mdl.getMaterials()
        method = self.mdl.getMethod()
        have_method = False

        methods = self.mdl.getLibPropertyMethods(material)
        for m in methods:
            if method == m[0]:
                have_method = True
            self.modelMethod.addItem(self.tr(m[0]), m[0], not m[1])

        # current method not in list, add it with warning
        if not have_method:
            self.modelMethod.addItem(self.tr(method), method, True)

        # update comboBoxMethod
        self.modelMethod.setItem(str_model=method)
        self.mdl.setMethod(method)

        self.updateReference(material, method)


    def updateReference(self, material, method):
        """
        update method list with material choice
        """

        self.comboBoxReference.hide()
        self.labelReference.hide()

        for nb in range(len(self.modelReference.getItems())):
            self.modelReference.delItem(0)

        reference = self.mdl.getReference()
        references = self.mdl.getAvailReferences(material, method)
        have_reference = False

        if references:
            for r in references:
                if reference == r:
                    have_reference = True
                self.modelReference.addItem(self.tr(r), r)
            self.comboBoxReference.show()
            self.labelReference.show()

            if not have_reference:
                reference = references[0]
            self.modelReference.setItem(str_model=reference)
            self.mdl.setReference(reference)


    @pyqtSlot(str)
    def slotPressure(self,  text):
        """
        Input PRESS.
        """
        if self.lineEditP0.validator().state == QValidator.Acceptable:
            p = from_qvariant(text, float)
            self.mdl.setPressure(p)


    @pyqtSlot(str)
    def slotTemperature(self,  text):
        """
        Input TEMPERATURE.
        """
        if self.lineEditT0.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTemperature(t)


    @pyqtSlot(str)
    def slotTempOxydant(self,  text):
        """
        Input oxydant TEMPERATURE.
        """
        if self.lineEditOxydant.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTempOxydant(t)


    @pyqtSlot(str)
    def slotTempFuel(self,  text):
        """
        Input fuel TEMPERATURE.
        """
        if self.lineEditFuel.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.mdl.setTempFuel(t)


    @pyqtSlot(str)
    def slotMassemol(self,  text):
        """
        Input Mass molar.
        """
        if self.lineEditMassMolar.validator().state == QValidator.Acceptable:
            m = from_qvariant(text, float)
            self.mdl.setMassemol(m)


    @pyqtSlot(str)
    def slotMaterial(self, text):
        """
        Method to call 'setMaterial'
        """
        choice = self.modelMaterial.dicoV2M[str(text)]
        old_choice = self.mdl.getMaterials()
        self.mdl.setMaterials(choice)
        self.updateMethod()
        self.updateTypeChoice(old_choice)


    @pyqtSlot(str)
    def slotMethod(self, text):
        """
        Method to call 'setMethod'
        """
        choice = self.modelMethod.dicoV2M[str(text)]
        self.mdl.setMethod(choice)
        self.updateReference(self.mdl.getMaterials(), choice)


    @pyqtSlot(str)
    def slotReference(self, text):
        """
        Method to call 'setReference'
        """
        choice = self.modelReference.dicoV2M[str(text)]
        self.mdl.setReference(choice)


    @pyqtSlot(str)
    def slotStateRho(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        self.__changeChoice(str(text), 'Rho', 'density')


    @pyqtSlot(str)
    def slotStateMu(self, text):
        """
        Method to call 'getState' with correct arguements for 'Mu'
        """
        self.__changeChoice(str(text), 'Mu', 'molecular_viscosity')


    @pyqtSlot(str)
    def slotStateSigma(self, text):
        """
        Method to call 'getState' with correct arguments for 'Sigma'
        """
        self.__changeChoice(str(text), 'Sigma', 'surface_tension')


    @pyqtSlot(str)
    def slotStateCp(self, text):
        """
        Method to call 'getState' with correct arguements for 'Cp'
        """
        self.__changeChoice(str(text), 'Cp', 'specific_heat')


    @pyqtSlot(str)
    def slotStateViscv0(self, text):
        """
        Method to call 'getState' with correct arguements for 'Viscv0'
        """
        self.__changeChoice(str(text), 'Viscv0', 'volume_viscosity')


    @pyqtSlot(str)
    def slotStateAl(self, text):
        """
        Method to call 'getState' with correct arguements for 'Al'
        """
        self.__changeChoice(str(text), 'Al', 'thermal_conductivity')


    @pyqtSlot(str)
    def slotStateDiff(self, text):
        """
        Method to set diffusion choice for the coefficient
        """
        choice = self.modelDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))

        if choice != 'user_law':
            self.pushButtonDiff.setEnabled(False)
            self.pushButtonDiff.setStyleSheet("background-color: None")
            self.pushButtonDiff.hide()
        else:
            self.pushButtonDiff.show()
            self.pushButtonDiff.setEnabled(True)
            name = self.mdl.m_sca.getScalarDiffusivityName(self.scalar)
            exp = self.mdl.m_sca.getDiffFormula(self.scalar, self.zone_id)
            if exp:
                self.pushButtonDiff.setStyleSheet("background-color: green")
                self.pushButtonDiff.setToolTip(exp)
            else:
                self.pushButtonDiff.setStyleSheet("background-color: red")

        self.mdl.m_sca.setScalarDiffusivityChoice(self.scalar, choice)


    @pyqtSlot(str)
    def slotNameDiff(self, text):
        """
        Method to set the variance scalar choosed
        """
        choice = self.modelNameDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))
        self.scalar = str(text)
        self.lineEditDiff.setText(str(self.mdl.m_sca.getScalarDiffusivityInitialValue(self.scalar)))

        mdl = self.mdl.m_sca.getScalarDiffusivityChoice(self.scalar)

        self.modelDiff.setItem(str_model=mdl)

        if  mdl!= 'user_law':
            self.pushButtonDiff.hide()
            self.pushButtonDiff.setEnabled(False)
            self.pushButtonDiff.setStyleSheet("background-color: None")
        else:
            self.pushButtonDiff.show()
            self.pushButtonDiff.setEnabled(True)
            name = self.mdl.m_sca.getScalarDiffusivityName(self.scalar)
            exp = self.mdl.m_sca.getDiffFormula(self.scalar, self.zone_id)
            if exp:
                self.pushButtonDiff.setStyleSheet("background-color: green")
                self.pushButtonDiff.setToolTip(exp)
            else:
                self.pushButtonDiff.setStyleSheet("background-color: red")


    def __changeChoice(self, text, sym, tag):
        """
        Input variable state
        """
        __model  = getattr(self, "model"      + sym)
        __line   = getattr(self, "lineEdit"   + sym)
        __combo  = getattr(self, "comboBox"   + sym)
        __label  = getattr(self, "label"      + sym)
        __button = getattr(self, "pushButton" + sym)
        __labelu = getattr(self, "labelUnit"  + sym)
        __labelv = getattr(self, "labelVar"  + sym)

        choice = __model.dicoV2M[text]
        log.debug("__changeChoice -> %s, %s" % (text, choice))

        if choice != 'user_law':
            __button.hide()
            __button.setEnabled(False)
            __button.setStyleSheet("background-color: None")
        else:
            __button.setEnabled(True)
            __button.show()
            exp = None
            if sym == "Rho":
                exp = self.mdl.getFormula('density', self.zone_id)
            elif sym == "Mu":
                exp = self.mdl.getFormula('molecular_viscosity', self.zone_id)
            elif sym == "Cp":
                exp = self.mdl.getFormula('specific_heat', self.zone_id)
            elif sym == "Viscv0":
                exp = self.mdl.getFormula('volume_viscosity', self.zone_id)
            elif sym == "Al":
                exp = self.mdl.getFormula('thermal_conductivity', self.zone_id)
            elif sym == "Diff":
                name = self.mdl.m_sca.getScalarDiffusivityName(self.scalar)
                exp = self.mdl.m_sca.getDiffFormula(self.scalar, self.zone_id)

            if exp:
                __button.setStyleSheet("background-color: green")
                __button.setToolTip(exp)
            else:
                __button.setStyleSheet("background-color: red")
        if choice == 'thermal_law':
            __line.hide()
            __label.hide()
            __labelu.hide()
            __labelv.hide()
        else:
            __line.show()
            __label.show()
            __labelu.show()
            __labelv.show()

        if self.is_main_zone:
            self.mdl.setPropertyMode(tag, choice)


    @pyqtSlot(str)
    def slotRho(self, text):
        """
        Update the density
        """
        if self.lineEditRho.validator().state == QValidator.Acceptable:
            rho = from_qvariant(text, float)
            self.mdl.setInitialValueDensity(rho)

    @pyqtSlot(str)
    def slotRho1(self, text):
        """
        Update the density of fluid 1 for VoF module
        """
        if self.lineEditRho1.validator().state == QValidator.Acceptable:
            rho = from_qvariant(text, float)
            self.mdl.setVofValueDensity(0, rho)

    @pyqtSlot(str)
    def slotRho2(self, text):
        """
        Update the density of fluid 2 for VoF module
        """
        if self.lineEditRho2.validator().state == QValidator.Acceptable:
            rho = from_qvariant(text, float)
            self.mdl.setVofValueDensity(1, rho)

    @pyqtSlot(str)
    def slotMu(self, text):
        """
        Update the molecular viscosity
        """
        if self.lineEditMu.validator().state == QValidator.Acceptable:
            mu = from_qvariant(text, float)
            self.mdl.setInitialValueViscosity(mu)


    @pyqtSlot(str)
    def slotMu1(self, text):
        """
        Update the molecular viscosity
        """
        if self.lineEditMu1.validator().state == QValidator.Acceptable:
            mu = from_qvariant(text, float)
            self.mdl.setVofValueViscosity(0, mu)


    @pyqtSlot(str)
    def slotMu2(self, text):
        """
        Update the molecular viscosity
        """
        if self.lineEditMu2.validator().state == QValidator.Acceptable:
            mu = from_qvariant(text, float)
            self.mdl.setVofValueViscosity(1, mu)


    @pyqtSlot(str)
    def slotSigma(self, text):
        """
        Update the surface tension
        """
        if self.lineEditSigma.validator().state == QValidator.Acceptable:
            sigma = from_qvariant(text, float)
            self.mdl.setInitialValueSurfaceTension(sigma)


    @pyqtSlot(str)
    def slotCp(self, text):
        """
        Update the specific heat
        """
        if self.lineEditCp.validator().state == QValidator.Acceptable:
            cp = from_qvariant(text, float)
            self.mdl.setInitialValueHeat(cp)


    @pyqtSlot(str)
    def slotViscv0(self, text):
        """
        Update the volumic viscosity
        """
        if self.lineEditViscv0.validator().state == QValidator.Acceptable:
            viscv0 = from_qvariant(text, float)
            self.mdl.setInitialValueVolumeViscosity(viscv0)


    @pyqtSlot(str)
    def slotAl(self, text):
        """
        Update the thermal conductivity
        """
        if self.lineEditAl.validator().state == QValidator.Acceptable:
            al = from_qvariant(text, float)
            self.mdl.setInitialValueCond(al)


    @pyqtSlot(str)
    def slotDiftl0(self, text):
        """
        Update the thermal conductivity
        """
        if self.lineEditDiftl0.validator().state == QValidator.Acceptable:
            diftl0 = from_qvariant(text, float)
            self.mdl.setInitialValueDyn(diftl0)


    @pyqtSlot(str)
    def slotDiff(self, text):
        """
        Update the thermal conductivity
        """
        if self.lineEditDiff.validator().state == QValidator.Acceptable:
            diff = from_qvariant(text, float)
            self.mdl.m_sca.setScalarDiffusivityInitialValue(self.scalar, diff)


    @pyqtSlot()
    def slotFormulaRho(self):
        """
        User formula for density
        """
        exp, req, sca, symbols_rho = \
            self.mdl.getFormulaPropertyComponents('density', self.zone_id)

        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarName()
        mdl = self.m_th.getThermalScalarModel()
        if mdl == "off":
            exa = FluidCharacteristicsView.density_wo
        elif mdl == "temperature_celsius":
            TempInContext = "("+s+" + 273.15)"
            exa = FluidCharacteristicsView.density.replace("temperature", TempInContext)
        elif mdl == "enthalpy":
            exa = FluidCharacteristicsView.density_h
        else:
            exa = FluidCharacteristicsView.density

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = 'density',
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_rho,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('density', str(result), self.zone_id)
            self.pushButtonRho.setToolTip(result)
            self.pushButtonRho.setStyleSheet("background-color: green")


    @pyqtSlot()
    def slotFormulaMu(self):
        """
        User formula for molecular viscosity
        """

        exp, req, sca, symbols_mu = \
                self.mdl.getFormulaPropertyComponents('molecular_viscosity', self.zone_id)

        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarName()
        mdl = self.m_th.getThermalScalarModel()
        if mdl == "off":
            exa = FluidCharacteristicsView.molecular_viscosity_wo
        elif mdl == "temperature_celsius":
            TempInContext = "("+s+" + 273.15)"
            exa = FluidCharacteristicsView.molecular_viscosity.replace("temperature", TempInContext)
        elif mdl == "enthalpy":
            exa = FluidCharacteristicsView.molecular_viscosity_h
        else:
            exa = FluidCharacteristicsView.molecular_viscosity

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = 'molecular_viscosity',
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_mu,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMu -> %s" % str(result))
            self.mdl.setFormula('molecular_viscosity', str(result), self.zone_id)
            self.pushButtonMu.setToolTip(result)
            self.pushButtonMu.setStyleSheet("background-color: green")


    @pyqtSlot()
    def slotFormulaCp(self):
        """
        User formula for specific heat
        """
        exp, req, sca, symbols_cp = \
                self.mdl.getFormulaPropertyComponents('specific_heat', self.zone_id)

        exa = FluidCharacteristicsView.specific_heat

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = 'specific_heat',
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_cp,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('specific_heat', str(result), self.zone_id)
            self.pushButtonCp.setToolTip(result)
            self.pushButtonCp.setStyleSheet("background-color: green")


    @pyqtSlot()
    def slotFormulaViscv0(self):
        """
        User formula for volumic viscosity
        """
        exp, req, sca, symbols_viscv0 = \
                self.mdl.getFormulaPropertyComponents('volume_viscosity', self.zone_id)

        exa = FluidCharacteristicsView.volume_viscosity

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = 'volume_viscosity',
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_viscv0,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaViscv0 -> %s" % str(result))
            self.mdl.setFormula('volume_viscosity', str(result), self.zone_id)
            self.pushButtonViscv0.setToolTip(result)
            self.pushButtonViscv0.setStyleSheet("background-color: green")


    @pyqtSlot()
    def slotFormulaAl(self):
        """
        User formula for thermal conductivity
        """
        exp, req, sca, symbols_al = \
                self.mdl.getFormulaPropertyComponents('thermal_conductivity', self.zone_id)

        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarName()
        mdl = self.m_th.getThermalScalarModel()
        if mdl == "temperature_celsius":
            TempInContext = "("+s+" + 273.15)"
            exa = FluidCharacteristicsView.thermal_conductivity.replace("temperature", TempInContext)
        elif mdl == "enthalpy":
            exa = FluidCharacteristicsView.thermal_conductivity_h
        else:
            exa = FluidCharacteristicsView.thermal_conductivity

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = 'thermal_conductivity',
                                expression    = exp,
                                required      = req,
                                symbols       = symbols_al,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaAl -> %s" % str(result))
            self.mdl.setFormula('thermal_conductivity', str(result), self.zone_id)
            self.pushButtonAl.setToolTip(result)
            self.pushButtonAl.setStyleSheet("background-color: green")


    @pyqtSlot()
    def slotFormulaDiff(self):
        """
        User formula for the diffusion coefficient
        """
        exp, req, sca, sym = \
                self.mdl.getFormulaDiffComponents(self.scalar, self.zone_id)

        exa = ''

        dname = self.mdl.m_sca.getScalarDiffusivityName(self.scalar)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'vol',
                                zone_name     = self.zone_name,
                                variable_name = dname,
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = sca,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDiff -> %s" % str(result))
            self.mdl.m_sca.setDiffFormula(self.scalar, str(result), self.zone_id)
            self.pushButtonDiff.setToolTip(result)
            self.pushButtonDiff.setStyleSheet("background-color: green")


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
