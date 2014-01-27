# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
import sys
if sys.version_info[0] == 2:
    import sip
    sip.setapi('QString', 2)

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel, setGreenColor
from Pages.FluidCharacteristicsForm import Ui_FluidCharacteristicsForm
from Pages.FluidCharacteristicsModel import FluidCharacteristicsModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.ReferenceValuesModel import ReferenceValuesModel
from Pages.CompressibleModel import CompressibleModel
from Pages.CoalCombustionModel import CoalCombustionModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.QMeiEditorView import QMeiEditorView

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
rho = 1.293*(273.15 / TempK);


# density for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

rho1 = 1.25051;
rho2 = 1.7832;
A = (Y1 / rho1) + (Y2 /rho2);
rho = 1.0 / A;

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

if ( TempK > 0 && TempK < 555) {
mu = mu_ref * (T0+CST / TempK+CST) * (TempK/T0)^(3./2.);
} else {
mu = -999.0;
}

"""
    specific_heat="""# specific heat for mixtures of gases
#
# Y1 -> mass fraction of component 1
# Y2 -> mass fraction of component 2

Cp1 = 520.3;
Cp2 = 1040.0;
cp = Y1 * Cp1 + Y2 *Cp2;
"""
    volumic_viscosity="""# volumic_viscosity
"""
    thermal_conductivity="""# oxygen
lambda = 6.2e-5 * TempK + 8.1e-3;

# nitrogen
lambda = 6.784141e-5 * TempK + 5.564317e-3;

# hydrogen
lambda = 4.431e-4 * TempK + 5.334e-2;

"""
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_FluidCharacteristicsForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.case.undoStopGlobal()

        self.mdl = FluidCharacteristicsModel(self.case)

        if CompressibleModel(self.case).getCompressibleModel() != 'off':
            list = [('density', 'Rho'),
                    ('molecular_viscosity', 'Mu'),
                    ('specific_heat', 'Cp'),
                    ('thermal_conductivity', 'Al'),
                    ('volumic_viscosity', 'Viscv0'),
                    ('dynamic_diffusion', 'Diftl0')]
        elif CoalCombustionModel(self.case).getCoalCombustionModel() != 'off':
            list = [('density', 'Rho'),
                    ('molecular_viscosity', 'Mu'),
                    ('specific_heat', 'Cp'),
                    ('dynamic_diffusion', 'Diftl0')]
        elif GasCombustionModel(self.case).getGasCombustionModel() != 'off':
            list = [('density', 'Rho'),
                    ('molecular_viscosity', 'Mu'),
                    ('specific_heat', 'Cp'),
                    ('dynamic_diffusion', 'Diftl0')]
        else:
            list = [('density', 'Rho'),
                    ('molecular_viscosity', 'Mu'),
                    ('specific_heat', 'Cp'),
                    ('thermal_conductivity', 'Al')]

        self.list_scalars = []
        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarLabel()
        if s:
            self.list_scalars.append((s, self.tr("Thermal scalar")))

        self.m_sca = DefineUserScalarsModel(self.case)
        for s in self.m_sca.getUserScalarLabelsList():
            self.list_scalars.append((s, self.tr("Additional scalar")))

        # Particular Widget initialization taking into account of "Calculation Features"
        mdl_atmo, mdl_joule, mdl_thermal, mdl_gas, mdl_coal, mdl_comp = self.mdl.getThermoPhysicalModel()

        # Combo models

        self.modelRho      = ComboModel(self.comboBoxRho, 3, 1)
        self.modelMu       = ComboModel(self.comboBoxMu, 3, 1)
        self.modelCp       = ComboModel(self.comboBoxCp, 3, 1)
        self.modelAl       = ComboModel(self.comboBoxAl, 3, 1)
        self.modelDiff     = ComboModel(self.comboBoxDiff, 2, 1)
        self.modelNameDiff = ComboModel(self.comboBoxNameDiff,1,1)
        self.modelViscv0   = ComboModel(self.comboBoxViscv0, 3, 1)
        self.modelDiftl0   = ComboModel(self.comboBoxDiftl0, 3, 1)

        self.modelRho.addItem(self.tr('constant'), 'constant')
        self.modelRho.addItem(self.tr('user law'), 'user_law')
        if mdl_atmo != 'off':
            self.modelRho.addItem(self.tr('defined in atphyv'), 'variable')
        elif mdl_joule == 'arc':
            self.modelRho.addItem(self.tr('defined in elphyv'), 'variable')
        else:
            self.modelRho.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')

        self.modelMu.addItem(self.tr('constant'), 'constant')
        self.modelMu.addItem(self.tr('user law'), 'user_law')
        if mdl_joule == 'arc':
            self.modelMu.addItem(self.tr('defined in elphyv'), 'variable')
        else:
            self.modelMu.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')

        self.modelCp.addItem(self.tr('constant'), 'constant')
        self.modelCp.addItem(self.tr('user law'), 'user_law')
        if mdl_joule == 'arc':
            self.modelCp.addItem(self.tr('defined in elphyv'), 'variable')
        else:
            self.modelCp.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')

        self.modelAl.addItem(self.tr('constant'), 'constant')
        self.modelAl.addItem(self.tr('user law'), 'user_law')
        if mdl_joule == 'arc':
            self.modelAl.addItem(self.tr('defined in elphyv'), 'variable')
        else:
            self.modelAl.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')

        self.modelDiff.addItem(self.tr('constant'), 'constant')
        self.modelDiff.addItem(self.tr('user law'), 'user_law')
        self.modelViscv0.addItem(self.tr('constant'), 'constant')
        self.modelViscv0.addItem(self.tr('user law'), 'user_law')
        self.modelViscv0.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')
        self.modelDiftl0.addItem(self.tr('constant'), 'constant')
        self.modelDiftl0.addItem(self.tr('user law'), 'user_law')
        self.modelDiftl0.addItem(self.tr('user subroutine (cs_user_physical_properties)'), 'variable')

        self.scalar = ""
        scalar_list = self.m_sca.getUserScalarLabelsList()
        for s in self.m_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)

        if scalar_list != []:
            self.scalar = scalar_list[0]
            for scalar in scalar_list:
                self.modelNameDiff.addItem(scalar)

        # Connections

        self.connect(self.comboBoxRho,      SIGNAL("activated(const QString&)"), self.slotStateRho)
        self.connect(self.comboBoxMu,       SIGNAL("activated(const QString&)"), self.slotStateMu)
        self.connect(self.comboBoxCp,       SIGNAL("activated(const QString&)"), self.slotStateCp)
        self.connect(self.comboBoxAl,       SIGNAL("activated(const QString&)"), self.slotStateAl)
        self.connect(self.comboBoxDiff,     SIGNAL("activated(const QString&)"), self.slotStateDiff)
        self.connect(self.comboBoxNameDiff, SIGNAL("activated(const QString&)"), self.slotNameDiff)
        self.connect(self.comboBoxViscv0,   SIGNAL("activated(const QString&)"), self.slotStateViscv0)
        self.connect(self.lineEditRho,      SIGNAL("textChanged(const QString &)"), self.slotRho)
        self.connect(self.lineEditMu,       SIGNAL("textChanged(const QString &)"), self.slotMu)
        self.connect(self.lineEditCp,       SIGNAL("textChanged(const QString &)"), self.slotCp)
        self.connect(self.lineEditAl,       SIGNAL("textChanged(const QString &)"), self.slotAl)
        self.connect(self.lineEditDiff,     SIGNAL("textChanged(const QString &)"), self.slotDiff)
        self.connect(self.lineEditDiftl0,   SIGNAL("textChanged(const QString &)"), self.slotDiftl0)
        self.connect(self.lineEditViscv0,   SIGNAL("textChanged(const QString &)"), self.slotViscv0)
        self.connect(self.pushButtonRho,    SIGNAL("clicked()"), self.slotFormulaRho)
        self.connect(self.pushButtonMu,     SIGNAL("clicked()"), self.slotFormulaMu)
        self.connect(self.pushButtonCp,     SIGNAL("clicked()"), self.slotFormulaCp)
        self.connect(self.pushButtonAl,     SIGNAL("clicked()"), self.slotFormulaAl)
        self.connect(self.pushButtonDiff,   SIGNAL("clicked()"), self.slotFormulaDiff)
        self.connect(self.pushButtonViscv0, SIGNAL("clicked()"), self.slotFormulaViscv0)

        # Validators

        validatorRho    = DoubleValidator(self.lineEditRho, min = 0.0)
        validatorMu     = DoubleValidator(self.lineEditMu, min = 0.0)
        validatorCp     = DoubleValidator(self.lineEditCp, min = 0.0)
        validatorAl     = DoubleValidator(self.lineEditAl, min = 0.0)
        validatorDiff   = DoubleValidator(self.lineEditDiff, min = 0.0)
        validatorViscv0 = DoubleValidator(self.lineEditViscv0, min = 0.0)
        validatorDiftl0 = DoubleValidator(self.lineEditDiftl0, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorMu.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)
        validatorDiff.setExclusiveMin(True)
        validatorDiftl0.setExclusiveMin(True)

        self.lineEditRho.setValidator(validatorRho)
        self.lineEditMu.setValidator(validatorMu)
        self.lineEditCp.setValidator(validatorCp)
        self.lineEditAl.setValidator(validatorAl)
        self.lineEditDiff.setValidator(validatorDiff)
        self.lineEditViscv0.setValidator(validatorViscv0)
        self.lineEditDiftl0.setValidator(validatorDiftl0)

        if scalar_list == []:
            self.groupBoxDiff.hide()
        else :
            self.groupBoxDiff.show()
            self.lineEditDiff.setText(str(self.m_sca.getScalarDiffusivityInitialValue(self.scalar)))

            diff_choice =  self.m_sca.getScalarDiffusivityChoice(self.scalar)
            self.modelDiff.setItem(str_model=diff_choice)
            self.modelNameDiff.setItem(str_model=str(self.scalar))
            if diff_choice  != 'user_law':
                self.pushButtonDiff.setEnabled(False)
                setGreenColor(self.pushButtonDiff, False)
            else:
                self.pushButtonDiff.setEnabled(True)
                setGreenColor(self.pushButtonDiff, True)

        #compressible
        self.groupBoxViscv0.hide()

        # combustion
        self.groupBoxDiftl0.hide()

        # Standard Widget initialization

        for tag, symbol in list:
            __model  = getattr(self, "model" + symbol)
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label" + symbol)
            if tag != 'dynamic_diffusion':
                c = self.mdl.getPropertyMode(tag)
                __model.setItem(str_model=c)
                if c == 'user_law':
                    __button.setEnabled(True)
                    __label.setText(self.tr("Reference value"))
                else:
                    __button.setEnabled(False)
                    __label.setText(self.tr("Reference value"))
            else:
                __label.setText(self.tr("Reference value"))

            self.mdl.getInitialValue(tag)
            __line.setText(str(self.mdl.getInitialValue(tag)))

        # no 'thermal_conductivity' if not Joule and not Thermal scalar and not
        if mdl_joule == 'off' and mdl_thermal == 'off' and mdl_atmo == 'off' and\
           CompressibleModel(self.case).getCompressibleModel() == 'off':
            self.groupBoxAl.hide()

        if mdl_gas != 'off' or mdl_coal != 'off':
            self.groupBoxDiftl0.show()

        for tag, symbol in list:
            __model  = getattr(self, "model" + symbol)
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label" + symbol)
            __combo  = getattr(self, "comboBox" + symbol)

            # Gas or coal combustion
            if mdl_gas != 'off' or mdl_coal != 'off':
                if tag == 'density':
                    __model.setItem(str_model='variable')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    self.mdl.setPropertyMode(tag, 'variable')
                    __label.setText(self.tr("Calculation by\n perfect gas law"))
                    __line.setText(str(""))
                    __line.setEnabled(False)
                elif tag == 'dynamic_diffusion':
                    __model.setItem(str_model='user_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                else:
                    __model.setItem(str_model='constant')
                    self.mdl.setPropertyMode(tag, 'constant')

            # Joule
            if mdl_joule == 'arc':
                __model.disableItem(str_model='constant')
                __model.disableItem(str_model='user_law')
                __model.disableItem(str_model='variable')
                __model.setItem(str_model='variable')
                __combo.setEnabled(False)
                __button.setEnabled(False)
                self.mdl.setPropertyMode(tag, 'variable')
            if mdl_joule == 'joule':
                __model.setItem(str_model='variable')
                __model.disableItem(str_model='constant')
                self.mdl.setPropertyMode(tag, 'variable')

            # Atmospheric Flows
            if mdl_atmo != 'off':
                if tag == 'density':
                    __model.disableItem(str_model='constant')
                    __model.disableItem(str_model='user_law')
                    __model.disableItem(str_model='variable')
                    __model.setItem(str_model='variable')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)

            # Compressible Flows
            if mdl_comp != 'off':
                if tag == 'density':
                    __model.setItem(str_model='variable')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    __combo.hide()
                    __button.hide()
                    self.mdl.setPropertyMode(tag, 'variable')
                    __line.setEnabled(True)
                self.groupBoxViscv0.hide()
                if tag == 'specific_heat':
                    __model.setItem(str_model='constant')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    self.mdl.setPropertyMode(tag, 'constant')
                    self.groupBoxCp.setTitle('Isobaric specific heat')

                if tag == 'volumic_viscosity':
                    __combo.setEnabled(True)
                    c = self.mdl.getPropertyMode(tag)
                    if c == 'user_law':
                        __button.setEnabled(True)
                    else:
                        __button.setEnabled(False)
                self.groupBoxViscv0.show()
            else:
                if tag == 'specific_heat':
                    self.groupBoxCp.setTitle('Specific heat')

        self.case.undoStartGlobal()


    @pyqtSignature("const QString &")
    def slotStateRho(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        self.__changeChoice(str(text), 'Rho', 'density')


    @pyqtSignature("const QString &")
    def slotStateMu(self, text):
        """
        Method to call 'getState' with correct arguements for 'Mu'
        """
        self.__changeChoice(str(text), 'Mu', 'molecular_viscosity')


    @pyqtSignature("const QString &")
    def slotStateCp(self, text):
        """
        Method to call 'getState' with correct arguements for 'Cp'
        """
        self.__changeChoice(str(text), 'Cp', 'specific_heat')


    @pyqtSignature("const QString &")
    def slotStateViscv0(self, text):
        """
        Method to call 'getState' with correct arguements for 'Viscv0'
        """
        self.__changeChoice(str(text), 'Viscv0', 'volumic_viscosity')


    @pyqtSignature("const QString &")
    def slotStateAl(self, text):
        """
        Method to call 'getState' with correct arguements for 'Al'
        """
        self.__changeChoice(str(text), 'Al', 'thermal_conductivity')


    @pyqtSignature("const QString &")
    def slotStateDiff(self, text):
        """
        Method to set diffusion choice for the coefficient
        """
        choice = self.modelDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))

        if choice != 'user_law':
            self.pushButtonDiff.setEnabled(False)
            setGreenColor(self.pushButtonDiff, False)
        else:
            self.pushButtonDiff.setEnabled(True)
            setGreenColor(self.pushButtonDiff, True)

        self.m_sca.setScalarDiffusivityChoice(self.scalar, choice)


    @pyqtSignature("const QString &")
    def slotNameDiff(self, text):
        """
        Method to set the variance scalar choosed
        """
        choice = self.modelNameDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))
        self.scalar = str(text)
        self.lineEditDiff.setText(str(self.m_sca.getScalarDiffusivityInitialValue(self.scalar)))

        self.modelDiff.setItem(str_model=self.m_sca.getScalarDiffusivityChoice(self.scalar))


    def __changeChoice(self, text, sym, tag):
        """
        Input variable state
        """
        __model  = getattr(self, "model"      + sym)
        __line   = getattr(self, "lineEdit"   + sym)
        __combo  = getattr(self, "comboBox"   + sym)
        __label  = getattr(self, "label"      + sym)
        __button = getattr(self, "pushButton" + sym)

        choice = __model.dicoV2M[text]
        log.debug("__changeChoice -> %s, %s" % (text, choice))

        if choice != 'user_law':
            __button.setEnabled(False)
            setGreenColor(__button, False)
        else:
            __button.setEnabled(True)
            setGreenColor(__button, True)

        self.mdl.setPropertyMode(tag, choice)


    @pyqtSignature("const QString &")
    def slotRho(self, text):
        """
        Update the density
        """
        if self.sender().validator().state == QValidator.Acceptable:
            rho = float(text)
            self.mdl.setInitialValueDensity(rho)


    @pyqtSignature("const QString &")
    def slotMu(self, text):
        """
        Update the molecular viscosity
        """
        if self.sender().validator().state == QValidator.Acceptable:
            mu = float(text)
            self.mdl.setInitialValueViscosity(mu)


    @pyqtSignature("const QString &")
    def slotCp(self, text):
        """
        Update the specific heat
        """
        if self.sender().validator().state == QValidator.Acceptable:
            cp = float(self.lineEditCp.text())
            self.mdl.setInitialValueHeat(cp)


    @pyqtSignature("const QString &")
    def slotViscv0(self, text):
        """
        Update the volumic viscosity
        """
        if self.sender().validator().state == QValidator.Acceptable:
            viscv0 = float(self.lineEditViscv0.text())
            self.mdl.setInitialValueVolumicViscosity(viscv0)


    @pyqtSignature("const QString &")
    def slotAl(self, text):
        """
        Update the thermal conductivity
        """
        if self.sender().validator().state == QValidator.Acceptable:
            al = float(text)
            self.mdl.setInitialValueCond(al)


    @pyqtSignature("const QString &")
    def slotDiftl0(self, text):
        """
        Update the thermal conductivity
        """
        if self.sender().validator().state == QValidator.Acceptable:
            diftl0 = float(text)
            self.mdl.setInitialValueDyn(diftl0)


    @pyqtSignature("const QString &")
    def slotDiff(self, text):
        """
        Update the thermal conductivity
        """
        if self.sender().validator().state == QValidator.Acceptable:
            diff = float(text)
            self.m_sca.setScalarDiffusivityInitialValue(self.scalar, diff)


    @pyqtSignature("")
    def slotFormulaRho(self):
        """
        User formula for density
        """
        exp = self.mdl.getFormula('density')
        req = [('rho', 'Density')]
        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarLabel()
        if s == "TempC":
            TempInContext = "("+s+" - 273.15)"
            exa = FluidCharacteristicsView.density.replace("TempK", TempInContext)
        else:
            exa = FluidCharacteristicsView.density

        setGreenColor(self.sender(), False)

        symbols_rho = []
        for s in self.list_scalars:
           symbols_rho.append(s)
        rho0_value = self.mdl.getInitialValueDensity()
        ref_pressure = ReferenceValuesModel(self.case).getPressure()
        symbols_rho.append(('rho0', 'Density (reference value) = ' + str(rho0_value)))
        symbols_rho.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_rho,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('density', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaMu(self):
        """
        User formula for molecular viscosity
        """
        exp = self.mdl.getFormula('molecular_viscosity')
        req = [('mu', 'Molecular Viscosity')]
        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarLabel()
        if s == "TempC":
            TempInContext = "("+s+" - 273.15)"
            exa = FluidCharacteristicsView.molecular_viscosity.replace("TempK", TempInContext)
        else:
            exa = FluidCharacteristicsView.molecular_viscosity

        symbols_mu = []
        for s in self.list_scalars:
           symbols_mu.append(s)
        mu0_value = self.mdl.getInitialValueViscosity()
        rho0_value = self.mdl.getInitialValueDensity()
        ref_pressure = ReferenceValuesModel(self.case).getPressure()
        symbols_mu.append(('mu0', 'Viscosity (reference value) = ' + str(mu0_value)))
        symbols_mu.append(('rho0', 'Density (reference value) = ' + str(rho0_value)))
        symbols_mu.append(('p0', 'Reference pressure = ' + str(ref_pressure)))
        symbols_mu.append(('rho', 'Density'))
        if CompressibleModel(self.case).getCompressibleModel() == 'on':
            symbols_mu.append(('T', 'Temperature'))
            ref_temperature = ReferenceValuesModel(self.case).getTemperature()
            symbols_mu.append(('t0', 'Reference temperature = '+str(ref_temperature)+' K'))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_mu,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMu -> %s" % str(result))
            self.mdl.setFormula('molecular_viscosity', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaCp(self):
        """
        User formula for specific heat
        """
        exp = self.mdl.getFormula('specific_heat')
        req = [('cp', 'Specific heat')]
        exa = FluidCharacteristicsView.specific_heat

        symbols_cp = []
        for s in self.list_scalars:
           symbols_cp.append(s)
        cp0_value = self.mdl.getInitialValueHeat()
        ref_pressure = ReferenceValuesModel(self.case).getPressure()
        symbols_cp.append(('cp0', 'Specific heat (reference value) = ' + str(cp0_value)))
        symbols_cp.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_cp,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('specific_heat', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaViscv0(self):
        """
        User formula for volumic viscosity
        """
        exp = self.mdl.getFormula('volumic_viscosity')
        req = [('viscv', 'Volumic viscosity')]
        exa = FluidCharacteristicsView.volumic_viscosity
        symbols_viscv0 = []
        for s in self.list_scalars:
           symbols_cp.append(s)
        viscv0_value = self.mdl.getInitialValueVolumicViscosity()
        ref_pressure = ReferenceValuesModel(self.case).getPressure()
        ref_temperature = ReferenceValuesModel(self.case).getTemperature()
        symbols_viscv0.append(('viscv0', 'Volumic viscosity (reference value) = '+str(viscv0_value)+' J/kg/K'))
        symbols_viscv0.append(('p0', 'Reference pressure = '+str(ref_pressure)+' Pa'))
        symbols_viscv0.append(('t0', 'Reference temperature = '+str(ref_temperature)+' K'))
        symbols_viscv0.append(('T', 'Temperature'))
        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_viscv0,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaViscv0 -> %s" % str(result))
            self.mdl.setFormula('volumic_viscosity', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaAl(self):
        """
        User formula for thermal conductivity
        """
        exp = self.mdl.getFormula('thermal_conductivity')
        req = [('lambda', 'Thermal conductivity')]
        self.m_th = ThermalScalarModel(self.case)
        s = self.m_th.getThermalScalarLabel()
        if s == "TempC":
            TempInContext = "("+s+" - 273.15)"
            exa = FluidCharacteristicsView.thermal_conductivity.replace("TempK", TempInContext)
        else:
            exa = FluidCharacteristicsView.thermal_conductivity

        symbols_al = []
        for s in self.list_scalars:
           symbols_al.append(s)
        lambda0_value = self.mdl.getInitialValueCond()
        ref_pressure = ReferenceValuesModel(self.case).getPressure()
        symbols_al.append(('lambda0', 'Thermal conductivity (reference value) = ' + str(lambda0_value)))
        symbols_al.append(('p0', 'Reference pressure = ' + str(ref_pressure)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbols_al,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaAl -> %s" % str(result))
            self.mdl.setFormula('thermal_conductivity', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaDiff(self):
        """
        User formula for the diffusion coefficient
        """
        name = self.m_sca.getScalarDiffusivityName(self.scalar)
        exp = self.m_sca.getDiffFormula(self.scalar)
        req = [(str(name), str(self.scalar)+'diffusion coefficient')]
        exa = ''
        sym = [('x','cell center coordinate'),
               ('y','cell center coordinate'),
               ('z','cell center coordinate'),]
        sym.append((str(self.scalar),str(self.scalar)))
        diff0_value = self.m_sca.getScalarDiffusivityInitialValue(self.scalar)
        sym.append((str(name)+'_ref', str(self.scalar)+' diffusion coefficient (reference value) = '+str(diff0_value)+' m^2/s'))
        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDiff -> %s" % str(result))
            self.m_sca.setDiffFormula(self.scalar, result)
            setGreenColor(self.pushButtonDiff, False)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
