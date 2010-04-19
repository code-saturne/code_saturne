# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
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

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel, setGreenColor
from FluidCharacteristicsForm import Ui_FluidCharacteristicsForm
from FluidCharacteristicsModel import FluidCharacteristicsModel
from DefineUserScalarsModel import DefineUserScalarsModel

from QMeiEditorView import QMeiEditorView

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

rho = 1.293 * (273.15 / Temp_K);

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
mu0 = 18.27e-6;

if ( Temp_K > 0 && Temp_K < 555) {
mu = mu0 * (T0+CST / Temp_K+CST) * (Temp_K/T0)^(3./2.);
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
    thermal_conductivity="""# oxygen
lambda = 6.2e-5 * Temp_K + 8.1e-3;

# nitrogen
lambda = 6.784141e-5 * Temp_K + 5.564317e-3;

# hydrogen
lambda = 4.431e-4 * Temp_K + 5.334e-2;

"""
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_FluidCharacteristicsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = FluidCharacteristicsModel(self.case)

        list = [('density', 'Rho'),
                ('molecular_viscosity', 'Mu'),
                ('specific_heat', 'Cp'),
                ('thermal_conductivity', 'Al')]

        self.list_scalars = []
        m_sca = DefineUserScalarsModel(self.case)
        s = m_sca.getThermalScalarLabel()
        if s:
            self.list_scalars.append((s, self.tr("Thermal scalar")))
        for s in m_sca.getUserScalarLabelsList():
            self.list_scalars.append((s, self.tr("Additional scalar")))

        # Combo models

        self.modelRho = ComboModel(self.comboBoxRho, 3, 1)
        self.modelMu  = ComboModel(self.comboBoxMu, 3, 1)
        self.modelCp  = ComboModel(self.comboBoxCp, 3, 1)
        self.modelAl  = ComboModel(self.comboBoxAl, 3, 1)

        self.modelRho.addItem(self.tr('constant'), 'constant')
        self.modelRho.addItem(self.tr('user law'), 'user_law')
        self.modelRho.addItem(self.tr('user subroutine (usphyv)'), 'variable')
        self.modelMu.addItem(self.tr('constant'), 'constant')
        self.modelMu.addItem(self.tr('user law'), 'user_law')
        self.modelMu.addItem(self.tr('user subroutine (usphyv)'), 'variable')
        self.modelCp.addItem(self.tr('constant'), 'constant')
        self.modelCp.addItem(self.tr('user law'), 'user_law')
        self.modelCp.addItem(self.tr('user subroutine (usphyv)'), 'variable')
        self.modelAl.addItem(self.tr('constant'), 'constant')
        self.modelAl.addItem(self.tr('user law'), 'user_law')
        self.modelAl.addItem(self.tr('user subroutine (usphyv)'), 'variable')

        # Connections

        self.connect(self.comboBoxRho, SIGNAL("activated(const QString&)"), self.slotStateRho)
        self.connect(self.comboBoxMu, SIGNAL("activated(const QString&)"), self.slotStateMu)
        self.connect(self.comboBoxCp, SIGNAL("activated(const QString&)"), self.slotStateCp)
        self.connect(self.comboBoxAl, SIGNAL("activated(const QString&)"), self.slotStateAl)
        self.connect(self.lineEditRho, SIGNAL("textChanged(const QString &)"), self.slotRho)
        self.connect(self.lineEditMu, SIGNAL("textChanged(const QString &)"), self.slotMu)
        self.connect(self.lineEditCp, SIGNAL("textChanged(const QString &)"), self.slotCp)
        self.connect(self.lineEditAl, SIGNAL("textChanged(const QString &)"), self.slotAl)
        self.connect(self.pushButtonRho, SIGNAL("clicked()"), self.slotFormulaRho)
        self.connect(self.pushButtonMu, SIGNAL("clicked()"), self.slotFormulaMu)
        self.connect(self.pushButtonCp, SIGNAL("clicked()"), self.slotFormulaCp)
        self.connect(self.pushButtonAl, SIGNAL("clicked()"), self.slotFormulaAl)

        # Validators

        validatorRho = DoubleValidator(self.lineEditRho, min = 0.0)
        validatorMu = DoubleValidator(self.lineEditMu, min = 0.0)
        validatorCp = DoubleValidator(self.lineEditCp, min = 0.0)
        validatorAl = DoubleValidator(self.lineEditAl, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorMu.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)

        self.lineEditRho.setValidator(validatorRho)
        self.lineEditMu.setValidator(validatorMu)
        self.lineEditCp.setValidator(validatorCp)
        self.lineEditAl.setValidator(validatorAl)

        # Standard Widget initialization

        for tag, symbol in list:
            __model  = getattr(self, "model" + symbol)
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label" + symbol)
            c = self.mdl.getPropertyMode(tag)
            __model.setItem(str_model=c)
            if c == 'user_law':
                __button.setEnabled(True)
                __label.setText(QString(self.tr("Reference value")))
            else:
                __button.setEnabled(False)
                __label.setText(QString(self.tr("Reference value")))
            self.mdl.getInitialValue(tag)
            __line.setText(QString(str(self.mdl.getInitialValue(tag))))

        # Particular Widget initialization taking into account of "Calculation Features"

        mdl_atmo, mdl_joule, mdl_thermal, mdl_gas, mdl_coal = self.mdl.getThermoPhysicalModel()

        # no 'thermal_conductivity' if not Joule and not Thermal scalar and not
        if mdl_joule == 'off' and mdl_thermal == 'off' and mdl_atmo == 'off':
            self.groupBoxAl.hide()

        for tag, symbol in list:
            __model  = getattr(self, "model" + symbol)
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __label  = getattr(self, "label" + symbol)
            __combo  = getattr(self, "comboBox" + symbol)

            # Gas or coal combustion
            if mdl_gas != 'off' or mdl_coal != 'off':
                if tag == 'density':
                    __model.setItem(str_model='user_law')
                    __combo.setEnabled(False)
                    __button.setEnabled(False)
                    self.mdl.setPropertyMode(tag, 'user_law')
                    __label.setText(QString(self.tr("Calculation by\n perfect gas law")))
                    __line.setText(QString(str("")))
                    __line.setEnabled(False)
                else:
                    __model.setItem(str_model='constant')
                    self.mdl.setPropertyMode(tag, 'constant')

            # Joule
            if mdl_joule != 'off':
                __model.setItem(str_model='user_law')
                __model.disableItem(str_model='constant')
                self.mdl.setPropertyMode(name, 'user_law')

            # Atmospheric Flows
            if mdl_atmo != 'off':
                if tag == 'density':
                    __model.disableItem(str_model='constant')


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
    def slotStateAl(self, text):
        """
        Method to call 'getState' with correct arguements for 'Al'
        """
        self.__changeChoice(str(text), 'Al', 'thermal_conductivity')


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
        rho, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setInitialValueDensity(rho)


    @pyqtSignature("const QString &")
    def slotMu(self, text):
        """
        Update the molecular viscosity
        """
        mu, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setInitialValueViscosity(mu)


    @pyqtSignature("const QString &")
    def slotCp(self, text):
        """
        Update the specific heat
        """
        cp, ok = self.lineEditCp.text().toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setInitialValueHeat(cp)


    @pyqtSignature("const QString &")
    def slotAl(self, text):
        """
        Update the thermal conductivity
        """
        al, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setInitialValueCond(al)


    @pyqtSignature("")
    def slotFormulaRho(self):
        """
        User formula for density
        """
        exp = self.mdl.getFormula('density')
        if not exp:
            exp = "rho ="
        req = [('rho', 'Density')]
        exa = FluidCharacteristicsView.density
        setGreenColor(self.sender(), False)

        dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = self.list_scalars,
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
        if not exp:
            exp = "mu ="
        req = [('mu', 'Molecular Viscosity')]
        exa = FluidCharacteristicsView.molecular_viscosity

        dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = self.list_scalars,
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
        if not exp:
            exp = "cp ="
        req = [('cp', 'Specific heat')]
        exa = FluidCharacteristicsView.specific_heat

        dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = self.list_scalars,
                                     examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('specific_heat', result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("")
    def slotFormulaAl(self):
        """
        User formula for thermal conductivity
        """
        exp = self.mdl.getFormula('thermal_conductivity')
        if not exp:
            exp = "lambda ="
        req = [('lambda', 'Thermal conductivity')]
        exa = FluidCharacteristicsView.thermal_conductivity

        dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = self.list_scalars,
                                     examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaRho -> %s" % str(result))
            self.mdl.setFormula('thermal_conductivity', result)
            setGreenColor(self.sender(), False)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
