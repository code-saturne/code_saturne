# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
This module contains the following classes:
- BoundaryConditionsTurbulenceInletView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BoundaryConditionsTurbulenceInletForm import Ui_BoundaryConditionsTurbulenceInletForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.Pages.TurbulenceNeptuneModel import TurbulenceModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsTurbulenceInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsTurbulenceInletView(QWidget, Ui_BoundaryConditionsTurbulenceInletForm) :
    """
    Boundary condition for turbulence
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsTurbulenceInletForm.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxTurbulence.activated[str].connect(self.__slotChoiceTurbulence)

        self.__modelTurbulence = ComboModel(self.comboBoxTurbulence, 3, 1)
        self.__modelTurbulence.addItem(self.tr("Calculation by hydraulic diameter"), 'hydraulic_diameter')
        self.__modelTurbulence.addItem(self.tr("Calculation by turbulent intensity"), 'turbulent_intensity')
        self.__modelTurbulence.addItem(self.tr("Calculation by formula"), 'formula')

        self.lineEditDiameter.textChanged[str].connect(self.__slotDiam)
        self.lineEditIntensity.textChanged[str].connect(self.__slotIntensity)
        self.lineEditDiameterIntens.textChanged[str].connect(self.__slotDiam)
        self.pushButtonTurb.clicked.connect(self.__slotTurbulenceFormula)

        validatorDiam = DoubleValidator(self.lineEditDiameter, min=0.)
        validatorDiam.setExclusiveMin(True)
        validatorIntensity = DoubleValidator(self.lineEditIntensity, min=0.)

        self.lineEditDiameter.setValidator(validatorDiam)
        self.lineEditDiameterIntens.setValidator(validatorDiam)
        self.lineEditIntensity.setValidator(validatorIntensity)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None
        self.__currentField = fieldId


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        turb_model =  TurbulenceModel(self.__case).getTurbulenceModel(self.__currentField)
        self.__modelTurbulence.enableItem(0)
        self.__modelTurbulence.enableItem(1)
        if turb_model != "none" and turb_model != 'mixing_length' and turb_model != 'tchen' and turb_model != 'r2-r12-tchen':
            turb_choice = boundary.getTurbulenceChoice(self.__currentField)
            self.__modelTurbulence.setItem(str_model=turb_choice)
            if turb_model == 'q2-q12' or turb_model == 'r2-q12':
                self.__modelTurbulence.disableItem(0)
                self.__modelTurbulence.disableItem(1)
            if turb_choice == "hydraulic_diameter":
                self.frameTurbDiameter.show()
                self.frameTurbIntensity.hide()
                d = boundary.getHydraulicDiameter(self.__currentField)
                self.lineEditDiameter.setText(str(d))
                self.pushButtonTurb.setEnabled(False)
            elif turb_choice == "turbulent_intensity":
                self.frameTurbIntensity.show()
                self.frameTurbDiameter.hide()
                i = boundary.getTurbulentIntensity(self.__currentField)
                d = boundary.getHydraulicDiameter(self.__currentField)
                self.lineEditIntensity.setText(str(i))
                self.lineEditDiameterIntens.setText(str(d))
                self.pushButtonTurb.setEnabled(False)
            elif turb_choice == "formula":
                self.frameTurbDiameter.hide()
                self.frameTurbIntensity.hide()
                self.pushButtonTurb.setEnabled(True)
                exp = self.__boundary.getTurbFormula(self.__currentField)
                if exp:
                    self.pushButtonTurb.setStyleSheet("background-color: green")
                    self.pushButtonTurb.setToolTip(exp)
                else:
                    self.pushButtonTurb.setStyleSheet("background-color: red")
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotChoiceTurbulence(self, text):
        """
        INPUT choice of method of calculation of the turbulence
        """
        turb_choice = self.__modelTurbulence.dicoV2M[str(text)]
        self.__boundary.setTurbulenceChoice(self.__currentField, turb_choice)

        self.frameTurbDiameter.hide()
        self.frameTurbIntensity.hide()
        self.pushButtonTurb.setEnabled(False)

        if turb_choice  == 'hydraulic_diameter':
            self.frameTurbDiameter.show()
            d = self.__boundary.getHydraulicDiameter(self.__currentField)
            self.lineEditDiameter.setText(str(d))
        elif turb_choice == 'turbulent_intensity':
            self.frameTurbIntensity.show()
            i = self.__boundary.getTurbulentIntensity(self.__currentField)
            self.lineEditIntensity.setText(str(i))
            d = self.__boundary.getHydraulicDiameter(self.__currentField)
            self.lineEditDiameterIntens.setText(str(d))
        elif turb_choice == "formula":
            self.pushButtonTurb.setEnabled(True)
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if exp:
                self.pushButtonTurb.setStyleSheet("background-color: green")
                self.pushButtonTurb.setToolTip(exp)
            else:
                self.pushButtonTurb.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def __slotDiam(self, text):
        """
        INPUT hydraulic diameter
        """
        if self.lineEditDiameter.validator().state == QValidator.Acceptable:
            diam = from_qvariant(text, float)
            self.__boundary.setHydraulicDiameter(self.__currentField, diam)


    @pyqtSlot(str)
    def __slotIntensity(self, text):
        """
        INPUT turbulent intensity
        """
        if self.lineEditIntensity.validator().state == QValidator.Acceptable:
            intens = from_qvariant(text, float)
            self.__boundary.setTurbulentIntensity(self.__currentField, intens)


    @pyqtSlot()
    def __slotTurbulenceFormula(self):
        """
        User formula for turbulence
        """
        turb_model =  TurbulenceModel(self.__case).getTurbulenceModel(self.__currentField)

        if turb_model in ('k-epsilon', 'k-epsilon_linear_production'):
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if not exp:
                exp = self.__boundary.getDefaultTurbFormula(turb_model)
            exa = """#example :
uref2 = 10.;
dh = 0.2;
re = sqrt(uref2)*dh*rho0/mu0;

if (re < 2000){
#     in this case u*^2 is directly calculated to not have a problem with
#     xlmbda=64/Re when Re->0

  ustar2 = 8.*mu0*sqrt(uref2)/rho0/dh;}

else if (re<4000){

  xlmbda = 0.021377 + 5.3115e-6*re;
  ustar2 = uref2*xlmbda/8.;}

else {

  xlmbda = 1/( 1.8*log(re)/log(10.)-1.64)^2;
  ustar2 = uref2*xlmbda/8.;}

cmu = 0.09;
kappa = 0.42;
k   = ustar2/sqrt(cmu);
eps = ustar2^1.5/(kappa*dh*0.1);"""
            req = [('k', "turbulent energy"),
            ('eps', "turbulent dissipation")]
            sym = [('x','face center coordinate'),
                   ('y','face center coordinate'),
                   ('z','face center coordinate'),
                   ('t','time'),
                   ('dt','time step'),
                   ('iter','number of time step')]
            dialog = QMeiEditorView(self,
                                    check_syntax = self.__case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.__boundary.setTurbFormula(self.__currentField, result)
                self.pushButtonTurb.setToolTip(result)
                self.pushButtonTurb.setStyleSheet("background-color: green")

        elif turb_model in ('rij-epsilon_ssg', 'rij-epsilon_ebrsm'):
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if not exp:
                exp = self.__boundary.getDefaultTurbFormula(turb_model)
            exa = """#exemple :
uref2 = 10.;
dh = 0.2;
re = sqrt(uref2)*dh*rho0/mu0;

if (re < 2000){
#     in this case u*^2 is directly calculated to not have a problem with
#     xlmbda=64/Re when Re->0

  ustar2 = 8.*mu0*sqrt(uref2)/rho0/dh;}

else if (re<4000){

  xlmbda = 0.021377 + 5.3115e-6*re;
  ustar2 = uref2*xlmbda/8.;}

else {

  xlmbda = 1/( 1.8*log(re)/log(10.)-1.64)^2;
  ustar2 = uref2*xlmbda/8.;}

cmu = 0.09;
kappa = 0.42;
k   = ustar2/sqrt(cmu);
eps = ustar2^1.5/(kappa*dh*0.1);
d2s3 = 2/3;
R11 = d2s3*k;
R22 = d2s3*k;
R33 = d2s3*k;
R12 = 0;
R13 = 0;
R23 = 0;
"""
            req = [('R11', "Reynolds stress R11"),
            ('R22', "Reynolds stress R22"),
            ('R33', "Reynolds stress R33"),
            ('R12', "Reynolds stress R12"),
            ('R13', "Reynolds stress R13"),
            ('R23', "Reynolds stress R23"),
            ('eps', "turbulent dissipation")]
            sym = [('x','face center coordinate'),
                   ('y','face center coordinate'),
                   ('z','face center coordinate'),
                   ('t','time'),
                   ('dt','time step'),
                   ('iter','number of time step')]
            dialog = QMeiEditorView(self,
                                    check_syntax = self.__case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.__boundary.setTurbFormula(self.__currentField, result)
                self.pushButtonTurb.setToolTip(result)
                self.pushButtonTurb.setStyleSheet("background-color: green")

        elif turb_model in ('tchen', 'q2-q12'):
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if not exp:
                exp = self.__boundary.getDefaultTurbFormula(turb_model)
            exa = """#example :
q2 = 5.e-05;
q12 = 0.0001;"""
            req = [('q2', "turbulent kinetic energy"),
            ('q12', "covariance")]
            sym = [('x','face center coordinate'),
                   ('y','face center coordinate'),
                   ('z','face center coordinate'),
                   ('t','time'),
                   ('dt','time step'),
                   ('iter','number of time step')]
            dialog = QMeiEditorView(self,
                                    check_syntax = self.__case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.__boundary.setTurbFormula(self.__currentField, result)
                self.pushButtonTurb.setToolTip(result)
                self.pushButtonTurb.setStyleSheet("background-color: green")

        elif turb_model in ('r2-q12'):
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if not exp:
                exp = self.__boundary.getDefaultTurbFormula(turb_model)
            exa = """#example :
R11 = 5e-05;
R22 = 5e-05;
R33 = 5e-05;
R12 = 5e-05;
R13 = 5e-05;
R23 = 5e-05;
q12 = 0.0001;"""
            req = [('R11', "Reynolds stress R11"),
            ('R22', "Reynolds stress R22"),
            ('R33', "Reynolds stress R33"),
            ('R12', "Reynolds stress R12"),
            ('R13', "Reynolds stress R13"),
            ('R23', "Reynolds stress R23"),
            ('q12', "covariance")]
            sym = [('x','face center coordinate'),
                   ('y','face center coordinate'),
                   ('z','face center coordinate'),
                   ('t','time'),
                   ('dt','time step'),
                   ('iter','number of time step')]
            dialog = QMeiEditorView(self,
                                    check_syntax = self.__case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.__boundary.setTurbFormula(self.__currentField, result)
                self.pushButtonTurb.setToolTip(result)
                self.pushButtonTurb.setStyleSheet("background-color: green")

        elif turb_model in ('r2-r12-tchen'):
            exp = self.__boundary.getTurbFormula(self.__currentField)
            if not exp:
                exp = self.__boundary.getDefaultTurbFormula(turb_model)
            exa = """#example :
R11 = 5e-05;
R22 = 5e-05;
R33 = 5e-05;
R12 = 5e-05;
R13 = 5e-05;
R23 = 5e-05;
R12-11 = 5e-05;
R12-22 = 5e-05;
R12-33 = 5e-05;
R12-12 = 5e-05;
R12-13 = 5e-05;
R12-23 = 5e-05;"""
            req = [('R11', "Reynolds stress R11"),
            ('R22', "Reynolds stress R22"),
            ('R33', "Reynolds stress R33"),
            ('R12', "Reynolds stress R12"),
            ('R13', "Reynolds stress R13"),
            ('R23', "Reynolds stress R23"),
            ('R12-11', "Reynolds stress R11"),
            ('R12-22', "Reynolds stress R22"),
            ('R12-33', "Reynolds stress R33"),
            ('R12-12', "Reynolds stress R12"),
            ('R12-13', "Reynolds stress R13"),
            ('R12-23', "Reynolds stress R23")]
            sym = [('x','face center coordinate'),
                   ('y','face center coordinate'),
                   ('z','face center coordinate'),
                   ('t','time'),
                   ('dt','time step'),
                   ('iter','number of time step')]
            dialog = QMeiEditorView(self,
                                    check_syntax = self.__case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotFormulaTurb -> %s" % str(result))
                self.__boundary.setTurbFormula(self.__currentField, result)
                self.pushButtonTurb.setToolTip(result)
                self.pushButtonTurb.setStyleSheet("background-color: green")


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
