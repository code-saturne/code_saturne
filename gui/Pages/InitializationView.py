# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
- InitializationView
"""

#-------------------------------------------------------------------------------
# Standard modules
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

from Pages.InitializationForm import Ui_InitializationForm

from Base.Toolbox import GuiParam
from Base.QtPage import IntValidator, DoubleValidator, ComboModel
from Pages.TurbulenceModel import TurbulenceModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.GasCombustionModel import GasCombustionModel
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.LocalizationModel import VolumicLocalizationModel, LocalizationModel
from Pages.InitializationModel import InitializationModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("InitializationView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class InitializationView(QWidget, Ui_InitializationForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_InitializationForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent

        self.init    = InitializationModel(self.case)
        self.turb    = TurbulenceModel(self.case)
        self.therm   = ThermalScalarModel(self.case)
        self.th_sca  = DefineUserScalarsModel(self.case)
        self.volzone = LocalizationModel('VolumicZone', self.case)


        # 0/ Read label names from XML file

        # Velocity
        [U, V, W] = self.init.getVelocityLabel()
        self.label_U.setText(QString(str(U)))
        self.label_V.setText(QString(str(V)))
        self.label_W.setText(QString(str(W)))

        # Thermal scalar
        namesca, unit = self.getThermalLabelAndUnit()
        self.labelThermal.setText(QString(str(namesca)))
        self.labelUnitThermal.setText(QString(str(unit)))
        self.th_sca_label = namesca

        # k-eps
        label_k   = self.init.getTurbulenceVariableLabel('turb_k')
        label_eps = self.init.getTurbulenceVariableLabel('turb_eps')
        self.labelKeps_k.setText(QString(str(label_k)))
        self.labelKeps_eps.setText(QString(str(label_eps)))

        # Rij
        for (coeff, attr) in [('R11', 'component_R11'),
                              ('R22', 'component_R22'),
                              ('R33', 'component_R33'),
                              ('R12', 'component_R12'),
                              ('R13', 'component_R13'),
                              ('R23', 'component_R23')]:
            tt = self.init.getTurbulenceVariableLabel(attr)
            line = getattr(self, "label" + coeff)
            line.setText(QString(str(tt)))

        tt = self.init.getTurbulenceVariableLabel('turb_eps')
        self.labelRijEps.setText(QString(str(tt)))

        # v2f
        labelv2f_k   = self.init.getTurbulenceVariableLabel('turb_k')
        labelv2f_eps = self.init.getTurbulenceVariableLabel('turb_eps')
        labelv2f_phi = self.init.getTurbulenceVariableLabel('turb_phi')
        labelv2f_fb  = self.init.getTurbulenceVariableLabel('turb_fb')
        self.labelv2f_k.setText(QString(str(labelv2f_k)))
        self.labelv2f_eps.setText(QString(str(labelv2f_eps)))
        self.labelv2f_phi.setText(QString(str(labelv2f_phi)))
        self.labelv2f_fb.setText(QString(str(labelv2f_fb)))

        # K-omega
        labelKomega_k     = self.init.getTurbulenceVariableLabel('turb_k')
        labelKomega_omega = self.init.getTurbulenceVariableLabel('turb_omega')
        self.labelKomega_k.setText(QString(str(labelKomega_k)))
        self.labelKomega_omega.setText(QString(str(labelKomega_omega)))

        # Spalart-Allmaras
        labelSA_nusa = self.init.getTurbulenceVariableLabel('turb_nusa')
        self.labelSA_nusa.setText(QString(str(labelSA_nusa)))

        # 1/ Combo box models

        self.modelZone = ComboModel(self.comboBoxZone, 1, 1)

        self.zone = ""
        zones = self.volzone.getZones()
        for zone in zones:
            if zone.getNature()['initialization'] == "on":
                label = zone.getLabel()
                name = str(zone.getCodeNumber())
                self.modelZone.addItem(self.tr(label), name)
                if label == "all_cells":
                    self.zone = name
                if not self.zone:
                    self.zone = name
        #    else:
        #        raise ValueError, "there is no initialization zone defined."

        self.modelZone.setItem(str_model = self.zone)

        self.modelTurbulence = ComboModel(self.comboBoxTurbulence, 3, 1)
        self.modelTurbulence.addItem(self.tr("Initialization by values for the selected zone"), 'values')
        self.modelTurbulence.addItem(self.tr("Initialization by reference velocity for all zones"), 'reference_velocity')
        self.modelTurbulence.addItem(self.tr("Initialization by reference velocity and reference length for all zones"), 'reference_velocity_length')

        # 2/ Connections

        self.connect(self.comboBoxZone,         SIGNAL("activated(const QString&)"),   self.slotZone)
        self.connect(self.lineEditU,            SIGNAL("textChanged(const QString&)"), self.slotU)
        self.connect(self.lineEditV,            SIGNAL("textChanged(const QString&)"), self.slotV)
        self.connect(self.lineEditW,            SIGNAL("textChanged(const QString&)"), self.slotW)
        self.connect(self.lineEditThermal,      SIGNAL("textChanged(const QString&)"), self.slotThermalValue)
        self.connect(self.comboBoxTurbulence,   SIGNAL("activated(const QString&)"),   self.slotChoice)
        self.connect(self.lineEditKeps_k,       SIGNAL("textChanged(const QString&)"), self.slotKeps_k)
        self.connect(self.lineEditKeps_eps,     SIGNAL("textChanged(const QString&)"), self.slotKeps_eps)
        self.connect(self.lineEditR11,          SIGNAL("textChanged(const QString&)"), self.slotR11)
        self.connect(self.lineEditR22,          SIGNAL("textChanged(const QString&)"), self.slotR22)
        self.connect(self.lineEditR33,          SIGNAL("textChanged(const QString&)"), self.slotR33)
        self.connect(self.lineEditR12,          SIGNAL("textChanged(const QString&)"), self.slotR12)
        self.connect(self.lineEditR13,          SIGNAL("textChanged(const QString&)"), self.slotR13)
        self.connect(self.lineEditR23,          SIGNAL("textChanged(const QString&)"), self.slotR23)
        self.connect(self.lineEditRijEps,       SIGNAL("textChanged(const QString&)"), self.slotRijEps)
        self.connect(self.lineEditv2f_k,        SIGNAL("textChanged(const QString&)"), self.slotv2f_k)
        self.connect(self.lineEditv2f_eps,      SIGNAL("textChanged(const QString&)"), self.slotv2f_eps)
        self.connect(self.lineEditv2f_phi,      SIGNAL("textChanged(const QString&)"), self.slotv2f_phi)
        self.connect(self.lineEditv2f_fb,       SIGNAL("textChanged(const QString&)"), self.slotv2f_fb)
        self.connect(self.lineEditKomega_k,     SIGNAL("textChanged(const QString&)"), self.slotKomega_k)
        self.connect(self.lineEditKomega_omega, SIGNAL("textChanged(const QString&)"), self.slotKomega_omega)
        self.connect(self.lineEditSA_nusa,      SIGNAL("textChanged(const QString&)"), self.slotSA_nusa)
        self.connect(self.lineEditRefVelocity,  SIGNAL("textChanged(const QString&)"), self.slotReferenceVelocity)
        self.connect(self.lineEditRefLength,    SIGNAL("textChanged(const QString&)"), self.slotReferenceVelocityAndLength)

        # 3/ Validator definitions

        validatorU            = DoubleValidator(self.lineEditU)
        validatorV            = DoubleValidator(self.lineEditV)
        validatorW            = DoubleValidator(self.lineEditW)
        validatorThermal      = DoubleValidator(self.lineEditThermal)
        validatorKeps_k       = DoubleValidator(self.lineEditKeps_k, min=0.)
        validatorKeps_eps     = DoubleValidator(self.lineEditKeps_eps, min=0.)
        validatorR11          = DoubleValidator(self.lineEditR11, min=0.)
        validatorR22          = DoubleValidator(self.lineEditR22, min=0.)
        validatorR33          = DoubleValidator(self.lineEditR33, min=0.)
        validatorR12          = DoubleValidator(self.lineEditR12, min=0.)
        validatorR13          = DoubleValidator(self.lineEditR13, min=0.)
        validatorR23          = DoubleValidator(self.lineEditR23, min=0.)
        validatorRijEps       = DoubleValidator(self.lineEditRijEps, min=0.)
        validatorv2f_k        = DoubleValidator(self.lineEditv2f_k, min=0.)
        validatorv2f_eps      = DoubleValidator(self.lineEditv2f_eps, min=0.)
        validatorv2f_phi      = DoubleValidator(self.lineEditv2f_phi, min=0.)
        validatorv2f_fb       = DoubleValidator(self.lineEditv2f_fb, min=0.)
        validatorKomega_k     = DoubleValidator(self.lineEditKomega_k, min=0.)
        validatorKomega_omega = DoubleValidator(self.lineEditKomega_omega, min=0.)
        validatorSA_nusa      = DoubleValidator(self.lineEditSA_nusa, min=0.)
        validatorRefVelocity  = DoubleValidator(self.lineEditRefVelocity, min=0.)
        validatorRefLength    = DoubleValidator(self.lineEditRefLength, min=0.)
        validatorRefLength.setExclusiveMin()

        self.lineEditU.setValidator(validatorU)
        self.lineEditV.setValidator(validatorV)
        self.lineEditW.setValidator(validatorW)
        self.lineEditThermal.setValidator(validatorThermal)
        self.lineEditKeps_k.setValidator(validatorKeps_k)
        self.lineEditKeps_eps.setValidator(validatorKeps_eps)
        self.lineEditR11.setValidator(validatorR11)
        self.lineEditR22.setValidator(validatorR22)
        self.lineEditR33.setValidator(validatorR33)
        self.lineEditR12.setValidator(validatorR12)
        self.lineEditR13.setValidator(validatorR13)
        self.lineEditR23.setValidator(validatorR23)
        self.lineEditRijEps.setValidator(validatorRijEps)
        self.lineEditv2f_k.setValidator(validatorv2f_k)
        self.lineEditv2f_eps.setValidator(validatorv2f_eps)
        self.lineEditv2f_phi.setValidator(validatorv2f_phi)
        self.lineEditv2f_fb.setValidator(validatorv2f_fb)
        self.lineEditKomega_k.setValidator(validatorKomega_k)
        self.lineEditKomega_omega.setValidator(validatorKomega_omega)
        self.lineEditSA_nusa.setValidator(validatorSA_nusa)
        self.lineEditRefVelocity.setValidator(validatorRefVelocity)
        self.lineEditRefLength.setValidator(validatorRefLength)

        # Initialize widget

        self.initializeVariables(self.zone)


    @pyqtSignature("const QString&")
    def slotZone(self, text):
        """
        INPUT label for choice of zone
        """
        self.zone = self.modelZone.dicoV2M[str(text)]
        self.initializeVariables(self.zone)


    @pyqtSignature("const QString&")
    def slotU(self, var):
        """
        """
        value, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setInitialVelocity(self.zone, 'velocity_U', value)


    @pyqtSignature("const QString&")
    def slotV(self, var):
        """
        """
        value, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setInitialVelocity(self.zone, 'velocity_V', value)


    @pyqtSignature("const QString&")
    def slotW(self, var):
        """
        """
        value, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setInitialVelocity(self.zone, 'velocity_W', value)


    @pyqtSignature("const QString&")
    def slotChoice(self, text):
        """
        INPUT choice of method of initialization
        """
        # Hide everything
        self.frameKeps.hide()
        self.frameRij.hide()
        self.framev2f.hide()
        self.frameKomega.hide()
        self.frameReference.hide()
        self.frameSA.hide()

        choice = self.modelTurbulence.dicoV2M[str(text)]
        log.debug("slotChoice choice =  %s "%str(choice))
        self.init.setInitialTurbulenceChoice(self.zone, choice)
        turb_model = self.turb.getTurbulenceModel()

        if choice == 'values':
            if turb_model in ('k-epsilon', 'k-epsilon-PL'):

                self.frameKeps.show()

                k   = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')
                self.lineEditKeps_k.setText(QString(str(k)))
                self.lineEditKeps_eps.setText(QString(str(eps)))

            elif turb_model in ('Rij-epsilon', 'Rij-SSG'):

                self.frameRij.show()

                R11 = self.init.getTurbulenceInitialValue(self.zone, 'component_R11')
                R22 = self.init.getTurbulenceInitialValue(self.zone, 'component_R22')
                R33 = self.init.getTurbulenceInitialValue(self.zone, 'component_R33')
                R12 = self.init.getTurbulenceInitialValue(self.zone, 'component_R12')
                R13 = self.init.getTurbulenceInitialValue(self.zone, 'component_R13')
                R23 = self.init.getTurbulenceInitialValue(self.zone, 'component_R23')
                eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')

                self.lineEditR11.setText(QString(R11))
                self.lineEditR22.setText(QString(R22))
                self.lineEditR33.setText(QString(R33))
                self.lineEditR12.setText(QString(R12))
                self.lineEditR13.setText(QString(R13))
                self.lineEditR23.setText(QString(R23))
                self.lineEditRijEps.setText(QString(eps))

            elif turb_model == 'v2f-phi':

                self.framev2f.show()

                k   = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')
                phi = self.init.getTurbulenceInitialValue(self.zone, 'turb_phi')
                fb  = self.init.getTurbulenceInitialValue(self.zone, 'turb_fb')

                self.lineEditv2f_k.setText(QString(k))
                self.lineEditv2f_eps.setText(QString(eps))
                self.lineEditv2f_phi.setText(QString(phi))
                self.lineEditv2f_fb.setText(QString(fb))

            elif turb_model == 'k-omega-SST':

                self.frameKomega.show()

                k     = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                omega = self.init.getTurbulenceInitialValue(self.zone, 'turb_omega')

                self.lineEditKomega_k.setText(QString(k))
                self.lineEditKomega_omega.setText(QString(omega))

            elif turb_model == 'Spalart-Allmaras':

                self.frameSA.show()

                nusa = self.init.getTurbulenceInitialValue(self.zone, 'turb_nusa')

                self.lineEditSA_nusa.setText(QString(nusa))

        elif choice == 'reference_velocity':

            self.frameReference.show()
            self.labelRefLength.hide()
            self.lineEditRefLength.hide()
            self.labelUnitRefLength.hide()

            v = self.init.getReferenceVelocity()
            self.lineEditRefVelocity.setText(QString(str(v)))

        elif choice == 'reference_velocity_length':

            self.frameReference.show()
            self.labelRefLength.show()
            self.lineEditRefLength.show()
            self.labelUnitRefLength.show()

            v, l = self.init.getReferenceVelocityAndLength()
            self.lineEditRefVelocity.setText(QString(str(v)))
            self.lineEditRefLength.setText(QString(str(l)))


    @pyqtSignature("const QString&")
    def slotKeps_k(self, var):
        """
        INPUT values for K-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_k', val)


    @pyqtSignature("const QString&")
    def slotKeps_eps(self, var):
        """
        INPUT values for K-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_eps', val)


    @pyqtSignature("const QString&")
    def slotR11(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R11', val)


    @pyqtSignature("const QString&")
    def slotR12(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R12', val)


    @pyqtSignature("const QString&")
    def slotR13(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R13', val)


    @pyqtSignature("const QString&")
    def slotR22(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R22', val)


    @pyqtSignature("const QString&")
    def slotR23(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R23', val)


    @pyqtSignature("const QString&")
    def slotR33(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'component_R33', val)


    @pyqtSignature("const QString&")
    def slotRijEps(self, var):
        """
        INPUT values for Rij-epsilon's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_eps', val)


    @pyqtSignature("const QString&")
    def slotv2f_k(self, var):
        """
        INPUT values for v2f's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_k', val)


    @pyqtSignature("const QString&")
    def slotv2f_eps(self, var):
        """
        INPUT values for v2f's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_eps', val)


    @pyqtSignature("const QString&")
    def slotv2f_phi(self, var):
        """
        INPUT values for v2f's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_phi', val)


    @pyqtSignature("const QString&")
    def slotv2f_fb(self, var):
        """
        INPUT values for v2f's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_fb', val)


    @pyqtSignature("const QString&")
    def slotKomega_k(self, var):
        """
        INPUT values for k-omega's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_k', val)


    @pyqtSignature("const QString&")
    def slotKomega_omega(self, var):
        """
        INPUT values for k-omega's initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_omega', val)


    @pyqtSignature("const QString&")
    def slotSA_nusa(self, var):
        """
        INPUT values for Spalart-Allmaras' initialization by values
        """
        val, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setTurbulenceInitialValue(self.zone, 'turb_nusa', val)


    @pyqtSignature("const QString&")
    def slotReferenceVelocity(self, var):
        """
        INPUT values for initialization by reference velocity or
        length and velocity of reference.
        """
        UREF, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setReferenceVelocity(UREF)


    @pyqtSignature("const QString&")
    def slotReferenceVelocityAndLength(self, var):
        """
        INPUT values for initialization by reference length.
        """
        long, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.init.setReferenceLength(long)

# FIXME: slotReferenceVelocityAndLength (long != -1e+12)
#            if long < 0.0:
#                if long != -1e+12:
#                    msg = self.tr("THE ONLY AUTHORIZED NEGATIVE VALUE IS: -1e+12")
#                    self.stbar.showMessage(msg, 2000)
#                    self.lineEditRefLength.setText(QString(""))
#                else:
#                    self.init.setReferenceLength(long)
#            else:
#                self.init.setReferenceLength(long)


    def getThermalLabelAndUnit(self):
        """
        Define the type of model is used.
        """
        model = self.therm.getThermalScalarModel()

        if model != 'off':
            th_sca_label = self.therm.getThermalScalarLabel()
            if model == "temperature_celsius":
                unit = "<sup>o</sup>C"
            elif model == "temperature_kelvin":
                unit = "Kelvin"
            elif model == "enthalpy":
                unit = "J/kg"
        else:
            th_sca_label = ''
            unit = None

        # Value to use in slotThermalValue
        self.th_sca_label = th_sca_label

        return th_sca_label, unit


    @pyqtSignature("const QString&")
    def slotThermalValue(self, var):
        """
        Input the initial value of thermal scalar
        """
        v, ok = var.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            min = self.th_sca.getScalarMinValue(self.th_sca_label)
            max = self.th_sca.getScalarMaxValue(self.th_sca_label)
            from Pages.DefineUserScalarsView import StandardItemModelScalars
            if StandardItemModelScalars(self.parent, self.th_sca,
                self.zone).checkInitMinMax(self.th_sca_label, v, min, max):
                self.th_sca.setScalarInitialValue(self.zone, self.th_sca_label, v)


    def initializeVariables(self, zone):
        """
        Initialize variables when a new volumic zone is choosen
        """
        # Velocity initialization

        [U, V, W] = self.init.getInitialVelocity(zone)
        self.lineEditU.setText(QString(str(U)))
        self.lineEditV.setText(QString(str(V)))
        self.lineEditW.setText(QString(str(W)))

        # Initialisation of Turbulence

        turb_model = self.turb.getTurbulenceModel()

        if turb_model not in ('k-epsilon',
                              'k-epsilon-PL',
                              'Rij-epsilon',
                              'Rij-SSG',
                              'v2f-phi',
                              'k-omega-SST',
                              'Spalart-Allmaras'):
            self.groupBoxTurbulence.hide()
        else:
            self.groupBoxTurbulence.show()

            # Select the initialization's choice

            choice = self.init.getInitialTurbulenceChoice(zone)
            self.modelTurbulence.setItem(str_model = choice)
            txt = self.comboBoxTurbulence.currentText()
            self.slotChoice(txt)

            if choice == 'reference_velocity':
                v = self.init.getReferenceVelocity()
                self.lineEditRefVelocity.setText(QString(str(v)))

            elif choice == 'reference_velocity_length':
                v,l = self.init.getReferenceVelocityAndLength()
                self.lineEditRefVelocity.setText(QString(str(v)))
                self.lineEditRefLength.setText(QString(str(l)))

            elif choice == 'values':

                if turb_model in ('k-epsilon', 'k-epsilon-PL'):

                    k   = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                    eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')
                    self.lineEditKeps_k.setText(QString(str(k)))
                    self.lineEditKeps_eps.setText(QString(str(eps)))

                elif turb_model in ('Rij-epsilon', 'Rij-epsilon-SSG'):

                    R11 = self.init.getTurbulenceInitialValue(self.zone, 'component_R11')
                    R22 = self.init.getTurbulenceInitialValue(self.zone, 'component_R22')
                    R33 = self.init.getTurbulenceInitialValue(self.zone, 'component_R33')
                    R12 = self.init.getTurbulenceInitialValue(self.zone, 'component_R12')
                    R13 = self.init.getTurbulenceInitialValue(self.zone, 'component_R13')
                    R23 = self.init.getTurbulenceInitialValue(self.zone, 'component_R23')
                    eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')

                    self.lineEditR11.setText(QString(str(R11)))
                    self.lineEditR22.setText(QString(str(R22)))
                    self.lineEditR33.setText(QString(str(R33)))
                    self.lineEditR12.setText(QString(str(R12)))
                    self.lineEditR13.setText(QString(str(R13)))
                    self.lineEditR23.setText(QString(str(R23)))
                    self.lineEditRijEps.setText(QString(str(eps)))

                elif turb_model == 'v2f-phi':

                    k   = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                    eps = self.init.getTurbulenceInitialValue(self.zone, 'turb_eps')
                    phi = self.init.getTurbulenceInitialValue(self.zone, 'turb_phi')
                    fb  = self.init.getTurbulenceInitialValue(self.zone, 'turb_fb')

                    self.lineEditv2f_k.setText(QString(str(k)))
                    self.lineEditv2f_eps.setText(QString(str(eps)))
                    self.lineEditv2f_phi.setText(QString(str(phi)))
                    self.lineEditv2f_fb.setText(QString(str(fb)))

                elif turb_model == 'k-omega-SST':

                    k     = self.init.getTurbulenceInitialValue(self.zone, 'turb_k')
                    omega = self.init.getTurbulenceInitialValue(self.zone, 'turb_omega')

                    self.lineEditKomega_k.setText(QString(str(k)))
                    self.lineEditKomega_omega.setText(QString(str(omega)))

                elif turb_model == 'Spalart-Allmaras':

                    nusa = self.init.getTurbulenceInitialValue(self.zone, 'turb_nusa')

                    self.lineEditSA_nusa.setText(QString(str(nusa)))


        # Initialisation of Model Variables if thermal model is selectionned

        name, unit = self.getThermalLabelAndUnit()

        if name == "":
            self.groupBoxThermal.hide()
        else:
            self.groupBoxThermal.show()

            ini = self.th_sca.getScalarInitialValue(zone, name)
            self.lineEditThermal.setText(QString(str(ini)))
            self.labelUnitThermal.setText(QString(str(unit)))


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
