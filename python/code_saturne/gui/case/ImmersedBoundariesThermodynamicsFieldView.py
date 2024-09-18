# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module define the 'ThermodynamicsField' page.
This module contains the following classes:
- ImmersedBoundariesThermodynamicsFieldView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
import logging

#-------------------------------------------------------------------------------
# EOS
#-------------------------------------------------------------------------------

from code_saturne.model.EosWrapper import eosWrapper

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
from code_saturne.gui.case.ImmersedBoundariesThermodynamicsField import Ui_ImmersedBoundariesThermodynamicsField
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#from code_saturne.model.EosWrapper import eosWrapper
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ThermodynamicsImmersedFieldView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# ImmersedBoundariesThermodynamicsFieldView class
#-------------------------------------------------------------------------------

class ImmersedBoundariesThermodynamicsFieldView(QWidget, Ui_ImmersedBoundariesThermodynamicsField):

    def __init__(self, parent = None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesThermodynamicsField.__init__(self)
        self.setupUi(self)

        self.parent    = None
        self.notebook  = None

        self.eos = eosWrapper()


    def setup(self, case, ibm, current_obj):
        self.case = case
        self.ibm = ibm
        self.current_obj = current_obj
        self.case.undoStopGlobal()

        self.is_FSI = self.ibm.getObjectFSI(self.current_obj)
        size_modelDensity = 2
        if self.is_FSI == "on":
            size_modelDensity = 1

        # Combo models
        self.modelDensity             = ComboModel(self.comboBoxDensity, size_modelDensity, 1)
        #self.modelStiffness           = ComboModel(self.comboBoxStiffness, 2, 1)
        self.modelStiffness           = ComboModel(self.comboBoxStiffness, 1, 1)
        #self.modelDamping             = ComboModel(self.comboBoxDamping, 2, 1)
        self.modelDamping             = ComboModel(self.comboBoxDamping, 1, 1)
        self.modelSpecificHeat        = ComboModel(self.comboBoxSpecificHeat, 2, 1)
        self.modelThermalConductivity = ComboModel(self.comboBoxThermalConductivity, 2, 1)

        self.modelDensity.addItem(self.tr('constant'), 'constant')
        if self.is_FSI == "off":
            self.modelDensity.addItem(self.tr('user law'), 'user_law')
        self.modelStiffness.addItem(self.tr('constant'), 'constant')
        #self.modelStiffness.addItem(self.tr('user law'), 'user law')
        self.modelDamping.addItem(self.tr('constant'), 'constant')
        #self.modelDamping.addItem(self.tr('user law'), 'user law')
        self.modelSpecificHeat.addItem(self.tr('constant'), 'constant')
        self.modelSpecificHeat.addItem(self.tr('user law'), 'user_law')
        self.modelThermalConductivity.addItem(self.tr('constant'), 'constant')
        self.modelThermalConductivity.addItem(self.tr('user law'), 'user_law')

        # Validators
        validatorRho = DoubleValidator(self.lineEditDensity, min = 0.0)
        validatorMass = DoubleValidator(self.lineEditMass, min = 0.0)
        validatorSti = DoubleValidator(self.lineEditStiffness, min = 0.0)
        validatorDam = DoubleValidator(self.lineEditDamping, min = 0.0)
        validatorCp = DoubleValidator(self.lineEditSpecificHeat, min = 0.0)
        validatorAl = DoubleValidator(self.lineEditThermalConductivity, min = 0.0)

        validatorRho.setExclusiveMin(True)
        validatorMass.setExclusiveMin(True)
        validatorCp.setExclusiveMin(True)
        validatorAl.setExclusiveMin(True)

        self.lineEditDensity.setValidator(validatorRho)
        self.lineEditMass.setValidator(validatorMass)
        self.lineEditStiffness.setValidator(validatorSti)
        self.lineEditDamping.setValidator(validatorDam)
        self.lineEditSpecificHeat.setValidator(validatorCp)
        self.lineEditThermalConductivity.setValidator(validatorAl)

        self.lineEditDensity.textChanged[str].connect(self.slotObjDensity)
        self.lineEditMass.textChanged[str].connect(self.slotObjMass)
        self.lineEditStiffness.textChanged[str].connect(self.slotObjStiffness)
        self.lineEditDamping.textChanged[str].connect(self.slotObjDamping)
        self.lineEditSpecificHeat.textChanged[str].connect(self.slotObjSpecificHeat)
        self.lineEditThermalConductivity.textChanged[str].connect(self.slotObjThermalConductivity)

        #radioButton
        self.radioButtonDensity.clicked.connect(self.__slotChoiceModeling)
        self.radioButtonMass.clicked.connect(self.__slotChoiceModeling)

        #Push Button for user law
        self.pushButtonDensity.clicked.connect(self.slotFormulaRho)
        self.pushButtonInertia.clicked.connect(self.slotFormulaInertia)
        self.pushButtonSpecificHeat.clicked.connect(self.slotFormulaCp)
        self.pushButtonThermalConductivity.clicked.connect(self.slotFormulaAl)
        self.pushButtonImposed.clicked.connect(self.slotPorousVelocityFormula)

        self.comboBoxes = {}
        self.comboBoxes['Rho']     = self.comboBoxDensity
        self.comboBoxes['Cp']      = self.comboBoxSpecificHeat
        self.comboBoxes['Al']      = self.comboBoxThermalConductivity
        self.comboBoxes['Sti']     = self.comboBoxStiffness
        self.comboBoxes['Dam']     = self.comboBoxDamping

        for k in self.comboBoxes.keys():
            self.comboBoxes[k].activated[str].connect(getattr(self,"slotState"+k))

        self.update()
        self.case.undoStartGlobal()


    def update(self):

        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            return

        is_CHT = self.ibm.getObjectCHT(self.current_obj)
        is_FSI = self.ibm.getObjectFSI(self.current_obj)

        if (is_CHT == "off" and is_FSI == "off"):
            self.groupBoxPhysicalProp.hide()

        if self.ibm.getObjectMoving(self.current_obj) == 'imposed':
            self.groupBoxImposed.show()

            exp = self.ibm.getObjectImposedMovingFormula(self.current_obj-1)
            if exp:
                self.pushButtonImposed.setToolTip(exp)
                self.pushButtonImposed.setStyleSheet("background-color: green")
            else:
                self.pushButtonImposed.setStyleSheet("background-color: red")

        else:
            self.groupBoxImposed.hide()

        # Need cp and lambda only if CHT activated
        if is_CHT == 'off':
            self.groupBoxSpecificHeat.hide()
            self.groupBoxThermalConductivity.hide()

            #Activate density only for FSI
            if (is_FSI == 'on'):
                self.groupBoxModeling.show()
                self.groupBoxStiffness.show()
                self.groupBoxDamping.show()
            else:
                self.groupBoxModeling.hide()
                self.groupBoxStiffness.hide()
                self.groupBoxDamping.hide()
        else:
            self.groupBoxModeling.show()
            if (is_FSI == "off"):
                self.ibm.setObjectModelingMode(self.current_obj, "density")
                self.radioButtonDensity.setEnabled(False)
                self.radioButtonMass.setEnabled(False)
                self.groupBoxStiffness.hide()
                self.groupBoxDamping.hide()
            else:
                self.groupBoxStiffness.show()
                self.groupBoxDamping.show()

            self.groupBoxSpecificHeat.show()
            self.groupBoxThermalConductivity.show()

        #Modeling
        modeling_mode = self.ibm.getObjectModelingMode(self.current_obj)
        if modeling_mode in ["density"]:
            self.groupBoxDensity.show()
            self.groupBoxMass.hide()
            self.radioButtonMass.setChecked(False)
            self.radioButtonDensity.setChecked(True)
        elif modeling_mode in ["mass"]:
            self.groupBoxMass.show()
            self.groupBoxDensity.hide()
            self.radioButtonDensity.setChecked(False)
            self.radioButtonMass.setChecked(True)
            self.lineEditMass.setText(str(self.ibm.getObjectMass(self.current_obj)))
            exp = self.ibm.getObjectInertiaFormula(self.current_obj-1)
            if exp:
                self.pushButtonInertia.setStyleSheet("background-color: green")
                self.pushButtonInertia.setToolTip(exp)
            else:
                self.pushButtonInertia.setStyleSheet("background-color: red")

        list = [('density', 'Density'),
                ('stiffness', 'Stiffness'),
                ('damping', 'Damping'),
                ('specific_heat', 'SpecificHeat'),
                ('thermal_conductivity', 'ThermalConductivity')]
        for tag, symbol in list :
            __line   = getattr(self, "lineEdit" + symbol)
            __button = getattr(self, "pushButton" + symbol)
            __model  = getattr(self, "model"      + symbol)

            mode = self.ibm.getObjectPropertyMode(self.current_obj, tag)
            if (is_FSI == 'on' and tag in ["density"]):
                #rho_sol becomes constant rho_fsi
                mode = "constant"
                self.ibm.setObjectPropertyMode(self.current_obj, mode, tag)
                __model.setItem(str_model=mode)

            if mode == "constant" :
                __button.setEnabled(False)
                __line.setEnabled(True)

                # Set correct values for each slot
                if tag == 'density':
                    __line.setText(str(self.ibm.getObjectDensity(self.current_obj)))
                elif tag == 'stiffness':
                    __line.setText(str(self.ibm.getObjectStiffness(self.current_obj)))
                elif tag == 'damping':
                    __line.setText(str(self.ibm.getObjectDamping(self.current_obj)))
                elif tag == 'specific_heat':
                    __line.setText(str(self.ibm.getObjectSpecificHeat(self.current_obj)))
                elif tag == 'thermal_conductivity':
                    __line.setText(str(self.ibm.getObjectThermalConductivity(self.current_obj)))

                __button.setStyleSheet("background-color: None")

            elif mode == "user_law" :
                __button.setEnabled(True)
                __line.setEnabled(False)

                exp = ''
                if tag == 'density':
                    exp = self.ibm.getObjectRhoFormula(self.current_obj-1)
                elif tag == 'specific_heat':
                    exp = self.ibm.getObjectCpFormula(self.current_obj-1)
                elif tag == 'thermal_conductivity':
                    exp = self.ibm.getObjectAlFormula(self.current_obj-1)

                if exp:
                    __button.setStyleSheet("background-color: green")
                    __button.setToolTip(exp)
                else:
                    __button.setStyleSheet("background-color: red")

            else :
                __button.setEnabled(False)
                __line.setEnabled(False)
                __line.setText("")
                __button.setStyleSheet("background-color: None")

    def __changeChoice(self, text, sym, tag):

        current_obj = self.current_obj

        __model  = getattr(self, "model"      + sym)
        __line   = getattr(self, "lineEdit"   + sym)
        __combo  = getattr(self, "comboBox"   + sym)
        __button = getattr(self, "pushButton" + sym)

        __button.setEnabled(False)
        __line.setEnabled(False)
        __line.setText("")
        __button.setStyleSheet("background-color: None")

        mode = __model.dicoV2M[str(text)]
        self.ibm.setObjectPropertyMode(current_obj, mode, tag)

        if mode == 'constant':
            __button.setEnabled(False)
            __line.setEnabled(True)
            __button.setStyleSheet("background-color: None")

            #__line.setText(str(self.ibm.defaultValues()['object_'+tag]))
            if tag == 'density':
                __line.setText(str(self.ibm.getObjectDensity(self.current_obj)))
            elif tag == 'specific_heat':
                __line.setText(str(self.ibm.getObjectSpecificHeat(self.current_obj)))
            elif tag == 'thermal_conductivity':
                __line.setText(str(self.ibm.getObjectThermalConductivity(self.current_obj)))
            elif tag == 'stiffness':
                __line.setText(str(self.ibm.getObjectStiffness(self.current_obj)))
            elif tag == 'damping':
                __line.setText(str(self.ibm.getObjectDamping(self.current_obj)))

        elif mode == 'user_law':
            __button.setEnabled(True)
            __line.setEnabled(False)
            __line.setText("")

            exp = ''
            if tag == 'density':
                exp = self.ibm.getObjectRhoFormula(self.current_obj-1)
            elif tag == 'specific_heat':
                exp = self.ibm.getObjectCpFormula(self.current_obj-1)
            elif tag == 'thermal_conductivity':
                exp = self.ibm.getObjectAlFormula(self.current_obj-1)

            if exp:
                __button.setStyleSheet("background-color: green")
                __button.setToolTip(exp)
            else:
                __button.setStyleSheet("background-color: red")


    @pyqtSlot(str)
    def slotObjDensity(self, text):
        if self.lineEditDensity.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectDensity(self.current_obj, val)


    @pyqtSlot(str)
    def slotObjMass(self, text):
        if self.lineEditMass.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectMass(self.current_obj, val)


    @pyqtSlot(str)
    def slotObjStiffness(self, text):
        if self.lineEditStiffness.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectStiffness(self.current_obj, val)


    @pyqtSlot(str)
    def slotObjDamping(self, text):
        if self.lineEditDamping.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectDamping(self.current_obj, val)


    @pyqtSlot(str)
    def slotObjSpecificHeat(self, text):
        if self.lineEditSpecificHeat.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectSpecificHeat(self.current_obj, val)

    @pyqtSlot(str)
    def slotObjThermalConductivity(self, text):
        if self.lineEditThermalConductivity.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectThermalConductivity(self.current_obj, val)


    @pyqtSlot(str)
    def slotStateRho(self, text):
        """
        Method to call 'getState' with correct arguements for 'Rho'
        """
        self.__changeChoice(str(text), 'Density', 'density')


    @pyqtSlot(str)
    def slotStateSti(self, text):
        """
        Method to call 'getState' with correct arguements for 'Sti'
        """
        self.__changeChoice(str(text), 'Stiffness', 'stiffness')

    @pyqtSlot(str)
    def slotStateDam(self, text):
        """
        Method to call 'getState' with correct arguements for 'Dam'
        """
        self.__changeChoice(str(text), 'Damping', 'damping')


    @pyqtSlot(str)
    def slotStateCp(self, text):
        """
        Method to call 'getState' with correct arguements for 'Cp'
        """
        self.__changeChoice(str(text), 'SpecificHeat', 'specific_heat')


    @pyqtSlot(str)
    def slotStateAl(self, text):
        """
        Method to call 'getState' with correct arguements for 'Al'
        """
        self.__changeChoice(str(text), 'ThermalConductivity', 'thermal_conductivity')


    @pyqtSlot()
    def slotPorousVelocityFormula(self):
        """
        Explicit formula for the velocity of the porous (solid) media
        """

        objId = self.current_obj

        exp, req, sym = self.ibm.getImposedMovingFormulaComponents(objId-1)

        exa = (
            "porous_velocity[0] = 0.;\n"
            "porous_velocity[1] = 0.;\n"
            "porous_velocity[2] = 0.;"
        )

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = name,
                                variable_name = 'porous_velocity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousVelocity -> %s" % str(result))
            self.pushButtonImposed.setStyleSheet("background-color: green")
            self.pushButtonImposed.setToolTip(exp)
            self.ibm.setObjectImposedMovingFormula(objId-1, result)


    @pyqtSlot()
    def slotFormulaRho(self):
        """
        User formula for density of the porous object
        """
        objId = self.current_obj

        exp, req, sym  = self.ibm.getFormulaRho(objId-1)
        exa = """# Density of air
rho = 1.293 * (273.15 / temperature);"""

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = name,
                                variable_name = 'porous_density',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousRho -> %s" % str(result))
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(exp)
            self.ibm.setObjectRhoFormula(objId-1, result)


    @pyqtSlot()
    def slotFormulaCp(self):
        """
        User formula for specific heat of the porous object
        """
        objId = self.current_obj

        exp, req, sym  = self.ibm.getFormulaCp(objId-1)
        exa = """cp = 4000.;"""

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = name,
                                variable_name = 'porous_specific_heat',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousCp -> %s" % str(result))
            self.pushButtonSpecificHeat.setStyleSheet("background-color: green")
            self.pushButtonSpecificHeat.setToolTip(exp)
            self.ibm.setObjectCpFormula(objId-1, result)


    @pyqtSlot()
    def slotFormulaAl(self):
        """
        User formula for specific heat of the porous object
        """
        objId = self.current_obj

        exp, req, sym  = self.ibm.getFormulaAl(objId-1)
        exa = """lambda = 1.e-5;"""

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = name,
                                variable_name = 'porous_thermal_conductivity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousAl -> %s" % str(result))
            self.pushButtonThermalConductivity.setStyleSheet("background-color: green")
            self.pushButtonThermalConductivity.setToolTip(exp)
            self.ibm.setObjectAlFormula(objId-1, result)


    @pyqtSlot()
    def slotFormulaInertia(self):
        """
        User formula for inertia of the porous object
        """
        objId = self.current_obj

        exp, req, sym  = self.ibm.getFormulaInertia(objId-1)
        exa = (
            "# Homogeneous sphere centered at the origin"
            " with mass m=5kg and radius R=0.1m\n\n"
            "m=5.;\n"
            "R=0.1;\n"
            "i11 = 0.4*m*R*R;\n"
            "i22 = 0.4*m*R*R;\n"
            "i33 = 0.4*m*R*R;\n"
            "i12 = 0.;\n"
            "i13 = 0.;\n"
            "i23 = 0.;\n"
        )

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_fsi',
                                zone_name     = name,
                                variable_name = 'porous_inertia',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousInertia -> %s" % str(result))
            self.pushButtonThermalConductivity.setStyleSheet("background-color: green")
            self.pushButtonThermalConductivity.setToolTip(exp)
            self.ibm.setObjectInertiaFormula(objId-1, result)



    @pyqtSlot()
    def __slotChoiceModeling(self):

        if self.radioButtonDensity.isChecked():
            self.ibm.setObjectModelingMode(self.current_obj,'density')
            self.groupBoxDensity.show()
            self.groupBoxMass.hide()
        elif self.radioButtonMass.isChecked():
            self.ibm.setObjectModelingMode(self.current_obj,'mass')
            self.groupBoxMass.show()
            self.groupBoxDensity.hide()

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
