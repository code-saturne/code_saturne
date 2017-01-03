# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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
- BoundaryConditionsVelocityMappedInletView
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

from code_saturne.Pages.BoundaryConditionsVelocityMappedInletForm import Ui_BoundaryConditionsVelocityMappedInletForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsVelocityMappedInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsVelocityMappedInletView(QWidget, Ui_BoundaryConditionsVelocityMappedInletForm):
    """
    Boundary condition for velocity in inlet, without particular physics.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsVelocityMappedInletForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        # Connections
        self.comboBoxVelocity.activated[str].connect(self.__slotChoiceVelocity)
        self.lineEditVelocity.textChanged[str].connect(self.__slotVelocityValue)

        self.lineEditTranslationX.textChanged[str].connect(self.__slotDirX)
        self.lineEditTranslationY.textChanged[str].connect(self.__slotDirY)
        self.lineEditTranslationZ.textChanged[str].connect(self.__slotDirZ)

        # Combo models
        self.modelVelocity = ComboModel(self.comboBoxVelocity, 6, 1)
        self.modelVelocity.addItem(self.tr("norm"), 'norm')
        self.modelVelocity.addItem(self.tr("mass flow rate"), 'flow1')
        self.modelVelocity.addItem(self.tr("volumic flow rate"), 'flow2')
        self.modelVelocity.addItem(self.tr("norm (user law)"), 'norm_formula')
        self.modelVelocity.addItem(self.tr("mass flow rate (user law)"), 'flow1_formula')
        self.modelVelocity.addItem(self.tr("volumic flow rate (user law)"), 'flow2_formula')

        # Validators
        validatorVelocity = DoubleValidator(self.lineEditVelocity)
        validatorX = DoubleValidator(self.lineEditTranslationX)
        validatorY = DoubleValidator(self.lineEditTranslationY)
        validatorZ = DoubleValidator(self.lineEditTranslationZ)

        # Apply validators
        self.lineEditVelocity.setValidator(validatorVelocity)
        self.lineEditTranslationX.setValidator(validatorX)
        self.lineEditTranslationY.setValidator(validatorY)
        self.lineEditTranslationZ.setValidator(validatorZ)

        self.pushButtonVelocityFormula.clicked.connect(self.__slotVelocityFormula)

        self.__case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        # Initialize velocity
        choice = self.__boundary.getVelocityChoice()
        self.modelVelocity.setItem(str_model=choice)
        self.__updateLabel()

        if choice[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            self.lineEditVelocity.setEnabled(False)
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(str(v))

        # Initialize direction
        v = self.__boundary.getTranslation('direction_x')
        self.lineEditTranslationX.setText(str(v))
        v = self.__boundary.getTranslation('direction_y')
        self.lineEditTranslationY.setText(str(v))
        v = self.__boundary.getTranslation('direction_z')
        self.lineEditTranslationZ.setText(str(v))

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(str)
    def __slotChoiceVelocity(self, text):
        """
        Private slot.

        Input the velocity boundary type choice (norm, ).

        @type text: C{QString}
        @param text: velocity boundary type choice.
        """
        c = self.modelVelocity.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setVelocityChoice(c)

        if c[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            exp = self.__boundary.getVelocity()
            if exp:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: green")
                self.pushButtonVelocityFormula.setToolTip(exp)
            else:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: red")
            self.lineEditVelocity.setEnabled(False)
            self.lineEditVelocity.setText("")
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.pushButtonVelocityFormula.setStyleSheet("background-color: None")
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(str(v))

        self.__updateLabel()


    def __updateLabel(self):
        """
        Update the unit for the velocity specification.
        """
        c = self.__boundary.getVelocityChoice()
        if c in ('norm', 'norm_formula'):
            self.labelUnitVelocity.setText(str('m/s'))
        elif c in ('flow1', 'flow1_formula'):
            self.labelUnitVelocity.setText(str('kg/s'))
        elif c in ('flow2', 'flow2_formula'):
            self.labelUnitVelocity.setText(str('m<sup>3</sup>/s'))


    @pyqtSlot(str)
    def __slotVelocityValue(self, text):
        """
        Private slot.

        New value associated to the velocity boundary type.

        @type text: C{QString}
        @param text: value
        """
        if self.lineEditVelocity.validator().state == QValidator.Acceptable:
            v = from_qvariant(text, float)
            self.__boundary.setVelocity(v)


    @pyqtSlot()
    def __slotVelocityFormula(self):
        """
        """
        exp = self.__boundary.getVelocity()
        c = self.__boundary.getVelocityChoice()
        if c == 'norm_formula':
            exa = "u_norm = 1.0;"
            req = [('u_norm', 'Norm of the velocity')]
        elif c == 'flow1_formula':
            exa = "q_m = 1.0;"
            req = [('q_m', 'mass flow rate')]
        elif c == 'flow2_formula':
            exa = "q_v = 1.0;"
            req = [('q_v', 'volumic flow rate')]

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        dialog = QMeiEditorView(self,
                                check_syntax = self.__case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.__boundary.setVelocity(str(result))
            self.pushButtonVelocityFormula.setStyleSheet("background-color: green")
            self.pushButtonVelocityFormula.setToolTip(result)


    @pyqtSlot(str)
    def __slotDirX(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.lineEditTranslationX.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setTranslation('direction_x', value)


    @pyqtSlot(str)
    def __slotDirY(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.lineEditTranslationY.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setTranslation('direction_y', value)


    @pyqtSlot(str)
    def __slotDirZ(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.lineEditTranslationZ.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setTranslation('direction_z', value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
