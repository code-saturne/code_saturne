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
- BoundaryConditionsVelocityInletView
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

from code_saturne.Pages.BoundaryConditionsVelocityInletForm import Ui_BoundaryConditionsVelocityInletForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.Pages.LocalizationModelNeptune import LocalizationModel, Zone
from code_saturne.Pages.BoundaryNeptune import Boundary

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsVelocityInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsVelocityInletView(QWidget, Ui_BoundaryConditionsVelocityInletForm) :
    """
    Boundary condition for velocity in inlet, without particular physics.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsVelocityInletForm.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxVelocity.activated[str].connect(self.__slotChoiceVelocity)
        self.lineEditVelocity.textChanged[str].connect(self.__slotVelocityValue)

        self.comboBoxDirection.activated[str].connect(self.__slotChoiceDirection)

        # Combo models
        self.modelVelocity = ComboModel(self.comboBoxVelocity, 4, 1)
        self.modelVelocity.addItem(self.tr("Norm"), 'norm')
        self.modelVelocity.addItem(self.tr("Mass flow rate"), 'flow1')
        self.modelVelocity.addItem(self.tr("Norm (user law)"), 'norm_formula')
        self.modelVelocity.addItem(self.tr("Mass flow rate (user law)"), 'flow1_formula')

        self.modelDirection = ComboModel(self.comboBoxDirection, 2, 1)
        self.modelDirection.addItem(self.tr("Normal direction to the inlet"), 'normal')
        self.modelDirection.addItem(self.tr("User profile"), 'formula')

        # Validators
        validatorVelocity = DoubleValidator(self.lineEditVelocity)

        # Apply validators
        self.lineEditVelocity.setValidator(validatorVelocity)

        self.pushButtonVelocityFormula.clicked.connect(self.__slotVelocityFormula)
        self.pushButtonDirectionFormula.clicked.connect(self.__slotDirectionFormula)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None
        self.__currentField = fieldId
        self.groupBoxCompressible.hide()
        self.groupBoxGasCombustion.hide()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        # Initialize velocity
        choice = self.__boundary.getVelocityChoice(self.__currentField)
        self.modelVelocity.setItem(str_model=choice)
        self.__updateLabel()

        if choice[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            self.lineEditVelocity.setEnabled(False)
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity(self.__currentField)
            self.lineEditVelocity.setText(str(v))

        # Initialize direction
        choice = self.__boundary.getDirectionChoice(self.__currentField)
        self.modelDirection.setItem(str_model=choice)
        text = self.modelDirection.dicoM2V[choice]
        if choice == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            self.frameDirectionCoordinates.hide()
        elif choice == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.frameDirectionCoordinates.hide()

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
        self.__boundary.setVelocityChoice(self.__currentField, c)

        if c[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            self.lineEditVelocity.setEnabled(False)
            self.lineEditVelocity.setText("")
            exp = self.__boundary.getVelocity(self.__currentField)
            if exp:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: green")
                self.pushButtonVelocityFormula.setToolTip(exp)
            else:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: red")
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.pushButtonVelocityFormula.setStyleSheet("background-color: None")
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity(self.__currentField)
            self.lineEditVelocity.setText(str(v))

        self.__updateLabel()


    def __updateLabel(self):
        """
        Update the unit for the velocity specification.
        """
        c = self.__boundary.getVelocityChoice(self.__currentField)
        if c in ('norm', 'norm_formula'):
            self.labelUnitVelocity.setText(str('m/s'))
        elif c in ('flow1', 'flow1_formula'):
            self.labelUnitVelocity.setText(str('kg/s'))


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
            self.__boundary.setVelocity(self.__currentField, v)


    @pyqtSlot()
    def __slotVelocityFormula(self):
        """
        """
        exp = self.__boundary.getVelocity(self.__currentField)
        c = self.__boundary.getVelocityChoice(self.__currentField)
        if c == 'norm_formula':
            req = [('u_norm', 'Norm of the velocity')]
            exa = "u_norm = 1.0;"
        elif c == 'flow1_formula':
            req = [('q_m', 'Mass flow rate')]
            exa = "q_m = 1.0;"

        if c == 'norm_formula':
            sym = [('x', "X face's gravity center"),
                   ('y', "Y face's gravity center"),
                   ('z', "Z face's gravity center"),
                   ('dt', 'time step'),
                   ('t', 'current time'),
                   ('iter', 'number of iteration')]

        elif c == 'flow1_formula':
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
            self.__boundary.setVelocity(self.__currentField, result)
            self.pushButtonVelocityFormula.setToolTip(result)
            self.pushButtonVelocityFormula.setStyleSheet("background-color: green")


    @pyqtSlot(str)
    def __slotChoiceDirection(self, text):
        """
        Input the direction type choice.
        """
        c = self.modelDirection.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setDirectionChoice(self.__currentField, c)

        if c == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            self.frameDirectionCoordinates.hide()
            exp = self.__boundary.getDirection(self.__currentField, 'direction_formula')
            if exp:
                self.pushButtonDirectionFormula.setStyleSheet("background-color: green")
                self.pushButtonDirectionFormula.setToolTip(exp)
            else:
                self.pushButtonDirectionFormula.setStyleSheet("background-color: red")
        elif c == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.pushButtonDirectionFormula.setStyleSheet("background-color: None")
            self.frameDirectionCoordinates.hide()


    @pyqtSlot()
    def __slotDirectionFormula(self):
        """
        """
        exp = self.__boundary.getDirection(self.__currentField, 'direction_formula')

        req = [('dir_x', 'Direction of the flow along X'),
               ('dir_y', 'Direction of the flow along Y'),
               ('dir_z', 'Direction of the flow along Z')]

        exa = "dir_x = 3.0;\ndir_y = 1.0;\ndir_z = 0.0;\n"

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
            log.debug("slotFormulaDirection -> %s" % str(result))
            self.__boundary.setDirection(self.__currentField, 'direction_formula', result)
            self.pushButtonDirectionFormula.setToolTip(result)
            self.pushButtonDirectionFormula.setStyleSheet("background-color: green")


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
