# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
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

#-------------------------------------------------------------------------------

"""
This module contains the following classes:
- BoundaryConditionsHydraulicHeadView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant

from code_saturne.Pages.BoundaryConditionsHydraulicHeadForm import \
     Ui_BoundaryConditionsHydraulicHeadForm
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsHydraulicHeadView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsHydraulicHeadView(QWidget, Ui_BoundaryConditionsHydraulicHeadForm):
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsHydraulicHeadForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None
        self.notebook = NotebookModel(self.__case)

        self.__case.undoStopGlobal()

        # Validators
        validatorHh   = DoubleValidator(self.lineEditValueHydraulicHead)
        validatorExHh = DoubleValidator(self.lineEditExHydraulicHead)

        # Apply validators
        self.lineEditValueHydraulicHead.setValidator(validatorHh)
        self.lineEditExHydraulicHead.setValidator(validatorHh)

        self.modelTypeHydraulic = ComboModel(self.comboBoxTypeHydraulicHead, 1, 1)
        self.modelTypeHydraulic.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelTypeHydraulic.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelTypeHydraulic.addItem(self.tr("Prescribed flux"), 'neumann')

        # Connections
        self.lineEditValueHydraulicHead.textChanged[str].connect(self.slotHydraulicHeadValue)
        self.lineEditExHydraulicHead.textChanged[str].connect(self.slotHydraulicHeadFlux)
        self.pushButtonHydraulicHead.clicked.connect(self.slotHydraulicHeadFormula)
        self.comboBoxTypeHydraulicHead.activated[str].connect(self.slotHydraulicHeadChoice)

        self.__case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        label = boundary.getLabel()
        self.nature  = boundary.getNature()
        self.__boundary = Boundary(self.nature, label, self.__case)
        self.initialize()


    def initialize(self):
        self.labelValueHydraulicHead.hide()
        self.labelExHydraulicHead.hide()
        self.lineEditValueHydraulicHead.hide()
        self.lineEditExHydraulicHead.hide()
        self.pushButtonHydraulicHead.setEnabled(False)
        self.pushButtonHydraulicHead.setStyleSheet("background-color: None")

        HydraulicChoice = self.__boundary.getHydraulicHeadChoice()
        self.modelTypeHydraulic.setItem(str_model = HydraulicChoice)
        if HydraulicChoice == 'dirichlet':
            self.labelValueHydraulicHead.show()
            self.lineEditValueHydraulicHead.show()
            h_head = self.__boundary.getHydraulicHeadValue()
            self.lineEditValueHydraulicHead.setText(str(h_head))
        elif HydraulicChoice == 'neumann':
            self.labelExHydraulicHead.show()
            self.lineEditExHydraulicHead.show()
            h_head = self.__boundary.getHydraulicHeadFlux()
            self.lineEditExHydraulicHead.setText(str(h_head))
        elif HydraulicChoice == 'dirichlet_formula':
            self.pushButtonHydraulicHead.setEnabled(True)

            exp = self.__boundary.getHydraulicHeadFormula()
            if exp:
                self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
                self.pushButtonHydraulicHead.setToolTip(exp)
            else:
                self.pushButtonHydraulicHead.setStyleSheet("background-color: red")

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(str)
    def slotHydraulicHeadValue(self, text):
        """
        INPUT hydraulic head value
        """
        if self.lineEditValueHydraulicHead.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.__boundary.setHydraulicHeadValue(t)


    @pyqtSlot(str)
    def slotHydraulicHeadFlux(self, text):
        """
        INPUT hydraulic head flux
        """
        if self.lineEditExHydraulicHead.validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.__boundary.setHydraulicHeadFlux(t)


    @pyqtSlot()
    def slotHydraulicHeadFormula(self):
        """
        """
        exp = self.__boundary.getHydraulicHeadFormula()
        exa = """#example: """
        req = [("H", "hydraulic head")]

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.__case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotHydraulicHeadFormula -> %s" % str(result))
            self.__boundary.setHydraulicHeadFormula(str(result))
            self.pushButtonHydraulicHead.setStyleSheet("background-color: green")
            self.pushButtonHydraulicHead.setToolTip(result)


    @pyqtSlot(str)
    def slotHydraulicHeadChoice(self, text):
        """
        INPUT label for choice of zone
        """
        HydraulicChoice = self.modelTypeHydraulic.dicoV2M[str(text)]
        self.__boundary.setHydraulicHeadChoice(HydraulicChoice)
        self.initialize()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
