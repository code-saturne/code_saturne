# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
- BoundaryConditionsElectricalView
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

from code_saturne.Pages.BoundaryConditionsElectricalForm import Ui_BoundaryConditionsElectricalForm
from code_saturne.model.ElectricalModel import ElectricalModel

from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.model.Boundary import Boundary
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsElectricalView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsElectricalView(QWidget, Ui_BoundaryConditionsElectricalForm):
    """
    Boundary condifition for the velocity part
    """
    def __init__(self, parent):
        """
        Constructor.
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsElectricalForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget.
        """
        self.__case = case
        self.__boundary = None
        self.__model = ElectricalModel(self.__case)
        self.species_list = []
        self.notebook = NotebookModel(self.__case)

        self.lineEditValuePotElec.textChanged[str].connect(self.slotPotElec)
        self.lineEditValuePotElecIm.textChanged[str].connect(self.slotPotElecIm)
        self.lineEditValueSpecies.textChanged[str].connect(self.slotSpecies)

        self.pushButtonPotVectorFormula.clicked.connect(self.slotPotVectorFormula)

        self.comboBoxTypePotElec.activated[str].connect(self.slotPotElecChoice)
        self.comboBoxTypePotElecIm.activated[str].connect(self.slotPotElecImChoice)
        self.comboBoxTypePotVector.activated[str].connect(self.slotPotVectorChoice)
        self.comboBoxSpecies.activated[str].connect(self.slotSpeciesChoice)
        self.comboBoxPotVector.activated[str].connect(self.slotPotVectorComponentChoice)

        ## Validators
        validatorPotElec      = DoubleValidator(self.lineEditValuePotElec)
        validatorPotElecIm    = DoubleValidator(self.lineEditValuePotElecIm)
        validatorSpecies      = DoubleValidator(self.lineEditValueSpecies, min=0.)

        self.lineEditValuePotElec.setValidator(validatorPotElec)
        self.lineEditValuePotElecIm.setValidator(validatorPotElecIm)
        self.lineEditValueSpecies.setValidator(validatorSpecies)


    def __setBoundary(self, boundary):
        """
        Set the current boundary
        """
        self.__boundary = boundary

        self.nature  = boundary.getNature()

        self.groupBoxPotElecIm.hide()
        self.groupBoxPotVector.hide()
        self.groupBoxMixture.hide()

        self.modelPotElec   = ComboModel(self.comboBoxTypePotElec, 1, 1)
        self.modelPotElecIm = ComboModel(self.comboBoxTypePotElecIm, 1, 1)
        self.modelPotVector = ComboModel(self.comboBoxTypePotVector, 1, 1)

        self.modelPotElec.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelPotElec.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelPotElec.addItem(self.tr("Prescribed flux"), 'neumann')
        self.modelPotElec.addItem(self.tr("Prescribed flux  (user law)"), 'neumann_formula')
        if self.__model.getScaling() == 'on':
            self.modelPotElec.addItem(self.tr("Implicit value (dpot)"), 'dirichlet_implicit')
        #TODO
        self.modelPotElec.disableItem(1)
        self.modelPotElec.disableItem(3)

        self.potElec = "elec_pot_r"
        self.modelPotElecLabel = ComboModel(self.comboBoxPotElec,1,1)
        self.modelPotElecLabel.addItem(self.tr(self.potElec),self.potElec)
        self.modelPotElecLabel.setItem(str_model = self.potElec)

        self.modelPotElecIm.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelPotElecIm.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelPotElecIm.addItem(self.tr("Prescribed flux"), 'neumann')
        self.modelPotElecIm.addItem(self.tr("Prescribed flux  (user law)"), 'neumann_formula')
        #TODO
        self.modelPotElecIm.disableItem(1)
        self.modelPotElecIm.disableItem(3)

        self.potElecIm = 'elec_pot_i'
        self.modelPotElecImLabel = ComboModel(self.comboBoxPotElecIm,1,1)
        self.modelPotElecImLabel.addItem(self.tr(self.potElecIm),self.potElecIm)
        self.modelPotElecImLabel.setItem(str_model = self.potElecIm)

        self.modelPotVector.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelPotVector.addItem(self.tr("Null flux"), 'neumann')
        self.modelPotVector.addItem(self.tr("Implicit flux"), 'neumann_implicit')
        self.modelPotVector.disableItem(0)

        self.potVect = 'vec_potential'
        self.modelPotVectLabel = ComboModel(self.comboBoxPotVector, 1, 1)
        self.modelPotVectLabel.addItem(self.tr('vec_potential'), 'vec_potential')
        self.modelPotVectLabel.setItem(str_model = self.potVect)

        if self.__model.getElectricalModel() == 'joule':
            if self.__model.getJouleModel() == 'three-phase' or \
               self.__model.getJouleModel() == 'three-phase+Transformer':
                self.groupBoxPotElecIm.show()
        elif self.__model.getElectricalModel() == 'arc':
            self.groupBoxPotVector.show()

            self.species = ""

            if self.nature == 'inlet':
                if self.__model.getGasNumber() > 1:
                    self.groupBoxMixture.show()
                    self.modelSpecies = ComboModel(self.comboBoxSpecies, 1, 1)
                    self.species_list = self.__model.getSpeciesLabelsList()
                    for species in self.species_list:
                        self.modelSpecies.addItem(self.tr(species), species)
                    self.species = self.species_list[0]
                    self.modelSpecies.setItem(str_model = self.species)

        self.initializeVariables()


    def initializeVariables(self):
        """
        Initialize widget
        """
        self.lineEditValuePotElec.hide()
        self.lineEditValuePotElecIm.hide()
        self.lineEditValuePotVector.hide()
        self.labelValuePotVector.hide()
        self.labelValuePotElec.hide()
        self.labelValuePotElecIm.hide()

        self.pushButtonPotVectorFormula.setEnabled(False)
        self.pushButtonPotVectorFormula.setStyleSheet("background-color: None")
        self.pushButtonPotElecFormula.setEnabled(False)
        self.pushButtonPotElecFormula.setStyleSheet("background-color: None")
        self.pushButtonPotElecImFormula.setEnabled(False)
        self.pushButtonPotElecImFormula.setStyleSheet("background-color: None")

        # Initialize electric potential
        self.potElec_type = self.__b.getElecScalarChoice(self.potElec)
        self.modelPotElec.setItem(str_model = self.potElec_type)

        if self.potElec_type == 'dirichlet' or self.potElec_type == 'neumann':
            self.lineEditValuePotElec.show()
            self.labelValuePotElec.show()
            v = self.__b.getElecScalarValue(self.potElec, self.potElec_type)
            self.lineEditValuePotElec.setText(str(v))

        # Initialize imaginary electric potential
        if self.__model.getElectricalModel() == 'joule':
            if self.__model.getJouleModel() == 'three-phase' or \
               self.__model.getJouleModel() == 'three-phase+Transformer':
                self.potElecIm_type = self.__b.getElecScalarChoice(self.potElecIm)
                self.modelPotElecIm.setItem(str_model = self.potElecIm_type)

                if self.potElecIm_type == 'dirichlet' or self.potElecIm_type == 'neumann':
                    self.lineEditValuePotElecIm.show()
                    self.labelValuePotElecIm.show()
                    v = self.__b.getElecScalarValue(self.potElecIm, self.potElecIm_type)
                    self.lineEditValuePotElecIm.setText(str(v))

        # Initialize potential vector
        if self.__model.getElectricalModel() == 'arc':
            self.potVec_type = self.__b.getPotentialVectorChoice(self.potVect)
            self.modelPotVector.setItem(str_model = self.potVec_type)

            if self.potVec_type == 'dirichlet_formula':
                self.pushButtonPotVectorFormula.setEnabled(True)
                exp = self.__b.getElecScalarFormula(self.potVect, self.potVec_type)
                if exp:
                    self.pushButtonPotVectorFormula.setStyleSheet("background-color: green")
                    self.pushButtonPotVectorFormula.setToolTip(exp)
                else:
                    self.pushButtonPotVectorFormula.setStyleSheet("background-color: red")

            # Initialize species
            if self.species :
                v = self.__b.getElecScalarValue(self.species, 'dirichlet')
                self.lineEditValueSpecies.setText(str(v))


    @pyqtSlot(str)
    def slotPotElecChoice(self, text):
        """
        INPUT choice for electric potential type
        """
        potElec_type = self.modelPotElec.dicoV2M[str(text)]
        self.__b.setElecScalarChoice(self.potElec, potElec_type)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotPotElecImChoice(self, text):
        """
        INPUT choice for imaginary electric potential type
        """
        potElecIm_type = self.modelPotElecIm.dicoV2M[str(text)]
        self.__b.setElecScalarChoice(self.potElecIm, potElecIm_type)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotPotVectorChoice(self, text):
        """
        INPUT choice for potential vector type
        """
        potVec_choice = self.modelPotVector.dicoV2M[str(text)]
        self.__b.setPotentialVectorChoice(self.potVect, potVec_choice)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotSpeciesChoice(self, text):
        """
        INPUT species choice
        """
        self.species = self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSlot(str)
    def slotPotVectorComponentChoice(self, text):
        """
        INPUT potential vector component choice
        """
        self.potVect = self.modelPotVectLabel.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSlot(str)
    def slotPotElec(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__b.setElecScalarValue(self.potElec, self.potElec_type, value)


    @pyqtSlot(str)
    def slotPotElecIm(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__b.setElecScalarValue(self.potElecIm, self.potElecIm_type, value)


    @pyqtSlot(str)
    def slotSpecies(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__b.setElecScalarValue(self.species, 'dirichlet', value)


    @pyqtSlot()
    def slotPotVectorFormula(self):
        """
        """
        exp = self.__b.getElecScalarFormula(self.potVect, self.potVec_type)
        exa = """#example: """

        if not exp:
            exp = self.potVect + "[0] = 0;\n" + \
                  self.potVect + "[1] = 0;\n" + \
                  self.potVect + "[2] = 0;\n"

        req = [(self.potVect + "[0]", 'vector potential X'),
               (self.potVect + "[1]", 'vector potential Y'),
               (self.potVect + "[2]", 'vector potential Z')]

        sym = [('x', "X cell's gravity center"),
               ('y', "Y cell's gravity center"),
               ('z', "Z cell's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotPotVectorFormula -> %s" % str(result))
            self.__b.setElecScalarFormula(self.potVect, self.potVec_type, str(result))
            self.pushButtonPotVectorFormula.setToolTip(result)
            self.pushButtonPotVectorFormula.setStyleSheet("background-color: green")


    def showWidget(self, b):
        """
        Show the widget.
        """
        self.__b = b
        if self.__model.getElectricalModel() != 'off':
            label = b.getLabel()
            nature = "joule_" + b.getNature()
            self.__b = Boundary(nature, label, self.__case)
            self.__setBoundary(b)

            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all.
        """
        self.hide()


    def tr(self, text):
        """
        Translation.
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
