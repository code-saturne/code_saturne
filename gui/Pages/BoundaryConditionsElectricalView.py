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
This module contains the following classes:
- BoundaryConditionsElectricalView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

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

from Pages.BoundaryConditionsElectricalForm import Ui_BoundaryConditionsElectricalForm
from Pages.ElectricalModel import ElectricalModel

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel, setGreenColor
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.QMeiEditorView import QMeiEditorView
from Pages.Boundary import Boundary

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

        self.connect(self.lineEditValuePotElec,   SIGNAL("textChanged(const QString &)"), self.slotPotElec)
        self.connect(self.lineEditValuePotElecIm, SIGNAL("textChanged(const QString &)"), self.slotPotElecIm)
        self.connect(self.lineEditValueSpecies,   SIGNAL("textChanged(const QString &)"), self.slotSpecies)

        self.connect(self.pushButtonPotVectorFormula, SIGNAL("clicked()"), self.slotPotVectorFormula)

        self.connect(self.comboBoxTypePotElec,   SIGNAL("activated(const QString&)"), self.slotPotElecChoice)
        self.connect(self.comboBoxTypePotElecIm, SIGNAL("activated(const QString&)"), self.slotPotElecImChoice)
        self.connect(self.comboBoxTypePotVector, SIGNAL("activated(const QString&)"), self.slotPotVectorChoice)
        self.connect(self.comboBoxSpecies,       SIGNAL("activated(const QString&)"), self.slotSpeciesChoice)

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

        self.potElec = self.__model.getScalarLabel('PotElecReal')
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

        self.potElecIm = self.__model.getScalarLabel('POT_EL_I')
        self.modelPotElecImLabel = ComboModel(self.comboBoxPotElecIm,1,1)
        self.modelPotElecImLabel.addItem(self.tr(self.potElecIm),self.potElecIm)
        self.modelPotElecImLabel.setItem(str_model = self.potElecIm)

        self.modelPotVector.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelPotVector.addItem(self.tr("Null flux"), 'neumann')
        self.modelPotVector.addItem(self.tr("Implicit flux"), 'neumann_implicit')
        self.modelPotVector.disableItem(0)

        self.potVect = self.__model.getScalarLabel('POT_VEC01')
        self.modelPotVectLabel = ComboModel(self.comboBoxPotVector,1,1)
        self.modelPotVectLabel.addItem(self.tr(self.potVect),self.potVect)
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
        setGreenColor(self.pushButtonPotVectorFormula, False)
        self.pushButtonPotElecFormula.setEnabled(False)
        setGreenColor(self.pushButtonPotElecFormula, False)
        self.pushButtonPotElecImFormula.setEnabled(False)
        setGreenColor(self.pushButtonPotElecImFormula, False)

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
                setGreenColor(self.pushButtonPotVectorFormula, True)

            # Initialize species
            if self.species :
                v = self.__b.getElecScalarValue(self.species, 'dirichlet')
                self.lineEditValueSpecies.setText(str(v))


    @pyqtSignature("const QString&")
    def slotPotElecChoice(self, text):
        """
        INPUT choice for electric potential type
        """
        potElec_type = self.modelPotElec.dicoV2M[str(text)]
        self.__b.setElecScalarChoice(self.potElec, potElec_type)
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotPotElecImChoice(self, text):
        """
        INPUT choice for imaginary electric potential type
        """
        potElecIm_type = self.modelPotElecIm.dicoV2M[str(text)]
        self.__b.setElecScalarChoice(self.potElecIm, potElecIm_type)
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotPotVectorChoice(self, text):
        """
        INPUT choice for potential vector type
        """
        potVec_choice = self.modelPotVector.dicoV2M[str(text)]
        self.__b.setPotentialVectorChoice(self.potVect, potVec_choice)
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotSpeciesChoice(self, text):
        """
        INPUT species choice
        """
        self.species = self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotPotElec(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__b.setElecScalarValue(self.potElec, self.potElec_type, value)


    @pyqtSignature("const QString&")
    def slotPotElecIm(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__b.setElecScalarValue(self.potElecIm, self.potElecIm_type, value)


    @pyqtSignature("const QString&")
    def slotSpecies(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__b.setElecScalarValue(self.species, 'dirichlet', value)


    @pyqtSignature("")
    def slotPotVectorFormula(self):
        """
        """
        exp = self.__b.getElecScalarFormula(self.potVect, self.potVec_type)
        exa = """#example: """

        if not exp:
            exp = "Ax = 0;\nAy = 0;\nAz = 0;"
        req = [('Ax', 'vector potential X'),
               ('Ay', 'vector potential Y'),
               ('Az', 'vector potential Z')]

        sym = [('x', "X cell's gravity center"),
               ('y', "Y cell's gravity center"),
               ('z', "Z cell's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        dialog = QMeiEditorView(self,expression = exp,
                                 required   = req,
                                 symbols    = sym,
                                 examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotPotVectorFormula -> %s" % str(result))
            self.__b.setElecScalarFormula(self.potVect, self.potVec_type, result)
            setGreenColor(self.pushButtonPotVectorFormula, False)


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

