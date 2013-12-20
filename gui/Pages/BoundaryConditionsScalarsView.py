# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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
- ScalarsBoundariesView
"""

#-------------------------------------------------------------------------------
# Standard modules
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

from BoundaryConditionsScalarsForm import Ui_BoundaryConditionsScalarsForm

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel, setGreenColor
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.DefineUserScalarsModel import DefineUserScalarsModel
from Pages.ThermalScalarModel import ThermalScalarModel
from Pages.QMeiEditorView import QMeiEditorView
from Pages.Boundary import Boundary
from Pages.CompressibleModel import CompressibleModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsScalarsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsScalarsView(QWidget, Ui_BoundaryConditionsScalarsForm):
    """
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsScalarsForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        self.connect(self.lineEditValueThermal, SIGNAL("textChanged(const QString &)"), self.slotValueThermal)
        self.connect(self.lineEditValueSpecies, SIGNAL("textChanged(const QString &)"), self.slotValueSpecies)
        self.connect(self.lineEditValueMeteo,   SIGNAL("textChanged(const QString &)"), self.slotValueMeteo)
        self.connect(self.lineEditExThermal,    SIGNAL("textChanged(const QString &)"), self.slotExThermal)
        self.connect(self.lineEditExSpecies,    SIGNAL("textChanged(const QString &)"), self.slotExSpecies)
        self.connect(self.lineEditExMeteo,      SIGNAL("textChanged(const QString &)"), self.slotExMeteo)

        self.connect(self.pushButtonThermal,    SIGNAL("clicked()"), self.slotThermalFormula)
        self.connect(self.pushButtonSpecies,    SIGNAL("clicked()"), self.slotSpeciesFormula)
        self.connect(self.pushButtonMeteo,      SIGNAL("clicked()"), self.slotMeteoFormula)
        self.connect(self.comboBoxThermal,      SIGNAL("activated(const QString&)"), self.slotThermalChoice)
        self.connect(self.comboBoxTypeThermal,  SIGNAL("activated(const QString&)"), self.slotThermalTypeChoice)
        self.connect(self.comboBoxSpecies,      SIGNAL("activated(const QString&)"), self.slotSpeciesChoice)
        self.connect(self.comboBoxTypeSpecies,  SIGNAL("activated(const QString&)"), self.slotSpeciesTypeChoice)
        self.connect(self.comboBoxMeteo,        SIGNAL("activated(const QString&)"), self.slotMeteoChoice)
        self.connect(self.comboBoxTypeMeteo,    SIGNAL("activated(const QString&)"), self.slotMeteoTypeChoice)

        ## Validators
        validatorValueThermal = DoubleValidator(self.lineEditValueThermal)
        validatorValueSpecies = DoubleValidator(self.lineEditValueSpecies)
        validatorValueMeteo   = DoubleValidator(self.lineEditValueMeteo)
        validatorExThermal    = DoubleValidator(self.lineEditExThermal)
        validatorExSpecies    = DoubleValidator(self.lineEditExSpecies)
        validatorExMeteo      = DoubleValidator(self.lineEditExMeteo)

        self.lineEditValueThermal.setValidator(validatorValueThermal)
        self.lineEditValueSpecies.setValidator(validatorValueSpecies)
        self.lineEditValueMeteo.setValidator(validatorValueMeteo)
        self.lineEditExThermal.setValidator(validatorExThermal)
        self.lineEditExSpecies.setValidator(validatorExSpecies)
        self.lineEditExMeteo.setValidator(validatorExMeteo)

        self.__case.undoStartGlobal()


    def __setBoundary(self, boundary):
        """
        Set the current boundary
        """
        self.__boundary = boundary

        self.nature  = boundary.getNature()
        self.therm   = ThermalScalarModel(self.__case)
        self.sca_mo  = DefineUserScalarsModel(self.__case)
        self.comp    = CompressibleModel(self.__case)

        self.modelTypeThermal = ComboModel(self.comboBoxTypeThermal, 1, 1)
        self.modelTypeSpecies = ComboModel(self.comboBoxTypeSpecies, 1, 1)
        self.modelTypeMeteo   = ComboModel(self.comboBoxTypeMeteo, 1, 1)

        self.modelTypeThermal.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelTypeSpecies.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelTypeMeteo.addItem(  self.tr("Prescribed value"), 'dirichlet')

        self.modelTypeThermal.addItem(self.tr("Prescribed value  (user law)"), 'dirichlet_formula')
        self.modelTypeSpecies.addItem(self.tr("Prescribed value (user law)"), 'dirichlet_formula')
        self.modelTypeMeteo.addItem(  self.tr("Prescribed value (user law)"), 'dirichlet_formula')

        if self.nature == 'outlet':
            self.modelTypeThermal.addItem(self.tr("Prescribed flux"), 'neumann')
            self.modelTypeSpecies.addItem(self.tr("Prescribed flux"), 'neumann')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed flux"), 'neumann')
        elif self.nature == 'wall':
            self.modelTypeThermal.addItem(self.tr("Prescribed flux"), 'neumann')
            self.modelTypeSpecies.addItem(self.tr("Prescribed flux"), 'neumann')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed flux"), 'neumann')
            self.modelTypeThermal.addItem(self.tr("Prescribed flux (user law)"), 'neumann_formula')
            self.modelTypeSpecies.addItem(self.tr("Prescribed flux (user law)"), 'neumann_formula')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed flux (user law)"), 'neumann_formula')
            self.modelTypeThermal.addItem(self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeSpecies.addItem(self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeMeteo.addItem(  self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeThermal.addItem(self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')
            self.modelTypeSpecies.addItem(self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')
            self.modelTypeMeteo.addItem(  self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')

        self.species = ""
        self.species_list = self.sca_mo.getUserScalarLabelsList()
        for s in self.sca_mo.getScalarsVarianceList():
            if s in self.species_list:
                self.species_list.remove(s)

        self.species = ""
        if self.species_list != []:
            self.groupBoxSpecies.show()
            self.modelSpecies = ComboModel(self.comboBoxSpecies, 1, 1)
            for species in self.species_list:
                self.modelSpecies.addItem(self.tr(species), species)
            self.species = self.species_list[0]
            self.modelSpecies.setItem(str_model = self.species)
        else:
            self.groupBoxSpecies.hide()

        self.meteo_list = ""
        self.meteo_list = self.sca_mo.getMeteoScalarsList()

        self.groupBoxMeteo.hide()

        if (self.meteo_list and (self.nature == 'inlet' or self.nature == 'outlet')):
            label = self.__boundary.getLabel()
            nature = "meteo_" + self.nature
            bb = Boundary(nature, label, self.__case)

            if bb.getMeteoDataStatus() == 'off':
                self.groupBoxMeteo.show()
                self.modelMeteo = ComboModel(self.comboBoxMeteo, 1, 1)
                for m in self.meteo_list:
                    self.modelMeteo.addItem(self.tr(m), m)
                self.meteo = self.meteo_list[0]
                self.modelMeteo.setItem(str_model = self.meteo)

        self.model_th = self.therm.getThermalScalarModel()
        if self.model_th != 'off' and self.comp.getCompressibleModel() == 'off':
            self.groupBoxThermal.show()
            self.modelThermal = ComboModel(self.comboBoxThermal,1,1)
            self.thermal = self.therm.getThermalScalarLabel()
            self.modelThermal.addItem(self.tr(self.thermal),self.thermal)
            self.modelThermal.setItem(str_model = self.thermal)
        else:
            self.groupBoxThermal.hide()

        self.initializeVariables()


    def initializeVariables(self):
        """
        Initialize widget
        """
        # Initalize exchange coef
        self.lineEditExThermal.hide()
        self.labelExThermal.hide()
        self.lineEditExSpecies.hide()
        self.labelExSpecies.hide()
        self.lineEditExMeteo.hide()
        self.labelExMeteo.hide()

        # Initalize thermal
        self.lineEditValueThermal.hide()
        self.labelValueThermal.hide()
        self.pushButtonThermal.setEnabled(False)
        setGreenColor(self.pushButtonThermal, False)

        if self.model_th != 'off' and self.comp.getCompressibleModel() == 'off':
            self.thermal_type = self.__boundary.getScalarChoice(self.thermal)
            self.modelTypeThermal.setItem(str_model = self.thermal_type)
            self.labelValueThermal.setText('Value')
            self.groupBoxThermal.setTitle('Thermal')

            if self.thermal_type in ('dirichlet', 'exchange_coefficient', 'neumann'):
                self.labelValueThermal.show()
                self.lineEditValueThermal.show()

                if self.thermal_type == 'exchange_coefficient':
                    self.lineEditExThermal.show()
                    self.labelExThermal.show()
                    v = self.__boundary.getScalarValue(self.thermal, 'dirichlet')
                    w = self.__boundary.getScalarValue(self.thermal, 'exchange_coefficient')
                    self.lineEditValueThermal.setText(str(v))
                    self.lineEditExThermal.setText(str(w))
                else:
                    v = self.__boundary.getScalarValue(self.thermal, self.thermal_type)
                    self.lineEditValueThermal.setText(str(v))

                if self.thermal_type == 'neumann':
                    self.labelValueThermal.setText('Flux')
                    if self.nature == 'outlet':
                        self.groupBoxThermal.setTitle('Thermal for backflow')

            elif self.thermal_type in ('exchange_coefficient_formula', 'dirichlet_formula', 'neumann_formula'):
                self.pushButtonThermal.setEnabled(True)
                setGreenColor(self.pushButtonThermal, True)

        # Initalize species
        self.labelValueSpecies.hide()
        self.lineEditValueSpecies.hide()
        self.pushButtonSpecies.setEnabled(False)
        setGreenColor(self.pushButtonSpecies, False)

        if self.species_list != None and self.species_list != []:
            self.species_type = self.__boundary.getScalarChoice(self.species)
            self.modelTypeSpecies.setItem(str_model = self.species_type)
            self.labelValueSpecies.setText('Value')
            self.groupBoxSpecies.setTitle('Species')

            if self.species_type in ('dirichlet', 'exchange_coefficient', 'neumann'):
                self.labelValueSpecies.show()
                self.lineEditValueSpecies.show()

                if self.species_type == 'exchange_coefficient':
                    self.lineEditExSpecies.show()
                    self.labelExSpecies.show()
                    v = self.__boundary.getScalarValue(self.species, 'dirichlet')
                    w = self.__boundary.getScalarValue(self.species, 'exchange_coefficient')
                    self.lineEditValueSpecies.setText(str(v))
                    self.lineEditExSpecies.setText(str(w))
                else:
                    v = self.__boundary.getScalarValue(self.species, self.species_type)
                    self.lineEditValueSpecies.setText(str(v))

                if self.species_type == 'neumann':
                    self.labelValueSpecies.setText('Flux')
                    if self.nature == 'outlet':
                        self.groupBoxSpecies.setTitle('Species for backflow')

            elif self.species_type in ('exchange_coefficient_formula', 'dirichlet_formula', 'neumann_formula'):
                self.pushButtonSpecies.setEnabled(True)
                setGreenColor(self.pushButtonSpecies, True)

        # Initalize meteo
        self.labelValueMeteo.hide()
        self.lineEditValueMeteo.hide()
        self.pushButtonMeteo.setEnabled(False)
        setGreenColor(self.pushButtonMeteo, False)

        if (self.meteo_list and (self.nature == 'inlet' or self.nature == 'outlet')):
            label = self.__boundary.getLabel()
            nature = "meteo_" + self.nature
            bb = Boundary(nature, label, self.__case)

            if bb.getMeteoDataStatus() == 'off':
                self.meteo_type = self.__boundary.getScalarChoice(self.meteo)
                self.modelTypeMeteo.setItem(str_model = self.meteo_type)
                self.labelValueMeteo.setText('Value')
                self.groupBoxMeteo.setTitle('Meteo')

                if self.meteo_type in ('dirichlet', 'exchange_coefficient', 'neumann'):
                    self.labelValueMeteo.show()
                    self.lineEditValueMeteo.show()

                    if self.meteo_type == 'exchange_coefficient':
                        self.lineEditExMeteo.show()
                        self.labelExMeteo.show()
                        v = self.__boundary.getScalarValue(self.meteo, 'dirichlet')
                        w = self.__boundary.getScalarValue(self.meteo, 'exchange_coefficient')
                        self.lineEditValueMeteo.setText(str(v))
                        self.lineEditExMeteo.setText(str(w))
                    else:
                        v = self.__boundary.getScalarValue(self.meteo, self.meteo_type)
                        self.lineEditValueMeteo.setText(str(v))

                if self.meteo_type == 'neumann':
                    self.labelValueMeteo.setText('Flux')
                    if self.nature == 'outlet':
                        self.groupBoxMeteo.setTitle('Meteo for backflow')

                if self.meteo_type in ('exchange_coefficient_formula', 'dirichlet_formula', 'neumann_formula'):
                    self.pushButtonMeteo.setEnabled(True)
                    setGreenColor(self.pushButtonMeteo, True)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        if DefineUserScalarsModel(self.__case).getScalarLabelsList() or\
           DefineUserScalarsModel(self.__case).getMeteoScalarsList() or\
           DefineUserScalarsModel(self.__case).getThermalScalarLabelsList():
            self.__setBoundary(boundary)
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSignature("const QString&")
    def slotThermalChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.thermal = self.modelThermal.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotThermalTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.thermal_type = self.modelTypeThermal.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.thermal, self.thermal_type)
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotSpeciesChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.species = self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotSpeciesTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.species_type = self.modelTypeSpecies.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.species, self.species_type)
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotMeteoChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.meteo = self.modelMeteo.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSignature("const QString&")
    def slotMeteoTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.meteo_type= self.modelTypeMeteo.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.meteo, self.meteo_type)
        self.initializeVariables()


    @pyqtSignature("")
    def slotThermalFormula(self):
        """
        """
        exp = self.__boundary.getScalarFormula(self.thermal, self.thermal_type)
        exa = """#example: """
        if self.thermal_type == 'dirichlet_formula':
            req = [(self.thermal, str(self.thermal))]
        elif self.thermal_type == 'neumann_formula':
            req = [("flux", "flux")]
        elif self.thermal_type == 'exchange_coefficient_formula':
            req = [(self.thermal, str(self.thermal)),("hc", "heat coefficient")]

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
            log.debug("slotThermalFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.thermal, self.thermal_type, result)
            setGreenColor(self.pushButtonThermal, False)


    @pyqtSignature("")
    def slotSpeciesFormula(self):
        """
        """
        exp = self.__boundary.getScalarFormula(self.species, self.species_type)
        exa = """#example: """
        if self.species_type == 'dirichlet_formula':
            req = [(self.species, str(self.species))]
        elif self.species_type == 'neumann_formula':
            req = [("flux", "flux")]
        elif self.species_type == 'exchange_coefficient_formula':
            req = [(self.species, str(self.species)),("hc", "heat coefficient")]

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
            log.debug("slotSpeciesFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.species, self.species_type, result)
            setGreenColor(self.pushButtonSpecies, False)


    @pyqtSignature("")
    def slotMeteoFormula(self):
        """
        """
        exp = self.__boundary.getScalarFormula(self.meteo, self.meteo_type)
        exa = """#example: """
        if self.meteo_type == 'dirichlet_formula':
            req = [(self.meteo, str(self.meteo))]
        elif self.meteo_type == 'neumann_formula':
            req = [("flux", "flux")]
        elif self.meteo_type == 'exchange_coefficient_formula':
            req = [(self.meteo, str(self.meteo)),("hc", "heat coefficient")]

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
            log.debug("slotMeteoFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.meteo, self.meteo_type, result)
            setGreenColor(self.pushButtonMeteo, False)


    @pyqtSignature("const QString&")
    def slotValueThermal(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            if self.thermal_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.thermal, self.thermal_type, value)
            elif self.thermal_type == 'exchange_coefficient':
                self.__boundary.setScalarValue(self.thermal, 'dirichlet', value)


    @pyqtSignature("const QString&")
    def slotValueSpecies(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            if self.species_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.species, self.species_type, value)
            elif self.species_type == 'exchange_coefficient' :
                self.__boundary.setScalarValue(self.species, 'dirichlet', value)


    @pyqtSignature("const QString&")
    def slotValueMeteo(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            if self.meteo_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.meteo, self.meteo_type, value)
            elif self.meteo_type == 'exchange_coefficient':
                self.__boundary.setScalarValue(self.meteo, 'dirichlet', value)


    @pyqtSignature("const QString&")
    def slotExThermal(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__boundary.setScalarValue(self.thermal, 'exchange_coefficient', value)


    @pyqtSignature("const QString&")
    def slotExSpecies(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__boundary.setScalarValue(self.species, 'exchange_coefficient', value)


    @pyqtSignature("const QString&")
    def slotExMeteo(self, var):
        """
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = float(var)
            self.__boundary.setScalarValue(self.meteo, 'exchange_coefficient', value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
