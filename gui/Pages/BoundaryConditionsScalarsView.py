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
- ScalarsBoundariesView
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
from code_saturne.Base.QtPage  import DoubleValidator, ComboModel, from_qvariant

from code_saturne.Pages.BoundaryConditionsScalarsForm import Ui_BoundaryConditionsScalarsForm
from code_saturne.model.LocalizationModel             import LocalizationModel, Zone
from code_saturne.model.DefineUserScalarsModel        import DefineUserScalarsModel
from code_saturne.model.ThermalScalarModel            import ThermalScalarModel
from code_saturne.Pages.QMeiEditorView                import QMeiEditorView
from code_saturne.model.Boundary                      import Boundary
from code_saturne.model.CompressibleModel             import CompressibleModel
from code_saturne.model.AtmosphericFlowsModel         import AtmosphericFlowsModel
from code_saturne.model.NotebookModel                 import NotebookModel

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
        self.notebook = NotebookModel(self.__case)

        self.lineEditValueThermal.textChanged[str].connect(self.slotValueThermal)
        self.lineEditValueSpecies.textChanged[str].connect(self.slotValueSpecies)
        self.lineEditValueMeteo.textChanged[str].connect(self.slotValueMeteo)
        self.lineEditExThermal.textChanged[str].connect(self.slotExThermal)
        self.lineEditExSpecies.textChanged[str].connect(self.slotExSpecies)
        self.lineEditExMeteo.textChanged[str].connect(self.slotExMeteo)

        self.pushButtonThermal.clicked.connect(self.slotThermalFormula)
        self.pushButtonSpecies.clicked.connect(self.slotSpeciesFormula)
        self.pushButtonMeteo.clicked.connect(self.slotMeteoFormula)
        self.comboBoxThermal.activated[str].connect(self.slotThermalChoice)
        self.comboBoxTypeThermal.activated[str].connect(self.slotThermalTypeChoice)
        self.comboBoxSpecies.activated[str].connect(self.slotSpeciesChoice)
        self.comboBoxTypeSpecies.activated[str].connect(self.slotSpeciesTypeChoice)
        self.comboBoxMeteo.activated[str].connect(self.slotMeteoChoice)
        self.comboBoxTypeMeteo.activated[str].connect(self.slotMeteoTypeChoice)

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
        self.atm     = AtmosphericFlowsModel(self.__case)

        self.modelTypeThermal = ComboModel(self.comboBoxTypeThermal, 1, 1)
        self.modelTypeSpecies = ComboModel(self.comboBoxTypeSpecies, 1, 1)
        self.modelTypeMeteo   = ComboModel(self.comboBoxTypeMeteo, 1, 1)

        self.modelTypeThermal.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelTypeSpecies.addItem(self.tr("Prescribed value"), 'dirichlet')
        self.modelTypeMeteo.addItem(  self.tr("Prescribed value"), 'dirichlet')

        self.modelTypeThermal.addItem(self.tr("Prescribed value (user law)"), 'dirichlet_formula')
        self.modelTypeSpecies.addItem(self.tr("Prescribed value (user law)"), 'dirichlet_formula')
        self.modelTypeMeteo.addItem(  self.tr("Prescribed value (user law)"), 'dirichlet_formula')

        if self.nature == 'outlet':
            self.modelTypeThermal.addItem(self.tr("Prescribed (outgoing) flux"), 'neumann')
            self.modelTypeSpecies.addItem(self.tr("Prescribed (outgoing) flux"), 'neumann')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed (outgoing) flux"), 'neumann')
        elif self.nature == 'wall':
            self.modelTypeThermal.addItem(self.tr("Prescribed (outgoing) flux"), 'neumann')
            self.modelTypeSpecies.addItem(self.tr("Prescribed (outgoing) flux"), 'neumann')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed (outgoing) flux"), 'neumann')
            self.modelTypeThermal.addItem(self.tr("Prescribed (outgoing) flux (user law)"), 'neumann_formula')
            self.modelTypeSpecies.addItem(self.tr("Prescribed (outgoing) flux (user law)"), 'neumann_formula')
            self.modelTypeMeteo.addItem(  self.tr("Prescribed (outgoing) flux (user law)"), 'neumann_formula')
            self.modelTypeThermal.addItem(self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeSpecies.addItem(self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeMeteo.addItem(  self.tr("Exchange coefficient"), 'exchange_coefficient')
            self.modelTypeThermal.addItem(self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')
            self.modelTypeSpecies.addItem(self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')
            self.modelTypeMeteo.addItem(  self.tr("Exchange coefficient (user law)"), 'exchange_coefficient_formula')
        elif self.nature == 'groundwater':
            self.modelTypeSpecies.addItem(self.tr("Prescribed (outgoing) flux"), 'neumann')

        self.species = ""
        self.species_list = self.sca_mo.getUserScalarNameList()
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

        self.model_th = self.therm.getThermalScalarModel()
        if self.model_th != 'off' and self.comp.getCompressibleModel() == 'off':
            self.groupBoxThermal.show()
            self.modelThermal = ComboModel(self.comboBoxThermal,1,1)
            self.thermal = self.therm.getThermalScalarName()
            self.modelThermal.addItem(self.tr(self.thermal),self.thermal)
            self.modelThermal.setItem(str_model = self.thermal)
        else:
            self.groupBoxThermal.hide()

        self.meteo_list = ""
        self.meteo_list = self.sca_mo.getMeteoScalarsNameList()

        self.groupBoxMeteo.hide()

        if (self.atm.getAtmosphericFlowsModel() != "off" and self.nature == 'wall'):
            self.modelMeteo = ComboModel(self.comboBoxMeteo, 1, 1)
            if len(self.meteo_list) > 0:
                self.groupBoxMeteo.show()
                for m in self.meteo_list:
                    self.modelMeteo.addItem(self.tr(m), m)
                self.meteo = self.meteo_list[0]
                self.modelMeteo.setItem(str_model = self.meteo)

        if (self.atm.getAtmosphericFlowsModel() != "off" and \
           (self.nature == 'inlet' or self.nature == 'outlet')):
            label = self.__boundary.getLabel()
            nature = "meteo_" + self.nature
            bb = Boundary(nature, label, self.__case)

            if bb.getMeteoDataStatus() == 'off':
                self.groupBoxMeteo.hide()
                self.groupBoxThermal.show()
                self.modelMeteo = ComboModel(self.comboBoxMeteo, 1, 1)
                if len(self.meteo_list) > 0:
                    self.groupBoxMeteo.show()
                    for m in self.meteo_list:
                        self.modelMeteo.addItem(self.tr(m), m)
                    self.meteo = self.meteo_list[0]
                    self.modelMeteo.setItem(str_model = self.meteo)
            else:
                self.groupBoxMeteo.hide()
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
        self.pushButtonThermal.setStyleSheet("background-color: None")

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
                exp = self.__boundary.getScalarFormula(self.thermal, self.thermal_type)
                if exp:
                    self.pushButtonThermal.setStyleSheet("background-color: green")
                    self.pushButtonThermal.setToolTip(exp)
                else:
                    self.pushButtonThermal.setStyleSheet("background-color: red")

        # Initalize species
        self.labelValueSpecies.hide()
        self.lineEditValueSpecies.hide()
        self.pushButtonSpecies.setEnabled(False)
        self.pushButtonSpecies.setStyleSheet("background-color: None")

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
                    if self.nature == 'groundwater':
                        self.labelValueSpecies.setText('Velocity')
                        self.labelExSpecies.setText('Concentration')
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
                exp = self.__boundary.getScalarFormula(self.species, self.species_type)
                if exp:
                    self.pushButtonSpecies.setStyleSheet("background-color: green")
                    self.pushButtonSpecies.setToolTip(exp)
                else:
                    self.pushButtonSpecies.setStyleSheet("background-color: red")

            if self.nature == 'groundwater':
                self.groupBoxSpecies.setTitle('Transport equation')

        # Initalize meteo
        self.labelValueMeteo.hide()
        self.lineEditValueMeteo.hide()
        self.pushButtonMeteo.setEnabled(False)
        self.pushButtonMeteo.setStyleSheet("background-color: None")

        if (self.meteo_list):
            label = self.__boundary.getLabel()
            if self.nature != 'wall':
                nature = "meteo_" + self.nature
            else:
                nature = self.nature
            bb = Boundary(nature, label, self.__case)

            if self.nature == 'wall' or bb.getMeteoDataStatus() == 'off':
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
                    exp = self.__boundary.getScalarFormula(self.meteo, self.meteo_type)
                    if exp:
                        self.pushButtonMeteo.setStyleSheet("background-color: green")
                        self.pushButtonMeteo.setToolTip(exp)
                    else:
                        self.pushButtonMeteo.setStyleSheet("background-color: red")


    def showWidget(self, boundary):
        """
        Show the widget
        """
        if DefineUserScalarsModel(self.__case).getScalarNameList() or\
           DefineUserScalarsModel(self.__case).getMeteoScalarsNameList() or\
           DefineUserScalarsModel(self.__case).getThermalScalarName():
            self.__setBoundary(boundary)
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(str)
    def slotThermalChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.thermal = self.modelThermal.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSlot(str)
    def slotThermalTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.thermal_type = self.modelTypeThermal.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.thermal, self.thermal_type)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotSpeciesChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.species = self.modelSpecies.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSlot(str)
    def slotSpeciesTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.species_type = self.modelTypeSpecies.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.species, self.species_type)
        self.initializeVariables()


    @pyqtSlot(str)
    def slotMeteoChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.meteo = self.modelMeteo.dicoV2M[str(text)]
        self.initializeVariables()


    @pyqtSlot(str)
    def slotMeteoTypeChoice(self, text):
        """
        INPUT label for choice of zone
        """
        self.meteo_type= self.modelTypeMeteo.dicoV2M[str(text)]
        self.__boundary.setScalarChoice(self.meteo, self.meteo_type)
        self.initializeVariables()


    @pyqtSlot()
    def slotThermalFormula(self):
        """
        """
        name = self.thermal
        exp = self.__boundary.getScalarFormula(self.thermal, self.thermal_type)
        exa = """#example: """
        if self.thermal_type == 'dirichlet_formula':
            req = [(name, str(name))]
        elif self.thermal_type == 'neumann_formula':
            req = [("flux", "flux")]
        elif self.thermal_type == 'exchange_coefficient_formula':
            req = [(name, str(name)),("hc", "heat coefficient")]

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
            log.debug("slotThermalFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.thermal, self.thermal_type, str(result))
            self.pushButtonThermal.setStyleSheet("background-color: green")
            self.pushButtonThermal.setToolTip(exp)


    @pyqtSlot()
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
            log.debug("slotSpeciesFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.species, self.species_type, str(result))
            self.pushButtonSpecies.setStyleSheet("background-color: green")
            self.pushButtonSpecies.setToolTip(exp)


    @pyqtSlot()
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
            req = [(self.meteo, str(self.meteo)),
                   ("hc", "heat coefficient")]

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
            log.debug("slotMeteoFormula -> %s" % str(result))
            self.__boundary.setScalarFormula(self.meteo, self.meteo_type, str(result))
            self.pushButtonMeteo.setStyleSheet("background-color: green")
            self.pushButtonMeteo.setToolTip(exp)


    @pyqtSlot(str)
    def slotValueThermal(self, var):
        """
        """
        if self.lineEditValueThermal.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            if self.thermal_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.thermal, self.thermal_type, value)
            elif self.thermal_type == 'exchange_coefficient':
                self.__boundary.setScalarValue(self.thermal, 'dirichlet', value)


    @pyqtSlot(str)
    def slotValueSpecies(self, var):
        """
        """
        if self.lineEditValueSpecies.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            if self.species_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.species, self.species_type, value)
            elif self.species_type == 'exchange_coefficient' :
                self.__boundary.setScalarValue(self.species, 'dirichlet', value)


    @pyqtSlot(str)
    def slotValueMeteo(self, var):
        """
        """
        if self.lineEditValueMeteo.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            if self.meteo_type in ('dirichlet', 'neumann'):
                self.__boundary.setScalarValue(self.meteo, self.meteo_type, value)
            elif self.meteo_type == 'exchange_coefficient':
                self.__boundary.setScalarValue(self.meteo, 'dirichlet', value)


    @pyqtSlot(str)
    def slotExThermal(self, var):
        """
        """
        if self.lineEditExThermal.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__boundary.setScalarValue(self.thermal, 'exchange_coefficient', value)


    @pyqtSlot(str)
    def slotExSpecies(self, var):
        """
        """
        if self.lineEditExSpecies.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__boundary.setScalarValue(self.species, 'exchange_coefficient', value)


    @pyqtSlot(str)
    def slotExMeteo(self, var):
        """
        """
        if self.lineEditExMeteo.validator().state == QValidator.Acceptable:
            value = from_qvariant(var, float)
            self.__boundary.setScalarValue(self.meteo, 'exchange_coefficient', value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
