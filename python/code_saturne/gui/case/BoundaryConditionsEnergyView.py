# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
- BoundaryConditionsEnergyView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsEnergy import Ui_BoundaryConditionsEnergy

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.model.MainFieldsModel import MainFieldsModel
from code_saturne.model.ThermodynamicsModel import ThermodynamicsModel
from code_saturne.model.ConjugateHeatTransferModel import ConjugateHeatTransferModel
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsEnergyView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Internal dictionary for units
#-------------------------------------------------------------------------------
bc_energy_units = {'dirichlet':'J/kg',
                   'dirichlet_formula':'J/kg',
                   'flux':'W/m2',
                   'flux_formula':'W/m2',
                   'temperature':'K',
                   'temperature_formula':'K',
                   'timp_K':'K',
                   'timp_K_formula':'K'}

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsEnergyView(QWidget, Ui_BoundaryConditionsEnergy) :
    """
    Boundary condition for energy
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsEnergy.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxEnergy.activated[str].connect(self.__slotChoiceEnergy)

        self.__modelEnergy = ComboModel(self.comboBoxEnergy, 1, 1)

        self.lineEditEnergy.textChanged[str].connect(self.__slotEnergy)

        validatorEner = DoubleValidator(self.lineEditEnergy)

        self.lineEditEnergy.setValidator(validatorEner)

        self.lineEditSyrthes.editingFinished.connect(self.__slotEnergySyrthes)

        # MEG formula
        self.pushButtonEnergy.clicked.connect(self.slotThermalFormula)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None
        self.__currentField = fieldId
        self.cht_model = ConjugateHeatTransferModel(self.case)
        self.notebook  = NotebookModel(self.case)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        # By default hide the syrthes instance line
        self.groupBoxSyrthes.hide()

        self.__modelEnergy.flushItems()

        _nature = self.__boundary.getNature()

        if _nature in ['inlet', 'outlet', 'wall']:
            self.__modelEnergy.addItem(self.tr("Imposed enthalpy"), 'dirichlet')
            self.__modelEnergy.addItem(self.tr("Imposed enthalpy (user law)"), 'dirichlet_formula')
            self.__modelEnergy.addItem(self.tr("Thermal flux"), 'flux')
            self.__modelEnergy.addItem(self.tr("Thermal flux (user law)"), 'flux_formula')
            if _nature == "wall":
                self.__modelEnergy.addItem(self.tr("Temperature"), 'temperature')
                self.__modelEnergy.addItem(self.tr("Temperature (user law)"), 'temperature_formula')
            else:
                self.__modelEnergy.addItem(self.tr("Temperature"), 'timp_K')
                self.__modelEnergy.addItem(self.tr("Temperature (user law)"), 'timp_K_formula')
            self.__modelEnergy.addItem(self.tr("Saturation enthalpy"), 'hsat_P')
            self.__modelEnergy.addItem(self.tr("SYRTHES coupling"), "syrthes_coupling")

        if _nature in ['inlet', 'outlet']:
            if MainFieldsModel(self.case).getFieldFromId(self.__currentField).enthalpy_model != 'off' :
                # No Syrthes coupling for inlet/outlet
                self.__modelEnergy.disableItem(str_model='syrthes_coupling')
                self.groupBoxSyrthes.hide()

                energychoice = self.__boundary.getEnthalpyChoice(self.__currentField)
                self.__modelEnergy.setItem(str_model=energychoice)

                self.__setSubWidgetsStatus(energychoice)
                if energychoice in ['dirichlet', 'flux', 'timp_K']:
                    val = self.__boundary.getEnthalpy(self.__currentField)
                    self.lineEditEnergy.setText(str(val))
                else:
                    self.pushButtonEnergy.setEnabled(True)

                self.show()
                _s = not ThermodynamicsModel(self.case).getMaterials(self.__currentField) == 'user_material'
                for c in ['hsat_P']:
                    self.__modelEnergy.setItemStatus(enable=_s, str_model=c)
            else :
                self.hideWidget()
        elif self.__boundary.getNature() == "wall" :
            if MainFieldsModel(self.case).hasEnthalpyResolvedField() :
                for sm in ['dirichlet', 'dirichlet_formula', 'hsat_P']:
                    self.__modelEnergy.setItemStatus(enable=False, str_model=sm)
                energychoice = self.__boundary.getEnthalpyChoice("none")

                if energychoice in ['dirichlet', 'dirichlet_formula', 'hsat_P']:
                    energychoice = 'flux'

                self.__modelEnergy.setItem(str_model=energychoice)

                self.__setSubWidgetsStatus(energychoice)

                val = self.__boundary.getEnthalpy("none")

                if energychoice == 'syrthes_coupling':
                    syrCompleter = \
                    QCompleter(self.cht_model.getSyrthesInstancesList())
                    self.lineEditSyrthes.setCompleter(syrCompleter)

                    self.lineEditSyrthes.setText(str(val))

                else:
                    if 'formula' != energychoice[-7:]:
                        self.lineEditEnergy.setText(str(val))

                self.show()
            else :
                self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()

    def __setSubWidgetsStatus(self, energyChoice):
        """
        Set hide/show status of subwidgets
        """

        # Update visibility status for different widgets
        _show_syr_box = True if energyChoice == 'syrthes_coupling' else False
        _need_formula = True if "_formula" in energyChoice else False
        _show_val = \
        True if (not _show_syr_box and energyChoice != "hsat_P") else False

        self.groupBoxSyrthes.setVisible(_show_syr_box)


        self.groupBoxEnergyValue.setVisible(_show_val)
        self.lineEditEnergy.setEnabled(not _need_formula)
#        self.labelValue.setVisible(_show_val)
#        self.pushButtonEnergy.setVisible(_show_button)
        self.pushButtonEnergy.setEnabled(_need_formula)
        if _show_val:
            self.labelEnergy.setText(bc_energy_units.get(energyChoice,''))


    @pyqtSlot(str)
    def __slotChoiceEnergy(self, text):
        """
        INPUT choice of method of calculation of the energy
        """
        energy_choice = self.__modelEnergy.dicoV2M[str(text)]
        old_choice = self.__boundary.getEnthalpyChoice(self.__currentField)
        # If we switch from coupling to BC condition, delete all previous
        # data
        if old_choice != energy_choice and old_choice == "syrthes_coupling":
            syrthes_name = self.__boundary.getEnthalpy(self.__currentField)
            bnd_label = self.__boundary.getLabel()
            self.cht_model.deleteSyrthesCoupling(syrthes_name, bnd_label)
            self.__boundary.setEnthalpy(self.__currentField, 0.)

        # Handle formula <-> constant conversion
        if old_choice != energy_choice:
           if old_choice[-7:] == 'formula' and energy_choice[-7:] != 'formula':
               self.__boundary.setEnthalpy(self.__currentField, "0.")
               self.pushButtonEnergy.setStyleSheet("background-color: grey")


        self.__boundary.setEnthalpyChoice(self.__currentField, energy_choice)

        if energy_choice != "hsat_P":
            val = self.__boundary.getEnthalpy(self.__currentField)
            if energy_choice == "syrthes_coupling":
                if val not in self.cht_model.getSyrthesInstancesList():
                    val = ""
                    self.__boundary.setEnthalpy(self.__currentField, val)
                self.lineEditSyrthes.setText(str(val))
            else:
                if energy_choice[-7:] != "formula":
                    self.lineEditEnergy.setText(str(val))
                else:
                    self.lineEditEnergy.setText(str(''))
                    if val not in ['', None]:
                        self.pushButtonEnergy.setToolTip(val)
                        self.pushButtonEnergy.setStyleSheet("background-color: green")
                    else:
                        self.pushButtonEnergy.setStyleSheet("background-color: red")
                self.labelEnergy.setText(bc_energy_units.get(energy_choice, ''))

        self.__setSubWidgetsStatus(energy_choice)


    @pyqtSlot(str)
    def __slotEnergy(self, text):
        """
        INPUT energy value
        """
        if self.lineEditEnergy.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setEnthalpy(self.__currentField, value)


    @pyqtSlot()
    def __slotEnergySyrthes(self):
        """
        Input syrthes instance name
        """

        value = str(self.lineEditSyrthes.text())
        if value:
            value_p = val = self.__boundary.getEnthalpy("none")
            bnd_label = self.__boundary.getLabel()
            if value != value_p:
                self.cht_model.deleteSyrthesCoupling(value_p, bnd_label)


            if value not in self.cht_model.getSyrthesInstancesList():
                self.cht_model.addSyrthesCoupling(value)

            self.cht_model.addBoundaryLabel(value, bnd_label)

            self.__boundary.setEnthalpy(self.__currentField, value)


    @pyqtSlot()
    def slotThermalFormula(self):
        """
        """
        energy_choice = self.__boundary.getEnthalpyChoice(self.__currentField)
        exp = self.__boundary.getEnthalpy(self.__currentField)

        if energy_choice == "dirichlet_formula":
            req = [('enthalpy', 'Specific enthalpy')]
        elif energy_choice == "flux_formula":
            req = [('flux', 'Heat flux [W/m2]')]
        elif energy_choice in ["temperature_formula", "timp_K_formula"]:
            req = [('temperature', 'Temperature [K]')]

        exa = self.__boundary.getDefaultEnthalpyFormula(energy_choice)

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration'),
               ('surface', 'Boundary zone surface')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent        = self,
                                function_type = "bnd",
                                zone_name     = self.__boundary._label,
                                variable_name = "enthalpy_"+str(self.__currentField),
                                expression    = str(exp),
                                required      = req,
                                symbols       = sym,
                                condition     = energy_choice,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotThermalFormula -> %s" % str(result))
            self.__boundary.setEnthalpy(self.__currentField, result)
            self.pushButtonEnergy.setToolTip(result)
            self.pushButtonEnergy.setStyleSheet("background-color: green")


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
