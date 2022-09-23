# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module defines the values of reference.

This module contains the following classes and function:
- FluidStructureInteractionAdvancedOptionsView
- StandardItemModel
- Coupling
- FluidStructureInteractionView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtPage import DoubleValidator, IntValidator
from code_saturne.gui.base.QtPage import from_qvariant
from code_saturne.gui.case.FluidStructureInteractionForm  import Ui_FluidStructureInteractionForm
from code_saturne.model.FluidStructureInteractionModel import FluidStructureInteractionModel
from code_saturne.gui.case.FluidStructureInteractionAdvancedOptionsDialogForm import \
Ui_FluidStructureInteractionAdvancedOptionsDialogForm

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FluidStructureInteractionView")

#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------

displacement_prediction_alpha          = 'displacement_prediction_alpha'
displacement_prediction_beta           = 'displacement_prediction_beta'
stress_prediction_alpha                = 'stress_prediction_alpha'
monitor_point_synchronisation          = 'monitor_point_synchronisation'

#-------------------------------------------------------------------------------
# Advanced dialog
#-------------------------------------------------------------------------------

class FluidStructureInteractionAdvancedOptionsView(QDialog,
                        Ui_FluidStructureInteractionAdvancedOptionsDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        # Init base classes
        QDialog.__init__(self, parent)
        Ui_FluidStructureInteractionAdvancedOptionsDialogForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        title = self.tr("Displacements prediction:")
        self.setWindowTitle(title)

        self.__default = default
        self.__result  = default.copy()
        self.__setValidator()
        self.__setInitialValues()

        self.case.undoStartGlobal()


    def __setValidator(self):
        """
        Set the validator
        """
        validator = DoubleValidator(self.lineEditDisplacementAlpha,
                                    min=0.0)
        self.lineEditDisplacementAlpha.setValidator(validator)

        validator = DoubleValidator(self.lineEditDisplacementBeta,
                                    min=0.0)
        self.lineEditDisplacementBeta.setValidator(validator)

        validator = DoubleValidator(self.lineEditStressAlpha, min=0.0)
        self.lineEditStressAlpha.setValidator(validator)


    def __setInitialValues(self):
        """
        Set the initial values for the 4 widgets
        """
        # Read from default
        displacementAlpha = str(self.__default[displacement_prediction_alpha])
        displacementBeta  = str(self.__default[displacement_prediction_beta ])
        stressAlpha       = str(self.__default[stress_prediction_alpha      ])

        isSynchronizationOn = self.__default[monitor_point_synchronisation] == 'on'

        # Update Widget
        self.lineEditDisplacementAlpha.setText(displacementAlpha)
        self.lineEditDisplacementBeta.setText(displacementBeta)
        self.lineEditStressAlpha.setText(stressAlpha)
        self.checkBoxSynchronization.setChecked(isSynchronizationOn)


    def get_result(self):
        """
        Method to get the result
        """
        return self.__result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        # Read value from widget
        displacementAlpha = float(self.lineEditDisplacementAlpha.text())
        displacementBeta  = float(self.lineEditDisplacementBeta.text())
        stressAlpha       = float(self.lineEditStressAlpha.text())

        if self.checkBoxSynchronization.isChecked():
            synchronization = 'on'
        else:
            synchronization = 'off'

        # Set result attributes
        self.__result[displacement_prediction_alpha] = displacementAlpha
        self.__result[displacement_prediction_beta ] = displacementBeta
        self.__result[stress_prediction_alpha      ] = stressAlpha
        self.__result[monitor_point_synchronisation] = synchronization

        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class FluidStructureInteractionView(QWidget, Ui_FluidStructureInteractionForm):
    """
    Main class.
    """
    def __init__(self, parent=None):
        """
        Constructor
        """
        # Init base classes
        QWidget.__init__(self, parent)

        Ui_FluidStructureInteractionForm.__init__(self)
        self.setupUi(self)

    def setup(self, case):

        self.case = case
        self.case.undoStopGlobal()
        self.__model = FluidStructureInteractionModel(case)

        self.__defineConnection()
        self.__addValidators()
        self.__setInitialValues()

        self.case.undoStartGlobal()


    def __defineConnection(self):
        """
        Define connection for widgets that do not depend on the boundary
        """
        self.lineEditNALIMX.textChanged[str].connect(self.__slotNalimx)

        self.lineEditEPALIM.textChanged[str].connect(self.__slotEpalim)
        self.pushButtonAdvanced.clicked.connect(self.__slotAdvanced)


    def __addValidators(self):
        """
        Add the validator for NALIMX and EPALIM
        """
        validatorNALIMX = IntValidator(self.lineEditNALIMX, min=1)
        self.lineEditNALIMX.setValidator(validatorNALIMX)

        validatorEPALIM = DoubleValidator(self.lineEditEPALIM, min=0.0)
        validatorEPALIM.setExclusiveMin(True)
        self.lineEditEPALIM.setValidator(validatorEPALIM)


    def __setInitialValues(self):
        """
        Set Widget initial values that do not depend on the boundary
        """
        nalimx = self.__model.getMaxIterations()
        self.lineEditNALIMX.setText(str(nalimx))
        epalim = self.__model.getPrecision()
        self.lineEditEPALIM.setText(str(epalim))


    @pyqtSlot(str)
    def __slotNalimx(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            nalimx = from_qvariant(text, int)
            self.__model.setMaxIterations(nalimx)


    @pyqtSlot(str)
    def __slotEpalim(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            epalim = from_qvariant(text, float)
            self.__model.setPrecision(epalim)


    @pyqtSlot()
    def __slotAdvanced(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        # Set the default value
        default = {}
        default[displacement_prediction_alpha] = self.__model.getDisplacementPredictionAlpha()
        default[displacement_prediction_beta ] = self.__model.getDisplacementPredictionBeta()
        default[stress_prediction_alpha      ] = self.__model.getStressPredictionAlpha()
        default[monitor_point_synchronisation] = \
                            self.__model.getMonitorPointSynchronisation()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        # run the dialog
        dialog = FluidStructureInteractionAdvancedOptionsView(self, self.case, default)
        if dialog.exec_():
            # Set the model with the dialog results
            result = dialog.get_result()
            log.debug("slotAdvanced -> %s" % str(result))
            self.__model.setDisplacementPredictionAlpha(result[displacement_prediction_alpha])
            self.__model.setDisplacementPredictionBeta(result[displacement_prediction_beta])
            self.__model.setStressPredictionAlpha(result[stress_prediction_alpha])
            self.__model.setMonitorPointSynchronisation(result[monitor_point_synchronisation])


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
