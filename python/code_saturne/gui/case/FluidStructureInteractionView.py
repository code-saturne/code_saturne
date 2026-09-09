# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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

# -------------------------------------------------------------------------------

"""
This module defines the values of reference.

This module contains the following classes and function:
- FluidStructureInteractionAdvancedOptionsView
- StandardItemModel
- Coupling
- FluidStructureInteractionView
"""

# -------------------------------------------------------------------------------
# Library modules import
# -------------------------------------------------------------------------------

import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtPage import DoubleValidator, IntValidator
from code_saturne.gui.base.QtPage import from_qvariant
from code_saturne.gui.case.FluidStructureInteractionForm import (
    Ui_FluidStructureInteractionForm,
)
from code_saturne.model.FluidStructureInteractionModel import (
    FluidStructureInteractionModel,
)
from code_saturne.gui.case.FluidStructureInteractionAdvancedOptionsDialogForm import (
    Ui_FluidStructureInteractionAdvancedOptionsDialogForm,
)

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

if QT_API == "PYQT6":
    from code_saturne.gui.case import resources_pages_rc

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FluidStructureInteractionView")

# -------------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------------

displacement_prediction_alpha = "displacement_prediction_alpha"
displacement_prediction_beta = "displacement_prediction_beta"
stress_prediction_alpha = "stress_prediction_alpha"
structure_time_plot = "monitor_point_synchronisation"

DISPLACEMENT_ACCELERATION_METHODS = {
    0: "none",
    1: "relaxation",
    2: "aitken",
}

DISPLACEMENT_ACCELERATION_INDICES = {
    value: key for key, value in DISPLACEMENT_ACCELERATION_METHODS.items()
}


DISPLACEMENT_PREDICTION_METHODS = {
    0: "none",
    1: "explicit_euler",
    2: "adams_bashforth",
    3: "user",
}

DISPLACEMENT_PREDICTION_INDICES = {
    value: key for key, value in DISPLACEMENT_PREDICTION_METHODS.items()
}


DISPLACEMENT_PREDICTION_COEFFICIENTS = {
    "none": (0.0, 0.0),
    "explicit_euler": (1.0, 0.0),
    "adams_bashforth": (1.0, 0.5),
}

# -------------------------------------------------------------------------------
# Advanced dialog
# -------------------------------------------------------------------------------


class FluidStructureInteractionAdvancedOptionsView(
    QDialog, Ui_FluidStructureInteractionAdvancedOptionsDialogForm
):
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
        self.__result = default.copy()
        self.__setValidator()
        self.__setInitialValues()

        self.case.undoStartGlobal()

    def __setValidator(self):
        """
        Set the validator
        """
        validator = DoubleValidator(self.lineEditDisplacementAlpha, min=0.0)
        self.lineEditDisplacementAlpha.setValidator(validator)

        validator = DoubleValidator(self.lineEditDisplacementBeta, min=0.0)
        self.lineEditDisplacementBeta.setValidator(validator)

        validator = DoubleValidator(self.lineEditStressAlpha, min=0.0)
        self.lineEditStressAlpha.setValidator(validator)

    def __setInitialValues(self):
        """
        Set the initial values for the 4 widgets
        """
        # Read from default
        displacementAlpha = str(self.__default[displacement_prediction_alpha])
        displacementBeta = str(self.__default[displacement_prediction_beta])
        stressAlpha = str(self.__default[stress_prediction_alpha])

        # Update Widget
        self.lineEditDisplacementAlpha.setText(displacementAlpha)
        self.lineEditDisplacementBeta.setText(displacementBeta)
        self.lineEditStressAlpha.setText(stressAlpha)

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
        displacementBeta = float(self.lineEditDisplacementBeta.text())
        stressAlpha = float(self.lineEditStressAlpha.text())

        # Set result attributes
        self.__result[displacement_prediction_alpha] = displacementAlpha
        self.__result[displacement_prediction_beta] = displacementBeta
        self.__result[stress_prediction_alpha] = stressAlpha

        QDialog.accept(self)

    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        QDialog.reject(self)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------


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

        self.checkBoxStructureTimePlot.stateChanged.connect(self.slotStructureTimePlot)

        self.spinBox_ast_log.valueChanged.connect(self.__slotVerbosityCA)
        self.spinBox_ast_viz.valueChanged.connect(self.__slotVisualizationCA)

        self.comboBoxDisplacementAcceleration.currentIndexChanged[int].connect(
            self.__slotDisplacementAccelerationMethodCA
        )

        self.lineEditDisplacementRelaxation.textChanged[str].connect(
            self.__slotDisplacementRelaxationCoefficientCA
        )

        self.comboBoxDisplacementPrediction.currentIndexChanged[int].connect(
            self.__slotDisplacementPredictionMethodCA
        )

        self.lineEditDisplacementPredictionAlpha.textChanged[str].connect(
            self.__slotDisplacementPredictionAlphaCA
        )

        self.lineEditDisplacementPredictionBeta.textChanged[str].connect(
            self.__slotDisplacementPredictionBetaCA
        )

    def __addValidators(self):
        """
        Add the validator for NALIMX and EPALIM
        """
        validatorNALIMX = IntValidator(self.lineEditNALIMX, min=1)
        self.lineEditNALIMX.setValidator(validatorNALIMX)

        validatorEPALIM = DoubleValidator(self.lineEditEPALIM, min=0.0)
        validatorEPALIM.setExclusiveMin(True)
        self.lineEditEPALIM.setValidator(validatorEPALIM)

        # Fixed displacement relaxation coefficient.
        validatorRelaxation = DoubleValidator(
            self.lineEditDisplacementRelaxation, min=0.0, max=1.0
        )
        self.lineEditDisplacementRelaxation.setValidator(validatorRelaxation)

        # User-defined displacement prediction coefficients.
        validatorAlpha = DoubleValidator(self.lineEditDisplacementPredictionAlpha)
        self.lineEditDisplacementPredictionAlpha.setValidator(validatorAlpha)

        validatorBeta = DoubleValidator(self.lineEditDisplacementPredictionBeta)
        self.lineEditDisplacementPredictionBeta.setValidator(validatorBeta)

    def __setInitialValues(self):
        """
        Initialise all widgets from the model without triggering
        user-editing slots.
        """

        # ----------------------------------------------------------
        # General coupling parameters
        # ----------------------------------------------------------

        nalimx = self.__model.getMaxIterations()

        self.lineEditNALIMX.blockSignals(True)
        self.lineEditNALIMX.setText(str(nalimx))
        self.lineEditNALIMX.blockSignals(False)

        epalim = self.__model.getPrecision()

        self.lineEditEPALIM.blockSignals(True)
        self.lineEditEPALIM.setText(str(epalim))
        self.lineEditEPALIM.blockSignals(False)

        structure_time_plot = self.__model.getInternalStructuresTimePlot()

        self.checkBoxStructureTimePlot.blockSignals(True)
        self.checkBoxStructureTimePlot.setChecked(structure_time_plot == "on")
        self.checkBoxStructureTimePlot.blockSignals(False)

        ast_log = self.__model.getVerbosityCA()

        self.spinBox_ast_log.blockSignals(True)
        self.spinBox_ast_log.setValue(ast_log)
        self.spinBox_ast_log.blockSignals(False)

        ast_vis = self.__model.getVisualizationCA()

        self.spinBox_ast_viz.blockSignals(True)
        self.spinBox_ast_viz.setValue(ast_vis)
        self.spinBox_ast_viz.blockSignals(False)

        # ----------------------------------------------------------
        # Displacement acceleration
        # ----------------------------------------------------------

        acceleration_method = self.__model.getDisplacementAccelerationMethodCA()

        acceleration_index = DISPLACEMENT_ACCELERATION_INDICES.get(
            acceleration_method, 0
        )

        self.comboBoxDisplacementAcceleration.blockSignals(True)

        try:
            self.comboBoxDisplacementAcceleration.setCurrentIndex(acceleration_index)
        finally:
            self.comboBoxDisplacementAcceleration.blockSignals(False)

        if acceleration_method == "relaxation":
            relaxation = self.__model.getDisplacementRelaxationCoefficientCA()
        else:
            relaxation = 1.0

        self.__model.setDisplacementRelaxationCoefficientCA(relaxation)

        self.lineEditDisplacementRelaxation.blockSignals(True)

        try:
            self.lineEditDisplacementRelaxation.setText(str(float(relaxation)))
        finally:
            self.lineEditDisplacementRelaxation.blockSignals(False)

        self.__updateDisplacementAccelerationWidgetsCA(acceleration_method)

        # ----------------------------------------------------------
        # Displacement prediction
        # ----------------------------------------------------------

        prediction_method = self.__model.getDisplacementPredictionMethodCA()

        prediction_index = DISPLACEMENT_PREDICTION_INDICES.get(prediction_method, 0)

        self.comboBoxDisplacementPrediction.blockSignals(True)

        try:
            self.comboBoxDisplacementPrediction.setCurrentIndex(prediction_index)
        finally:
            self.comboBoxDisplacementPrediction.blockSignals(False)

        if prediction_method == "user":
            alpha = self.__model.getDisplacementPredictionAlphaCA()
            beta = self.__model.getDisplacementPredictionBetaCA()

            self.__setPredictionCoefficientValuesCA(alpha, beta)

            # Do not overwrite user-defined values.
            update_values = False

        else:
            # Values will be imposed by the selected method.
            update_values = True

        # Mandatory even if the initial combo-box index is already 0.
        self.__updateDisplacementPredictionWidgetsCA(
            prediction_method, update_values=update_values
        )

    @Slot(str)
    def __slotNalimx(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        if self.sender().validator().state == QValidator.State.Acceptable:
            nalimx = int(text)
            self.__model.setMaxIterations(nalimx)

    @Slot(str)
    def __slotEpalim(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        if self.sender().validator().state == QValidator.State.Acceptable:
            epalim = from_qvariant(text, float)
            self.__model.setPrecision(float(epalim))

    @Slot(int)
    def slotStructureTimePlot(self, val):

        if val == 0:
            self.__model.setInternalStructuresTimePlot("off")
        else:
            self.__model.setInternalStructuresTimePlot("on")

    @Slot(int)
    def __slotVerbosityCA(self, text):
        """
        Set value for code_aster logging
        """
        self.__model.setVerbosityCA(int(text))
        log.debug("__slotVerbosityCA-> %s" % text)

    @Slot(int)
    def __slotVisualizationCA(self, text):
        """
        Set value for code_aster logging
        """
        self.__model.setVisualizationCA(int(text))
        log.debug("__slotVisualizationCA-> %s" % text)

    def __updateDisplacementAccelerationWidgetsCA(self, method):
        """
        Update the availability of the relaxation coefficient.
        """
        editable = method == "relaxation"

        widgets = (
            self.labelDisplacementRelaxation,
            self.lineEditDisplacementRelaxation,
        )

        for widget in widgets:
            widget.setEnabled(editable)

        self.lineEditDisplacementRelaxation.setReadOnly(not editable)

    @Slot(int)
    def __slotDisplacementAccelerationMethodCA(self, index):
        """
        Update the displacement acceleration method.
        """
        method = DISPLACEMENT_ACCELERATION_METHODS.get(index, "none")

        self.__model.setDisplacementAccelerationMethodCA(method)

        if method == "relaxation":
            value = self.__model.getDisplacementRelaxationCoefficientCA()
        else:
            value = 0.0

        self.__model.setDisplacementRelaxationCoefficientCA(value)

        self.lineEditDisplacementRelaxation.blockSignals(True)

        try:
            self.lineEditDisplacementRelaxation.setText(str(float(value)))
        finally:
            self.lineEditDisplacementRelaxation.blockSignals(False)

        self.__updateDisplacementAccelerationWidgetsCA(method)

    @Slot(str)
    def __slotDisplacementRelaxationCoefficientCA(self, text):
        """
        Update the fixed displacement relaxation coefficient.
        """
        if not self.lineEditDisplacementRelaxation.isEnabled():
            return

        validator = self.lineEditDisplacementRelaxation.validator()

        if validator.state == QValidator.State.Acceptable:
            value = from_qvariant(text, float)

            self.__model.setDisplacementRelaxationCoefficientCA(float(value))

    def __setPredictionCoefficientValuesCA(self, alpha, beta):
        """
        Update alpha and beta without triggering their slots.
        """
        alpha_widget = self.lineEditDisplacementPredictionAlpha
        beta_widget = self.lineEditDisplacementPredictionBeta

        alpha_widget.blockSignals(True)
        beta_widget.blockSignals(True)

        try:
            alpha_widget.setText(str(float(alpha)))
            beta_widget.setText(str(float(beta)))
        finally:
            alpha_widget.blockSignals(False)
            beta_widget.blockSignals(False)

    def __updateDisplacementPredictionWidgetsCA(self, method, update_values=True):
        """
        Update displacement prediction coefficients.

        For predefined methods, alpha and beta are imposed and
        both fields are disabled.

        For the User method, both fields are enabled and editable.
        """
        is_user = method == "user"

        if update_values and not is_user:
            alpha, beta = DISPLACEMENT_PREDICTION_COEFFICIENTS.get(method, (0.0, 0.0))

            self.__setPredictionCoefficientValuesCA(alpha, beta)

            self.__model.setDisplacementPredictionAlphaCA(float(alpha))
            self.__model.setDisplacementPredictionBetaCA(float(beta))

        labels = (
            self.labelDisplacementPredictionAlpha,
            self.labelDisplacementPredictionBeta,
        )

        editors = (
            self.lineEditDisplacementPredictionAlpha,
            self.lineEditDisplacementPredictionBeta,
        )

        for label in labels:
            label.setEnabled(is_user)

        for editor in editors:
            editor.setEnabled(is_user)
            editor.setReadOnly(not is_user)

    @Slot(int)
    def __slotDisplacementPredictionMethodCA(self, index):
        """
        Update the displacement prediction method.
        """
        method = DISPLACEMENT_PREDICTION_METHODS.get(index, "none")

        self.__model.setDisplacementPredictionMethodCA(method)

        if method == "user":
            alpha = self.__model.getDisplacementPredictionAlphaCA()
            beta = self.__model.getDisplacementPredictionBetaCA()

            self.__setPredictionCoefficientValuesCA(alpha, beta)

            self.__updateDisplacementPredictionWidgetsCA(method, update_values=False)

        else:
            self.__updateDisplacementPredictionWidgetsCA(method, update_values=True)

    @Slot(str)
    def __slotDisplacementPredictionAlphaCA(self, text):
        """
        Update the user-defined displacement prediction alpha.
        """
        method = self.__model.getDisplacementPredictionMethodCA()

        if method != "user":
            return

        validator = self.lineEditDisplacementPredictionAlpha.validator()

        if validator.state == QValidator.State.Acceptable:
            value = from_qvariant(text, float)

            self.__model.setDisplacementPredictionAlphaCA(float(value))

    @Slot(str)
    def __slotDisplacementPredictionBetaCA(self, text):
        """
        Update the user-defined displacement prediction beta.
        """
        method = self.__model.getDisplacementPredictionMethodCA()

        if method != "user":
            return

        validator = self.lineEditDisplacementPredictionBeta.validator()

        if validator.state == QValidator.State.Acceptable:
            value = from_qvariant(text, float)

            self.__model.setDisplacementPredictionBetaCA(float(value))

    @Slot()
    def __slotAdvanced(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        # Set the default value
        default = {}
        default[displacement_prediction_alpha] = (
            self.__model.getDisplacementPredictionAlpha()
        )
        default[displacement_prediction_beta] = (
            self.__model.getDisplacementPredictionBeta()
        )
        default[stress_prediction_alpha] = self.__model.getStressPredictionAlpha()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        # run the dialog
        dialog = FluidStructureInteractionAdvancedOptionsView(self, self.case, default)
        if dialog.exec():
            # Set the model with the dialog results
            result = dialog.get_result()
            log.debug("slotAdvanced -> %s" % str(result))
            self.__model.setDisplacementPredictionAlpha(
                result[displacement_prediction_alpha]
            )
            self.__model.setDisplacementPredictionBeta(
                result[displacement_prediction_beta]
            )
            self.__model.setStressPredictionAlpha(result[stress_prediction_alpha])


# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
