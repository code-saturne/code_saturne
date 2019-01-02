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
This module defines the values of reference.

This module contains the following classes and function:
- FluidStructureInteractionAdvancedOptionsView
- StandardItemModel
- Coupling
- LineEditCoupling
- FormulaCoupling
- CheckBoxCoupling
- CouplingManager
- FluidStructureInteractionView
"""

#-------------------------------------------------------------------------------
# Library modules import
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

from code_saturne.Base.QtPage                          import DoubleValidator, IntValidator
from code_saturne.Base.QtPage                          import to_qvariant, from_qvariant
from code_saturne.Pages.FluidStructureInteractionForm  import Ui_FluidStructureInteractionForm
from code_saturne.Pages.FluidStructureInteractionModel import FluidStructureInteractionModel
from code_saturne.Pages.LocalizationModel              import LocalizationModel
from code_saturne.Pages.Boundary                       import Boundary
from code_saturne.Pages.FluidStructureInteractionAdvancedOptionsDialogForm import \
Ui_FluidStructureInteractionAdvancedOptionsDialogForm

from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.NotebookModel import NotebookModel

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


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModel(QStandardItemModel):
    """
    StandardItemModel for table view
    """

    def __init__(self):
        """
        StandarItemModel for fluid structure interaction tableView
        """
        QStandardItemModel.__init__(self)

        # Define header
        self.headers = [self.tr("Structure number"),
                        self.tr("Label"),
                        self.tr("Location")]
        self.setColumnCount(len(self.headers))

        # Set attributes
        self.__data    = []


    def data(self, index, role):
        """
        Called when table need to read the data
        """
        # return value only for Qt.DisplayRole
        if index.isValid() and role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            return to_qvariant(self.__data[row][col])

        return to_qvariant()


    def flags(self, index):
        """
        Define which column is editable/selectable/enable.
        """
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        """
        Return the header column data.
        """
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        """
        Set a the data when table changed
        """
        raise Exception('Cannot edit column')


    def addItem(self, zone):
        """
        Add an element in the table view.
        """
        index = len(self.__data)
        line = [index + 1, zone.getLabel(), zone.getLocalization()]

        # Add the line and increment line number
        self.__data.append(line)
        row = self.rowCount()
        self.setRowCount(row + 1)


    def getLabel(self, index):
        """
        return the label
        """
        row = index.row()
        [index, label, localization] = self.__data[row]
        return label


#-------------------------------------------------------------------------------
# Coupling base class
#-------------------------------------------------------------------------------

class Coupling:
    """
    Coupling is the base class to manage all widgets whose value depend on the
    boundary.

    It provides getBoundaryDefinedValue/setBoundaryDefinedValue methods which
    gets/sets the value of a boundary attribute. Getter and setter are specified
    in the constructor

    It also automatically enables/disables the widget when boundary is present
    or not

    Derived classes can override onBoundarySet method to specify the initial
    value of the widget base on the boundary
    """

    def __init__(self, widget, getterStr, setterStr ):
        """
        Constructor. getterStr and setterStr are string
        """
        self.__widget    = widget
        self.__getterStr = getterStr
        self.__setterStr = setterStr
        self.__boundary  = None

        # As no boundary is selected disable the widget
        widget.setEnabled(False)


    def getWidget(self):
        """
        Return the widget
        """
        return self.__widget


    def setBoundary(self, boundary):
        """
        Set the current boundary
        """
        self.__boundary = boundary

        #Enable widget
        self.__widget.setEnabled(True)

        # call onBoundarySet for derived class
        self.onBoundarySet()


    def onBoundarySet(self):
        """
        Called when boundary is set. Nothing by default
        """
        pass


    def getBoundaryDefinedValue(self):
        """
        Return the value of the boundary using the getter function
        """
        return getattr(self.__boundary, self.__getterStr)()


    def setBoundaryDefinedValue(self, value):
        """
        Set the value of the boundary using the setter function
        """
        getattr(self.__boundary, self.__setterStr)(value)


#-------------------------------------------------------------------------------
# LineEdit Coupling class
#-------------------------------------------------------------------------------

class LineEditCoupling(Coupling):
    """
    LineEdit that depends on a boundary
    """

    def __init__(self, lineEdit, getter, setter):
        """
        Constructor
        """
        Coupling.__init__(self, lineEdit, getter, setter)

        # Add validator.
        validator = DoubleValidator(lineEdit)
        lineEdit.setValidator(validator)
        lineEdit.textChanged[str].connect(self.__slotTextChanged)


    def onBoundarySet(self):
        """
        Called when boundary is set. Update lineEdit text
        """
        value  = self.getBoundaryDefinedValue()
        self.getWidget().setText(str(value))


    # NOTE: using a decorated slot to connect to a signal is usually recommended,
    # as it is slightly faster and uses less memory, but is causes a crash
    # with PyQt5 (not PyQt4) when connecting to a signal from another class.

    # @pyqtSlot(str)
    def __slotTextChanged(self, text):
        """
        Update the model
        """
        self.setBoundaryDefinedValue(text)

#-------------------------------------------------------------------------------
# Formula Coupling class
#-------------------------------------------------------------------------------

class FormulaCoupling(Coupling):
    """
    Formula button that depend on a boundary
    """

    def __init__(self, button, getter, setter, default, required, symbols, examples):
        """
        Constructor
        """
        Coupling.__init__(self, button, getter, setter)

        self.__default  = default
        self.__required = required
        self.__examples = examples
        self.__symbols  = symbols

        button.clicked.connect(self.__slotFormula)


    def onBoundarySet(self):
        """
        Called when boundary is set.
        """
        # call getter to create default value if needed
        self.getBoundaryDefinedValue()


    # NOTE: as above, do not use decorator to avoid crash in PyQt5.

    # @pyqtSlot(bool)
    def __slotFormula(self, checked):
        """
        Run formula editor.
        """
        # Read current expression
        exp = self.getBoundaryDefinedValue()

        if not exp:
            exp = self.__default

        # run the editor
        dialog = QMeiEditorView(self.getWidget(),
                                check_syntax = None,
                                expression = exp,
                                required   = self.__required,
                                symbols    = self.__symbols,
                                examples   = self.__examples)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("FormulaCoupling -> %s" % str(result))
            self.setBoundaryDefinedValue(result)

#-------------------------------------------------------------------------------
# CheckBoxCouplings Coupling class
#-------------------------------------------------------------------------------

class CheckBoxCoupling(Coupling):
    """
    CheckBox that depend on a boundary
    """

    def __init__(self, checkBox, getter, setter):
        """
        Constructor
        """
        Coupling.__init__(self, checkBox, getter, setter)
        checkBox.stateChanged[int].connect(self.__slotStateChanged)


    def onBoundarySet(self):
        """
        Called when boundary is set.
        """
        value = self.getBoundaryDefinedValue()

        if value == "on":
            state = Qt.Checked
        else:
            state = Qt.Unchecked
        self.getWidget().setCheckState(state)

    # NOTE: as above, do not use decorator to avoid crash in PyQt5.

    # @pyqtSlot(int)
    def __slotStateChanged(self, state):
        """
        Called when checkbox state changed
        """
        value = "off"
        if state == Qt.Checked:
            value = "on"
        self.setBoundaryDefinedValue(value)

#-------------------------------------------------------------------------------
# CouplingManager class
#-------------------------------------------------------------------------------

class CouplingManager:
    """
    Manage and initialize coupling derived objects
    """

    def __init__(self, mainView, case,
                 internalTableView, internalTableModel,
                 externalTableView, externalTableModel):
        """
        Constructor
        """
        self.case               = case
        self.case.undoStopGlobal()
        self.__internalTableView = internalTableView
        self.__externalTableView = externalTableView
        self.__internalTableModel = internalTableModel
        self.__externalTableModel = externalTableModel
        self.__internalCouplings = []
        self.__externalCouplings = []

        # Init widgets
        self.__initLineEditCouplings(mainView)
        self.__initFormulaCouplings (mainView)
        self.__initCheckBoxCouplings(mainView)
        self.case.undoStartGlobal()


    def __initLineEditCouplings(self, mainView):
        """
        Initialize the creation of LineEditCoupling
        """
        couplings = []
        couplings.append(LineEditCoupling(mainView.lineEditInitialDisplacementX,
                                          "getInitialDisplacementX",
                                          "setInitialDisplacementX"))
        couplings.append(LineEditCoupling(mainView.lineEditInitialDisplacementY,
                                          "getInitialDisplacementY",
                                          "setInitialDisplacementY"))
        couplings.append(LineEditCoupling(mainView.lineEditInitialDisplacementZ,
                                          "getInitialDisplacementZ",
                                          "setInitialDisplacementZ"))

        couplings.append(LineEditCoupling(mainView.lineEditEquilibriumDisplacementX,
                                          "getEquilibriumDisplacementX",
                                          "setEquilibriumDisplacementX"))

        couplings.append(LineEditCoupling(mainView.lineEditEquilibriumDisplacementY,
                                          "getEquilibriumDisplacementY",
                                          "setEquilibriumDisplacementY"))
        couplings.append(LineEditCoupling(mainView.lineEditEquilibriumDisplacementZ,
                                          "getEquilibriumDisplacementZ",
                                          "setEquilibriumDisplacementZ"))

        couplings.append(LineEditCoupling(mainView.lineEditInitialVelocityX,
                                          "getInitialVelocityX",
                                          "setInitialVelocityX"))
        couplings.append(LineEditCoupling(mainView.lineEditInitialVelocityY,
                                          "getInitialVelocityY",
                                          "setInitialVelocityY"))
        couplings.append(LineEditCoupling(mainView.lineEditInitialVelocityZ,
                                          "getInitialVelocityZ",
                                          "setInitialVelocityZ"))
        self.__internalCouplings.extend(couplings)


    def __initCheckBoxCouplings(self, mainView):
        """
        Initialize the creation of the checkbox coupling
        """
        couplings = []
        couplings.append(CheckBoxCoupling(mainView.checkBoxDDLX,
                                          "getDDLX", "setDDLX"))
        couplings.append(CheckBoxCoupling(mainView.checkBoxDDLY,
                                          "getDDLY", "setDDLY"))
        couplings.append(CheckBoxCoupling(mainView.checkBoxDDLZ,
                                          "getDDLZ", "setDDLZ"))
        self.__externalCouplings.extend(couplings)


    def __initFormulaCouplings(self, mainView):
        """
        Initialize the creation of the formula button
        """
        default = "%(t)s11 = ;"
        defaultRequired = [('%(t)s11', '%(n)s matrix of the structure (1,1)'),
                           ('%(t)s22', '%(n)s matrix of the structure (2,2)'),
                           ('%(t)s33', '%(n)s matrix of the structure (3,3)'),
                           ('%(t)s12', '%(n)s matrix of the structure (1,2)'),
                           ('%(t)s13', '%(n)s matrix of the structure (1,3)'),
                           ('%(t)s23', '%(n)s matrix of the structure (2,3)'),
                           ('%(t)s21', '%(n)s matrix of the structure (2,1)'),
                           ('%(t)s31', '%(n)s matrix of the structure (3,1)'),
                           ('%(t)s32', '%(n)s matrix of the structure (3,2)')]
        symbols = [('dt', 'time step'),
                   ('t', 'current time'),
                   ('nbIter', 'number of iteration')]

        # Add notebook symbols
        self.notebook = NotebookModel(self.case)
        for (nme, val) in self.notebook.getNotebookList():
            symbols.append((nme, 'value (notebook) = ' + str(val)))

        m_default = default % {'t':'m'}
        c_default = default % {'t':'c'}
        k_default = default % {'t':'k'}

        m_default_required = []
        c_default_required = []
        k_default_required = []
        for v, s in defaultRequired:
            m_default_required.append((v % {'t':'m'}, s % {'n':'mass'}))
            c_default_required.append((v % {'t':'c'}, s % {'n':'damping'}))
            k_default_required.append((v % {'t':'k'}, s % {'n':'stiffness'}))

        m_examples = """# Mass of the structure: 5 kg
#
m11 = 5;\nm22 = 5;\nm33 = 5;\nm12 = 0;\nm13 = 0;\nm23 = 0;\nm21 = 0;\nm31 = 0;\nm32 = 0;
"""
        c_examples = """# Damping of the structure: 3 kg.s
#
c11 = 3;\nc22 = 3;\nc33 = 3;\nc12 = 0;\nc13 = 0;\nc23 = 0;\nc21 = 0;\nc31 = 0;\nc32 = 0;
"""
        k_examples = """# Stiffness of the structure: 2 N/m
#
k11 = 2;\nk22 = 2;\nk33 = 2;\nk12 = 0;\nk13 = 0;\nk23 = 0;\nk21 = 0;\nk31 = 0;\nk32 = 0;
"""

        couplings = []
        couplings.append(FormulaCoupling(mainView.pushButtonMassMatrix,
                                         "getMassMatrix", "setMassMatrix",
                                         m_default, m_default_required,
                                         symbols, m_examples))

        couplings.append(FormulaCoupling(mainView.pushButtonDampingMatrix,
                                         "getDampingMatrix", "setDampingMatrix",
                                         c_default, c_default_required,
                                         symbols, c_examples))

        couplings.append(FormulaCoupling(mainView.pushButtonStiffnessMatrix,
                                         "getStiffnessMatrix", "setStiffnessMatrix",
                                         k_default, k_default_required,
                                         symbols, k_examples))

        defaultFluidForce  = "fx = "
        requiredFluidForce = [('fx', 'force applied to the structure along X'),
                              ('fy', 'force applied to the structure along Y'),
                              ('fz', 'force applied to the structure along Z')]
        symbolsFluidForce = symbols[:];
        symbolsFluidForce.append(('fluid_fx', 'force of flow along X'))
        symbolsFluidForce.append(('fluid_fy', 'force of flow along Y'))
        symbolsFluidForce.append(('fluid_fz', 'force of flow along Z'))

        examplesFluidForce = """# The fluid force is zero in the Y direction.
#
fx = fluid_fx;\nfy = 0;\nfz = fluid_fz;"""
        couplings.append( FormulaCoupling(mainView.pushButtonFluidForce,
                                          "getFluidForceMatrix",
                                          "setFluidForceMatrix",
                                          defaultFluidForce,
                                          requiredFluidForce,
                                          symbolsFluidForce,
                                          examplesFluidForce))
        self.__internalCouplings.extend(couplings)


    # NOTE: as above, do not use decorator to avoid crash in PyQt5.

    # @pyqtSlot(QItemSelection, QItemSelection)
    def slotInternalSelectionChanged(self, selected, deselected):
        """
        Called when internal tableView selection changed
        """
        self.__selectionChanged(self.__internalTableView,
                                self.__internalTableModel,
                                self.__internalCouplings, selected)

    # NOTE: as above, do not use decorator to avoid crash in PyQt5.

    # @pyqtSlot(QItemSelection, QItemSelection)
    def slotExternalSelectionChanged(self, selected, deselected):
        """
        Called when external tableView selection changed
        """
        self.__selectionChanged(self.__externalTableView,
                                self.__externalTableModel,
                                self.__externalCouplings, selected)


    def __selectionChanged(self, tableView, tableModel, couplings, selected):
        """
        Called when a tableView selection changed
        """
        # Get Boundary
        index = tableView.currentIndex()
        label = tableModel.getLabel(index)

        boundary = Boundary("coupling_mobile_boundary", label, self.case)

        # Set boundary for coupling
        for coupling in couplings:
            coupling.setBoundary(boundary)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class FluidStructureInteractionView(QWidget, Ui_FluidStructureInteractionForm):
    """
    Main class.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        # Init base classes
        QWidget.__init__(self, parent)

        Ui_FluidStructureInteractionForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.__model = FluidStructureInteractionModel(case)

        self.__defineConnection()
        self.__addValidators()
        self.__setInitialValues()

        # Use localization model for column 0, 1, 3
        modelLocalization  = LocalizationModel("BoundaryZone", case)

        # Store modelLocalization as attribut to avoid garbage collector to clean it
        self.__modelLocalization = modelLocalization

        # Initialize the internal and external TableViewItemModel
        self.__internalTableModel = self.__createTableViewItemModel(modelLocalization,
                                                                    'internal_coupling')
        self.__externalTableModel = self.__createTableViewItemModel(modelLocalization,
                                                                    'external_coupling')

        # Coupling Manager
        couplingManager = CouplingManager(self, case,
                                          self.tableInternalCoupling,
                                          self.__internalTableModel,
                                          self.tableExternalCoupling,
                                          self.__externalTableModel)
        # Avoid garbage collector to delete couplingManager
        self.__couplingManager = couplingManager

        # Initialize internal / external table view
        self.__initTableView(self.tableInternalCoupling,
                             self.__internalTableModel,
                             couplingManager.slotInternalSelectionChanged)

        self.__initTableView(self.tableExternalCoupling,
                             self.__externalTableModel,
                             couplingManager.slotExternalSelectionChanged)
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


    def __createTableViewItemModel(self, modelLocalization, filterALE):
        """
        Create the table view item model
        """
        tableViewItemModel = StandardItemModel()

        # Populate QTableView model
        for zone in modelLocalization.getZones():
            boundary = Boundary("mobile_boundary", zone.getLabel(), self.case)
            if boundary.getALEChoice() == filterALE:
                tableViewItemModel.addItem(zone)
        return tableViewItemModel


    def __initTableView(self, tableView, tableViewItemModel, slotSelectionChanged):
        """
        Initialize the main table view
        """
        # Set the model
        tableView.setModel(tableViewItemModel)

        # set the column size
        if QT_API == "PYQT4":
            tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            tableView.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            tableView.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        # Connect slot when selection changed
        tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        selectionModel = QItemSelectionModel(tableViewItemModel, tableView)
        tableView.setSelectionModel(selectionModel)
        selectionModel.selectionChanged.connect(slotSelectionChanged)


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


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
