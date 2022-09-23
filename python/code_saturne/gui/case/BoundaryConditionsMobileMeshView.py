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
This module contains the following classes:
- LineEditCoupling
- FormulaCoupling
- CouplingManager
- BoundaryConditionsMobileMeshView
"""

#-------------------------------------------------------------------------------
# Standard modules
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

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel
from code_saturne.gui.base.QtPage import DoubleValidator, IntValidator

from code_saturne.gui.case.BoundaryConditionsMobileMeshForm import Ui_BoundaryConditionsMobileMeshForm
from code_saturne.model.MobileMeshModel import MobileMeshModel
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsMobileMeshView")
log.setLevel(GuiParam.DEBUG)

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

    def __init__(self, widget, getterStr, setterStr):
        """
        Constructor. getterStr and setterStr are string
        """
        self.__widget    = widget
        self.__getterStr = getterStr
        self.__setterStr = setterStr
        self.__boundary  = None


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


    def getBoundaryName(self):
        """
        Return the name of the boundary using the getter function
        """
        return getattr(self.__boundary, "getZoneName")()


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
# Coupling Formula class
#-------------------------------------------------------------------------------

class FormulaCoupling(Coupling):
    """
    Formula button that depend on a boundary
    """

    def __init__(self, button, parent, getter, setter,
                 default, required, symbols, examples):
        """
        Constructor
        """
        super(FormulaCoupling, self).__init__(button, getter, setter)

        self.parent     = parent
        self.__default  = default
        self.__required = required
        self.__examples = examples
        self.__symbols  = symbols

        self.object_type = ''
        if getter == 'getMassMatrix':
            self.object_type = 'mass_matrix'
        elif getter == 'getStiffnessMatrix':
            self.object_type = 'stiffness_matrix'
        elif getter == 'getDampingMatrix':
            self.object_type = 'damping_matrix'
        elif getter == 'getFluidForceMatrix':
            self.object_type = 'fluid_force'

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
        name = str(self.getBoundaryName())
        exp = self.getBoundaryDefinedValue()

        if not exp:
            exp = self.__default

        # run the editor
        dialog = QMegEditorView(self.parent,
                                function_type = 'fsi',
                                zone_name = name,
                                variable_name = self.object_type,
                                expression = exp,
                                required = self.__required,
                                symbols = self.__symbols,
                                examples = self.__examples)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("FormulaCoupling -> %s" % str(result))
            self.setBoundaryDefinedValue(result)

#-------------------------------------------------------------------------------
# CouplingManager class
#-------------------------------------------------------------------------------

class CouplingManager:
    """
    Manage and initialize coupling derived objects
    """

    def __init__(self, mainView, case):
        """
        Constructor
        """
        self.case = case
        self.case.undoStopGlobal()
        self.__internalCouplings = []

        # Init widgets
        self.__initLineEditCouplings(mainView)
        self.__initFormulaCouplings (mainView)
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
                   ('iter', 'current iteration')]

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
                                         mainView,
                                         "getMassMatrix", "setMassMatrix",
                                         m_default, m_default_required,
                                         symbols, m_examples))

        couplings.append(FormulaCoupling(mainView.pushButtonDampingMatrix,
                                         mainView,
                                         "getDampingMatrix", "setDampingMatrix",
                                         c_default, c_default_required,
                                         symbols, c_examples))

        couplings.append(FormulaCoupling(mainView.pushButtonStiffnessMatrix,
                                         mainView,
                                         "getStiffnessMatrix", "setStiffnessMatrix",
                                         k_default, k_default_required,
                                         symbols, k_examples))

        defaultFluidForce  = "fx = 0;"
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
        couplings.append(FormulaCoupling(mainView.pushButtonFluidForce,
                                         mainView,
                                         "getFluidForceMatrix",
                                         "setFluidForceMatrix",
                                         defaultFluidForce,
                                         requiredFluidForce,
                                         symbolsFluidForce,
                                         examplesFluidForce))
        self.__internalCouplings.extend(couplings)


    def setBoundary(self, boundary):
        """
        Assign a boundary.
        """
        # Set boundary for coupling
        for coupling in self.__internalCouplings:
            coupling.setBoundary(boundary)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsMobileMeshView(QWidget,
                                       Ui_BoundaryConditionsMobileMeshForm):
    """
    Boundary condifition for mobil mesh (ALE and/or Fluid-interaction)
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsMobileMeshForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None

        self.case.undoStopGlobal()

        self.__model = MobileMeshModel(self.case)
        self.notebook = NotebookModel(self.case)

        self.__comboModel = ComboModel(self.comboMobilBoundary, 6, 1)
        self.__comboModel.addItem(self.tr("Fixed boundary"), "fixed_boundary")
        self.__comboModel.addItem(self.tr("Sliding boundary"), "sliding_boundary")
        self.__comboModel.addItem(self.tr("Internal coupling"), "internal_coupling")
        self.__comboModel.addItem(self.tr("External coupling"), "external_coupling")
        self.__comboModel.addItem(self.tr("Fixed velocity"), "fixed_velocity")
        self.__comboModel.addItem(self.tr("Fixed displacement"), "fixed_displacement")
        self.comboMobilBoundary.activated[str].connect(self.__slotCombo)
        self.pushButtonMobilBoundary.clicked.connect(self.__slotFormula)

        # Coupling Manager
        self.__couplingManager = CouplingManager(self, case)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def __slotFormula(self):
        """
        Run formula editor.
        """
        exp = self.__boundary.getALEFormula()
        c = self.__boundary.getALEChoice();

        if c == "fixed_velocity":
            if not exp:
                exp = 'mesh_velocity[0] = 0.;\nmesh_velocity[1] = 0.;\nmesh_velocity[2] = 0.;'
            req = [('mesh_velocity[0]', 'Fixed velocity of the mesh'),
                   ('mesh_velocity[1]', 'Fixed velocity of the mesh'),
                   ('mesh_velocity[2]', 'Fixed velocity of the mesh')]
            exa = 'mesh_velocity[0] = 0.;\nmesh_velocity[1] = 0.;\nmesh_velocity[2] = 1.;'
        elif c == "fixed_displacement":
            if not exp:
                exp = 'mesh_displacement[0] = 0.;\nmesh_displacement[1] = 0.;\nmesh_displacement[2] = 0.;'
            req = [('mesh_displacement[0]', 'Fixed displacement of the mesh'),
                   ('mesh_displacement[1]', 'Fixed displacement of the mesh'),
                   ('mesh_displacement[2]', 'Fixed displacement of the mesh')]
            exa = 'mesh_displacement[0] = 0.;\nmesh_displacement[1] = 0.;\nmesh_displacement[2] = 1.;'

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration'),
               ('surface', 'Boundary zone surface')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent = self,
                                function_type = 'bnd',
                                zone_name     = self.__boundary._label,
                                variable_name = 'mesh_velocity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                condition     = c,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMobileMeshBoundary -> %s" % str(result))
            self.__boundary.setFormula(str(result))
            self.pushButtonMobilBoundary.setStyleSheet("background-color: green")
            self.pushButtonMobilBoundary.setToolTip(result)


    @pyqtSlot(str)
    def __slotCombo(self, text):
        """
        Called when the combobox changed.
        """
        modelData = self.__comboModel.dicoV2M[str(text)]

        if modelData == self.__boundary.getALEChoice():
            return

        self.__boundary.setALEChoice(modelData)
        exp = self.__boundary.getALEFormula()

        # Hide/Show formula button.
        # Formula is always reset when changing values, so set
        # color to red.
        self.update_view(modelData)
        if exp:
            self.pushButtonMobilBoundary.setStyleSheet("background-color: red")
            self.pushButtonMobilBoundary.setToolTip(exp)
        else:
            self.pushButtonMobilBoundary.setStyleSheet("background-color: red")


    def update_view(self, modelData):
        """
        Show the widgets matching the model
        """
        self.__comboModel.setItem(str_model=modelData)
        if modelData in ["fixed_velocity", "fixed_displacement"]:
            self.pushButtonMobilBoundary.show()
        else:
            self.pushButtonMobilBoundary.hide()
        if modelData == "internal_coupling":
            self.groupBoxStructureVelPos.show()
            self.groupBoxStructureCharacteristics.show()
            self.groupBoxForceApplied.show()
        else:
            self.groupBoxStructureVelPos.hide()
            self.groupBoxStructureCharacteristics.hide()
            self.groupBoxForceApplied.hide()


    def showWidget(self, b):
        """
        Show the widget
        """
        if self.__model.getMethod() != "off":
            boundary = Boundary("coupling_mobile_boundary",
                                b.getLabel(), self.case)
            self.__boundary = boundary
            modelData = b.getALEChoice()
            if modelData == 'internal_coupling':
                self.__couplingManager.setBoundary(boundary)
            self.update_view(modelData)
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
