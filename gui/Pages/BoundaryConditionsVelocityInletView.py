# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module contains the following classes:
- BoundaryConditionsVelocityInletView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

try:
    import mei
    _have_mei = True
except:
    _have_mei = False

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.BoundaryConditionsVelocityInletForm import Ui_BoundaryConditionsVelocityInletForm

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel, setGreenColor
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary
if _have_mei:
    from QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsVelocityInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsVelocityInletView(QWidget, Ui_BoundaryConditionsVelocityInletForm):
    """
    Boundary condition for velocity in inlet, without particular physics.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsVelocityInletForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        # Connections
        self.connect(self.comboBoxVelocity, SIGNAL("activated(const QString&)"), self.__slotChoiceVelocity)
        self.connect(self.lineEditVelocity, SIGNAL("textChanged(const QString &)"), self.__slotVelocityValue)

        self.connect(self.comboBoxDirection, SIGNAL("activated(const QString&)"), self.__slotChoiceDirection)
        self.connect(self.lineEditDirectionX, SIGNAL("textChanged(const QString &)"), self.__slotDirX)
        self.connect(self.lineEditDirectionY, SIGNAL("textChanged(const QString &)"), self.__slotDirY)
        self.connect(self.lineEditDirectionZ, SIGNAL("textChanged(const QString &)"), self.__slotDirZ)

        # Combo models
        self.modelVelocity = ComboModel(self.comboBoxVelocity, 6, 1)
        self.modelVelocity.addItem(self.tr("norm"), 'norm')
        self.modelVelocity.addItem(self.tr("mass flow rate"), 'flow1')
        self.modelVelocity.addItem(self.tr("volumic flow rate"), 'flow2')
        self.modelVelocity.addItem(self.tr("norm (user law)"), 'norm_formula')
        self.modelVelocity.addItem(self.tr("mass flow rate (user law)"), 'flow1_formula')
        self.modelVelocity.addItem(self.tr("volumic flow rate (user law)"), 'flow2_formula')

        self.modelDirection = ComboModel(self.comboBoxDirection, 3, 1)
        self.modelDirection.addItem(self.tr("normal direction to the inlet"), 'normal')
        self.modelDirection.addItem(self.tr("specified coordinates"), 'coordinates')
        self.modelDirection.addItem(self.tr("user profile"), 'formula')

        # Validators
        validatorVelocity = DoubleValidator(self.lineEditVelocity)
        validatorX = DoubleValidator(self.lineEditDirectionX)
        validatorY = DoubleValidator(self.lineEditDirectionY)
        validatorZ = DoubleValidator(self.lineEditDirectionZ)

        # Apply validators
        self.lineEditVelocity.setValidator(validatorVelocity)
        self.lineEditDirectionX.setValidator(validatorX)
        self.lineEditDirectionY.setValidator(validatorY)
        self.lineEditDirectionZ.setValidator(validatorZ)

        if _have_mei:
            self.connect(self.pushButtonVelocityFormula, SIGNAL("clicked()"), self.__slotVelocityFormula)
            self.connect(self.pushButtonDirectionFormula, SIGNAL("clicked()"), self.__slotDirectionFormula)
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.pushButtonDirectionFormula.setEnabled(False)
            self.modelVelocity.disableItem(str_model="norm_formula")
            self.modelVelocity.disableItem(str_model="flow1_formula")
            self.modelVelocity.disableItem(str_model="flow2_formula")
            self.modelDirection.disableItem(str_model="formula")


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        # Initialize velocity
        choice = self.__boundary.getVelocityChoice()
        self.modelVelocity.setItem(str_model=choice)
        self.__updateLabel()

        if not _have_mei:
            if self.__boundary.getVelocityChoice()[-7:] == "formula":
                c = self.__boundary.defaultValues()['velocityChoice']
                self.__boundary.setVelocityChoice(c)
            if self.__boundary.getDirectionChoice() == "formula":
                c = self.__boundary.defaultValues()['directionChoice']
                self.__boundary.setDirectionChoice(c)

        if choice[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            self.lineEditVelocity.setEnabled(False)
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(QString(str(v)))

        # Initialize direction
        choice = self.__boundary.getDirectionChoice()
        self.modelDirection.setItem(str_model=choice)
        text = self.modelDirection.dicoM2V[choice]
        if choice == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            self.frameDirectionCoordinates.hide()
        elif choice == "coordinates":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.frameDirectionCoordinates.show()
            v = self.__boundary.getDirection('direction_x')
            self.lineEditDirectionX.setText(QString(str(v)))
            v = self.__boundary.getDirection('direction_y')
            self.lineEditDirectionY.setText(QString(str(v)))
            v = self.__boundary.getDirection('direction_z')
            self.lineEditDirectionZ.setText(QString(str(v)))
        elif choice == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.frameDirectionCoordinates.hide()

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSignature("const QString&")
    def __slotChoiceVelocity(self, text):
        """
        Private slot.

        Input the velocity boundary type choice (norm, ).

        @type text: C{QString}
        @param text: velocity boundary type choice.
        """
        c = self.modelVelocity.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setVelocityChoice(c)

        if c[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            setGreenColor(self.pushButtonVelocityFormula, True)
            self.lineEditVelocity.setEnabled(False)
            self.lineEditVelocity.setText(QString(""))
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            setGreenColor(self.pushButtonVelocityFormula, False)
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(QString(str(v)))

        self.__updateLabel()


    def __updateLabel(self):
        """
        Update the unit for the velocity specification.
        """
        c = self.__boundary.getVelocityChoice()
        if c in ('norm', 'norm_formula'):
            self.labelUnitVelocity.setText(QString(str('m/s')))
        elif c in ('flow1', 'flow1_formula'):
            self.labelUnitVelocity.setText(QString(str('kg/s')))
        elif c in ('flow2', 'flow2_formula'):
            self.labelUnitVelocity.setText(QString(str('m<sup>3</sup>/s')))


    @pyqtSignature("const QString&")
    def __slotVelocityValue(self, text):
        """
        Private slot.

        New value associated to the velocity boundary type.

        @type text: C{QString}
        @param text: value
        """
        v, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setVelocity(v)


    @pyqtSignature("")
    def __slotVelocityFormula(self):
        """
        """
        exp = self.__boundary.getVelocity()
        c = self.__boundary.getVelocityChoice()
        req = [('u_norm', 'Norm of the velocity')]
        if c == 'norm_formula':
            exa = "u_norm = 1.0;"
        elif c == 'flow1_formula':
            exa = "q_m = 1.0;"
        elif c == 'flow2_formula':
            exa = "q_v = 1.0;"

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        dialog = QMeiEditorView(self, expression = exp,
                                      required   = req,
                                      symbols    = sym,
                                      examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.__boundary.setVelocity(result)
            setGreenColor(self.pushButtonVelocityFormula, False)


    @pyqtSignature("const QString&")
    def __slotChoiceDirection(self, text):
        """
        Input the direction type choice.
        """
        c = self.modelDirection.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setDirectionChoice(c)

        if c == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            setGreenColor(self.pushButtonDirectionFormula, True)
            self.frameDirectionCoordinates.hide()
        elif c == "coordinates":
            self.pushButtonDirectionFormula.setEnabled(False)
            setGreenColor(self.pushButtonDirectionFormula, False)
            self.frameDirectionCoordinates.show()
            v = self.__boundary.getDirection('direction_x')
            self.lineEditDirectionX.setText(QString(str(v)))
            v = self.__boundary.getDirection('direction_y')
            self.lineEditDirectionY.setText(QString(str(v)))
            v = self.__boundary.getDirection('direction_z')
            self.lineEditDirectionZ.setText(QString(str(v)))
        elif c == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            setGreenColor(self.pushButtonDirectionFormula, False)
            self.frameDirectionCoordinates.hide()


    @pyqtSignature("const QString&")
    def __slotDirX(self, text):
        """
        INPUT value into direction of inlet flow
        """
        value, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setDirection('direction_x', value)


    @pyqtSignature("const QString&")
    def __slotDirY(self, text):
        """
        INPUT value into direction of inlet flow
        """
        value, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setDirection('direction_y', value)


    @pyqtSignature("const QString&")
    def __slotDirZ(self, text):
        """
        INPUT value into direction of inlet flow
        """
        value, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setDirection('direction_z', value)


    @pyqtSignature("")
    def __slotDirectionFormula(self):
        """
        """
        exp = self.__boundary.getDirection('direction_formula')

        req = [('dir_x', 'Direction of the flow along X'),
               ('dir_y', 'Direction of the flow along Y'),
               ('dir_z', 'Direction of the flow along Z')]

        exa = "dir_x = 3.0;\ndir_y = 1.0;\ndir_z = 0.0;\n"

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        dialog = QMeiEditorView(self,expression = exp,
                                     required   = req,
                                     symbols    = sym,
                                     examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDirection -> %s" % str(result))
            self.__boundary.setDirection('direction_formula', result)
            setGreenColor(self.pushButtonDirectionFormula, False)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
