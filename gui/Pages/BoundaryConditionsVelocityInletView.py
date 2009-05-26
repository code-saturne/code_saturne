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

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.BoundaryConditionsVelocityInletForm import Ui_BoundaryConditionsVelocityInletForm

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary

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
        Setup teh widget
        """
        self.__case = case
        self.__boundary = None

        # Connections
        self.connect(self.comboBoxVelocity, SIGNAL("activated(const QString&)"), self.__slotChoiceVelocity)
        self.connect(self.lineEditVelocity, SIGNAL("textChanged(const QString &)"), self.__slotVelocityValue)

        self.connect(self.lineEditXVelocity, SIGNAL("textChanged(const QString &)"), self.__slotDirX)
        self.connect(self.lineEditYVelocity, SIGNAL("textChanged(const QString &)"), self.__slotDirY)
        self.connect(self.lineEditZVelocity, SIGNAL("textChanged(const QString &)"), self.__slotDirZ)

        # Combo models
        self.modelVelocity = ComboModel(self.comboBoxVelocity, 6, 1)
        self.modelVelocity.addItem(self.tr("Velocity"), 'norm')
        self.modelVelocity.addItem(self.tr("Mass flow rate"), 'flow1')
        self.modelVelocity.addItem(self.tr("Volumic flow rate"), 'flow2')
        self.modelVelocity.addItem(self.tr("Velocity and direction"), 'norm+direction')
        self.modelVelocity.addItem(self.tr("Mass flow rate and direction"), 'flow1+direction')
        self.modelVelocity.addItem(self.tr("Volumic flow rate and direction"), 'flow2+direction')

        # Validators
        validatorVelocity = DoubleValidator(self.lineEditVelocity)
        validatorX = DoubleValidator(self.lineEditXVelocity)
        validatorY = DoubleValidator(self.lineEditYVelocity)
        validatorZ = DoubleValidator(self.lineEditZVelocity)
        
        # Apply validators
        self.lineEditVelocity.setValidator(validatorVelocity)
        self.lineEditXVelocity.setValidator(validatorX)
        self.lineEditYVelocity.setValidator(validatorY)
        self.lineEditZVelocity.setValidator(validatorZ)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        choice = self.__boundary.getVelocityChoice()
        self.modelVelocity.setItem(str_model=choice)
        text = self.modelVelocity.dicoM2V[choice]
        self.__slotChoiceVelocity(QString(text))
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

        # update the value associated to the velocity boundary type choice
        v = self.__boundary.getVelocity()
        self.lineEditVelocity.setText(QString(str(v)))

        cc = string.split(c, '+')[0]
        if cc  == 'norm':
            self.labelUnitVelocity.setText(QString(str('m/s')))
        elif cc == 'flow1':
            self.labelUnitVelocity.setText(QString(str('kg/s')))
        elif cc == 'flow2':
            self.labelUnitVelocity.setText(QString(str('m<sup>3</sup>/s')))

        if c[-9:] == 'direction':
            x = self.__boundary.getDirection('direction_x')
            y = self.__boundary.getDirection('direction_y')
            z = self.__boundary.getDirection('direction_z')
            self.lineEditXVelocity.setText(QString(str(x)))
            self.lineEditYVelocity.setText(QString(str(y)))
            self.lineEditZVelocity.setText(QString(str(z)))
            self.frameVelocity.show()
        else:
            self.__boundary.deleteDirectionNodes()
            self.frameVelocity.hide()


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


    def tr(self, text):
        """
        Translation
        """
        return text


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
