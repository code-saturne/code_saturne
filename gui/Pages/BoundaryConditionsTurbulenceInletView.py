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
- BoundaryConditionsTurbulenceInletView
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

from Pages.BoundaryConditionsTurbulenceInletForm import Ui_BoundaryConditionsTurbulenceInletForm
from TurbulenceModel import TurbulenceModel

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsTurbulenceInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsTurbulenceInletView(QWidget, Ui_BoundaryConditionsTurbulenceInletForm):
    """
    Boundary condition for turbulence
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsTurbulenceInletForm.__init__(self)
        self.setupUi(self)

    
    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.connect(self.comboBoxTurbulence, SIGNAL("activated(const QString&)"), self.__slotChoiceTurbulence)

        self.__modelTurbulence = ComboModel(self.comboBoxTurbulence, 2, 1)
        self.__modelTurbulence.addItem(self.tr("Calculation by hydraulic diameter"), 'hydraulic_diameter')
        self.__modelTurbulence.addItem(self.tr("Calculation by turbulent intensity"), 'turbulent_intensity')
        
        self.connect(self.lineEditDiameter, SIGNAL("textChanged(const QString &)"), self.__slotDiam)
        self.connect(self.lineEditIntensity, SIGNAL("textChanged(const QString &)"), self.__slotIntensity)
        self.connect(self.lineEditDiameterIntens, SIGNAL("textChanged(const QString &)"), self.__slotDiam)

        validatorDiam = DoubleValidator(self.lineEditDiameter, min=0.)
        validatorDiam.setExclusiveMin(True)
        validatorIntensity = DoubleValidator(self.lineEditIntensity, min=0.)

        self.lineEditDiameter.setValidator(validatorDiam)
        self.lineEditDiameterIntens.setValidator(validatorDiam)
        self.lineEditIntensity.setValidator(validatorIntensity)


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        if TurbulenceModel(self.__case).getTurbulenceVariable():
            turb_choice = boundary.getTurbulenceChoice()
            self.__modelTurbulence.setItem(str_model=turb_choice)
            if turb_choice == "hydraulic_diameter":
                self.frameTurbDiameter.show()
                self.frameTurbIntensity.hide()
                d = boundary.getHydraulicDiameter()
                self.lineEditDiameter.setText(QString(str(d)))
            elif turb_choice == "turbulent_intensity":
                self.frameTurbIntensity.show()
                self.frameTurbDiameter.hide()
                i = boundary.getTurbulentIntensity()
                d = boundary.getHydraulicDiameter()
                self.lineEditIntensity.setText(QString(str(i)))
                self.lineEditDiameterIntens.setText(QString(str(d)))
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSignature("const QString&")
    def __slotChoiceTurbulence(self, text):
        """
        INPUT choice of method of calculation of the turbulence
        """
        turb_choice = self.__modelTurbulence.dicoV2M[str(text)]
        self.__boundary.setTurbulenceChoice(turb_choice)

        self.frameTurbDiameter.hide()
        self.frameTurbIntensity.hide()

        if turb_choice  == 'hydraulic_diameter':
            self.frameTurbDiameter.show()
            d = self.__boundary.getHydraulicDiameter()
            self.lineEditDiameter.setText(QString(str(d)))
        elif turb_choice == 'turbulent_intensity':
            self.frameTurbIntensity.show()
            i = self.__boundary.getTurbulentIntensity()
            self.lineEditIntensity.setText(QString(str(i)))
            d = self.__boundary.getHydraulicDiameter()
            self.lineEditDiameterIntens.setText(QString(str(d)))


    @pyqtSignature("const QString&")
    def __slotDiam(self, text):
        """
        INPUT hydraulic diameter
        """
        diam, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setHydraulicDiameter(diam)


    @pyqtSignature("const QString&")
    def __slotIntensity(self, text):
        """
        INPUT turbulent intensity
        """
        intens, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setTurbulentIntensity(intens)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
