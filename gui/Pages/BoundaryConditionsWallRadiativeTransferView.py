# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
- StandardItemModelBoundaries
- StandardItemModelScalars
- RadiativeBoundariesView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.BoundaryConditionsWallRadiativeTransferForm import Ui_BoundaryConditionsWallRadiativeTransferForm
from Pages.ThermalRadiationModel import ThermalRadiationModel

from Base.Toolbox import GuiParam
from Base.QtPage import IntValidator, DoubleValidator, ComboModel, setGreenColor
#from Pages.RadiativeBoundariesModel import RadiativeBoundariesModel
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsWallRadiativeTransferView")

#-------------------------------------------------------------------------------
# StandarItemModel class to display scalars properties
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):
    def __init__(self, bdModel):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Wall radiative\ncaracteristics"),
                        self.tr("Value"),
                        self.tr("Unit")]
        self.setColumnCount(len(self.headers))
        self.bdModel = bdModel
        self.liste = self.getListVariablesForCondition()
        self.setRowCount(len(self.liste))
        log.debug("StandardItemModelScalars.__init__  liste = %s " % str(self.liste))

        self.dataScalars = {}
        self.dataScalars["EPSP"]  = bdModel.getEmissivity()
        self.dataScalars["XLAMP"] = bdModel.getThermalConductivity()
        self.dataScalars["EPAP"]  = bdModel.getThickness()
        self.dataScalars["TEXTP"] = bdModel.getExternalTemperatureProfile()
        self.dataScalars["TINTP"] = bdModel.getInternalTemperatureProfile()
        self.dataScalars["FLUX"]  = bdModel.getFlux()


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            if index.column() == 0:
                return QVariant(self.liste[index.row()][1])
            elif index.column() == 1:
                key = self.liste[index.row()][3]
                return QVariant(self.dataScalars[key])
            elif index.column() == 2:
                return QVariant(self.liste[index.row()][2])
        if role == Qt.ToolTipRole:
            kword = self.liste[index.row()][3]
            return QVariant(self.tr("Code_Saturne keyword: " + kword))
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        if index.column() == 1:
            row = index.row()
            key = self.liste[row][3]
            tag = self.liste[row][4]
            val, ok = value.toDouble()
            self.bdModel.setValRay(val, tag)
            self.dataScalars[key] = val
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def deleteAll(self):
        self.dataScalars = []
        self.setRowCount(0)


    def insertItem(self):
        self.dataScalars.append()
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataScalars[row]


    def getListVariablesForCondition(self):
        """
        Get list of variables for condition choosed
        """
        cond = self.bdModel.getRadiativeChoice()

        if cond == 'itpimp':
            liste = [(0, self.tr("Emissivite"), '',  'EPSP',  'emissivity'),
                     (1, self.tr("Initial temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
        if cond == 'ipgrno':
            liste = [(0, self.tr("Emissivity"), '',  'EPSP',  'emissivity'),
                     (1, self.tr("Conductivity"), 'W/m/K', 'XLAMP', 'thermal_conductivity'),
                     (2, self.tr("Thickness"), 'm', 'EPAP' , 'thickness'),
                     (3, self.tr("Profile of external temperature"), 'K', 'TEXTP', 'external_temperature_profile'),
                     (4, self.tr("Profile of internal temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
##        if cond == 'iprefl':
##            list = [(0, self.xlamp,t.XLAMP, 'W/m/K', 'XLAMP'),
##                    (1, self.epap, t.EPAP,  'm', 'EPAP'),
##                    (2, self.textp,t.TEXTP, 'K', 'TEXTP'),
##                    (3, self.tintp,t.TINTP, 'K', 'TINTP')]
##            self.f43 = Tix.Frame(self.f4, relief=FLAT)
##            self.f43.pack(side=TOP, fill=X, pady=10)
##            frad = self.f43
        if cond == 'ifgrno':
            liste = [(0, self.tr("Emissivity"),'', 'EPSP', 'emissivity'),
                     (1, self.tr("Flux of conduction"), 'W/m2', 'FLUX',  'flux'),
                     (2, self.tr("Inital temperature"), 'K', 'TINTP', 'internal_temperature_profile')]
##        if cond == 'ifrefl':
##            list = [(0, self.flux, t.FLUX, 'W/m2', 'FLUX'),
##                    (1, self.tintp, t.TINTP, 'K', 'TINTP')]
##            self.f45 = Tix.Frame(self.f4, relief=FLAT)
##            self.f45.pack(side=TOP, fill=X, pady=10)
##            frad = self.f45
        return liste

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsWallRadiativeTransferView(QWidget,
                                                  Ui_BoundaryConditionsWallRadiativeTransferForm):
    """
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsWallRadiativeTransferForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        # Create the Page layout.

        # Combo
        self.modelRadiative = ComboModel(self.comboBoxRadiative,3,1)
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and profile of fixed internal temperature"), 'itpimp')
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and profile of fixed external temperature"), 'ipgrno')
##         self.modelRadiative.addItem(self.tr("Paroi reflechissante\n"\
##                                               " + profil de temperature externe impose"), 'iprefl')
        self.modelRadiative.addItem(self.tr("Gray or black wall\n"\
                                            " and flux of fixed conduction"), 'ifgrno')
##         self.modelRadiative.addItem(self.tr("Paroi reflechissante\n"\
##                                               " + flux de conduction impose en paroi"), 'ifrefl')

        # Validator
        validatorZone = IntValidator(self.lineEditZone, min=0)
        validatorZone.setExclusiveMin(True)
        self.lineEditZone.setValidator(validatorZone)

        # Connections
        self.connect(self.comboBoxRadiative,
                     SIGNAL("activated(const QString&)"),
                     self.slotRadiativeChoice)
        self.connect(self.lineEditZone,
                     SIGNAL("textChanged(const QString &)"),
                     self.slotZone)


    def showWidget(self, b):
        """
        Show the widget
        """
        if ThermalRadiationModel(self.__case).getRadiativeModel() != "off":
            label = b.getLabel()
            self.__boundary = Boundary('radiative_wall', label, self.__case)
            choice = self.__boundary.getRadiativeChoice()
            self.modelRadiative.setItem(str_model=choice)

            if hasattr(self, "modelScalars"):
                del self.modelScalars
            self.modelScalars = StandardItemModelScalars(self.__boundary)
            self.tableViewScalars.setModel(self.modelScalars)

            self.nb_zone = self.__boundary.getOutputRadiativeZone()
            self.lineEditZone.setText(QString(str(self.nb_zone)))

            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all the widget
        """
        self.hide()


    @pyqtSignature("const QString&")
    def slotZone(self, text):
        nb_zone, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.__boundary.setOutputRadiativeZone(nb_zone)
            return nb_zone


    @pyqtSignature("const QString&")
    def slotRadiativeChoice(self, text):
        cond = self.modelRadiative.dicoV2M[str(text)]
        log.debug("slotRadiativeChoice cond = %s "%cond)
        self.__boundary.setRadiativeChoice(cond)
        self.modelScalars.deleteAll()
        self.modelScalars = StandardItemModelScalars(self.__boundary)
        self.tableViewScalars.setModel(self.modelScalars)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
