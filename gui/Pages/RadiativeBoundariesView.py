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

from RadiativeBoundariesForm import Ui_RadiativeBoundariesForm

from Base.Toolbox import GuiParam
import Base.QtPage as QtPage
#from Pages.RadiativeBoundariesModel import RadiativeBoundariesModel
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("RadiativeBoundariesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class to display boundaries in a QTreeView 
#-------------------------------------------------------------------------------

class StandardItemModelBoundaries(QStandardItemModel):

    def __init__(self):
        QStandardItemModel.__init__(self) 
        self.headers = [self.tr("Label"), self.tr("Zone"),
                        self.tr("Nature"), self.tr("Localization")]
        self.setColumnCount(len(self.headers))
        self.dataBoundary = []

    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            return QVariant(self.dataBoundary[row][col])
        return QVariant()

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
    
    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()
    
    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

    def insertItem(self, label, codeNumber, var_nature, local):
        line = [label, codeNumber, var_nature, local]
        self.dataBoundary.append(line) 
        row = self.rowCount()
        self.setRowCount(row+1)

    def getItem(self, row):
        return self.dataBoundary[row]

#-------------------------------------------------------------------------------
# StandarItemModel class to display scalars properties
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):

    def __init__(self, bdModel, liste):
        QStandardItemModel.__init__(self) 
        self.headers = [self.tr("Scalar name"), self.tr("Scalar value"), self.tr("Unit")]
        self.setColumnCount(len(self.headers))
        self.bdModel = bdModel
        self.liste = liste
        self.setRowCount(len(self.liste))
        log.debug("StandardItemModelScalars.__init__  liste = %s "%str(liste))

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

    def insertItem(self):
        self.dataScalars.append() 
        row = self.rowCount()
        self.setRowCount(row+1)

    def getItem(self, row):
        return self.dataScalars[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class RadiativeBoundariesView(QWidget, Ui_RadiativeBoundariesForm):
    """
    """
    
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_RadiativeBoundariesForm.__init__(self)
        self.setupUi(self)
        
        self.case = case

        # Create the Page layout.

        # Model and QTreeView for Boundaries
        self.modelBoundaries = StandardItemModelBoundaries()
        self.treeViewBoundaries.setModel(self.modelBoundaries)
        self.treeViewBoundaries.setAlternatingRowColors(True)
        self.treeViewBoundaries.setSelectionBehavior(QAbstractItemView.SelectRows)

        # Combo
        self.modelRadiative = QtPage.ComboModel(self.comboBoxRadiative,3,1)
        self.modelRadiative.addItem(self.tr("Paroi grise ou noire\n"\
                                            " + profil de temperature interne impose"), 'itpimp')
        self.modelRadiative.addItem(self.tr("Paroi grise ou noire\n"\
                                            " + profil de temperature externe impose"), 'ipgrno')
##         self.modelRadiative.addItem(self.tr("Paroi reflechissante\n"\
##                                               " + profil de temperature externe impose"), 'iprefl')
        self.modelRadiative.addItem(self.tr("Paroi grise ou noire\n"\
                                            " + flux de conduction impose en paroi"), 'ifgrno')
##         self.modelRadiative.addItem(self.tr("Paroi reflechissante\n"\
##                                               " + flux de conduction impose en paroi"), 'ifrefl')

        # Validator
        validatorZone = QtPage.IntValidator(self.lineEditZone, min=0)
        validatorZone.setExclusiveMin(True)
        self.lineEditZone.setValidator(validatorZone)

        # Connections
        self.connect(self.treeViewBoundaries, SIGNAL("clicked(const QModelIndex &)"), self.slotSelectBoundary)
        self.connect(self.comboBoxRadiative,  SIGNAL("activated(const QString&)"), self.slotRadiativeChoice)
        self.connect(self.lineEditZone, SIGNAL("textChanged(const QString &)"), self.slotZone)

        d = LocalizationModel('BoundaryZone', self.case)
        for zone in d.getZones():
            nature = zone.getNature()
            if nature == 'wall':
                label = zone.getLabel()
                bdModel = Boundary("radiative_wall", label, self.case)
                bdModel.getRadiativeChoice()
                codeNumber = zone.getCodeNumber()
                local = zone.getLocalization()
                self.modelBoundaries.insertItem(label, codeNumber, nature, local)

        self.frameRadiative.hide()

    def setVariablesForCondition(self, cond):
        """
        Put variables for the type of condition choised.
        """
        liste = self.getListVariablesForCondition(cond)
        self.setRadiativeWall(liste)


    def getListVariablesForCondition(self, cond):
        """
        Get list of variables for condition choosed
        """
        if cond == 'itpimp':
            liste = [(0, self.tr("Emissivite"), '',  'EPSP',  'emissivity'), 
                     (1, self.tr("Température initiale"), 'K', 'TINTP', 'internal_temperature_profile')]
        if cond == 'ipgrno':
            liste = [(0, self.tr("Emissivite"), '',  'EPSP',  'emissivity'), 
                     (1, self.tr("Conductivité"), 'W/m/K', 'XLAMP', 'thermal_conductivity'), 
                     (2, self.tr("Thickness"), 'm', 'EPAP' , 'thickness'),
                     (3, self.tr("Température externe"), 'K', 'TEXTP', 'external_temperature_profile'),
                     (4, self.tr("Température initiale"), 'K', 'TINTP', 'internal_temperature_profile')]
##        if cond == 'iprefl':
##            list = [(0, self.xlamp,t.XLAMP, 'W/m/K', 'XLAMP'), 
##                    (1, self.epap, t.EPAP,  'm', 'EPAP'), 
##                    (2, self.textp,t.TEXTP, 'K', 'TEXTP'), 
##                    (3, self.tintp,t.TINTP, 'K', 'TINTP')]
##            self.f43 = Tix.Frame(self.f4, relief=FLAT)
##            self.f43.pack(side=TOP, fill=X, pady=10)
##            frad = self.f43
        if cond == 'ifgrno':
            liste = [(0, self.tr("Emissivite"),'', 'EPSP', 'emissivity'),   
                     (1, self.tr("Flux de conduction"), 'W/m2', 'FLUX',  'flux'),
                     (2, self.tr("Température initiale"), 'K', 'TINTP', 'internal_temperature_profile')]
##        if cond == 'ifrefl':
##            list = [(0, self.flux, t.FLUX, 'W/m2', 'FLUX'), 
##                    (1, self.tintp, t.TINTP, 'K', 'TINTP')]
##            self.f45 = Tix.Frame(self.f4, relief=FLAT)
##            self.f45.pack(side=TOP, fill=X, pady=10)
##            frad = self.f45
        return liste
        

    def setRadiativeWall(self, liste):
        """
        Create the Page layout.
        """
        bdModel = Boundary("radiative_wall", self.label, self.case)

        if hasattr(self, "modelScalars"): del self.modelScalars
        self.modelScalars = StandardItemModelScalars(bdModel, liste)
        self.tableViewScalars.setModel(self.modelScalars)
        self.tableViewScalars.setAlternatingRowColors(True)
        self.tableViewScalars.setSelectionBehavior(QAbstractItemView.SelectRows)
        
        # number of zone
        self.nb_zone = bdModel.getOutputRadiativeZone()
        self.lineEditZone.setText(QString(str(self.nb_zone)))
        
        
    @pyqtSignature("const QModelIndex&")
    def slotSelectBoundary(self, index):
        log.debug("slotSelectBoundary")
        label, codeNumber, nature, local = self.modelBoundaries.getItem(index.row())
        self.label = label

##         self.mdl = Boundary('radiative_wall', label, self.case)
        from Pages.ThermalRadiationModel import ThermalRadiationModel
##         if ThermalRadiationModel.getRadiativeModel == "off":
##         else:
        self.putTypeofCondition()
        self.frameRadiative.show()


    @pyqtSignature("const QString&")
    def slotZone(self, text):
        nb_zone, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            bdModel = Boundary("radiative_wall", self.label, self.case)
            bdModel.setOutputRadiativeZone(nb_zone)
            return nb_zone


    def putTypeofCondition(self):
        """
        Put variables for the type of condition choised.
        """
        bdModel = Boundary("radiative_wall", self.label, self.case)
        cond = bdModel.getRadiativeChoice()
        self.modelRadiative.setItem(str_model=cond)
        self.setVariablesForCondition(cond)        


    @pyqtSignature("const QString&")
    def slotRadiativeChoice(self, text):
        cond = self.modelRadiative.dicoV2M[str(text)]
        log.debug("slotRadiativeChoice cond = %s "%cond)
        cindex = self.treeViewBoundaries.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            label, zone, nature, local = self.modelBoundaries.getItem(row)
            bdModel = Boundary("radiative_wall", label, self.case)
            bdModel.setRadiativeChoice(cond)
            self.setVariablesForCondition(cond)
    
        
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