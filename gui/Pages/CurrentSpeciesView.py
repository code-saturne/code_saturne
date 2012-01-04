# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
This module contains the following classes and function:
- CurrentSpeciesView
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

from Base.Toolbox import GuiParam
from Pages.CurrentSpeciesForm import Ui_CurrentSpeciesForm
from Pages.CurrentSpeciesModel import CurrentSpeciesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CurrentSpeciesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class for enthalpy-Temperature tabutation
#-------------------------------------------------------------------------------

class StandardItemModelEnthalpyTemp(QStandardItemModel):

    def __init__(self, mdl):
        """
        Constructor.
        """
        QStandardItemModel.__init__(self)
        self.setColumnCount(2)
        self.setRowCount(3)
        self.modelSpecies = mdl
        self.species = self.modelSpecies.getSpecies()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        # Display
        if role == Qt.DisplayRole:
            if index.column() == 0:
                if index.row() == 0:
                    return QVariant("points number of tabulation")
                elif index.row() == 1:
                    return QVariant("Min temperature")
                elif index.row() == 2:
                    return QVariant("Max temperature")
                else:
                    return QVariant()
            elif index.column() == 1:
                if index.row() == 0:
                    text = str(self.species.getEnthalpyTempTabNb())
                    return QVariant(text)
                elif index.row() == 1:
                    text = str(self.species.getMinTempTab()) + " K"
                    return QVariant(text)
                elif index.row() == 2:
                    text = str(self.species.getMaxTempTab()) + " K"
                    return QVariant(text)
                else:
                    return QVariant()
        return QVariant()


    def headerData(self, section, orientation, role):
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# StandarItemModel class for composition
#-------------------------------------------------------------------------------

class StandardItemModelComposition(QStandardItemModel):

    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.modelSpecies = mdl
        self.species = self.modelSpecies.getSpecies()
        nbspecies = self.species.getCurrentSpeciesNb()
        self.setColumnCount(nbspecies)
        nbelem = self.species.getElementarySpeciesList()
        self.setRowCount(len(nbelem))


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            elemSpeciesName = self.species.getElementarySpeciesList()
            composition = self.species.getCurrentSpeciesCompositionList()[row]
            text=str(composition[col])
            return QVariant(text)
        return QVariant()


    def headerData(self, section, orientation, role):
        currentSpeciesName = self.species.getCurrentSpeciesList()
        currentElementarySpecies = self.species.getElementarySpeciesList()
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            text = str(currentSpeciesName[section])
            return QVariant(text)
        if orientation == Qt.Vertical and role == Qt.DisplayRole:
            text = str(currentElementarySpecies[section])
            return QVariant(text)
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# StandarItemModel class for elementary species molar mass
#-------------------------------------------------------------------------------

class StandardItemModelMolarMass(QStandardItemModel):
    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)
        self.modelSpecies = mdl
        self.species = self.modelSpecies.getSpecies()
        self.setColumnCount(2)
        nbspecies = self.species.getElementarySpeciesNb()
        self.setRowCount(nbspecies)


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            row = index.row()
            elemSpeciesName = self.species.getElementarySpeciesList()
            molarMass = self.species.getElementarySpeciesMolarMassesList()
            if index.column() == 0:
                text=str(elemSpeciesName[row])
                return QVariant(text)
            elif index.column() == 1:
                text = str(molarMass[row]) + " kg/mol"
                return QVariant(text)
        return QVariant()


    def headerData(self, section, orientation, role):
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class CurrentSpeciesView(QWidget, Ui_CurrentSpeciesForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_CurrentSpeciesForm.__init__(self)
        self.setupUi(self)

        self.case = case

        # widgets layout

        self.modelSpecies = CurrentSpeciesModel(self.case)

        self.modelEnthalpyTemp = StandardItemModelEnthalpyTemp(self.modelSpecies)
        self.modelComposition = StandardItemModelComposition(self.modelSpecies)
        self.modelMolarMass = StandardItemModelMolarMass(self.modelSpecies)

        self.tableViewEnthalpyTemp.setModel(self.modelEnthalpyTemp)
        self.tableViewComposition.setModel(self.modelComposition)
        self.tableViewMolarMass.setModel(self.modelMolarMass)

        self.tableViewEnthalpyTemp.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
        self.tableViewComposition.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        self.tableViewMolarMass.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)

        #col = self.tableViewEnthalpyTemp.model().rowCount() -1
        self.tableViewEnthalpyTemp.horizontalHeader().setResizeMode(1, QHeaderView.Stretch)
        #col = self.tableViewMolarMass.model().rowCount() -1
        self.tableViewMolarMass.horizontalHeader().setResizeMode(1, QHeaderView.Stretch)

        self.tableViewEnthalpyTemp.horizontalHeader().hide()
        self.tableViewMolarMass.horizontalHeader().hide()

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------


if __name__ == "__main__":
    pass


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
