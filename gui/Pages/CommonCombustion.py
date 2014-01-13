# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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

import os, re, sys
import string

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

class Species:

    def __init__(self):
        """ constructor """

        self.currentSpeciesNb                 = 0
        self.enthalpyTempTabNb                = 0
        self.currentSpeciesList               = []
        self.minTempTab                       = 0
        self.maxTempTab                       = 0
        self.elementarySpeciesNb              = 0
        self.elementarySpeciesMolarMassesList = []
        self.elementarySpeciesList            = []
        self.currentSpeciesCompositionList    = []

        self.Comment = {}
        self.Comment["currentSpeciesNb"]   = "Nb especes courantes"
        self.Comment["enthalpyTempTabNb"]  = "Nb de points de tabulation ENTH-TEMP"
        self.Comment["minTempTab"]         = "Tmin"
        self.Comment["maxTempTab"]         = "Tmax"
        self.Comment["elementarySpeciesNb"]= "Nb especes elementaires / Composition C,H,O,N"


    def getDefault(self) :
        self.currentSpeciesNb                 = 8
        self.enthalpyTempTabNb                = 8
        self.currentSpeciesList               = ['CH4','C2H4','CO','O2','CO2','H2O','N2','C(S)']
        self.minTempTab                       = 300.
        self.maxTempTab                       = 2400.
        self.elementarySpeciesNb              = 4
        self.elementarySpeciesMolarMassesList = [0.012, 0.001, 0.016, 0.014]
        self.elementarySpeciesList            = ['C','H','O','N']
        self.currentSpeciesCompositionList    = [[1,2,1,0,1,0,0,1],[4,4,0,0,0,2,0,0],[0,0,1,2,2,1,0,0],[0,0,0,0,0,0,2,0]]
        return


    def getCurrentSpeciesNb(self) :
        return self.currentSpeciesNb


    def getEnthalpyTempTabNb(self) :
        return self.enthalpyTempTabNb


    def setCurrentSpeciesNb(self, value) :
        self.currentSpeciesNb = value


    def setEnthalpyTempTabNb(self, value) :
        self.enthalpyTempTabNb = value


    def getCurrentSpeciesList(self):
        return self.currentSpeciesList


    def getElementarySpeciesList(self):
        return self.elementarySpeciesList


    def getMinTempTab(self):
        return self.minTempTab


    def getMaxTempTab(self):
        return self.maxTempTab


    def getElementarySpeciesNb(self):
        return self.elementarySpeciesNb


    def getElementarySpeciesMolarMassesList(self):
        return self.elementarySpeciesMolarMassesList


    def getCurrentSpeciesCompositionList(self):
        return self.currentSpeciesCompositionList


    def setCurrentSpeciesList(self, value):
        self.currentSpeciesList = value


    def setElementarySpeciesList(self, value):
        self.elementarySpeciesList = value


    def setMinTempTab(self, value):
        self.minTempTab = value


    def setMaxTempTab(self, value):
        self.maxTempTab = value


    def setElementarySpeciesNb(self, value):
        self.elementarySpeciesNb = value


    def setElementarySpeciesMolarMassesList(self, value):
        self.elementarySpeciesMolarMassesList = value


    def setCurrentSpeciesCompositionList(self, value):
        self.currentSpeciesCompositionList = value


    def addElementarySpecie(self, name, molarMass):
        print("Not yet implemented")


    def cancelElementarySpecie(self, number):
        print("Not yet implemented")


    def addCurrentSpecie(self, name, composition):
        print("Not yet implemented")


    def cancelElementarySpecie(self, number):
        print("Not yet implemented")

