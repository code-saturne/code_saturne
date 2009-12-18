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

import os, re, sys
import string

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.XMLvariables import Model
import Pages.CommonCombustion as CommonCombustion

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Coal(Model):
    def __init__(self):
        self.classesNumber = 1
        self.initDiameterClasses = [0.]
        self.CDryComposition = 0
        self.HDryComposition = 0
        self.ODryComposition = 0
        self.Dry = 0
        self.PCIValue = 0
        self.thermalCapacity = 0
        self.density = 0
        self.CDryCompositionCoke = 0
        self.HDryCompositionCoke = 0
        self.ODryCompositionCoke = 0
        self.PCICokeValue = 0
        self.ashesRatio = 0
        self.ashesFormingEnthalpy = 0
        self.ashesThermalCapacity = 0
        self.humidity = 0
        self.IY1CH = 0
        self.Y1CH = 0
        self.IY2CH = 0
        self.Y2CH = 0
        self.A1CH = 0
        self.A2CH = 0
        self.E1CH = 0
        self.E2CH = 0
        self.AHETCH_O2 = 0
        self.EHETCH_O2 = 0
        self.IOCHET_O2 = 0
        self.AHETCH_CO2 = 0
        self.EHETCH_CO2 = 0
        self.IOCHET_CO2 = 0


    def setClassesNumber(self, value):
        self.classesNumber = value


    def setInitDiameterClasses(self, value):
        self.initDiameterClasses = value


    def addInitDiameterClasses(self, value):
        self.initDiameterClasses.append(value)
        self.classesNumber += 1
        self.isLowerOrEqual(self.classesNumber, 10)


    def updateInitDiameterClasses(self, number, value):
        self.initDiameterClasses[number-1] = value


    def cancelInitDiameterClasses(self, number):
        self.classesNumber -= 1
        self.initDiameterClasses.pop(number-1)


    def setCDryComposition(self, value):
        self.CDryComposition = value


    def setODryComposition(self, value):
        self.ODryComposition = value


    def setHDryComposition(self, value):
        self.HDryComposition = value


    def setDry(self, value):
        self.Dry = value


    def setPCIValue(self, value):
        self.PCIValue = value


    def setThermalCapacity(self, value):
        self.thermalCapacity = value


    def setDensity(self, value):
        self.density = value


    def setCDryCompositionCoke(self, value):
        self.CDryCompositionCoke = value


    def setHDryCompositionCoke(self, value):
        self.HDryCompositionCoke = value


    def setODryCompositionCoke(self, value):
        self.ODryCompositionCoke = value


    def setPCICokeValue(self, value):
        self.PCICokeValue = value


    def setAshesRatio(self, value):
        self.ashesRatio = value


    def setAshesFormingEnthalpy(self, value):
        self.ashesFormingEnthalpy = value


    def setAshesThermalCapacity(self, value):
        self.ashesThermalCapacity = value


    def setHumidity(self, value):
        self.humidity = value


    def setIY1CH(self, value):
        self.IY1CH = value


    def setY1CH(self, value):
        self.Y1CH = value


    def setIY2CH(self, value):
        self.IY2CH = value


    def setY2CH(self, value):
        self.Y2CH = value


    def setA1CH(self, value):
        self.A1CH = value


    def setA2CH(self, value):
        self.A2CH = value


    def setE1CH(self, value):
        self.E1CH = value


    def setE2CH(self, value):
        self.E2CH = value


    def setAHETCH_O2(self, value):
        self.AHETCH_O2 = value


    def setEHETCH_O2(self, value):
        self.EHETCH_O2 = value


    def setIOCHET_O2(self, value):
        self.IOCHET_O2 = value


    def setAHETCH_CO2(self, value):
        self.AHETCH_CO2 = value


    def setEHETCH_CO2(self, value):
        self.EHETCH_CO2 = value


    def setIOCHET_CO2(self, value):
        self.IOCHET_CO2 = value


    def getClassesNumber(self):
        return self.classesNumber


    def getInitDiameterClasses(self):
        return self.initDiameterClasses


    def getCDryComposition(self):
        return self.CDryComposition


    def getODryComposition(self):
        return self.ODryComposition


    def getHDryComposition(self):
        return self.HDryComposition


    def getDry(self):
        return self.Dry


    def getPCIValue(self):
        return self.PCIValue


    def getThermalCapacity(self):
        return self.thermalCapacity


    def getDensity(self):
        return self.density


    def getHumidity(self):
        return self.humidity


    def getCDryCompositionCoke(self):
        return self.CDryCompositionCoke


    def getHDryCompositionCoke(self):
        return self.HDryCompositionCoke


    def getODryCompositionCoke(self):
        return self.ODryCompositionCoke


    def getPCICokeValue(self):
        return self.PCICokeValue


    def getAshesRatio(self):
        return self.ashesRatio


    def getAshesFormingEnthalpy(self):
        return self.ashesFormingEnthalpy


    def getAshesThermalCapacity(self):
        return self.ashesThermalCapacity


    def getIY1CH(self):
        return self.IY1CH


    def getY1CH(self):
        return self.Y1CH


    def getIY2CH(self):
        return self.IY2CH


    def getY2CH(self):
        return self.Y2CH


    def getA1CH(self):
        return self.A1CH


    def getA2CH(self):
        return self.A2CH


    def getE1CH(self):
        return self.E1CH


    def getE2CH(self):
        return self.E2CH


    def getAHETCH_O2(self):
        return self.AHETCH_O2


    def getEHETCH_O2(self):
        return self.EHETCH_O2


    def getIOCHET_O2(self):
        return self.IOCHET_O2


    def getAHETCH_CO2(self):
        return self.AHETCH_CO2


    def getEHETCH_CO2(self):
        return self.EHETCH_CO2


    def getIOCHET_CO2(self):
        return self.IOCHET_CO2

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Coals(Model):

    def __init__(self):
        """ constructor """
        #
        # "CARACTERISTIQUES CHARBON" part
        self.Comment = {}
        self.coalNumber                         = 0
        self.Comment["coalNumber"]              = "Nb de charbons"
        self.classesNumberList                  = []
        self.Comment["classesNumberList"]       = "Nb de classes"
        self.initDiameterClassesList            = []
        self.Comment["initDiameterClassesList"] = "Diametre initial de la classe (m)"
        self.CDryCompositionList                = []
        self.Comment["CDryCompositionList"]     = "Composition C (%) sur sec"
        self.HDryCompositionList                = []
        self.Comment["HDryCompositionList"]     = "Composition H (%) sur sec"
        self.ODryCompositionList                = []
        self.Comment["ODryCompositionList"]     = "Composition O (%) sur sec"
        self.DryList                            = []
        self.Comment["DryList"]                 = "PCI (J/kg) sur pur (0) sur sec 1"
        self.PCIValueList                       = []
        self.thermalCapacityList                = []
        self.Comment["thermalCapacityList"]     = "CP moyen du charbon (J/kg/K)"
        self.densityList                        = []
        self.Comment["densityList"]             = "Masse volumique (kg/m3)"
        self.humidityList                       = []
        self.Comment["humidityList"]            = "Humidite du charbon(%)"
        self.CDryCompositionCokeList            = []
        self.Comment["CDryCompositionCokeList"] = "Composition C (%) sur sec"
        self.HDryCompositionCokeList            = []
        self.Comment["HDryCompositionCokeList"] = "Composition H (%) sur sec"
        self.ODryCompositionCokeList            = []
        self.Comment["ODryCompositionCokeList"] = "Composition O (%) sur sec"
        self.PCICokeValueList                   = []
        self.Comment["PCICokeValueList"]        = "PCI sur sec (J/kg)"
        self.ashesRatioList                     = []
        self.Comment["ashesRatioList"]          = "Taux de cendres en masse (%)"
        self.ashesFormingEnthalpyList           = []
        self.Comment["ashesFormingEnthalpyList"]= "Enthalpie de formation des cendres (J/kg)"
        self.ashesThermalCapacityList           = []
        self.Comment["ashesThermalCapacityList"]= "CP des cendres (J/kg/K)"
        #
        # "DEVOLATILISATION" part
        self.IY1CHList                          = []
        self.Comment["Y1CHList"]                = "Coef. stoe. (adim) Y1 si = 0 calcul automatique"
        self.Y1CHList                           = []
        self.Comment["Y2CHList"]                = "Coef. stoe. (adim) Y2 si = 0 calcul automatique"
        self.IY2CHList                          = []
        self.Y2CHList                           = []
        self.A1CHList                           = []
        self.Comment["A1CHList"]                = "Facteur pre-exponentielle A1 (s-1)"
        self.A2CHList                           = []
        self.Comment["A2CHList"]                = "Facteur pre-exponentielle A2 (s-1)"
        self.E1CHList                           = []
        self.Comment["E1CHList"]                = "Energie d'activation E1 (J/mol)"
        self.E2CHList                           = []
        self.Comment["E2CHList"]                = "Energie d'activation E2 (J/mol)"
        #
        # "COMBUSTION HETEROGENE" part
        self.AHETCH_O2List                      = []
        self.Comment["AHETCH_O2List"]           = "Constante pre-exponentielle (kg/m2/s/atm)"
        self.EHETCH_O2List                      = []
        self.Comment["EHETCH_O2List"]           = "Energie d'activation (kcal/mol)"
        self.IOCHET_O2List                      = []
        self.Comment["IOCHET_O2List"]           = "Ordre de la reaction 0.5 si = 0 1 si = 1"
        self.AHETCH_CO2List                     = []
        self.Comment["AHETCH_CO2List"]          = "Constante pre-exponentielle (kg/m2/s/atm)"
        self.EHETCH_CO2List                     = []
        self.Comment["EHETCH_CO2List"]          = "Energie d'activation (kcal/mol)"
        self.IOCHET_CO2List                     = []
        self.Comment["IOCHET_CO2List"]          = "Ordre de la reaction 0.5 si = 0 1 si = 1"


    def getNumber(self):
        return self.coalNumber


    #
    # "CARACTERISTIQUES CHARBON" part

    def deleteCoal(self, number):
        self.coalNumber -= 1
        nb = number - 1
        self.classesNumberList.pop(nb)
        self.initDiameterClassesList.pop(nb)
        self.CDryCompositionList.pop(nb)
        self.HDryCompositionList.pop(nb)
        self.ODryCompositionList.pop(nb)
        self.DryList.pop(nb)
        self.PCIValueList.pop(nb)
        self.thermalCapacityList.pop(nb)
        self.densityList.pop(nb)
        self.humidityList.pop(nb)
        self.CDryCompositionCokeList.pop(nb)
        self.HDryCompositionCokeList.pop(nb)
        self.ODryCompositionCokeList.pop(nb)
        self.PCICokeValueList.pop(nb)
        self.ashesRatioList.pop(nb)
        self.ashesFormingEnthalpyList.pop(nb)
        self.ashesThermalCapacityList.pop(nb)
        self.IY1CHList.pop(nb)
        self.Y1CHList.pop(nb)
        self.IY2CHList.pop(nb)
        self.Y2CHList.pop(nb)
        self.A1CHList.pop(nb)
        self.A2CHList.pop(nb)
        self.E1CHList.pop(nb)
        self.E2CHList.pop(nb)
        self.AHETCH_O2List.pop(nb)
        self.EHETCH_O2List.pop(nb)
        self.IOCHET_O2List.pop(nb)
        self.AHETCH_CO2List.pop(nb)
        self.EHETCH_CO2List.pop(nb)
        self.IOCHET_CO2List.pop(nb)


    def addCoal(self, coal):
        self.coalNumber += 1
        self.isLowerOrEqual(self.coalNumber, 3)
        self.classesNumberList.append(coal.getClassesNumber())
        self.initDiameterClassesList.append(coal.getInitDiameterClasses())
        self.CDryCompositionList.append(coal.getCDryComposition())
        self.HDryCompositionList.append(coal.getHDryComposition())
        self.ODryCompositionList.append(coal.getODryComposition())
        self.DryList.append(coal.getDry())
        self.PCIValueList.append(coal.getPCIValue())
        self.thermalCapacityList.append(coal.getThermalCapacity())
        self.densityList.append(coal.getDensity())
        self.humidityList.append(coal.getHumidity())
        self.CDryCompositionCokeList.append(coal.getCDryCompositionCoke())
        self.HDryCompositionCokeList.append(coal.getHDryCompositionCoke())
        self.ODryCompositionCokeList.append(coal.getODryCompositionCoke())
        self.PCICokeValueList.append(coal.getPCICokeValue())
        self.ashesRatioList.append(coal.getAshesRatio())
        self.ashesFormingEnthalpyList.append(coal.getAshesFormingEnthalpy())
        self.ashesThermalCapacityList.append(coal.getAshesThermalCapacity())
        self.IY1CHList.append(coal.getIY1CH())
        self.Y1CHList.append(coal.getY1CH())
        self.IY2CHList.append(coal.getIY2CH())
        self.Y2CHList.append(coal.getY2CH())
        self.A1CHList.append(coal.getA1CH())
        self.A2CHList.append(coal.getA2CH())
        self.E1CHList.append(coal.getE1CH())
        self.E2CHList.append(coal.getE2CH())
        self.AHETCH_O2List.append(coal.getAHETCH_O2())
        self.EHETCH_O2List.append(coal.getEHETCH_O2())
        self.IOCHET_O2List.append(coal.getIOCHET_O2())
        self.AHETCH_CO2List.append(coal.getAHETCH_CO2())
        self.EHETCH_CO2List.append(coal.getEHETCH_CO2())
        self.IOCHET_CO2List.append(coal.getIOCHET_CO2())


    def updateCoal(self,number, coal) :
        nb = number - 1
        self.classesNumberList[nb]        = coal.getClassesNumber()
        self.initDiameterClassesList[nb]  = coal.getInitDiameterClasses()
        self.CDryCompositionList[nb]      = coal.getCDryComposition()
        self.HDryCompositionList[nb]      = coal.getHDryComposition()
        self.ODryCompositionList[nb]      = coal.getODryComposition()
        self.DryList[nb]                  = coal.getDry()
        self.PCIValueList[nb]             = coal.getPCIValue()
        self.thermalCapacityList[nb]      = coal.getThermalCapacity()
        self.densityList[nb]              = coal.getDensity()
        self.humidityList[nb]             = coal.getHumidity()
        self.CDryCompositionCokeList[nb]  = coal.getCDryCompositionCoke()
        self.HDryCompositionCokeList[nb]  = coal.getHDryCompositionCoke()
        self.ODryCompositionCokeList[nb]  = coal.getODryCompositionCoke()
        self.PCICokeValueList[nb]         = coal.getPCICokeValue()
        self.ashesRatioList[nb]           = coal.getAshesRatio()
        self.ashesFormingEnthalpyList[nb] = coal.getAshesFormingEnthalpy()
        self.ashesThermalCapacityList[nb] = coal.getAshesThermalCapacity()
        self.IY1CHList[nb]                = coal.getIY1CH()
        self.Y1CHList[nb]                 = coal.getY1CH()
        self.IY2CHList[nb]                = coal.getIY2CH()
        self.Y2CHList[nb]                 = coal.getY2CH()
        self.A1CHList[nb]                 = coal.getA1CH()
        self.A2CHList[nb]                 = coal.getA2CH()
        self.E1CHList[nb]                 = coal.getE1CH()
        self.E2CHList[nb]                 = coal.getE2CH()
        self.AHETCH_O2List[nb]            = coal.getAHETCH_O2()
        self.EHETCH_O2List[nb]            = coal.getEHETCH_O2()
        self.IOCHET_O2List[nb]            = coal.getIOCHET_O2()
        self.AHETCH_CO2List[nb]           = coal.getAHETCH_CO2()
        self.EHETCH_CO2List[nb]           = coal.getEHETCH_CO2()
        self.IOCHET_CO2List[nb]           = coal.getIOCHET_CO2()


    def getCoal(self, number):
        coal = Coal()
        nb = number - 1
        coal.setClassesNumber(self.classesNumberList[nb])
        coal.setInitDiameterClasses(self.initDiameterClassesList[nb])
        coal.setCDryComposition(self.CDryCompositionList[nb])
        coal.setHDryComposition(self.HDryCompositionList[nb])
        coal.setODryComposition(self.ODryCompositionList[nb])
        coal.setDry(self.DryList[nb])
        coal.setPCIValue(self.PCIValueList[nb])
        coal.setThermalCapacity(self.thermalCapacityList[nb])
        coal.setDensity(self.densityList[nb])
        coal.setHumidity(self.humidityList[nb])
        coal.setCDryCompositionCoke(self.CDryCompositionCokeList[nb])
        coal.setHDryCompositionCoke(self.HDryCompositionCokeList[nb])
        coal.setODryCompositionCoke(self.ODryCompositionCokeList[nb])
        coal.setPCICokeValue(self.PCICokeValueList[nb])
        coal.setAshesRatio(self.ashesRatioList[nb])
        coal.setAshesFormingEnthalpy(self.ashesFormingEnthalpyList[nb])
        coal.setAshesThermalCapacity(self.ashesThermalCapacityList[nb])
        coal.setIY1CH(self.IY1CHList[nb])
        coal.setY1CH(self.Y1CHList[nb])
        coal.setIY2CH(self.IY2CHList[nb])
        coal.setY2CH(self.Y2CHList[nb])
        coal.setA1CH(self.A1CHList[nb])
        coal.setA2CH(self.A2CHList[nb])
        coal.setE1CH(self.E1CHList[nb])
        coal.setE2CH(self.E2CHList[nb])
        coal.setAHETCH_O2(self.AHETCH_O2List[nb])
        coal.setEHETCH_O2(self.EHETCH_O2List[nb])
        coal.setIOCHET_O2(self.IOCHET_O2List[nb])
        coal.setAHETCH_CO2(self.AHETCH_CO2List[nb])
        coal.setEHETCH_CO2(self.EHETCH_CO2List[nb])
        coal.setIOCHET_CO2(self.IOCHET_CO2List[nb])
        return coal


    def getCoalNumber(self):
        return self.coalNumber


    def getClassesNumberList(self):
        return self.classesNumberList


    def getInitDiameterClassesList(self):
        return self.initDiameterClassesList


    def getCDryCompositionList(self):
        return self.CDryCompositionList


    def getHDryCompositionList(self):
        return self.HDryCompositionList


    def getODryCompositionList(self):
        return self.ODryCompositionList


    def getDryList(self):
        return self.DryList


    def getPCIValueList(self):
        return self.PCIValueList


    def getThermalCapacityList(self):
        return self.thermalCapacityList


    def getDensityList(self):
        return self.densityList


    def getHumidityList(self):
        return self.humidityList


    def getCDryCompositionCokeList(self):
        return self.CDryCompositionCokeList


    def getHDryCompositionCokeList(self):
        return self.HDryCompositionCokeList


    def getODryCompositionCokeList(self):
        return self.ODryCompositionCokeList


    def getPCICokeValueList(self):
        return self.PCICokeValueList


    def getAshesRatioList(self):
        return self.ashesRatioList


    def getAshesFormingEnthalpyList(self):
        return self.ashesFormingEnthalpyList


    def getAshesThermalCapacityList(self):
        return self.ashesThermalCapacityList


    def setCoalNumber(self, value):
        self.coalNumber = value


    def setClassesNumberList(self, value):
        self.classesNumberList = value


    def setInitDiameterClassesList(self, value):
        self.initDiameterClassesList = value


    def setCDryCompositionList(self, value):
        self.CDryCompositionList = value


    def setHDryCompositionList(self, value):
        self.HDryCompositionList = value


    def setODryCompositionList(self, value):
        self.ODryCompositionList = value


    def setDryList(self, value):
        self.DryList = value


    def setPCIValueList(self, value):
        self.PCIValueList = value


    def setThermalCapacityList(self, value):
        self.thermalCapacityList = value


    def setDensityList(self, value):
        self.densityList = value


    def setHumidityList(self, value):
        self.humidityList = value


    def setCDryCompositionCokeList(self, value):
        self.CDryCompositionCokeList = value


    def setHDryCompositionCokeList(self, value):
        self.HDryCompositionCokeList = value


    def setODryCompositionCokeList(self, value):
        self.ODryCompositionCokeList = value


    def setPCICokeValueList(self, value):
        self.PCICokeValueList = value


    def setAshesRatioList(self, value):
        self.ashesRatioList = value


    def setAshesFormingEnthalpyList(self, value):
        self.ashesFormingEnthalpyList = value


    def setAshesThermalCapacityList(self, value):
        self.ashesThermalCapacityList = value

    #
    # "DEVOLATISATION" part
    def getIY1CHList(self):
        return self.IY1CHList


    def getY1CHList(self):
        return self.Y1CHList


    def getIY2CHList(self):
        return self.IY2CHList


    def getY2CHList(self):
        return self.Y2CHList


    def getA1CHList(self):
        return self.A1CHList


    def getA2CHList(self):
        return self.A2CHList


    def getE1CHList(self):
        return self.E1CHList


    def getE2CHList(self):
        return self.E2CHList


    def setIY1CHList(self, value):
        self.IY1CHList = value


    def setY1CHList(self, value):
        self.Y1CHList = value


    def setIY2CHList(self, value):
        self.IY2CHList  = value


    def setY2CHList(self, value):
        self.Y2CHList = value


    def setA1CHList(self, value):
        self.A1CHList = value


    def setA2CHList(self, value):
        self.A2CHList = value


    def setE1CHList(self, value):
        self.E1CHList = value


    def setE2CHList(self, value):
        self.E2CHList = value
    #
    # "COMBUSTION HETEROGENE" part
    def getAHETCH_O2List(self):
        return self.AHETCH_O2List


    def getEHETCH_O2List(self):
        return self.EHETCH_O2List


    def getIOCHET_O2List(self):
        return self.IOCHET_O2List


    def getAHETCH_CO2List(self):
        return self.AHETCH_CO2List


    def getEHETCH_CO2List(self):
        return self.EHETCH_CO2List


    def getIOCHET_CO2List(self):
        return self.IOCHET_CO2List


    def setAHETCH_O2List(self, value):
        self.AHETCH_O2List = value


    def setEHETCH_O2List(self, value):
        self.EHETCH_O2List = value


    def setIOCHET_O2List(self, value):
        self.IOCHET_O2List = value


    def setAHETCH_CO2List(self, value):
        self.AHETCH_CO2List = value


    def setEHETCH_CO2List(self, value):
        self.EHETCH_CO2List = value


    def setIOCHET_CO2List(self, value):
        self.IOCHET_CO2List = value

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class RadiativTransfer:
    def __init__(self):
        """ constructor """
        #
        # "RAYONNEMENT" part
        self.Comment = {}
        #self.radiativTransfer                 = 0
        #self.Comment["radiativTransfer"]   = "Rayonnement : (0 : pas de rayonnement ; 1 : constant donne ci-dessous ; 2 : par Modak)"
        self.absorptionCoeff                  = 0
        self.Comment["absorptionCoeff"]    = "Coeff absorption pour le melange gazeux constant"

    #
    # "RAYONNEMENT" part
    #def getRadiativTransfer(self) :
    #    return self.radiativTransfer


    def getAbsorptionCoeff(self):
        return self.absorptionCoeff


    #def setRadiativTransfer(self, value) :
    #    self.radiativTransfer = value


    def setAbsorptionCoeff(self, value):
        self.absorptionCoeff = value

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Oxydant:
    """
    Oxydant compisition in O2,N2,H2O,CO2
    """
    def __init__(self):
        self.oxydantNumber = 1
        self.O2  = 0.25
        self.N2  = 0.75
        self.H2O = 0
        self.CO2 = 0


    def setOxydantNumber(self, value):
        self.oxydantNumber = value


    def setO2(self, value):
        self.O2 = value


    def setN2(self, value):
        self.N2 = value


    def setH2O(self, value):
        self.H2O = value


    def setCO2(self, value):
        self.CO2 = value


    def getOxydantNumber(self):
        return self.oxydantNumber


    def getO2(self):
        return self.O2


    def getN2(self):
        return self.N2


    def getH2O(self):
        return self.H2O


    def getCO2(self):
        return self.CO2

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Oxydants(Model):
    def __init__(self):
        """ constructor """
        self.Comment = {}
        self.oxydantNumber = 0
        self.Comment["oxydantNumber"] = ""
        self.oxydantNumberList = []
        self.Comment["oxydantNumberList"] = ""
        self.O2List  = []
        self.Comment["O2List"] = ""
        self.N2List  = []
        self.Comment["N2List"] = ""
        self.H2OList = []
        self.Comment["H2OList"] = ""
        self.CO2List = []
        self.Comment["CO2List"] = ""


    def deleteOxydant(self, number):
        self.oxydantNumber -= 1
        nb = number - 1
        self.oxydantNumberList.pop(nb)
        self.O2List.pop(nb)
        self.N2List.pop(nb)
        self.H2OList.pop(nb)
        self.CO2List.pop(nb)


    def addOxydant(self, oxy):
        self.oxydantNumber += 1
        self.isLowerOrEqual(self.oxydantNumber, 3)
        self.oxydantNumberList.append(oxy.getOxydantNumber())
        self.O2List.append(oxy.getO2())
        self.N2List.append(oxy.getN2())
        self.H2OList.append(oxy.getH2O())
        self.CO2List.append(oxy.getCO2())


    def updateOxydant(self, number, oxy) :
        nb = number - 1
        self.oxydantNumberList[nb] = oxy.getOxydantNumber()
        self.O2List[nb] = oxy.getO2()
        self.N2List[nb] = oxy.getN2()
        self.H2OList[nb] = oxy.getH2O()
        self.CO2List[nb] = oxy.getCO2()


    def getOxydant(self, number):
        oxy = Oxydant()
        nb = number - 1
        oxy.setOxydantNumber(self.oxydantNumberList[nb])
        oxy.setO2(self.O2List[nb])
        oxy.setN2(self.N2List[nb])
        oxy.setH2O(self.H2OList[nb])
        oxy.setCO2(self.CO2List[nb])
        return oxy


    def getNumber(self):
        return self.oxydantNumber


    def setNumber(self, value):
        self.oxydantNumber = value
        self.oxydantNumberList = []
        for i in range(0, value):
            self.oxydantNumberList.append(i)


    def getOxydantNumberList(self):
        return self.oxydantNumberList


    def setOxydantNumberList(self, value):
        self.oxydantNumberList = value
        self.oxydantNumber = len(value)


    def getO2List(self):
        return self.O2List


    def getN2List(self):
        return self.N2List


    def getH2OList(self):
        return self.H2OList


    def getCO2List(self):
        return self.CO2List


    def setO2List(self, value):
        self.O2List = value


    def setN2List(self, value):
        self.N2List = value


    def setH2OList(self, value):
        self.H2OList = value


    def setCO2List(self, value):
        self.CO2List = value

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class CoalThermoChemistryModel:
    def __init__(self, fileName, case):
        """ constructor """
        self.case = case

        # A remplacer par le bon repertoire
        self.fileName = fileName
        self.varComment = {}

        self.species = CommonCombustion.Species()
        self.radiativTransfer = RadiativTransfer()
        self.coals = Coals()
        self.oxydants = Oxydants()

        # Load file
        if not self.load():

            # create default coal
            coal = Coal()
            self.coals.addCoal(coal)

            # create default species
            self.species.getDefault()

            # create default oxydant
            oxy = Oxydant()
            self.oxydants.addOxydant(oxy)


    def getCoals(self):
        return self.coals


    def getSpecies(self):
        return self.species


    def getOxydants(self):
        return self.oxydants


    def load(self):
        """ read thermochimestry file"""
        #FIXME bug to obtain case_path
        filePath = self.case['data_path']+"/"+self.fileName
        try :
            ThermoChFile = open(filePath, "r")
        except :
            return 0

        # "THERMOCHIMIE" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^THERMOCHIMIE")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath
            msg = msg + "\nNo string 'THERMOCHIMIE' in the first line\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        value = self.__readIntValue(line, "currentSpeciesNb", filePath)
        self.species.setCurrentSpeciesNb(value)

        line = ThermoChFile.readline()
        value = self.__readIntValue(line, "enthalpyTempTabNb", filePath)
        self.species.setEnthalpyTempTabNb(value)

        # "ESPECE COURANTE" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^ESPECES COURANTES")
        if not RegExpr.match(line):
            msg =  "Reading file error: " + filePath
            msg = msg + "\nNo string 'ESPECES COURANTES' in the first line\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        try:
            self.species.setCurrentSpeciesList(line.split())
        except:
            msg = "Reading file error: " + filePath + " bad list of current species"
            raise ValueError, msg

        line = ThermoChFile.readline()
        value = self.__readFloatValue(line, "MinTempTab", filePath)
        self.species.setMinTempTab(value)

        line = ThermoChFile.readline()
        value = self.__readFloatValue(line, "MaxTempTab", filePath)
        self.species.setMaxTempTab(value)

        line = ThermoChFile.readline()
        value = self.__readIntValue(line, "ElementarySpeciesNb", filePath)
        self.species.setElementarySpeciesNb(value)

        ElementarySpeciesMolarMassesList = []
        CurrentSpeciesCompositionList    = []
        for elemSpec in range(self.species.getElementarySpeciesNb()):
            try:
                line = ThermoChFile.readline()
                valuesList = line.split()
                ElementarySpeciesMolarMassesList.append(float(valuesList[0]))
                composition = []
                for currentSpec in range(self.species.getCurrentSpeciesNb()):
                    composition.append(int(valuesList[1 + currentSpec]))
                CurrentSpeciesCompositionList.append(composition)

            except:
                msg = "Reading file error: " + filePath
                msg = msg + "\nElementarySpeciesMolarMassesList or CurrentSpeciesCompositionList reading\n"
                raise ValueError, msg
        self.species.setElementarySpeciesList(["C","H","O","N"])
        self.species.setElementarySpeciesMolarMassesList(ElementarySpeciesMolarMassesList)
        self.species.setCurrentSpeciesCompositionList(CurrentSpeciesCompositionList)

        # "RAYONNEMENT" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^RAYONNEMENT")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + "\nNo string 'RAYONNEMENT'\n"
            raise ValueError, msg

        #line = ThermoChFile.readline()
        #value = self.__readIntValue(line, "RadiativTransfer", filePath)
        #self.radiativTransfer.setRadiativTransfer(value)

        line = ThermoChFile.readline()
        value = self.__readFloatValue(line, "AbsorptionCoeff", filePath)
        self.radiativTransfer.setAbsorptionCoeff(value)

        # "CARACTERISTIQUES CHARBONS" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^CARACTERISTIQUES CHARBONS")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + "\nNo string 'CARACTERISTIQUES CHARBONS'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        value = self.__readIntValue(line, "CoalNumber", filePath)
        self.coals.setCoalNumber(value)

        line = ThermoChFile.readline()
        values = self.__readIntCoalValues(line, "ClassesNumberList", filePath)
        self.coals.setClassesNumberList(values)

        line = ThermoChFile.readline()
        try:
            strList = line.split()
            InitDiameterClassesList = []
            pt = 0
            for coal in range(self.coals.coalNumber):
                InitDiameter = []
                for Class in range(self.coals.classesNumberList[coal]):
                    InitDiameter.append(float(strList[pt+Class]))
                pt += self.coals.classesNumberList[coal]
                InitDiameterClassesList.append(InitDiameter)

            self.coals.setInitDiameterClassesList(InitDiameterClassesList)

        except:
            msg = "Reading file error: " + filePath + "\nbad value for InitDiameterClassesList\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "CDryCompositionList", filePath)
        self.coals.setCDryCompositionList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "HDryCompositionList", filePath)
        self.coals.setHDryCompositionList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "ODryCompositionList", filePath)
        self.coals.setODryCompositionList(values)

        line = ThermoChFile.readline()
        intValues, floatValues = self.__readIntFloatCoalValues(line, "PCIValue", filePath)
        self.coals.setDryList(intValues)
        self.coals.setPCIValueList(floatValues)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "thermalCapacity", filePath)
        self.coals.setThermalCapacityList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "density", filePath)
        self.coals.setDensityList(values)

        # "Coke" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^Coke")
        if not RegExpr.match(line):
            msg = "_reading file error: " + filePath + "\nNo string 'Coke'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "CDryCompositionCokeList", filePath)
        self.coals.setCDryCompositionCokeList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "HDryCompositionCokeList", filePath)
        self.coals.setHDryCompositionCokeList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "ODryCompositionCokeList", filePath)
        self.coals.setODryCompositionCokeList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "PCICokeValue", filePath)
        self.coals.setPCICokeValueList(values)

        # "Cendres" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^Cendres")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + "\nNo string 'Cendres'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "ashesRatioList", filePath)
        self.coals.setAshesRatioList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "AshesFormingEnthalpy", filePath)
        self.coals.setAshesFormingEnthalpyList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "AshesThermalCapacity", filePath)
        self.coals.setAshesThermalCapacityList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "humidity", filePath)
        self.coals.setHumidityList(values)

        # "Devolatilisation (Kobayashi)" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("^Parametres de devolatilisation")
        if not RegExpr.match(line):
            msg =  "Reading file error: " + filePath + \
                   "\nNo string 'Parametres de devolatilisation (modele de Kobayashi)'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        intValues, floatValues = self.__readIntFloatCoalValues(line, "IY1CH", filePath)
        self.coals.setIY1CHList(intValues)
        self.coals.setY1CHList(floatValues)

        line = ThermoChFile.readline()
        intValues, floatValues = self.__readIntFloatCoalValues(line, "IY2CH", filePath)
        self.coals.setIY2CHList(intValues)
        self.coals.setY2CHList(floatValues)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "A1CH", filePath)
        self.coals.setA1CHList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "A2CH", filePath)
        self.coals.setA2CHList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "E1CH", filePath)
        self.coals.setE1CHList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "E2CH", filePath)
        self.coals.setE2CHList(values)

        # "Devolatisation (Kobayashi)" part
        line = ThermoChFile.readline()
        RegExpr = re.compile("Parametres de combustion heterogene")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + \
                   "\nNo string '^Parametres de combustion heterogene O2 (modele a sphere retrecissante)'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "AHETCH_O2", filePath)
        self.coals.setAHETCH_O2List(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "EHETCH_O2", filePath)
        self.coals.setEHETCH_O2List(values)

        line = ThermoChFile.readline()
        values = self.__readIntCoalValues(line, "IOCHET_O2", filePath)
        self.coals.setIOCHET_O2List(values)

        line = ThermoChFile.readline()
        RegExpr = re.compile("^Parametres de combustion heterogene")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + \
                   "\nNo string 'Parametres de combustion heterogene CO2 (modele a sphere retrecissante)'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "AHETCH_CO2", filePath)
        self.coals.setAHETCH_CO2List(values)

        line = ThermoChFile.readline()
        values = self.__readFloatCoalValues(line, "EHETCH_CO2", filePath)
        self.coals.setEHETCH_CO2List(values)

        line = ThermoChFile.readline()
        values = self.__readIntCoalValues(line, "IOCHET_CO2", filePath)
        self.coals.setIOCHET_CO2List(values)

        line = ThermoChFile.readline()
        RegExpr = re.compile("^CARACTERISTIQUES OXYDANTS")
        if not RegExpr.match(line):
            msg = "Reading file error: " + filePath + \
                   "\nNo string 'CARACTERISTIQUES OXYDANTS (O2,N2,H2O,CO2)'\n"
            raise ValueError, msg

        line = ThermoChFile.readline()
        value = self.__readIntValue(line, "oxydantNumber", filePath)
        self.oxydants.setNumber(value)

        line = ThermoChFile.readline()
        values = self.__readFloatOxydantsValues(line, "O2 for oxydants", filePath)
        self.oxydants.setO2List(values)

        line = ThermoChFile.readline()
        values = self.__readFloatOxydantsValues(line, "N2 for oxydants", filePath)
        self.oxydants.setN2List(values)

        line = ThermoChFile.readline()
        values = self.__readFloatOxydantsValues(line, "H2O for oxydants", filePath)
        self.oxydants.setH2OList(values)

        line = ThermoChFile.readline()
        values = self.__readFloatOxydantsValues(line, "CO2 for oxydants", filePath)
        self.oxydants.setCO2List(values)

        ThermoChFile.close()

        return 1


    def save(self):
        """ write and save thermochimestry file"""
        #
        # Openning
        filePath = self.case['data_path']+"/"+self.fileName
        try :
            ThermoChFile = open(filePath, "w")
        except :
            msg = "Openning file error: " + filePath
            raise ValueError, msg
        #
        # writing
        ThermoChFile.write("THERMOCHIMIE\n")
        chain = str(self.species.getCurrentSpeciesNb())
        comment = self.createComment(len(chain),self.species.Comment["currentSpeciesNb"])
        ThermoChFile.write(chain + comment + "\n")

        chain = str(self.species.getEnthalpyTempTabNb())
        comment = self.createComment(len(chain),self.species.Comment["enthalpyTempTabNb"])
        ThermoChFile.write(chain + comment + "\n")
        #
        ThermoChFile.write("ESPECES COURANTES\n")
        ThermoChFile.write(self.__preWriteList(self.species.getCurrentSpeciesList())+"\n")

        chain = str(self.species.getMinTempTab())
        comment = self.createComment(len(chain),self.species.Comment["minTempTab"])
        ThermoChFile.write(chain + comment + "\n")

        chain = str(self.species.getMaxTempTab())
        comment = self.createComment(len(chain),self.species.Comment["maxTempTab"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.species.getElementarySpeciesNb())
        comment = self.createComment(len(chain),self.species.Comment["elementarySpeciesNb"])
        ThermoChFile.write(chain + comment +"\n")

        list2write=[]
        ind = 0
        for i in self.species.getElementarySpeciesMolarMassesList():
            listTMP = []
            listTMP.append(i)
            for j in self.species.getCurrentSpeciesCompositionList()[ind]:
                listTMP.append(j)
            ind += 1
            list2write.append(listTMP)
        ThermoChFile.write(self.__preWriteList(list2write))
        #
        ThermoChFile.write("RAYONNEMENT\n")

        #chain = str(self.radiativTransfer.getRadiativTransfer())
        #comment = self.createComment(len(chain),self.radiativTransfer.Comment["radiativTransfer"])
        #ThermoChFile.write(chain + comment +"\n")

        chain = str(self.radiativTransfer.getAbsorptionCoeff())
        comment = self.createComment(len(chain),self.radiativTransfer.Comment["absorptionCoeff"])
        ThermoChFile.write(chain + comment +"\n")
        #
        ThermoChFile.write("CARACTERISTIQUES CHARBONS\n")

        chain = str(self.coals.getCoalNumber())
        comment = self.createComment(len(chain),self.coals.Comment["coalNumber"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getClassesNumberList()))
        comment = self.createComment(len(chain),self.coals.Comment["classesNumberList"])
        ThermoChFile.write(chain + comment +"\n")

        listTMP = []
        for i in self.coals.getInitDiameterClassesList():
            for j in i:
                listTMP.append(j)
        chain = str(self.__preWriteList(listTMP))
        comment = self.createComment(len(chain),self.coals.Comment["initDiameterClassesList"])
        ThermoChFile.write(chain + comment + "\n")

        chain = str(self.__preWriteList(self.coals.getCDryCompositionList()))
        comment = self.createComment(len(chain),self.coals.Comment["CDryCompositionList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getHDryCompositionList()))
        comment = self.createComment(len(chain),self.coals.Comment["HDryCompositionList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getODryCompositionList()))
        comment = self.createComment(len(chain),self.coals.Comment["ODryCompositionList"])
        ThermoChFile.write(chain + comment +"\n")

        list2write=[]

        for coal in range(self.coals.coalNumber):
            list2write.append(self.coals.getDryList()[coal])
            list2write.append(self.coals.getPCIValueList()[coal])

        chain = str(self.__preWriteList(list2write))
        comment = self.createComment(len(chain),self.coals.Comment["DryList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getThermalCapacityList()))
        comment = self.createComment(len(chain),self.coals.Comment["thermalCapacityList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getDensityList()))
        comment = self.createComment(len(chain),self.coals.Comment["densityList"])
        ThermoChFile.write(chain + comment +"\n")
        #
        ThermoChFile.write("Coke\n")

        chain = str(self.__preWriteList(self.coals.getCDryCompositionCokeList()))
        comment = self.createComment(len(chain),self.coals.Comment["CDryCompositionCokeList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getHDryCompositionCokeList()))
        comment = self.createComment(len(chain),self.coals.Comment["HDryCompositionCokeList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getODryCompositionCokeList()))
        comment = self.createComment(len(chain),self.coals.Comment["ODryCompositionCokeList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getPCICokeValueList()))
        comment = self.createComment(len(chain),self.coals.Comment["PCICokeValueList"])
        ThermoChFile.write(chain + comment +"\n")
        #
        ThermoChFile.write("Cendres\n")

        chain = str(self.__preWriteList(self.coals.getAshesRatioList()))
        comment = self.createComment(len(chain),self.coals.Comment["ashesRatioList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getAshesFormingEnthalpyList()))
        comment = self.createComment(len(chain),self.coals.Comment["ashesFormingEnthalpyList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getAshesThermalCapacityList()))
        comment = self.createComment(len(chain),self.coals.Comment["ashesThermalCapacityList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = str(self.__preWriteList(self.coals.getHumidityList()))
        comment = self.createComment(len(chain),self.coals.Comment["humidityList"])
        ThermoChFile.write(chain + comment +"\n")

        #
        ThermoChFile.write("Parametres de devolatilisation (modele de Kobayashi)\n")

        list2write=[]
        ind = 0
        for i in self.coals.getIY1CHList():
            j = self.coals.getY1CHList()[ind]
            list2write.append(i)
            list2write.append(j)
            ind += 1
        chain = self.__preWriteList(list2write)
        comment = self.createComment(len(chain),self.coals.Comment["Y1CHList"])
        ThermoChFile.write(chain + comment +"\n")

        list2write=[]

        ind = 0
        for i in self.coals.getIY2CHList():
            j = self.coals.getY2CHList()[ind]
            list2write.append(i)
            list2write.append(j)
            ind += 1

        chain = self.__preWriteList(list2write)
        comment = self.createComment(len(chain),self.coals.Comment["Y2CHList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getA1CHList())
        comment = self.createComment(len(chain),self.coals.Comment["A1CHList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getA2CHList())
        comment = self.createComment(len(chain),self.coals.Comment["A2CHList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getE1CHList())
        comment = self.createComment(len(chain),self.coals.Comment["E1CHList"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getE2CHList())
        comment = self.createComment(len(chain),self.coals.Comment["E2CHList"])
        ThermoChFile.write(chain + comment +"\n")

        ThermoChFile.write("Parametres de combustion heterogene O2 (modele a sphere retrecissante)\n")

        chain = self.__preWriteList(self.coals.getAHETCH_O2List())
        comment = self.createComment(len(chain),self.coals.Comment["AHETCH_O2List"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getEHETCH_O2List())
        comment = self.createComment(len(chain),self.coals.Comment["EHETCH_O2List"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getIOCHET_O2List())
        comment = self.createComment(len(chain),self.coals.Comment["IOCHET_O2List"])
        ThermoChFile.write(chain + comment +"\n")

        ThermoChFile.write("Parametres de combustion heterogene CO2 (modele a sphere retrecissante)\n")

        chain = self.__preWriteList(self.coals.getAHETCH_CO2List())
        comment = self.createComment(len(chain),self.coals.Comment["AHETCH_CO2List"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getEHETCH_CO2List())
        comment = self.createComment(len(chain),self.coals.Comment["EHETCH_CO2List"])
        ThermoChFile.write(chain + comment +"\n")

        chain = self.__preWriteList(self.coals.getIOCHET_CO2List())
        comment = self.createComment(len(chain),self.coals.Comment["IOCHET_CO2List"])
        ThermoChFile.write(chain + comment +"\n")

        ThermoChFile.write("CARACTERISTIQUES OXYDANTS (O2,N2,H2O,CO2)\n")

        chain = str(self.oxydants.getNumber())
        ThermoChFile.write(chain + "\n")

        chain = self.__preWriteList(self.oxydants.getO2List())
        ThermoChFile.write(chain + "\n")

        chain = self.__preWriteList(self.oxydants.getN2List())
        ThermoChFile.write(chain + "\n")

        chain = self.__preWriteList(self.oxydants.getH2OList())
        ThermoChFile.write(chain + "\n")

        chain = self.__preWriteList(self.oxydants.getCO2List())
        ThermoChFile.write(chain + "\n")

        ThermoChFile.close()


    def createComment(self, length, comment):

        nb = length
        blancString = ""
        while nb < 60 :
            blancString += " "
            nb += 1
        return blancString+comment


    def __readIntValue(self, line, valueName, filePath):

        try:
            value = int(line.split()[0])
        except:
            msg = "Reading file error: " + filePath + "\n"\
                  "bad value for " + valueName
            raise ValueError, msg
        return value


    def __readFloatValue(self, line, valueName, filePath):

        try:
            value = float(line.split()[0])
        except:
            msg = "Reading file error: " + filePath + "\n" \
                  "bad value for " + valueName
            raise ValueError, msg
        return value


    def __readFloatCoalValues(self, line, valueName, filePath):
        """ read coal values"""
        try:
            strList = line.split()
            value = []
            for coal in range(self.coals.coalNumber):
                value.append(float(strList[coal]))
        except:
            msg = "Reading file error: " + filePath + "\n"\
                  "bad value for "+ valueName
            raise ValueError, msg
        return value


    def __readIntCoalValues(self, line, valueName, filePath):
        """ read coal values"""
        try:
            strList = line.split()
            value = []
            for coal in range(self.coals.coalNumber):
                value.append(int(strList[coal]))
        except:
            msg =  "Reading file error: " + filePath + " bad value for "+ valueName
            raise ValueError, msg
        return value


    def __readIntFloatCoalValues(self, line, valueName, filePath):
        """ read coal values"""
        try:
            strList = line.split()
            intValue = []
            floatValue = []
            for coal in range(0,2*self.coals.coalNumber,2):
                intValue.append(int(strList[coal]))
                floatValue.append(float(strList[coal+1]))
        except:
            msg = "Reading file error: " + filePath + " bad value for "+ valueName
            raise ValueError, msg
        return intValue, floatValue


    def __preWriteList(self, List):
        """build a line with elements list"""
        string = ""
        for i in List:
            if type(i) == type([]) :
                stringTMP = self.__preWriteList(i)
                string += stringTMP + "\n"
            else:
                if string != "":
                    string += " " + str(i) + " "
                else:
                    string += str(i)
        return string


    def __readFloatOxydantsValues(self, line, valueName, filePath):
        """ read oxydants values"""
        try:
            strList = line.split()
            values = []
            for oxy in range(self.oxydants.getNumber()):
                values.append(float(strList[oxy]))
        except:
            msg = "Reading file error: " + filePath + "\nbad value for " + valueName
            raise ValueError, msg
        return values

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
