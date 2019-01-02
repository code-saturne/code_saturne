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

from code_saturne.Pages.LocalizationModel import BoundaryZone, Zone
from code_saturne.Pages.BoundaryNeptune import *

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class LocalizationModel(object):
    """
    Abstract class
    """
    def __new__(cls, typeZone, case):
        """
        """
        if typeZone == 'BoundaryZone':
            return BoundaryZone.__new__(BoundaryLocalizationModel, case)
        else:
            raise ValueError("Unknown type zone")


    def __init__(self, typeZone, case):
        """
        """
        self.case = case
        self._initModel()
        self._typeZone = typeZone


    def getLocalizationsZonesList(self):
        """
        Return list of localizations used by zones
        """
        zones = self.getZones()
        locals = []
        for zone in zones:
            locals.append(zone.getLocalization())

        return locals


    def getLabelsZonesList(self):
        """
        Return list of labels used by zones
        """
        zones = self.getZones()
        labels = []
        for zone in zones:
            labels.append(zone.getLabel())

        return labels


    def getCodeNumbersList(self):
        """
        Return list of code numbers used
        """
        zones = self.getZones()
        codes = []
        for zone in zones:
            codes.append(zone.getCodeNumber())

        return codes


    def getZones(self):
        """
        Return zones list after XML file reading (virtual method)
        """
        return []


    def getMaxCodeNumber(self):
        """
        Return maximum of code number's values to put on name
        """
        zones = self.getZones()
        codeNumber = 0
        for zone in zones:
            codeNumber = max(codeNumber, zone.getCodeNumber())

        return codeNumber


    def getMaxNumberNature(self, nature):
        """
        Return maximum of nature number's values to put on name
        """
        Model().isInList(nature, self._natureList)


    def renameLabel(self, label, newLabel):
        """
        """
        if label == newLabel: return

        labels = self.getLabelsZonesList()

        Model().isInList(label, labels)
        Model().isNotInList(newLabel, labels)


    def renameName(self, codeNumber, newCodeNumber):
        """
        """
        if codeNumber == newCodeNumber: return

        labels = self.getLabelsZonesList()

        Model().isInt(newCodeNumber)


    def replaceLocalization(self, local, newLocal):
        """
        """
        if local == newLocal: return

        Model().isNotInList(newLocal, self.getLocalizationsZonesList())


    def setLocalization(self, label, localization):
        """
        Define a new localization for the current zone (zone.getLabel == label).
        """
        labels = self.getLabelsZonesList()
        Model().isInList(label, labels)


    def setNature(self, label, nature):
        """
        Define a new nature number for the current zone (zone.getLabel == label)
        """
        # Set nature: nothing here, see other setNature reimplementation methods
        pass


    def addZone(self, newZone = None):
        """
        Add a new zone. Management of default values.
        """
        if newZone == None:
            newZone = Zone(self._typeZone, case = self.case)

        zones = self.getZones()

        # Set localization

        if newZone.getLocalization() == newZone.defaultValues()['localization']:
            oldLocalization = ""
            newLocalization = ""
            for zone in zones:
                oldLocalization = zone.getLocalization()
                if oldLocalization != "all[]":
                    if newLocalization == "":
                        newLocalization = oldLocalization
                    else:
                        newLocalization += " or " + oldLocalization

            if newLocalization == "":
                newLocalization = newZone.defaultValues()['localization']
            if newLocalization != "all[]":
                newLocalization = "not (" + newLocalization + ")"

            newZone.setLocalization(newLocalization)
        else:
            # No control on localization is availiable
            pass

        # Set code number

        if newZone.getCodeNumber() == newZone.defaultValues()['codeNumber']:
            newZone.setCodeNumber(self.getMaxCodeNumber() + 1)
        else:
            codes = []
            for zone in zones:
                codes.append(zone.getCodeNumber())
            Model().isNotInList(newZone.getCodeNumber(), codes)

        # Set label: search a free label (Type of newLabel is: "default_n")

        newLabel = newZone.getLabel()
        if newLabel == newZone.defaultValues()['label']:
            code = 1
            inLabels = 1
            while (inLabels):
                inLabels = 0
                for zone in zones:
                    if newLabel+str(code) == zone.getLabel():
                        inLabels = 1
                        break
                code += 1
            newLabel = newLabel + str(code-1)

            newZone.setLabel(newLabel)
        else:
            labels = []
            for zone in zones:
                labels.append(zone.getLabel())
            Model().isNotInList(newLabel, labels)

        # Set nature: nothing here, see other addZone reimplementation methods

        return newZone


    def replaceZone(self, old_zone, new_zone):
        """
        Replace a zone by another in the XML file
        """
        newLabel = new_zone.getLabel()
        if newLabel == new_zone.defaultValues()['label']:
            newLabel = old_zone.getLabel()
        self.renameLabel(old_zone.getLabel(), newLabel)

        newCodeNumber = new_zone.getCodeNumber()
        if newCodeNumber == new_zone.defaultValues()['codeNumber']:
            newCodeNumber = old_zone.getCodeNumber()
        self.renameName(old_zone.getCodeNumber(), newCodeNumber)

        newLocal = new_zone.getLocalization()
        if newLocal == new_zone.defaultValues()['localization']:
            newLocal = old_zone.getLocalization()
        self.replaceLocalization(old_zone.getLocalization(), newLocal)

        return newLabel, newCodeNumber, newLocal


    def deleteZone(self, label):
        """
        Delete the current zone (zone.getLabel == label)
        """
        pass


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class BoundaryLocalizationModel(LocalizationModel):
    """
    """
    def __new__(cls, case):
        """
        Constructor
        """
        return object.__new__(cls)


    def _initModel(self):
        """
        Initialize mode
        """
        self.__XMLBoundaryConditionsNode = self.case.xmlInitNode('boundary_conditions')
        self.__natureList = Zone('BoundaryZone', case = self.case).getNatureList()


    @Variables.noUndo
    def getZones(self):
        """
        Get zones in the XML file
        """
        XMLZonesNodes = self.__XMLBoundaryConditionsNode.xmlGetChildNodeList('boundary', 'label', 'name', 'nature')
        #
        # XML file reading
        zones = []
        for node in XMLZonesNodes:
            label = str(node['label'])
            nature = str(node['nature'])
            codeNumber = int(node['name'])
            localization = str(node.xmlGetTextNode())
            zone = Zone('BoundaryZone', case = self.case, label = label, codeNumber = codeNumber, localization = localization, nature = nature)
            zones.append(zone)
        return zones


    @Variables.noUndo
    def getMaxNumberNature(self, nature):
        """
        Return maximum of nature number's values to put on name
        """

        XMLZonesNodes = self.__XMLBoundaryConditionsNode.xmlGetChildNodeList('boundary', 'label', 'name', 'nature')
        max = 0
        #
        # XML file reading
        zones = []
        for node in XMLZonesNodes:
            if node['nature'] == nature:
                max = max + 1
        return max


    @Variables.undoLocal
    def setLabel(self, label, newLabel):
        """
        Define a new label for the current zone (zone.getLabel == label)
        Update XML file
        """
        LocalizationModel.renameLabel(self, label, newLabel)
        #
        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlGetChildNode('boundary', 'name', 'nature', label = label )
        node['label'] = newLabel
        nature = node['nature']
        #
        # Update references
        XMLZonesNodes = self.__XMLBoundaryConditionsNode.xmlGetChildNodeList(nature, label = label)
        for node in XMLZonesNodes:
            node['label'] = newLabel


    @Variables.undoLocal
    def setLocalization(self, label, localization):
        """
        Define a new localization for the current zone (zone.getLabel == label)
        Update XML file
        """
        LocalizationModel.setLocalization(self, label, localization)
        #
        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlGetChildNode('boundary', 'name', 'nature', label = label)
        node.xmlSetTextNode(localization)


    @Variables.undoLocal
    def setCodeNumber(self, label, codeNumber):
        """
        Define a new code number for the current zone (zone.getLabel == label)
        Update XML file
        """
        LocalizationModel.setCodeNumber(self, label, codeNumber)
        #
        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlGetChildNode('boundary', 'name', 'nature', label = label)
        node['name'] = str(codeNumber)


    @Variables.undoGlobal
    def setNature(self, label, nature):
        """
        Define a new Nature for the current zone (zone.getLabel == label)
        Update XML file
        """
        LocalizationModel.setNature(self, label, nature)

        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlGetChildNode('boundary', 'name', 'nature', label = label)
        oldNature = node['nature']
        node['nature'] = str(nature)

        # Delete oldNature boundary
        Boundary(oldNature, label, self.case).delete()

        # Create nature boundary
        Boundary(nature, label, self.case)


    @Variables.undoGlobal
    def addZone(self, zone = None):
        """
        Add a new zone in the XML file
        """
        newZone = LocalizationModel.addZone(self, zone)

        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlInitNode('boundary',
                                                            label = newZone.getLabel(),
                                                            name = str(newZone.getCodeNumber()),
                                                            nature = newZone.getNature())
        node.xmlSetTextNode(newZone.getLocalization())

        # Create nature boundary
        Boundary(newZone.getNature(), newZone.getLabel(), self.case)

        return newZone


    @Variables.undoGlobal
    def replaceZone(self, old_zone, new_zone):
        """
        Replace a zone by another in the XML file
        """
        # if nature of zone change we delete old zone
        if (new_zone.getNature() != old_zone.getNature()):
            Boundary(old_zone.getNature(), old_zone.getLabel(), self.case).delete()

            newLabel, newCodeNumber, newLocal = LocalizationModel.replaceZone(self, old_zone, new_zone)

            newNature = new_zone.getNature()
            Model().isInList(newNature, self.__natureList)

            node = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary',
                                                               label = old_zone.getLabel())

            node['label'] = newLabel
            node['name'] = newCodeNumber
            node['nature'] = newNature
            node.xmlSetTextNode(newLocal)

            if (new_zone.getNature() != old_zone.getNature()):
                Boundary(new_zone.getNature(), new_zone.getLabel(), self.case)
        else:
            label      = old_zone.getLabel()
            codeNumber = old_zone.getCodeNumber()
            localis    = old_zone.getLocalization()

            newLabel = new_zone.getLabel()
            if label != newLabel:
                self.setLabel(label, newLabel)

            newCodeNumber = new_zone.getCodeNumber()
            if codeNumber != newCodeNumber:
                self.setCodeNumber(codeNumber, newCodeNumber)

            newLocal = new_zone.getLocalization()
            if localis != newLocal:
                self.setLocalization(newLabel, newLocal)


    @Variables.undoGlobal
    def deleteZone(self, label):
        """
        Delete a zone in the XML file
        """
        LocalizationModel.deleteZone(self, label)

        # Get Nature
        node = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary', 'name', 'nature', label = label)
        nature = node['nature']
        node.xmlRemoveNode()

        # Delete nature boundary
        Boundary(nature, label, self.case).delete()


    def getXMLBoundaryConditionsNode(self):

        return self.__XMLBoundaryConditionsNode


#-------------------------------------------------------------------------------
# DefineUsersScalars test case
#-------------------------------------------------------------------------------
class LocalizationTestCase(ModelTest):
    """
    """
    def checkLocalizationInstantiation(self):
        """Check whether the LocalizationModel class could be instantiated"""
        model = None
        model = LocalizationModel("BoundaryZone",self.case)
        assert model != None, 'Could not instantiate LocalizationModel'


    def checkAddAndDeleteZone(self):
        """Check whether the zone could be added and deleted"""
        model = LocalizationModel("BoundaryZone",self.case)
        model.addZone()
        doc = '''<boundary_conditions>
                         <boundary label="BC_1" name="1" nature="wall">
                                 all[]
                         </boundary>
                         <wall field_id="none" label="BC_1"/>
                 </boundary_conditions>'''
        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not add  zone'
        model.deleteZone(label='BC_1')
        doc = '''<boundary_conditions/>'''

        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not delete zone'


    def checkReplaceZone(self):
        """Check whether the zone could be replaced for boundary conditions."""
        model = LocalizationModel("BoundaryZone", self.case)
        zone1 = Zone("BoundaryZone", label='entre1', localization="porte", nature='inlet')
        model.addZone(zone1)
        zone4 = Zone("BoundaryZone", label='hublot', localization="2 et 3", nature='symmetry')
        model.replaceZone(zone1, zone4)
        doc = '''<boundary_conditions>
                         <boundary label="hublot" name="1" nature="symmetry">
                                 2 et 3
                         </boundary>
                 </boundary_conditions>'''

        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not replace zone in localizationModel for boundaries conditions'


    def checkSetLabelAndLocalization(self):
        """Check whether the zone could be replaced for boundary conditions."""
        model = LocalizationModel("BoundaryZone", self.case)
        zone1 = Zone("BoundaryZone", label='entre1', localization="porte", nature='inlet')
        model.addZone(zone1)
        model.setLabel('entre1','entry1')
        doc = '''<boundary_conditions>
                         <boundary label="entry1" name="1" nature="inlet">
                                 porte
                         </boundary>
                 </boundary_conditions>'''


        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not set Label'
        model.setLocalization('entry1','door')
        doc = '''<boundary_conditions>
                         <boundary label="entry1" name="1" nature="inlet">
                                 door
                         </boundary>
                 </boundary_conditions>'''
        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not set Localization'


    def checkSetNature(self):
        """Check whether the zone could be replaced for boundary conditions."""
        model = LocalizationModel("BoundaryZone", self.case)
        zone1 = Zone("BoundaryZone", label='entre1', localization="porte", nature='inlet')
        model.addZone(zone1)
        model.setNature('entre1','wall')
        doc = '''<boundary_conditions>
                         <boundary label="entre1" name="1" nature="wall">
                                 porte
                         </boundary>
                         <wall field_id="none" label="entre1"/>
                 </boundary_conditions>'''
        assert self.case.xmlGetNode('boundary_conditions') == self.xmlNodeFromString(doc),\
           'Could not set Nature'


def suite1():
    testSuite = unittest.makeSuite(LocalizationTestCase, "check")
    return testSuite


def runTest1():
    print("LocalizationTestCase")
    runner = unittest.TextTestRunner()
    runner.run(suite1())
