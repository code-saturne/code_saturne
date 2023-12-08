# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
- Zone
- BoundaryZone
- VolumicZone
- LocalizationModel
- VolumicLocalizationModel
- BoundaryLocalizationModel
- LocalizationVolumicTestCase
- LocalizationSurfacicTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest, types

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.XMLengine import *

try:
    from code_saturne.model.BoundaryNeptune import Boundary as BoundaryNCFD
except:
    pass
from code_saturne.model.Boundary import Boundary

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Zone(object):
    """
    Zone API
    """
    def __new__(cls, typeZone, case = None, label = None, codeNumber = None,
                localization = None, nature = None):
        """
        Factory
        """
        if typeZone == 'BoundaryZone':
            return BoundaryZone.__new__(BoundaryZone, label, codeNumber, localization, nature)
        elif typeZone == 'VolumicZone':
            return VolumicZone.__new__(VolumicZone, label, codeNumber, localization)
        else:
            raise ValueError("Unknown type zone")

    def __init__(self, typeZone, case = None, label = None, codeNumber = None,
                 localization = None, nature = None):
        """
        """
        self.case = case
        self._initNatureList(self.case)

        if label:
            self._label = label
        else:
            self._label = self.defaultValues()['label']
        if codeNumber:
            self._codeNumber = codeNumber
        else:
            self._codeNumber = self.defaultValues()['codeNumber']
        if localization:
            self._localization = localization
        else:
            self._localization = self.defaultValues()['localization']
        if nature:
            if typeZone == 'VolumicZone' and type(nature) == str:
                self._nature = self.defaultValues()['nature'].copy()
                for n in nature.split(':'):
                    self._nature[n] = "on"
            else:
                self._nature = nature
        else:
            self._nature = self.defaultValues()['nature']

    def _initNatureList(self, case):
        self._natureList = []
        self._natureDict = {}
        self.case = case

    def setLabel(self, text):
        if Model().isStr(text):
            self._label = text

    def getLabel(self):
        return self._label

    def setLocalization(self, text):
        if Model().isStr(text):
            self._localization = text

    def getLocalization(self):
        return self._localization

    def setCodeNumber(self, number):
        if Model().isPositiveInt(number):
            self._codeNumber = number

    def getCodeNumber(self):
        return self._codeNumber

    def setNature(self, text):
        if Model().isInList(text, self._natureList):
            self._nature = text

    def getNature(self):
        return self._nature

    def getNatureList(self):
        return self._natureList

    def isNatureActivated(self, text):
        return {"on": True, "off": False}[self._nature.get(text, "off")]

    def getModel2ViewDictionary(self):
        return self._natureDict

    def defaultValues(self):
        dico = {}
        dico['codeNumber'] = -1
        dico['localization'] = "all[]"
        return dico

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class BoundaryZone(Zone):
    """
    """
    def __new__(cls, label = None, codeNumber = None, localization = None, nature = None):
        """
        Constructor
        """
        return object.__new__(cls)


    def _initNatureList(self, case):

        self.case = case

        self._natureDict = {}
        self._natureDict['wall']              = self.tr("Wall")
        self._natureDict['inlet']             = self.tr("Inlet")
        self._natureDict['outlet']            = self.tr("Outlet")
        self._natureDict['symmetry']          = self.tr("Symmetry")

        if case != None:
            if self.case.module_name() == 'code_saturne':
                self._natureList = ['wall', 'inlet', 'outlet', 'symmetry',
                                    'free_inlet_outlet', 'imposed_p_outlet']
                self._natureDict['free_inlet_outlet'] = self.tr("Free inlet/outlet")
                self._natureDict['imposed_p_outlet'] = self.tr("Imposed P Outlet")
                from code_saturne.model.MobileMeshModel import MobileMeshModel
                if MobileMeshModel(self.case).getMethod() != "off":
                    self._natureDict['free_surface'] = self.tr("Free surface")
                    self._natureList = ['wall', 'inlet', 'outlet', 'symmetry',
                                        'free_inlet_outlet', 'free_surface',
                                        'imposed_p_outlet']
                del MobileMeshModel
                from code_saturne.model.GroundwaterModel import GroundwaterModel
                if GroundwaterModel(self.case).getGroundwaterModel() != "off":
                    self._natureDict['groundwater'] = self.tr("Groundwater flow")
                    self._natureList = ['wall', 'inlet', 'outlet', 'symmetry',
                                        'free_inlet_outlet', 'groundwater']
                del GroundwaterModel
                from code_saturne.model.LagrangianModel import LagrangianModel
            else:
                self._natureList = ['wall', 'inlet', 'outlet', 'symmetry']
        else:
            self._natureList = ['wall', 'inlet', 'outlet', 'symmetry',
                                'free_inlet_outlet',  'imposed_p_outlet']
            self._natureDict['free_inlet_outlet'] = self.tr("Free inlet/outlet")
            self._natureDict['imposed_p_outlet'] = self.tr("Imposed P Outlet")


    def defaultValues(self):
        """
        """
        dico = Zone.defaultValues(self)
        dico['label'] = 'BC_'
        dico['nature'] = self._natureList[0]
        return dico


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class VolumicZone(Zone):
    """
    """
    def __new__(cls, label = None, codeNumber = None, localization = None, nature = None):
        return object.__new__(cls)


    def _initNatureList(self, case):
        self._natureList = ['initialization']
        self.case = case

        self._natureDict = {}
        self._natureDict['initialization'] = self.tr("Initialization")
        self._natureDict['physical_properties'] = self.tr("Physical properties")
        self._natureList.append('physical_properties')

        if self.case.module_name() == 'code_saturne':
            from code_saturne.model.GroundwaterModel import GroundwaterModel
            if GroundwaterModel(self.case).getGroundwaterModel() != "groundwater":
                self._natureList.append('solid')
                self._natureDict['solid'] = self.tr("Solid")

                self._natureList.append('head_losses')
                self._natureList.append('porosity')
                self._natureList.append('momentum_source_term')

                self._natureDict['head_losses'] = self.tr("Head losses")
                self._natureDict['porosity'] = self.tr("Porosity")
                self._natureDict['momentum_source_term'] = self.tr("Momentum source\n term")
            else:
                self._natureList.append('momentum_source_term')
                self._natureDict['momentum_source_term'] = self.tr("Volumic source\n term")
            del GroundwaterModel

            from code_saturne.model.ThermalScalarModel import ThermalScalarModel
            if ThermalScalarModel(self.case).getThermalScalarModel() != 'off':
                self._natureList.append('thermal_source_term')
                self._natureDict['thermal_source_term']  = self.tr("Thermal source term")
            del ThermalScalarModel

            self.node_models = self.case.xmlGetNode('thermophysical_models')
            node_darcy = self.node_models.xmlGetNode('groundwater_model')
            if node_darcy:
                if node_darcy['model'] != 'off':
                    self._natureList.append('groundwater_law')
                    self._natureDict['groundwater_law']  = self.tr("Groundwater\n volumic law")

        elif self.case.module_name() == 'neptune_cfd':
            self._natureList.append('head_losses')
            self._natureList.append('porosity')

            self._natureDict['head_losses']          = self.tr("Head losses")
            self._natureDict['porosity']             = self.tr("Porosity")
            # Thermal source term for enthalpy
            # TODO see how we can remove the import of MainFieldsModel
            from code_saturne.model.MainFieldsModel import MainFieldsModel
            if len(MainFieldsModel(self.case).getFieldIdList()) > 0:
                if MainFieldsModel(self.case).getFieldFromId(1).enthalpy_model != "off":
                    self._natureList.append('thermal_source_term')
                    self._natureDict['thermal_source_term'] = self.tr("Thermal source term")
            del MainFieldsModel


        node = self.case.xmlGetNode('additional_scalars')
        number = len(node.xmlGetNodeList('variable', type='user'))
        if number > 0:
            self._natureList.append('scalar_source_term')
            self._natureDict['scalar_source_term']   = self.tr("Scalar source term")


    def defaultValues(self):
        dico = Zone.defaultValues(self)
        dico['label'] = 'Zone_'
        dico['nature'] = {}
        return dico

    def isNatureActivated(self, text):
        if text == "source_term":
            status = False
            for source_term_type in ["momentum", "thermal", "scalar"]:
                status = status or super().isNatureActivated(source_term_type + "_source_term")
            return status
        else:
            return super().isNatureActivated(text)

    def tr(self, text):
        """
        Translation
        """
        return text

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
        elif typeZone == 'VolumicZone':
            return VolumicZone.__new__(VolumicLocalizationModel, case)
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

    def getSortedZoneLabels(self):
        zone_labels = self.getLabelsZonesList()
        zone_ids = self.getCodeNumbersList()
        zone_ids = map(int, zone_ids)
        sorted_labels = [None for i in range(len(zone_labels))]
        for unsorted_id, sorted_id in enumerate(zone_ids):
            sorted_labels[sorted_id - 1] = zone_labels[unsorted_id]
        return sorted_labels

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

    def selectZone(self, value, criterium):
        """ Return first zone satisfying criterium """
        zones = self.getZones()
        for zone in zones:
            if criterium == "label":
                if zone.getLabel() == value:
                    return zone
            elif criterium == "codeNumber":
                if zone.getCodeNumber() == value:
                    return zone
            else:
                raise ValueError

    def addZone(self, newZone=None, checkPresence=True):
        """
        Add a new zone. Management of default values.
        """
        if newZone is None:
            newZone = Zone(self._typeZone, case=self.case)

        if newZone.getLocalization() == newZone.defaultValues()['localization'] or newZone.getLabel() == newZone.defaultValues()['label'] or checkPresence == True:
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
            # No control on localization is available
            pass

        # Set code number

        if newZone.getCodeNumber() == newZone.defaultValues()['codeNumber']:
            newZone.setCodeNumber(self.getMaxCodeNumber() + 1)
        elif checkPresence == True:
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
        elif checkPresence == True:
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


    def mergeZones(self, label, localization, lst):
        """
        Merge zones (zone.getLabel == label)
        lst : list of zone to merge
        """
        pass

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class VolumicLocalizationModel(LocalizationModel):
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
        XMLSolutionDomainNode = self.case.xmlInitNode('solution_domain')
        self.__XMLVolumicConditionsNode = XMLSolutionDomainNode.xmlInitNode('volumic_conditions')
        self.__natureOptions = Zone('VolumicZone', case = self.case).getNatureList()
        self._tagList = ['formula', 'head_loss']
        self.node_models = self.case.xmlGetNode('thermophysical_models')
        self.node_veloce = self.node_models.xmlGetNode('velocity_pressure')
        self.scalar_node = self.case.xmlGetNode('additional_scalars')
        self.losses_node = self.case.xmlGetNode('head_losses')


    @Variables.noUndo
    def getZones(self):
        """
        Get zones in the XML file
        """
        XMLZonesNodes = self.__XMLVolumicConditionsNode.xmlGetChildNodeList('zone', 'label', 'id')

        # XML file reading
        zones = []
        for node in XMLZonesNodes:
            label = str(node['label'])
            codeNumber = int(node['id'])
            localization = str(node.xmlGetTextNode())
            nature = self.getNature(label)
            zone = Zone('VolumicZone',
                        case=self.case,
                        label=label,
                        codeNumber=codeNumber,
                        localization=localization,
                        nature=nature)
            zones.append(zone)
        return zones

    @Variables.noUndo
    def getCodeNumberOfZoneLabel(self, label):
        """
        Get zones in the XML file
        """
        XMLZonesNodes = self.__XMLVolumicConditionsNode.xmlGetChildNodeList('zone', 'label', 'id')

        # XML file reading
        for node in XMLZonesNodes:
            if node['label'] == label:
                codeNumber = node['id']

        return codeNumber


    @Variables.undoLocal
    def setLocalization(self, label, localization):
        """
        Define a new localization for the current zone (zone.getLabel == label)
        Update XML file
        """
        node = self.__XMLVolumicConditionsNode.xmlGetChildNode('zone', 'id', label = label)
        node.xmlSetTextNode(localization)


    @Variables.noUndo
    def getCodeNumbersList(self, codeNumber=None):
        """
        Define a new code number for the current zone (zone.getLabel == label)
        Update XML file
        """
        XMLZonesNodesList = self.case.xmlGetNodeList('zone', 'label', 'id')
        codeList = []
        for node in XMLZonesNodesList:
            codeList.append(node['id'])
        return codeList


    @Variables.noUndo
    def getNature(self, label):
        """
        Define a new Nature for the current zone (zone.getLabel == label)
        Update XML file
        """
        node = self.__XMLVolumicConditionsNode.xmlGetChildNode('zone', 'id', label = label)
        nature = {}
        for option in self.__natureOptions:
            if option in node.xmlGetAttributeDictionary():
                if node[option] == 'on':
                    nature[option] = 'on'
                else:
                    nature[option] = 'off'
            else:
                nature[option] = 'off'
        return nature


    @Variables.undoGlobal
    def setNature(self, label, nature):
        """
        Define a new Nature for the current zone (zone.getLabel == label)
        Update XML file
        """
        node = self.__XMLVolumicConditionsNode.xmlGetChildNode('zone', 'id', label = label)
        oldNature = self.getNature(label)
        if oldNature != nature:
            for option in self.__natureOptions:
                if option not in list(nature.keys()):
                    nature[option] = 'off'
            for k,v in list(nature.items()):
                node[k] = v
                if node[k] == 'off':
                    del(node[k])


    @Variables.undoGlobal
    def addZone(self, zone=None, checkPresence=True):
        """
        Add a new zone in the XML file
        """
        newZone = LocalizationModel.addZone(self, zone, checkPresence)

        # XML file updating
        node = self.__XMLVolumicConditionsNode.xmlInitNode('zone',
                                                           label = newZone.getLabel(),
                                                           id = newZone.getCodeNumber())

        for k, v in list(newZone.getNature().items()):
            node[k] = v
            if node[k] == 'off':
                del(node[k])

        node.xmlSetTextNode(newZone.getLocalization())

        return newZone


    @Variables.undoGlobal
    def replaceZone(self, old_zone, new_zone):
        """
        Replace a zone by another in the XML file
        """
        newLabel, newCodeNumber, newLocal = LocalizationModel.replaceZone(self, old_zone, new_zone)

        node = self.__XMLVolumicConditionsNode.xmlGetNode('zone',
                                                          label = old_zone.getLabel())

        # if codeNumber is modified, we must modify zone in initialization of variables
        if old_zone.getCodeNumber() != newCodeNumber:
            node['id'] = newCodeNumber

        node['label'] = newLabel
        node.xmlSetTextNode(newLocal)
        for k, v in list(new_zone.getNature().items()):
            node[k] = v

        for k in node.xmlGetAttributeDictionary():
            if node[k] == 'off':
                del(node[k])

        # update data in the entire case
        lst = self.__natureOptions
        lst.append('initial_value')
        lst.append('head_losses')
        for tag in lst:
            for n in self.case.xmlGetNodeList(tag, zone=old_zone.getCodeNumber()):
                n['zone'] = newCodeNumber
            for n in self.case.xmlGetNodeList(tag, id=old_zone.getCodeNumber()):
                n['zone_id'] = newCodeNumber
            for n in self.case.xmlGetNodeList(tag, label=old_zone.getLabel()):
                n['label'] = newLabel


    @Variables.undoGlobal
    def deleteZone(self, label):
        """
        Delete one zone in the XML file
        """
        LocalizationModel.deleteZone(self, label)
        #
        # Delete node
        node = self.case.xmlGetNode('zone', label=label)
        if node:
            name = node['id']
            node.xmlRemoveNode()

            # Delete the other nodes for zone initializations
            n_d = self.case.xmlGetNodeWithAttrList('zone_id', zone_id=name)
            for n in n_d:
                n.xmlRemoveNode()

            self.renumberZones(name)


    @Variables.undoGlobal
    def mergeZones(self, label, localization, lst):
        """
        Merge zones in the XML file
        """
        LocalizationModel.mergeZones(self, label, localization, lst)

        node = self.__XMLVolumicConditionsNode.xmlGetNode('zone', 'id', label = label)
        node.xmlSetTextNode(localization)

        lst.reverse()

        for z in lst:
            n = self.__XMLVolumicConditionsNode.xmlGetNode('zone', 'label', id = str(z + 1))
            n.xmlRemoveNode()

            for tag in self._tagList:
                nodeList = self.case.xmlGetNodeList(tag, zone_id = str(z + 1))
                for node in nodeList:
                    node.xmlRemoveNode()

        self.renumberZones()


    @Variables.undoGlobal
    def renumberZones(self, start_id='0'):
        """
        Renumber zones in the XML file and update dependent id's
        """

        volCondNode = self.__XMLVolumicConditionsNode
        XMLZonesNodes = volCondNode.xmlGetChildNodeList('zone', 'label', 'id')

        n_d = self.case.xmlGetNodeWithAttrList('zone_id')

        count = 1

        l = []
        for node in XMLZonesNodes:
            l.append(int(node['id']))
        l.sort()

        for zoneid in l:
            if zoneid > int(start_id):
                node = volCondNode.xmlGetNode('zone', id=zoneid)

                new_id = str(count)
                old_id = node['id']
                node['id'] = new_id

                for n in n_d:
                    if n['zone_id'] == old_id:
                        n['zone_id'] = new_id

            count = count + 1


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
        #LocalizationModel._initModel(self)
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
        node = self.__XMLBoundaryConditionsNode.xmlGetChildNode('boundary', 'name', 'nature',
                                                                label = label)
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

        if self.case.module_name() == 'code_saturne':
            # Delete oldNature boundary
            Boundary(oldNature, label, self.case).delete()
            # Create nature boundary
            Boundary(nature, label, self.case)
        elif self.case.module_name() == 'neptune_cfd':
            # Delete oldNature boundary
            BoundaryNCFD(oldNature, label, self.case).delete()
            # Create nature boundary
            BoundaryNCFD(nature, label, self.case)


    @Variables.undoGlobal
    def addZone(self, zone=None, checkPresence=True):
        """
        Add a new zone in the XML file
        """
        newZone = LocalizationModel.addZone(self, zone, checkPresence)

        # XML file updating
        node = self.__XMLBoundaryConditionsNode.xmlInitNode('boundary',
                                                            label = newZone.getLabel(),
                                                            name = str(newZone.getCodeNumber()),
                                                            nature = newZone.getNature())
        node.xmlSetTextNode(newZone.getLocalization())

        # Create nature boundary
        if self.case.module_name() == 'code_saturne':
            Boundary(newZone.getNature(), newZone.getLabel(), self.case)
        elif self.case.module_name() == 'neptune_cfd':
            BoundaryNCFD(newZone.getNature(), newZone.getLabel(), self.case)

        return newZone


    @Variables.undoGlobal
    def replaceZone(self, old_zone, new_zone):
        """
        Replace a zone by another in the XML file
        """
        if (new_zone.getNature() != old_zone.getNature()):
            if self.case.module_name() == 'code_saturne':
                Boundary(old_zone.getNature(), old_zone.getLabel(), self.case).delete()
            elif self.case.module_name() == 'neptune_cfd':
                BoundaryNCFD(old_zone.getNature(), old_zone.getLabel(), self.case).delete()
            newLabel, newCodeNumber, newLocal = LocalizationModel.replaceZone(self, old_zone, new_zone)

            newNature = new_zone.getNature()
            Model().isInList(newNature, self.__natureList)

            node = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary',
                                                            label = old_zone.getLabel())

            node['label'] = newLabel
            node['name'] = newCodeNumber
            node['nature'] = newNature
            node.xmlSetTextNode(newLocal)

            if self.case.module_name() == 'code_saturne':
                Boundary(new_zone.getNature(), new_zone.getLabel(), self.case)
            elif self.case.module_name() == 'neptune_cfd':
                BoundaryNCFD(new_zone.getNature(), new_zone.getLabel(), self.case)
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
        node = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary', 'name', 'nature',
                                                           label = label)
        nature = node['nature']
        node.xmlRemoveNode()

        # Delete nature boundary
        if self.case.module_name() == 'code_saturne':
            Boundary(nature, label, self.case).delete()
        elif self.case.module_name() == 'neptune_cfd':
            BoundaryNCFD(nature, label, self.case).delete()

        self.renumberZones()


    @Variables.undoGlobal
    def mergeZones(self, label, localization, lst):
        """
        Merge zones in the XML file
        """
        LocalizationModel.mergeZones(self, label, localization, lst)

        node = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary', 'name', 'nature',
                                                           label = label)
        node.xmlSetTextNode(localization)

        lst.reverse()

        for z in lst:
            n = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary', 'nature', 'label',
                                                            name = z + 1)
            label = n['label']
            nature = n['nature']
            if self.case.module_name() == 'code_saturne':
                Boundary(nature, label, self.case).delete()
            elif self.case.module_name() == 'neptune_cfd':
                BoundaryNCFD(nature, label, self.case).delete()
            n.xmlRemoveNode()

        self.renumberZones()


    @Variables.undoGlobal
    def renumberZones(self):
        """
        Merge zones in the XML file
        """

        count = 1
        l = self.getCodeNumbersList()
        l.sort()
        for z in l:
            n = self.__XMLBoundaryConditionsNode.xmlGetNode('boundary', 'nature', 'label',
                                                            name = z)
            n['name'] = str(count)
            nature = n['nature']
            count = count + 1


#-------------------------------------------------------------------------------
# LocalizationModel test case for volumic zones
#-------------------------------------------------------------------------------

class LocalizationVolumicTestCase(ModelTest):
    """
    Unittest.
    """
    def checkLocalizationInstantiation(self):
        """Check whether the LocalizationModel class could be instantiated."""
        model = None
        model = LocalizationModel("VolumicZone", self.case)
        assert model != None, 'Could not instantiate LocalizationVolumicModel'


def suite1():
    testSuite = unittest.makeSuite(LocalizationVolumicTestCase, "check")
    return testSuite


def runTest1():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite1())

#-------------------------------------------------------------------------------
# LocalizationModel test case for boundary conditions
#-------------------------------------------------------------------------------

class LocalizationSurfacicTestCase(ModelTest):
    """
    Unittest.
    """
    def checkLocalizationSurfacicInstantiation(self):
        """
        Check whether the LocalizationModel class could be instantiated
        for boundary conditions.
        """
        model = None
        model = LocalizationModel("BoundaryZone", self.case)
        assert model != None, 'Could not instantiate LocalizationSurfacicModel'


    def checkAddAndDeleteZone(self):
        """Check whether the zone could be added and deleted for boundary conditions."""
        model = LocalizationModel("BoundaryZone", self.case)
        node = self.case.xmlGetNode('boundary_conditions')
        zone1 = Zone("BoundaryZone", label='entre1', localization="porte", nature='inlet')
        zone2 = Zone("BoundaryZone", label='entre2', localization="fenetre", nature='inlet')
        zone3 = Zone("BoundaryZone", label='plafond', localization="not porte", nature='wall')
        model.addZone(zone1)
        model.addZone(zone2)
        model.addZone(zone3)

        doc = '''<boundary_conditions>
                        <boundary label="entre1" name="1" nature="inlet">
                                porte
                        </boundary>
                        <inlet label="entre1">
                            <velocity_pressure choice="norm">
                            <norm>1</norm>
                            </velocity_pressure>
                            <turbulence choice="hydraulic_diameter">
                                <hydraulic_diameter>1</hydraulic_diameter>
                            </turbulence>
                        </inlet>
                        <boundary label="entre2" name="2" nature="inlet">
                                fenetre
                        </boundary>
                        <inlet label="entre2">
                            <velocity_pressure choice="norm">
                                <norm>1</norm>
                            </velocity_pressure>
                            <turbulence choice="hydraulic_diameter">
                                <hydraulic_diameter>1</hydraulic_diameter>
                            </turbulence>
                        </inlet>
                        <boundary label="plafond" name="3" nature="wall">
                                not porte
                        </boundary>
                        <wall label="plafond">
                            <velocity_pressure choice="off"/>
                        </wall>
                  </boundary_conditions>'''

        model.deleteZone("entre2")

        doc = '''<boundary_conditions>
                        <boundary label="entre1" name="1" nature="inlet">
                                porte
                        </boundary>
                        <inlet label="entre1">
                            <velocity_pressure choice="norm">
                            <norm>1</norm>
                            </velocity_pressure>
                            <turbulence choice="hydraulic_diameter">
                                <hydraulic_diameter>1</hydraulic_diameter>
                            </turbulence>
                        </inlet>
                        <boundary label="plafond" name="3" nature="wall">
                                not porte
                        </boundary>
                        <wall label="plafond">
                            <velocity_pressure choice="off"/>
                        </wall>
                  </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not delete zone in localizationModel for boundaries conditions'


    def checkReplaceZone(self):
        """Check whether the zone could be replaced for boundary conditions."""
        model = LocalizationModel("BoundaryZone", self.case)
        node = self.case.xmlGetNode('boundary_conditions')
        zone1 = Zone("BoundaryZone", label='entre1', localization="porte", nature='inlet')
        zone2 = Zone("BoundaryZone", label='entre2', localization="fenetre", nature='inlet')
        zone3 = Zone("BoundaryZone", label='plafond', localization="not porte", nature='wall')
        model.addZone(zone1)
        model.addZone(zone2)
        model.addZone(zone3)
        zone4 = Zone("BoundaryZone", label='hublot', localization="2 et 3", nature='symmetry')
        model.replaceZone(zone2, zone4)

        doc = '''<boundary_conditions>
                        <boundary label="entre1" name="1" nature="inlet">
                                porte
                        </boundary>
                        <inlet label="entre1">
                            <velocity_pressure choice="norm">
                            <norm>1</norm>
                            </velocity_pressure>
                            <turbulence choice="hydraulic_diameter">
                                <hydraulic_diameter>1</hydraulic_diameter>
                            </turbulence>
                        </inlet>
                        <boundary label="hublot" name="2" nature="symmetry">
                                2 et 3
                        </boundary>
                        <boundary label="plafond" name="3" nature="wall">
                                not porte
                        </boundary>
                        <wall label="plafond">
                            <velocity_pressure choice="off"/>
                        </wall>
                        <symmetry label="hublot"/>
                  </boundary_conditions>'''

        assert node == self.xmlNodeFromString(doc),\
           'Could not replace zone in localizationModel for boundaries conditions'


def suite2():
    testSuite = unittest.makeSuite(LocalizationSurfacicTestCase, "check")
    return testSuite


def runTest2():
    print(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite2())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
