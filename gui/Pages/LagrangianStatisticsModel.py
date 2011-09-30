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
This module defines the lagrangian two phase flow modelling management.

This module contains the following classes and function:
- LagrangianStatisticsModel
- LagrangianStatisticsTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest, logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Pages.LagrangianModel import LagrangianModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianStatisticsModel")
log.setLevel(Tool.GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# lagrangian model class
#-------------------------------------------------------------------------------

class LagrangianStatisticsModel(Model):
    """
    Manage the input/output markups in the xml doc about Lagrangian module.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_lagr = self.case.root().xmlInitNode('lagrangian', 'model')
        self.node_stat = self.node_lagr.xmlInitChildNode('statistics')


    def _defaultLagrangianStatisticsValues(self):
        """
        Return a dictionnary which contains default values.
        """
        default = {}
        default['restart'] = "off"
        default['statistics_groups_of_particles'] = 0
        default['volume_statistics'] = "off"
        default['iteration_start_volume'] = 1
        default['threshold_volume'] = 0.

        for v in self._defaultVariablesNamesVolume():
            default[v] = v

        default['boundary_statistics'] = "off"
        default['iteration_start_boundary'] = 1
        default['threshold_boundary'] = 0.

        for v in self.getVariablesNamesBoundary():
            default[v] = v

        default['monitoring_point'] = "off"
        default['listing_printing'] = "off"
        default['postprocessing_recording'] = "off"

        return default


    def _defaultVariablesNamesVolume(self):
        names = self.getVariablesNamesVolume()
        volume_names = []
        for name in names:
            if name == "statistical_weight":
                volume_names.append(name)
            else:
                volume_names.append("mean_" + name)
                volume_names.append("variance_" + name)
        return volume_names


    # not private, used in View
    def getVariablesNamesVolume(self):

        names = ["statistical_weight",
                 "velocity_U", "velocity_V", "velocity_W",
                 "volume_fraction", "resident_time", "temperature",
                 "diameter", "shrinking_core_diameter",
                 "raw_coal_mass_fraction", "char_mass_fraction" ]
        return names


    # not private, used in View
    def getVariablesNamesBoundary(self):
        names = ["impacts", "mass_flux",
                 "angle", "velocity", "coal_fouling"]
        return names


    def setRestartStatisticsStatus(self, status):
        """
        Update the restart status markup from the XML document.
        """
        self.isOnOff(status)
        node_restart = self.node_stat.xmlInitChildNode('restart', 'status')
        node_restart['status'] = status


    def getRestartStatisticsStatus(self):
        """
        Return status of restart file.
        """
        node_restart = self.node_stat.xmlInitChildNode('restart', 'status')
        status = node_restart['status']
        if not status:
            status = self._defaultLagrangianStatisticsValues()['restart']
            self.setRestartStatisticsStatus(status)
        return status


    def setGroupOfParticlesValue(self, value):
        """
        Update the value of group of particles.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_stat.xmlSetData('statistics_groups_of_particles', value)


    def getGroupOfParticlesValue(self):
        """
        Return the value of group of particles.
        """
        npart = self.node_stat.xmlGetInt('statistics_groups_of_particles')
        if npart == None:
            npart = self._defaultLagrangianStatisticsValues()['statistics_groups_of_particles']
            self.setGroupOfParticlesValue(npart)
        return npart


    # Volume functions
    # ----------------
    def setVolumeStatisticsStatus(self, status):
        """
        """
        self.isOnOff(status)
        self.node_volume['status'] = status


    def getVolumeStatisticsStatus(self):
        """
        """
        self.node_volume = self.node_stat.xmlInitChildNode('volume', 'status')
        status = self.node_volume['status']
        if not status:
            status = self._defaultLagrangianStatisticsValues()['volume_statistics']
            self.setVolumeStatisticsStatus(status)
        return status


    def setIterationStartVolume(self, value):
        """
        Update the iteration value for start of volume statistics calculation.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_volume.xmlSetData('iteration_start_volume', value)


    def getIterationStartVolume(self):
        """
        Return the iteration value for start of volume statistics calculation.
        """
        value = self.node_volume.xmlGetInt('iteration_start_volume')
        if value == None:
            value = self._defaultLagrangianStatisticsValues()['iteration_start_volume']
            self.setIterationStartVolume(value)
        return value


    def setThresholdValueVolume(self, value):
        """
        Update the limit statistical weight value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_volume.xmlSetData('threshold_volume', value)


    def getThresholdValueVolume(self):
        """
        Return the limit statistical weight value.
        """
        value = self.node_volume.xmlGetDouble('threshold_volume')
        if not value:
            value = self._defaultLagrangianStatisticsValues()['threshold_volume']
            self.setThresholdValueVolume(value)
        return value


    def getPropertyLabelFromNameVolume(self, name):
        node = self.node_volume.xmlInitChildNode('property', name=name)
        label = node['label']
        if not label:
            label = self._defaultLagrangianStatisticsValues()[name]
            self.setPropertyLabelFromNameVolume(label, label)
        return label


    def setPropertyLabelFromNameVolume(self, name, label):
        node = self.node_volume.xmlInitChildNode('property', name=name)
        node['label'] = label


    def getMonitoringStatusFromName(self, name):
        node = self.node_volume.xmlInitChildNode('property', name=name)
        node2 = node.xmlGetChildNode('monitoring_point', 'status')
        if not node2:
            return "on"
        else:
            return "off" # node2['status']


    def setMonitoringStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_volume.xmlInitChildNode('property', name=name)
        node2 = node.xmlInitChildNode('monitoring_point', 'status')
        if status == "on":
            node.xmlRemoveChild('monitoring_point')
        elif status == "off":
            node2['status'] = status


    # Boundary functions
    # ------------------
    def setBoundaryStatisticsStatus(self, status):
        """
        """
        self.isOnOff(status)
        self.node_boundary['status'] = status


    def getBoundaryStatisticsStatus(self):
        """
        """
        self.node_boundary = self.node_stat.xmlInitChildNode('boundary', 'status')
        status = self.node_boundary['status']
        if not status:
            status = self._defaultLagrangianStatisticsValues()['boundary_statistics']
            self.setBoundaryStatisticsStatus(status)
        return status


    def setIterationStartBoundary(self, value):
        """
        Update iteration value for start of boundary statistics calculation.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_boundary.xmlSetData('iteration_start_boundary', value)


    def getIterationStartBoundary(self):
        """
        Return the iteration value for start of boundary statistics calculation.
        """
        value = self.node_boundary.xmlGetInt('iteration_start_boundary')
        if value == None:
            value = self._defaultLagrangianStatisticsValues()['iteration_start_boundary']
            self.setIterationStartBoundary(value)
        return value


    def setThresholdValueBoundary(self, value):
        """
        Update the limit statistical weight value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_boundary.xmlSetData('threshold_boundary', value)


    def getThresholdValueBoundary(self):
        """
        Return the limit statistical weight value.
        """
        value = self.node_boundary.xmlGetDouble('threshold_boundary')
        if not value:
            value = self._defaultLagrangianStatisticsValues()['threshold_boundary']
            self.setThresholdValueBoundary(value)
        return value


    def getPropertyLabelFromNameBoundary(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        label = node['label']
        if not label:
            label = self._defaultLagrangianStatisticsValues()[name]
            self.setPropertyLabelFromNameBoundary(label, label)
        return label


    def setPropertyLabelFromNameBoundary(self, name, label):
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        node['label'] = label


    def getListingPrintingStatusFromName(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        node2 = node.xmlGetChildNode('listing_printing', 'status')
        if not node2:
            return "on"
        else:
            return "off" # node2['status']


    def setListingPrintingStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        node2 = node.xmlInitChildNode('listing_printing', 'status')
        if status == "on":
            node.xmlRemoveChild('listing_printing')
        elif status == "off":
            node2['status'] = status


    def getPostprocessingStatusFromName(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        node2 = node.xmlGetChildNode('postprocessing_recording', 'status')
        if not node2:
            return "on"
        else:
            return "off" # node2['status']


    def setPostprocessingStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_boundary.xmlInitChildNode('property', name=name) #, support="boundary")
        node2 = node.xmlInitChildNode('postprocessing_recording', 'status')
        if status == "on":
            node.xmlRemoveChild('postprocessing_recording')
        elif status == "off":
            node2['status'] = status

#-------------------------------------------------------------------------------
# LagrangianStatistics test case
#-------------------------------------------------------------------------------

class LagrangianStatisticsTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        from Base.XMLengine import Case
        from Base.XMLinitialize import XMLinit
        self.case = Case()
        XMLinit(self.case)


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case


    def checkLagrangianStatisticsInstantiation(self):
        """
        Check whether the LagrangianStatisticsModel class could be instantiated
        """
        model = None
        model = LagrangianStatisticsModel(self.case)

        assert model != None, 'Could not instantiate LagrangianStatisticsModel'


    def checkLagrangianStatisticsDefaultValues(self):
        """
        Check the default values
        """
        model = LagrangianStatisticsModel(self.case)
        doc = """"""

        assert model.node_output == self.xmlNodeFromString(doc),\
               'Could not get default values for model'


#    def checkGetandSetLagrangianStatisticsModel(self):
#        """
#        Check whether the LagrangianStatisticsModel could be set/get
#        """
#        mdl = LagrangianStatisticsModel(self.case)
#        lagrangianStatus = mdl.lagrangianStatus()
#        assert lagrangianStatus == ('off','one_way','two_way','frozen')
#        for name in lagrangianStatus:
#            mdl.setLagrangianStatisticsModel(name)
#            name2 = mdl.getLagrangianStatisticsModel()
#            assert name == name2 ,\
#                   'Could not use the get/setLagrangianStatisticsModel method for model name %s '%name


    def checkSetandGetRestart(self):
        """
        Check whether the restart method could be set and get
        """
        mdl = LagrangianStatisticsModel(self.case)
        # default
        status = mdl.getRestart()
        assert status == 'off' ,\
        'Could not get default values for restart status'
        mdl.setRestart('on')
        doc = """
        """

        assert mdl.node_output == self.xmlNodeFromString(doc) ,\
            'Could not set values for restart status'


def suite():
    testSuite = unittest.makeSuite(LagrangianStatisticsTestCase, "check")
    return testSuite


def runTest():
    print("LagrangianStatisticsTestCase TODO*********.")
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
