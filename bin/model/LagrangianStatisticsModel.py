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

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.LagrangianModel import LagrangianModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianStatisticsModel")
log.setLevel(GuiParam.DEBUG)

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
        default['iteration_start'] = 1
        default['iteration_steady_start'] = default['iteration_start']
        default['threshold'] = 0.

        for v in self._defaultVariablesNamesVolume():
            default[v] = v

        default['boundary_statistics'] = "off"

        for v in self.getVariablesNamesBoundary():
            default[v] = v

        default['listing_printing'] = "off"
        default['postprocessing_recording'] = "off"

        return default


    def _defaultVariablesNamesVolume(self):
        names = self.getVariablesNamesVolume()
        volume_names = []
        for name in names:
            if name == "Part_statis_weight":
                volume_names.append(name)
            else:
                volume_names.append(name)
                volume_names.append("var_" + name)
        return volume_names


    @Variables.noUndo
    def getVariablesNamesVolume(self):

        names = ["Part_vol_frac", "Part_velocity",
                 "Part_resid_time", "Part_statis_weight"]
        return names


    @Variables.noUndo
    def getVariablesNamesBoundary(self):
        names = ["Part_bndy_mass_flux","Part_impact_number",
                 "Part_impact_angle", "Part_impact_velocity"]
        if LagrangianModel(self.case).getCoalFouling() == 'on':
            names = ["Part_bndy_mass_flux"      ,"Part_impact_number",
                     "Part_impact_angle"        , "Part_impact_velocity",
                     "Part_fouled_impact_number", "Part_fouled_mass_flux",
                     "Part_fouled_diam"         , "Part_fouled_Xck"]
        return names


    @Variables.undoLocal
    def setRestartStatisticsStatus(self, status):
        """
        Update the restart status markup from the XML document.
        """
        self.isOnOff(status)
        node_restart = self.node_stat.xmlInitChildNode('restart', 'status')
        node_restart['status'] = status


    @Variables.noUndo
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


    @Variables.undoLocal
    def setGroupOfParticlesValue(self, value):
        """
        Update the value of group of particles.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_stat.xmlSetData('statistics_groups_of_particles', value)


    @Variables.noUndo
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
    @Variables.undoLocal
    def setVolumeStatisticsStatus(self, status):
        """
        """
        self.isOnOff(status)
        self.node_volume['status'] = status


    @Variables.noUndo
    def getVolumeStatisticsStatus(self):
        """
        """
        self.node_volume = self.node_stat.xmlInitChildNode('volume', 'status')
        status = self.node_volume['status']
        if not status:
            status = self._defaultLagrangianStatisticsValues()['volume_statistics']
            self.setVolumeStatisticsStatus(status)
        return status


    @Variables.undoLocal
    def setIterationStart(self, value):
        """
        Update the iteration value for start of statistics calculation.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        self.node_stat.xmlSetData('iteration_start', value)


    @Variables.noUndo
    def getIterationStart(self):
        """
        Return the iteration value for start of volume statistics calculation.
        """
        value = self.node_stat.xmlGetInt('iteration_start')
        if value == None:
            value = self._defaultLagrangianStatisticsValues()['iteration_start']
            self.setIterationStart(value)
        return value


    @Variables.undoLocal
    def setIterSteadyStart(self, value):
        """
        Update the iteration value for start of steady statistics calculation.
        """
        self.isInt(value)
        self.isGreaterOrEqual(value,0)
        self.node_stat.xmlSetData('iteration_steady_start', value)


    @Variables.noUndo
    def getIterSteadyStart(self):
        """
        Return the iteration value for start of steady statistics calculation.
        """
        value = self.node_stat.xmlGetInt('iteration_steady_start')
        if value == None:
            value = self._defaultLagrangianStatisticsValues()['iteration_steady_start']
            self.setIterSteadyStart(value)
        return value


    @Variables.undoLocal
    def setThresholdValue(self, value):
        """
        Update the limit statistical weight value.
        """
        self.isFloat(value)
        self.isGreaterOrEqual(value, 0)
        self.node_stat.xmlSetData('threshold', value)


    @Variables.noUndo
    def getThresholdValue(self):
        """
        Return the limit statistical weight value.
        """
        value = self.node_stat.xmlGetDouble('threshold')
        if not value:
            value = self._defaultLagrangianStatisticsValues()['threshold']
            self.setThresholdValue(value)
        return value


    @Variables.noUndo
    def getPostprocessingVolStatusFromName(self, name):
        node = self.node_volume.xmlInitChildNode('property', name=name)
        node2 = node.xmlGetChildNode('postprocessing_recording', 'status')
        if not node2:
            return "on"
        else:
            return "off"


    @Variables.undoLocal
    def setPostprocessingVolStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_volume.xmlInitChildNode('property', name=name)
        node2 = node.xmlInitChildNode('postprocessing_recording', 'status')
        if status == "on":
            node.xmlRemoveChild('postprocessing_recording')
        elif status == "off":
            node2['status'] = status


    # Boundary functions
    # ------------------
    @Variables.undoLocal
    def setBoundaryStatisticsStatus(self, status):
        """
        """
        self.isOnOff(status)
        self.node_boundary['status'] = status


    @Variables.noUndo
    def getBoundaryStatisticsStatus(self):
        """
        """
        self.node_boundary = self.node_stat.xmlInitChildNode('boundary', 'status')
        status = self.node_boundary['status']
        if not status:
            status = self._defaultLagrangianStatisticsValues()['boundary_statistics']
            self.setBoundaryStatisticsStatus(status)
        return status


    @Variables.noUndo
    def getPropertyLabelFromNameBoundary(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name)
        label = node['label']
        if not label:
            label = self._defaultLagrangianStatisticsValues()[name]
            self.setPropertyLabelFromNameBoundary(label, label)
        return label


    @Variables.undoLocal
    def setPropertyLabelFromNameBoundary(self, name, label):
        node = self.node_boundary.xmlInitChildNode('property', name=name)
        node['label'] = label


    @Variables.noUndo
    def getListingPrintingStatusFromName(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name)
        node2 = node.xmlGetChildNode('listing_printing', 'status')
        if not node2:
            return "on"
        else:
            return "off" # node2['status']


    @Variables.undoLocal
    def setListingPrintingStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_boundary.xmlInitChildNode('property', name=name)
        node2 = node.xmlInitChildNode('listing_printing', 'status')
        if status == "on":
            node.xmlRemoveChild('listing_printing')
        elif status == "off":
            node2['status'] = status


    @Variables.noUndo
    def getPostprocessingStatusFromName(self, name):
        node = self.node_boundary.xmlInitChildNode('property', name=name)
        node2 = node.xmlGetChildNode('postprocessing_recording', 'status')
        if not node2:
            return "on"
        else:
            return "off" # node2['status']


    @Variables.undoLocal
    def setPostprocessingStatusFromName(self, name, status):
        self.isOnOff(status)
        node = self.node_boundary.xmlInitChildNode('property', name=name)
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
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        self.case = Case()
        XMLinit(self.case).initialize()


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
