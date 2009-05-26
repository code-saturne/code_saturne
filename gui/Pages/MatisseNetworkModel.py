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
This module defines the storage system type.

This module contains the following classes and function:
- MatisseNetworkModel
- MatisseNetworwTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
import Pages.MatisseTypeModel as MatisseType
import Pages.MatisseGeomModel as MatisseGeom

#-------------------------------------------------------------------------------
# Matisse thermic class
#-------------------------------------------------------------------------------

class MatisseNetworkModel:
    def __init__(self, case):
        """
        Constructor.
        """
        self.case               = case
        self.node_matisse       = self.case.root().xmlInitChildNode('matisse')
        self.node_compute       = self.node_matisse.xmlInitChildNode('compute')
        self.node_phymodel      = self.node_compute.xmlInitChildNode('physical_model')

        #
        # Vars Filters
        self.node_map           = self.node_compute.xmlInitChildNode('map')
        self.node_hl            = self.node_map.xmlInitChildNode('headloss')
        self.node_network       = self.node_hl.xmlInitChildNode('network')
        self.node_line          = self.node_network.xmlInitChildNode('line')
        self.node_row           = self.node_network.xmlInitChildNode('row')

        self.list_line_area     = self.node_line.xmlGetNodeList('area','label') 
        self.list_row_area      = self.node_row.xmlGetNodeList('area','label') 


    def defaultMatisseNetworkValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}

        #
        # double
        modelgeom = MatisseGeom.MatisseGeomModel(self.case)
        hreso = modelgeom.getMatisseGeomDoubleVar('nchest')
        step = modelgeom.getMatisseGeomDoubleVar('epchel')
        
        default['nbcellreso'] = hreso
        default['nbcellplen'] = 0.
        default['hreso'] = hreso*step
        default['hplen'] = 0.

        return default


    def SetArea(self, areatype, num, label, bmin, bmax):
        """
        Add to the XML doc the new values.
        """
        if areatype == 'line':
            if label: self.list_line_area[num]['label'] = label
            if bmin: self.list_line_area[num].xmlSetData('min',bmin)
            if bmax: self.list_line_area[num].xmlSetData('max',bmax)
        elif areatype == 'row':
            if label: self.list_row_area[num]['label'] = label
            if bmin: self.list_row_area[num].xmlSetData('min',bmin)
            if bmax: self.list_row_area[num].xmlSetData('max',bmax)
        else :
            print areatype, ": Unknown area type"
            sys.exit(1)


    def NewArea(self, areatype, label, bmin, bmax):
        """
        Add to the XML doc the new values.
        """

        if areatype == 'line':
            node = self.node_line.xmlAddChild('area', label=label)
        elif areatype == 'row':
            node = self.node_row.xmlAddChild('area', label=label)
        else :
            print areatype + ": Unknown area type"
            sys.exit(1)
            
        node.xmlAddChild('min').xmlSetTextNode(bmin)
        node.xmlAddChild('max').xmlSetTextNode(bmax)

        if areatype == 'line':
            self.list_line_area.append(node)
        elif areatype == 'row':
            self.list_row_area.append(node)
        else :
            print areatype + ": Unknown area type"
            sys.exit(1)


    def GetArea(self, areatype, num):
        """
        Get area Values.
        """
        if areatype == 'line':
            label = self.list_line_area[num]['label']
            bmin = self.list_line_area[num].xmlGetString('min')
            bmax = self.list_line_area[num].xmlGetString('max')
        elif areatype == 'row':
            label = self.list_row_area[num]['label']
            bmin = self.list_row_area[num].xmlGetString('min')
            bmax = self.list_row_area[num].xmlGetString('max')
        elif areatype == 'height':
            label = self.list_height_area[num]['label']
            bmin = self.list_height_area[num].xmlGetString('min')
            bmax = self.list_height_area[num].xmlGetString('max')
        else :
            print "GetArea" + areatype + " : Unknown area type"
            sys.exit(1)

        return label, bmin, bmax


    def GetAreas(self, areatype):
        """
        Get areas Values.
        """
        llabel=[]
        lbmin=[]
        lbmax=[]

        if areatype == 'line':
            nodes_nb = len(self.list_line_area)
        elif areatype == 'row':
            nodes_nb = len(self.list_row_area)
        else :
            print "GetAreas '" + areatype +  "': Unknown area type"
            sys.exit(1)

        for i in range(0,nodes_nb):
            label, bmin, bmax = self.GetArea(areatype, i)
            llabel.append(label)
            lbmin.append(bmin)
            lbmax.append(bmax)

        return llabel, lbmin, lbmax

        
    def EraseArea(self, areatype, num):
        """
        Remove Area.
        """
        if areatype == 'line':
            node = self.list_line_area.pop(num)
        elif areatype == 'row':
            node = self.list_row_area.pop(num)
        else :
            print "EraseErea" + areatype, " : Unknown area type"
            sys.exit(1)

        node.xmlRemoveChild('min')
        node.xmlRemoveChild('max')
        node.xmlRemoveNode()


    def setMatisseNetworkVar(self, tag, val):
        """
        """
        if tag.find('reso') >= 0:
            tag = tag[:-1]
        self.node_phymodel.xmlSetData(tag, val)


    def getMatisseNetworkDoubleVar(self, tag):
        """
        """
        if tag.find('reso') >= 0:
            tagtmp = tag[:-1]
        else:
            tagtmp = tag
        val = self.node_phymodel.xmlGetDouble(tagtmp)

        if val == "" or val == None: 
            self.node_phymodel.xmlInitChildNode(tagtmp)
            val = self.defaultMatisseNetworkValues()[tagtmp]
            self.setMatisseNetworkVar(tag, val)

        return val

#-------------------------------------------------------------------------------
# MatisseNetwork Model Test Class
#-------------------------------------------------------------------------------

class MatisseNetworkModelTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        """This method is executed before all "check" methods."""
        from Base.XMLengine import Case, XMLDocument
        from Base.XMLinitialize import XMLinit
        Tool.GuiParam.lang = 'en'
        self.case = Case(None)
        XMLinit(self.case)
        self.doc = XMLDocument()

    def tearDown(self):
        """This method is executed after all "check" methods."""
        del self.case
        del self.doc

    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()

    def checkMatisseNetworkModelInstantiation(self):
        """Check whether the NOMModel class could be instantiated"""
        model = None
        model = MatisseNetworkModel(self.case)
        assert model != None, 'Could not instantiate MatisseNetworkModel'
        
def suite():
    testSuite = unittest.makeSuite(MatisseNetworkModelTestCase, "check")
    return testSuite

def runTest():
    print "MatisseNetworkModelTestCase - A FAIRE**************"
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------