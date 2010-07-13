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
- MatisseThermicModel
- MatisseThermicModelTestCase
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

class MatisseThermicModel:
    """
    Manage the input/output markups in the xml doc about matisse thermiclic load
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case               = case
        self.node_matisse       = self.case.root().xmlInitChildNode('matisse')
        self.node_compute       = self.node_matisse.xmlInitChildNode('compute')
        self.node_phymodel      = self.node_compute.xmlInitChildNode('physical_model')

        self.node_mdcnt         = self.node_phymodel.xmlInitChildNode('imdcnt','status')

        self.status = ('on',
                       'off')

        #
        # Vars Filters
        self.node_map           = self.node_compute.xmlInitChildNode('map')
        self.node_th_cp         = self.node_map.xmlInitChildNode('thermal_capacity')
        self.node_line          = self.node_th_cp.xmlInitChildNode('line')
        self.node_row           = self.node_th_cp.xmlInitChildNode('row')
        self.node_height        = self.node_th_cp.xmlInitChildNode('height')

        self.list_line_area     = self.node_line.xmlGetNodeList('area','label')
        self.list_line_min      = self.node_line.xmlGetNodeList('min')
        self.list_line_max      = self.node_line.xmlGetNodeList('max')
        self.list_line_value    = self.node_line.xmlGetNodeList('value')

        self.list_row_area      = self.node_row.xmlGetNodeList('area','label')
        self.list_row_min       = self.node_row.xmlGetNodeList('min')
        self.list_row_max       = self.node_row.xmlGetNodeList('max')
        self.list_row_value     = self.node_row.xmlGetNodeList('value')

        self.list_height_area   = self.node_height.xmlGetNodeList('area','label')
        self.list_height_min    = self.node_height.xmlGetNodeList('min')
        self.list_height_max    = self.node_height.xmlGetNodeList('max')
        self.list_height_value  = self.node_height.xmlGetNodeList('value')


    def defaultMatisseThermicValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}

        #
        # bool
        default['imdcnt']='off'

        #
        # double
        default['puicon'] = 1000.
        default['tinit'] = 20.
        default['tcrit'] = 80.
        default['emicon'] = 0.4
        default['emimur'] = 0.8
        default['hepcnt'] = 6.
        default['dhpcnt'] = 0.

        #
        # Vars Filters
        default['maplabel'] = 'default'
        default['maplmin']  = 0.
        default['maprmin']  = 0.
        default['maphmin']  = 0.
        default['maplmax']  = 0.
        default['maprmax']  = 0.
        default['maphmax']  = 0.
        default['maplval']  = 1.
        default['maprval']  = 1.
        default['maphval']  = 1.

        return default


    def SetArea(self, areatype, num, label, bmin, bmax, val):
        """
        Add to the XML doc the new values.
        """
        if areatype == 'line' :
            if label: self.list_line_area[num]['label'] = label
            if bmin: self.list_line_area[num].xmlSetData('min',bmin)
            if bmax: self.list_line_area[num].xmlSetData('max',bmax)
            if val: self.list_line_area[num].xmlSetData('value',val)
        elif areatype == 'row' :
            if label: self.list_row_area[num]['label'] = label
            if bmin: self.list_row_area[num].xmlSetData('min',bmin)
            if bmax: self.list_row_area[num].xmlSetData('max',bmax)
            if val: self.list_row_area[num].xmlSetData('value',val)
        elif areatype == 'height' :
            if label: self.list_height_area[num]['label'] = label
            if bmin: self.list_height_area[num].xmlSetData('min',bmin)
            if bmax: self.list_height_area[num].xmlSetData('max',bmax)
            if val: self.list_height_area[num].xmlSetData('value',val)
        else :
            print(areatype, " : Unknown area type")
            sys.exit(1)


    def NewArea(self, areatype, label, bmin, bmax, val):
        """
        Add to the XML doc the new values.
        """
        dlabel, dmin, dmax, dval = self.DefaultArea(areatype)

        if areatype == 'line' :
            if label :
                node = self.node_line.xmlAddChild('area', label=label)
            else :
                node = self.node_line.xmlAddChild('area', label=dlabel)
        elif areatype == 'row' :
            if label :
                node = self.node_row.xmlAddChild('area', label=label)
            else :
                node = self.node_row.xmlAddChild('area', label=dlabel)
        elif areatype == 'height' :
            if label :
                node = self.node_height.xmlAddChild('area', label=label)
            else :
                node = self.node_height.xmlAddChild('area', label=dlabel)
        else :
            print(areatype + " : Unknown area type")
            sys.exit(1)

        nmin = node.xmlAddChild('min')
        nmax = node.xmlAddChild('max')
        nval = node.xmlAddChild('value')

        if bmin != None and bmin != "":
            nmin.xmlSetTextNode(bmin)
        else:
            nmin.xmlSetTextNode(dmin)
        if bmax != None and bmax != "":
            nmax.xmlSetTextNode(bmax)
        else:
            nmax.xmlSetTextNode(dmax)
        if val != None and val != "":
            nval.xmlSetTextNode(val)
        else:
            nval.xmlSetTextNode(dval)

        if areatype == 'line' :
            self.list_line_area.append(node)
        elif areatype == 'row' :
            self.list_row_area.append(node)
        elif areatype == 'height' :
            self.list_height_area.append(node)
        else :
            print(areatype + " : Unknown area type")
            sys.exit(1)


    def DefaultArea(self,areatype):
        """
        Return default values of a area
        """
        modelgeom = MatisseGeom.MatisseGeomModel(self.case)
        dlabel = "default"
        dmin = 0
        dval = 1.0

        if areatype == 'line' :
            dmax = modelgeom.getMatisseGeomDoubleVar('nptran')
        elif areatype == 'row' :
            dmax = modelgeom.getMatisseGeomDoubleVar('nplgrs')
        elif areatype == 'height' :
            dmax = modelgeom.getMatisseGeomDoubleVar('nchest')
        else :
            print("DefaultArea" + areatype + " : Unknown area type")
            sys.exit(1)

        return dlabel, dmin, dmax, dval


    def GetArea(self, areatype, num):
        """
        Get area Values.
        """
        if areatype == 'line' :
            label = self.list_line_area[num]['label']
            bmin = self.list_line_area[num].xmlGetString('min')
            bmax = self.list_line_area[num].xmlGetString('max')
            val = self.list_line_area[num].xmlGetString('value')
        elif areatype == 'row' :
            label = self.list_row_area[num]['label']
            bmin = self.list_row_area[num].xmlGetString('min')
            bmax = self.list_row_area[num].xmlGetString('max')
            val = self.list_row_area[num].xmlGetString('value')
        elif areatype == 'height' :
            label = self.list_height_area[num]['label']
            bmin = self.list_height_area[num].xmlGetString('min')
            bmax = self.list_height_area[num].xmlGetString('max')
            val = self.list_height_area[num].xmlGetString('value')
        else :
            print("GetArea" + areatype + " : Unknown area type")
            sys.exit(1)

        return label, bmin, bmax, val


    def GetAreas(self, areatype):
        """
        Get areas Values.
        """
        llabel=[]
        lbmin=[]
        lbmax=[]
        lval=[]

        if areatype == 'line' :
            nodes_nb = len(self.list_line_area)
        elif areatype == 'row' :
            nodes_nb = len(self.list_row_area)
        elif areatype == 'height' :
            nodes_nb = len(self.list_height_area)
        else :
            print("GetAreas '" + areatype +  "' : Unknown area type")
            sys.exit(1)

        for i in range(0,nodes_nb) :
            label, bmin, bmax, val = self.GetArea(areatype, i)
            llabel.append(label)
            lbmin.append(bmin)
            lbmax.append(bmax)
            lval.append(val)

        return llabel, lbmin, lbmax, lval


    def EraseArea(self, areatype, num):
        """
        Remove Area.
        """
        if areatype == 'line' :
            node = self.list_line_area.pop(num)
        elif areatype == 'row' :
            node = self.list_row_area.pop(num)
        elif areatype == 'height' :
            node = self.list_height_area.pop(num)
        else :
            print("EraseErea" + areatype, " : Unknown area type")
            sys.exit(1)

        node.xmlRemoveChild('min')
        node.xmlRemoveChild('max')
        node.xmlRemoveChild('value')
        node.xmlRemoveNode()


    def setMatisseThermicVar(self, tag, val):
        """
        """
        self.node_phymodel.xmlSetData(tag, val)

        tinit = None
        tcrit = None

        if tag == 'tinit':
            tinit = val
        elif tag == 'tcrit':
            tcrit = val

        import Pages.MatisseModel as Matisse
        Matisse.MatisseThermUpdate(self.case,tinit,tcrit).compute()
        del Matisse


    def getMatisseThermicDoubleVar(self,tag):
        """
        """
        val = self.node_phymodel.xmlGetDouble(tag)

        if val == "" or val == None:
            self.node_phymodel.xmlInitChildNode(tag)
            val = self.defaultMatisseThermicValues()[tag]
            self.setMatisseThermicVar(tag, val)

        return val


    def setNatConvPanacheStatus(self, stat):
        """
        Input natural convection panache status
        """
        if stat not in self.status :
            sys.exit(1)

        self.node_mdcnt['status'] = stat


    def getNatConvPanacheStatus(self):
        """
        Return natural convection panache status
        """
        stat = self.node_mdcnt['status']

        if stat not in self.status :
            stat = self.defaultMatisseThermicValues()['imdcnt']
            self.setNatConvPanacheStatus(stat)
        return stat


#-------------------------------------------------------------------------------
# MatisseThermic Model Test class
#-------------------------------------------------------------------------------
class MatisseThermicModelTestCase(unittest.TestCase):
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

    def checkMatisseThermicModelInstantiation(self):
        """Check whether the MatisseThermicModel class could be instantiated"""
        model = None
        model = MatisseThermicModel(self.case)
        assert model != None, 'Could not instantiate MatisseThermicModel'


def suite():
    testSuite = unittest.makeSuite(MatisseThermicModelTestCase, "check")
    return testSuite

def runTest():
    print("MatisseThermicModelTestCase - TODO**************")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
