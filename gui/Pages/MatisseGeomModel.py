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
- MatisseGeomModel
- MatisseGeomTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import shutil, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *
import Base.Toolbox as Tool
import Pages.MatisseTypeModel as MatisseType

#-------------------------------------------------------------------------------
# Matisse Geom model class
#-------------------------------------------------------------------------------

class MatisseGeomModel:
    """
    Manage the input/output markups in the xml doc about Turbulence
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case                = case
        self.node_matisse        = self.case.root().xmlInitChildNode('matisse')
        self.node_mesh           = self.node_matisse.xmlInitChildNode('mesh')
        self.node_compute        = self.node_matisse.xmlInitChildNode('compute')
        self.node_geom_mesh      = self.node_mesh.xmlInitChildNode('geometry')
        self.node_geom_compute   = self.node_compute.xmlInitChildNode('geometry')

        self.mesh_geom_vars = ['nechrg','nergrs','neclrg','neciel','nergch',
                               'jeuchr','jeurcl','jeuclr','jeurch',
                               'hbdtoi']

        self.compute_geom_vars = ['epchem','nptran','nplgrs','nelgrs',
                                  'nchest','netran','epregi','hconve',
                                  'rconve','hchali','hcheva','hfttoi',
                                  'ptrres','frdtra','plgres','epchel',
                                  'dmcont']
        typeModel = MatisseType.MatisseTypeModel(self.case)
        self.matisseType = typeModel.getMatisseType()


    def defaultMatisseGeomValues(self):
        """
        Return in a dictionnary which contains default values
        """
        default = {}
        #
        # node geom_mesh
        default['nechrg'] = 1
        default['nergrs'] = 5
        default['neclrg'] = 5
        default['nergch'] = 1
        default['neciel'] = 20
        default['jeuchr'] = 0.2
        default['jeurcl'] = 2.
        default['jeuclr'] = 2.
        default['jeurch'] = 0.2
        default['hbdtoi'] = 15.

        #
        # node geom_compute
        default['epchem'] = 2.
        default['nptran'] = 1
        default['nplgrs'] = 24
        default['nelgrs'] = 1
        default['nchest'] = 24
        default['netran'] = 1
        default['epregi'] = 1.
        default['hconve'] = 15.
        default['rconve'] = 4.93
        default['hchali'] = 24.
        default['hcheva'] = 24.
        default['hfttoi'] = 20.
        default['ptrres'] = 0.8
        default['frdtra'] = 1.
        default['plgres'] = 0.7
        default['epchel'] = 0.333
        default['dmcont'] = 0.347

        return default


    def setMatisseGeomVar(self, tag, val):
        """
        """
        if tag in self.mesh_geom_vars :
            self.node_geom_mesh.xmlSetData(tag, val)
        elif tag in self.compute_geom_vars :
            self.node_geom_compute.xmlSetData(tag, val)
        else :
            print(tag + ": unknown parameter")
            sys.exit(1)
        self.updateMeshAndProbes()


    def getMatisseGeomIntVar(self,tag):
        """
        """
        if tag in self.mesh_geom_vars :
            val = self.node_geom_mesh.xmlGetInt(tag)
        elif tag in self.compute_geom_vars :
            val = self.node_geom_compute.xmlGetInt(tag)
        else :
            print(tag + ": unknown parameter")
            sys.exit(1)

        if val == "" or val == None:
            if tag in self.mesh_geom_vars :
                self.node_geom_mesh.xmlInitChildNode(tag)
            elif tag in self.compute_geom_vars :
                self.node_geom_compute.xmlInitChildNode(tag)

            val = self.defaultMatisseGeomValues()[tag]
            self.setMatisseGeomVar(tag, val)

        return val


    def getMatisseGeomDoubleVar(self,tag):
        """
        """
        if tag in self.mesh_geom_vars :
            val = self.node_geom_mesh.xmlGetDouble(tag)
        elif tag in self.compute_geom_vars :
            val = self.node_geom_compute.xmlGetDouble(tag)
        else :
            print(tag + ": unknown parameter")
            sys.exit(1)

        if val == "" or val == None:
            if tag in self.mesh_geom_vars :
                self.node_geom_mesh.xmlInitChildNode(tag)
            elif tag in self.compute_geom_vars :
                self.node_geom_compute.xmlInitChildNode(tag)

            val = self.defaultMatisseGeomValues()[tag]
            self.setMatisseGeomVar(tag, val)

        return val


    def updateMeshAndProbes(self):
        """
        """
        self.replaceText={}
        self.replaceText['PTTUB']   = str(self.getMatisseGeomDoubleVar('ptrres'))
        self.replaceText['NLIGN']   = str(self.getMatisseGeomIntVar('nptran'))
        self.replaceText['PLTUB']   = str(self.getMatisseGeomDoubleVar('plgres'))
        self.replaceText['NRANG']   = str(self.getMatisseGeomIntVar('nplgrs'))
        self.replaceText['PVTUB']   = str(self.getMatisseGeomDoubleVar('epchel'))
        self.replaceText['NZTUB']   = str(self.getMatisseGeomIntVar('nchest'))
        self.replaceText['ABAMON']  = str(self.getMatisseGeomDoubleVar('epregi'))
        self.replaceText['EPSCHEM'] = str(self.getMatisseGeomDoubleVar('epchem'))
        self.replaceText['HCHALIM'] = str(self.getMatisseGeomDoubleVar('hchali'))
        self.replaceText['HCHEVAC'] = str(self.getMatisseGeomDoubleVar('hcheva'))
        self.replaceText['NCELT']   = str(self.getMatisseGeomIntVar('netran'))
        self.replaceText['NCELL']   = str(self.getMatisseGeomIntVar('nelgrs'))
        self.replaceText['GAPGECH'] = str(self.getMatisseGeomDoubleVar('jeurch'))
        self.replaceText['NELGECH'] = str(self.getMatisseGeomIntVar('nergch'))
        self.replaceText['GAPREG']  = str(self.getMatisseGeomDoubleVar('jeuclr'))
        self.replaceText['NELREG']  = str(self.getMatisseGeomIntVar('neclrg'))
        self.replaceText['GAPGRE']  = str(self.getMatisseGeomDoubleVar('jeurcl'))
        self.replaceText['NELGRE']  = str(self.getMatisseGeomIntVar('nergrs'))
        self.replaceText['GAPCHEG'] = str(self.getMatisseGeomDoubleVar('jeuchr'))
        self.replaceText['NELCHEG'] = str(self.getMatisseGeomIntVar('nechrg'))
        self.replaceText['HCONVER'] = str(self.getMatisseGeomDoubleVar('hconve'))
        self.replaceText['RCONVER'] = str(self.getMatisseGeomDoubleVar('rconve'))
        self.replaceText['NCELCIEL']= str(self.getMatisseGeomIntVar('neciel'))
        self.replaceText['HTOITBAS']= str(self.getMatisseGeomDoubleVar('hbdtoi'))
        self.replaceText['HTOITHAU']= str(self.getMatisseGeomDoubleVar('hfttoi'))

        self.meshUpdate()
        self.probesUpdate()


    def meshUpdate(self):
        """
        Update mesh file
        """
        #
        # copy of default mesh
        print("Not available in the current version.")
        sys.exit(1)

        # TODO: To be adapated, once GUI and Kernel are mergerd
        defaultMeshDir = cs_home[:-1] + '/data/mati'

        if (self.matisseType == 'vault') or (self.matisseType == 'djw'):
            meshFile = 'vault.dat'
            geomFile = 'vault.geom'
            desFile  = 'vault.des'
        elif self.matisseType == 'emm':
            meshFile = 'emm.dat'
            geomFile = 'emm.geom'
            desFile  = 'emm.des'
        else:
            print(self.matisseType + ": matisseType unknown")
            sys.exit(1)

        geomFileInitName = defaultMeshDir + '/' + geomFile
        geomFileName = self.case['mesh_path'] + '/' + geomFile

        meshFileName = self.case['mesh_path'] + '/' + meshFile
        shutil.copyfile(defaultMeshDir + '/' + meshFile, meshFileName)

        geomFileInit = open(geomFileInitName, 'r')
        geomFile     = open(geomFileName, 'w')

        #
        # update mesh
        while 1:
            line = geomFileInit.readline()
            if line != '':
                for param in self.replaceText:
                    newLine = line.replace(param,self.replaceText[param])
                    line = newLine
                geomFile.write(newLine)
            else:
                break

        #
        # update node <solution_domain> in XML file
        node_preprocessor  = self.case.root().xmlGetNode('solution_domain')
        node_meshes_list= node_preprocessor.xmlInitChildNode('meshes_list')
        mesh_nodes      = node_meshes_list.xmlGetNodeList('mesh', 'name')

        if mesh_nodes != []:
            for i in range(len(mesh_nodes)):
                if i == 0:
                    mesh_nodes[i]['format'] = 'des'
                    mesh_nodes[i]['name']   = desFile
                else:
                    mesh_nodes[i].xmlRemoveNode()
        else :
            node_meshes_list.xmlInitChildNode('mesh', format='des',name=desFile)


    def probesUpdate(self):
        """
        Update probes
        """
        #
        # update probes
        cy = float(self.replaceText['PLTUB'])
        nr = int(self.replaceText['NRANG'])
        cz = float(self.replaceText['PVTUB'])
        nz = int(self.replaceText['NZTUB'])

        probes=[]

        probes.append(["1",[0.02, 0.25*cy*nr/11., 1*0.95*cz*nz/7.]])
        probes.append(["2",[0.02, 0.25*cy*nr/11., 2*0.95*cz*nz/7.]])
        probes.append(["3",[0.02, 0.25*cy*nr/11., 3*0.95*cz*nz/7.]])
        probes.append(["4",[0.02, 0.25*cy*nr/11., 4*0.95*cz*nz/7.]])
        probes.append(["5",[0.02, 0.25*cy*nr/11., 5*0.95*cz*nz/7.]])
        probes.append(["6",[0.02, 0.25*cy*nr/11., 6*0.95*cz*nz/7.]])
        probes.append(["7",[0.02, 0.25*cy*nr/11., 7*0.95*cz*nz/7.]])

        probes.append(["8",[0.02, 4.*cy*nr/11., 1*0.95*cz*nz/7.]])
        probes.append(["9",[0.02, 4.*cy*nr/11., 2*0.95*cz*nz/7.]])
        probes.append(["10",[0.02, 4.*cy*nr/11., 3*0.95*cz*nz/7.]])
        probes.append(["11",[0.02, 4.*cy*nr/11., 4*0.95*cz*nz/7.]])
        probes.append(["12",[0.02, 4.*cy*nr/11., 5*0.95*cz*nz/7.]])
        probes.append(["13",[0.02, 4.*cy*nr/11., 6*0.95*cz*nz/7.]])
        probes.append(["14",[0.02, 4.*cy*nr/11., 7*0.95*cz*nz/7.]])

        probes.append(["15",[0.02, 7.*cy*nr/11., 1*0.95*cz*nz/7.]])
        probes.append(["16",[0.02, 7.*cy*nr/11., 2*0.95*cz*nz/7.]])
        probes.append(["17",[0.02, 7.*cy*nr/11., 3*0.95*cz*nz/7.]])
        probes.append(["18",[0.02, 7.*cy*nr/11., 4*0.95*cz*nz/7.]])
        probes.append(["19",[0.02, 7.*cy*nr/11., 5*0.95*cz*nz/7.]])
        probes.append(["20",[0.02, 7.*cy*nr/11., 6*0.95*cz*nz/7.]])
        probes.append(["21",[0.02, 7.*cy*nr/11., 7*0.95*cz*nz/7.]])

        probes.append(["22",[0.02, 10.*cy*nr/11., 1*0.95*cz*nz/7.]])
        probes.append(["23",[0.02, 10.*cy*nr/11., 2*0.95*cz*nz/7.]])
        probes.append(["24",[0.02, 10.*cy*nr/11., 3*0.95*cz*nz/7.]])
        probes.append(["25",[0.02, 10.*cy*nr/11., 4*0.95*cz*nz/7.]])
        probes.append(["26",[0.02, 10.*cy*nr/11., 5*0.95*cz*nz/7.]])
        probes.append(["27",[0.02, 10.*cy*nr/11., 6*0.95*cz*nz/7.]])
        probes.append(["28",[0.02, 10.*cy*nr/11., 7*0.95*cz*nz/7.]])

        #
        # update XML
        node_control = self.case.root().xmlGetChildNode('analysis_control')
        node_output  = node_control.xmlGetChildNode('output')
        node_probes  = node_output.xmlGetNodeList('probe','name')

        import Pages.OutputControlModel as OutputControl
        output = OutputControl.OutputControlModel(self.case)
        for probe in probes :
            try:
                output.replaceMonitoringPointCoordinates(x = probe[1][0],
                                                         y = probe[1][1],
                                                         z = probe[1][2],
                                                         name = probe[0])
            except:
                output.addMonitoringPoint(x = probe[1][0],
                                          y = probe[1][1],
                                          z = probe[1][2])


#-------------------------------------------------------------------------------
# MatisseGeom Model test case
#-------------------------------------------------------------------------------


class MatisseGeomTestCase(unittest.TestCase):
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


    def checkMatisseGeomInstantiation(self):
        """
        Check whether the TurbulenceModel class could be instantiated
        """
        model = None
        model = MatisseGeomModel(self.case)
        assert model != None, 'Could not instantiate MatisseGeomModel'


def suite():
    testSuite = unittest.makeSuite(MatisseGeomTestCase, "check")
    return testSuite


def runTest():
    print("MatisseGeomTestCase - TODO**************")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
