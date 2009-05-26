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
"""


#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------


import sys, os, shutil, unittest
import filecmp, os.path


#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from Base.Common import *
from Base.Toolbox import *
import Pages.BodyForcesModel as BodyForcesModel


#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------


class MatisseInit :
    """
    This class is called only by XMLmodel for initializing:
     - Turbulence
     - Gravity
     - Initilization zones
     - Thermal Scalar
     - Additional scalars
     - Properties density , thermal_conductivity, molecular_viscosity
     - MatisseThermicModel
     - MatisseThermUpdate
     - MatisseGeomModel
    """
    def __init__(self, case):

        self.case = case

        #
        # save old geom files
        geomvault = self.case['mesh_path'] +'/vault.geom' 
        geomemm = self.case['mesh_path'] +'/emm.geom'
        try:
            shutil.copyfile(geomvault,geomvault+'~')
        except:
            try:
                shutil.copyfile(geomemm,geomemm+'~')
            except:
                pass
        
        #
        # turbulence
        import Pages.TurbulenceModel as Turbulence
        Turbulence.TurbulenceModel(self.case).setTurbulenceModel("off")
        del Turbulence
        
        #
        # gravity
        model_bodyForce = BodyForcesModel.BodyForcesModel(self.case)
        model_bodyForce.setGravity('gravity_z', -9.81)
        model_bodyForce.setHydrostaticPressure('on')

        #
        # Zones initialization
        from Pages.LocalizationModel import LocalizationModel, Zone
        model_zone = LocalizationModel('VolumicZone', self.case)
        zones =[["zone1","1","0"],["zone2","2","2"],["zone3","3","3"],["zone4","4","4"],
                ["zone5","5","6"],["zone6","6","7"],["zone7","7","8"],["zone8","8","9"],
                ["zone9","9","10"]]
        for i in range(len(zones)):
            if zones[i][0] not in model_zone.getLabelsZonesList():
                zone = Zone('VolumicZone',
                            label = zones[i][0],
                            codeNumber = zones[i][1],
                            localization = zones[i][2],
                            nature = 'initialization')
                model_zone.addZone(zone)

        #
        # Thermal Scalar
        import Pages.ThermalScalarModel as ThermalScalar
        thermal = ThermalScalar.ThermalScalarModel(self.case)
        thermal.setThermalModel('temperature_celsius')
        thermal.node_therm.xmlRemoveChild('initial_value')

        #
        # Additional scalars
        import Pages.DefineUserScalarsModel as DefineUserScalars
        temp = DefineUserScalars.DefineUserScalarsModel(self.case)

        t = PageText()
        for i in range(len(zones)):
            temp.setScalarValues('T_PeauCol', zones[i][1], 0.0, -1e+12, 1e+12, t.NO_VARIANCE)
        temp.setScalarDiffusivityChoice('T_PeauCol', 'variable')

        for i in range(len(zones)):
            temp.setScalarValues('T_PeauMur', zones[i][1], 0.0, -1e+12, 1e+12, t.NO_VARIANCE)
        temp.setScalarDiffusivityChoice('T_PeauMur', 'variable')

        temp.node_user_sca.xmlRemoveChild('initial_value')
        temp.node_user_sca.xmlRemoveChild('variance')
        
        del DefineUserScalars
        
        #
        # properties density , thermal_conductivity, molecular_viscosity
        import Pages.FluidCharacteristicsModel as FluidCharacteristics
        fluid = FluidCharacteristics.FluidCharacteristicsModel(self.case)
        fluid.setPropertyMode('density','variable')
        fluid.setPropertyMode('molecular_viscosity','variable')
        fluid.setPropertyMode('thermal_conductivity','variable')
        fluid.setPropertyMode('specific_heat','constant')
        del FluidCharacteristics

        import Pages.MatisseThermicModel as MatisseThermic
        model = MatisseThermic.MatisseThermicModel(self.case)
        tinit = model.getMatisseThermicDoubleVar('tinit')
        tcrit = model.getMatisseThermicDoubleVar('tcrit')
        MatisseThermUpdate(self.case, tinit, tcrit).compute()
        del MatisseThermic

        from Pages.MatisseTypeModel import MatisseTypeModel
        m = MatisseTypeModel(self.case)
        m.getMatisseType()

        from Pages.MatisseGeomModel import MatisseGeomModel
        m = MatisseGeomModel(self.case)
        m.updateMeshAndProbes()


class MatisseThermUpdate :
    def __init__(self , case, tinit, tcrit):

        self.case = case
        self.tinit = tinit
        self.tcrit = tcrit

    def compute(self):
        """
        Update XML
        """
        #
        # Additional Scalars (init)
        import Pages.DefineUserScalarsModel as DefineUserScalars
        mdl = DefineUserScalars.DefineUserScalarsModel(self.case)

        for node in mdl.node_user_sca.xmlGetNodeList('scalar'):
            label = node['label']
            if node['type'] == "thermal":
                if self.tinit : mdl.setScalarInitialValue("1",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("2",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("3",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("4",label,self.tinit)
                if self.tcrit : mdl.setScalarInitialValue("5",label,self.tcrit)
                if self.tinit : mdl.setScalarInitialValue("6",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("7",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("8",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("9",label,self.tinit)
            else:
                if self.tinit : mdl.setScalarInitialValue("1",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("2",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("3",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("4",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("5",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("6",label,self.tinit)
                if self.tcrit : mdl.setScalarInitialValue("7",label,self.tcrit)
                if self.tinit : mdl.setScalarInitialValue("8",label,self.tinit)
                if self.tinit : mdl.setScalarInitialValue("9",label,self.tinit)
        del DefineUserScalars


class MatisseMeshRunning :
    """
    """
    def __init__(self, case):
        """
        """
        self.case = case
        self.ok = True
        #
        # update node <solution_domain> in XML file 
        node_preprocessor  = self.case.root().xmlGetNode('solution_domain')
        node_meshes_list= node_preprocessor.xmlInitChildNode('meshes_list')
        mesh_nodes      = node_meshes_list.xmlGetNodeList('mesh', 'name')

        if len(mesh_nodes) == 1:
            try:
                desFile = mesh_nodes[0]['name'] 
                geomFile = desFile.split(".")[0]+'.geom'
                datFile = desFile.split(".")[0]+'.dat'
            except:
                print "MatisseMeshRunning: see the node mesh"
                sys.exit(1)

            #
            # Files
            newGeom = self.case['mesh_path'] + '/' + geomFile
            oldGeom = self.case['mesh_path'] + '/' + geomFile + '~'
            
            oldDes = self.case['mesh_path']  + '/' + desFile
            datFile = self.case['mesh_path'] + '/' + datFile

            #
            # Checks
            simail = False
            if not os.path.isfile(newGeom) or \
               not os.path.isfile(datFile) or \
               not os.path.isfile(oldDes)  :

                import Pages.MatisseGeomModel as MatisseGeom 
                MatisseGeom.MatisseGeomModel(self.case).updateMeshAndProbes()
                del MatisseGeom
                simail = True

            else :
                if not os.path.isfile(oldGeom) :
                    simail = True
                else:    
                    if not filecmp.cmp(oldGeom, newGeom):
                        simail = True

            #
            # Simail
            if simail:
                if os.path.isfile(oldDes) :
                    os.remove(oldDes)
                newDes = oldDes
                shutil.copyfile(newGeom, oldGeom)
                os.chdir(self.case['mesh_path'])
                cmd = "xsimail -b " + datFile + " | tee " + self.case['case_path'] + "/RESU/listsim"
                os.system(cmd)
                os.chdir(self.case['case_path'])
                if not os.path.isfile(newDes)  :
                    self.ok = False
            
        else :
            print "MatisseMeshRunning: see meshes_list"
            sys.exit(1)

#-------------------------------------------------------------------------------
# Page text class
#-------------------------------------------------------------------------------

class PageText:
    """
    Storage of all texts and messages for this page.
    """
    def __init__(self):
        if GuiParam.lang == 'fr':
            self.NO_VARIANCE = "no"
        else:
            self.NO_VARIANCE = "no"

            
#-------------------------------------------------------------------------------
# Matisse Init test class
#-------------------------------------------------------------------------------

class MatisseInitTestCase(unittest.TestCase):
    """
    """
##    def setUp(self):
##        """This method is executed before all "check" methods."""
##        from Base.XMLengine import Case, XMLDocument
##        from Base.XMLinitialize import XMLinit
##        GuiParam.lang = 'en'
##        self.case = Case(None)
##        XMLinit(self.case)
##        self.doc = XMLDocument()
##
##    def tearDown(self):
##        """This method is executed after all "check" methods."""
##        del self.case
##        del self.doc
##
##    def xmlNodeFromString(self, string):
##        """Private method to return a xml node from string"""
##        return self.doc.parseString(string).root()
##
##    def checkMatisseInstantiation(self):
##        """Check whether the NOMModel class could be instantiated"""
##        model = None
##        model = MatisseInit(self.case)
##        assert model != None, 'Could not instantiate MatisseInit'


#-------------------------------------------------------------------------------
# Matisse ThermUpdate test class
#-------------------------------------------------------------------------------

class MatisseThermUpdateTestCase(unittest.TestCase):
    """
    """
##    def setUp(self):
##        """This method is executed before all "check" methods."""
##        from Base.XMLengine import Case, XMLDocument
##        from Base.XMLinitialize import XMLinit
##        GuiParam.lang = 'en'
##        self.case = Case(None)
##        XMLinit(self.case)
##        self.doc = XMLDocument()
##
##    def tearDown(self):
##        """This method is executed after all "check" methods."""
##        del self.case
##        del self.doc
##
##    def xmlNodeFromString(self, string):
##        """Private method to return a xml node from string"""
##        return self.doc.parseString(string).root()
##
##    def checkMatisseThermUpdateInstantiation(self):
##        """Check whether the NOMModel class could be instantiated"""
##        model = None
##        model = MatisseThermUpdate(self.case, '0','0')
##        assert model != None, 'Could not instantiate MatisseThermUpdate'


def suite1():
    testSuite = unittest.makeSuite(MatisseInitTestCase, "check")
    return testSuite

def suite2():
    testSuite = unittest.makeSuite(MatisseThermUpdateTestCase, "check")
    return testSuite

def runTest():
    runner = unittest.TextTestRunner()

    print "MatisseInitTestCase: - A FAIRE************"
    runner.run(suite1())

    print "MatisseThermUpdateTestCase: - A FAIRE************"
    runner.run(suite2())        
        

        