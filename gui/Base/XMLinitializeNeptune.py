# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
This module defines the XML data model in which the user defines the physical
options of the treated case.

This module contains the following classe:
- XMLinit
- XMLinitTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.XMLvariablesNeptune import Variables
import code_saturne.Base.Toolbox
from code_saturne.Pages.LocalizationModel import Zone, LocalizationModel
from code_saturne.Pages.OutputControlModel import OutputControlModel

#-------------------------------------------------------------------------------
# Detection of EOS
#-------------------------------------------------------------------------------

EOS = 1
try:
   import eosAva
except:
   EOS = 0
else :
   import eosAva

#-------------------------------------------------------------------------------
# class XMLinit
#-------------------------------------------------------------------------------

class XMLinitNeptune(Variables):
    """
    This class initialize the XML contents of the case.
    """
    def __init__(self, case):
        """
        """
        self.case = case


    def initialize(self, prepro = False):
        """
        Verify that all Heading exist only once in the XMLDocument and create
        the missing heading.
        """
        msg = self.__initHeading(prepro)
        if msg:
            return msg

        if not prepro:
            self.__backwardCompatibility()


            # Initialization (order is important, see turbulenceModelsList method)

            self.XMLNodeAna      = self.case.xmlGetNode('analysis_control')
            self.XMLNodeAverage  = self.XMLNodeAna.xmlInitNode('time_averages')
            self.XMLThermo       = self.case.xmlGetNode('thermophysical_models')
            self.__XMLThermo     = self.XMLThermo.xmlInitNode('thermodynamics')
            self.XMLNodeVariable = self.XMLThermo.xmlInitNode('variables')
            OutputControlModel(self.case).addDefaultWriter()
            OutputControlModel(self.case).addDefaultMesh()
            self.XMLUserScalar   = self.case.xmlGetNode('additional_scalars')
            self.XMLUser         = self.XMLUserScalar.xmlInitNode('users')

            # First Volume Zone definition for all cells -> initialization

            zones = LocalizationModel("VolumicZone", self.case).getZones()
            iok = 0
            for zone in zones:
                if zone.getLabel() == 'all_cells':
                    iok = 1
            if iok == 0:
                zone = Zone("VolumicZone", case = self.case, label = 'all_cells', localization = 'all[]')
                LocalizationModel("VolumicZone", self.case).addZone(zone)
                zone = LocalizationModel("VolumicZone", self.case).getCodeNumberOfZoneLabel('all_cells')

            # If EOS is not avalaible, check if EOS is needed by the current set up.

            if EOS == 0:
                from code_saturne.Pages.MainFieldsModel import MainFieldsModel
                from code_saturne.Pages.ThermodynamicsModel import ThermodynamicsModel

                for fieldId in MainFieldsModel(self.case).getFieldIdList():
                    if ThermodynamicsModel(self.case).getMaterials(fieldId) != "user_material":
                        msg = "The current GUI does not found EOS, but this file of parameters has" \
                              " been generated with EOS. \n\n Please check the disponibility of "  \
                              "the prerequisite EOS."

            # Initialize fields
            if not self.XMLThermo.xmlGetNode('fields'):
                from code_saturne.Pages.MainFieldsModel import MainFieldsModel
                MainFieldsModel(self.case).setPredefinedFlow("None")
                del MainFieldsModel

        return msg


    def __initHeading(self, prepro):
        """
        Create if necessary headings from the root element of the case.
        """
        msg = ""
        tagList = ('solution_domain',
                   'analysis_control',
                   'calculation_management')
        if not prepro :
            tagList += ('thermophysical_models',
                        'additional_scalars',
                        'closure_modeling',
                        'boundary_conditions',
                        'numerical_parameters')

        for tag in tagList:
            nodeList = self.case.root().xmlInitChildNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. "\
                      "There is more than one occurence of the tag: \n\n" + tag + \
                      "\n\nThe application will finish. Sorry."

        for tag in tagList:
            nodeList = self.case.xmlInitNodeList(tag)

            if len(nodeList) > 1:
                msg = "There is an error with the use of the initHeading method. "\
                      "There is more than one occurence of the tag: \n\n" + tag + \
                      "\n\nThe application will finish. Sorry."

        return msg


    def __backwardCompatibility(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        if self.case.root()["solver_version"]:
            vers = self.case.root()["solver_version"]
            history = vers.split(";")
            cur_vers = history[len(history) - 1]
            if history[len(history) - 1] == self.case['package'].version:
                self.__backwardCompatibilityCurrentVersion()
            else:
                self.__backwardCompatibilityOldVersion(cur_vers)
                self.__backwardCompatibilityCurrentVersion()
                his = ""
                for v in history:
                    his = his + v + ";"
                his = his + self.case['package'].version
                self.case.root().xmlSetAttribute(solver_version = his)

        else:
            vers = self.case['package'].version
            self.case.root().xmlSetAttribute(solver_version = vers)

            # apply all backwardCompatibility we don't know when it was create
            self.__backwardCompatibilityOldVersion("-1")
            self.__backwardCompatibilityCurrentVersion()


    def __backwardCompatibilityOldVersion(self, from_vers):
        """
        Change XML in order to ensure backward compatibility for old version
        there is nothing to do for 2.1 to 2.2
        """
        if from_vers == "-1":
            self.__backwardCompatibilityFrom_2_0()
            self.__backwardCompatibilityFrom_2_2()
        elif from_vers == "2.0":
            self.__backwardCompatibilityFrom_2_0()
            self.__backwardCompatibilityFrom_2_2()
        elif from_vers == "2.2":
            self.__backwardCompatibilityFrom_2_2()
        elif from_vers <= "5.0":
            self.__backwardCompatibilityFrom_4_2()
            self.__backwardCompatibilityFrom_5_0()


    def __backwardCompatibilityFrom_2_0(self):
        """
        Change XML in order to ensure backward compatibility from 2.0 to 2.1
        """
        # time averages
        idTA = 0
        for Node in self.case.xmlGetNodeList('time_average'):
            idTA = idTA + 1
            if Node['name'] == None:
                Node['name'] = "TimeAverage_" + str(idTA)

        # Profiles
        for node in self.case.xmlGetNodeList('profile'):
            if node:
                n = node.xmlGetNode("output_type")
                if n == None:
                    freq = node.xmlGetInt("output_frequency")
                    if freq == -1:
                        node.xmlSetData('output_type', "end")
                    else:
                        node.xmlSetData('output_type', "frequency")

        # Users
        nodeInit = self.case.xmlGetNode('users')
        if nodeInit:
            nodeList = nodeInit.xmlGetNodeList('variable')
            for i in range(len(nodeList)):
                varNode = nodeList[i]
                if varNode['support'] == None:
                    varNode['support'] = "cells"

        # <enthalpy_model field_id="2" model="user_function"/>
        for node in self.case.xmlGetNodeList('enthalpy_model'):
            if node['model'] == "user_function":
                node['model'] = "no_source_term"

        nfield = self.case.xmlGetNode('fields')
        if nfield != None:
            nfield1 = nfield.xmlGetNode('field', field_id="1")
            nh = nfield1.xmlGetNode('hresolution')
            if nh['status'] == "off":
                nodeSurf = self.case.xmlGetNode('property', name='surface_tension')
                if nodeSurf == None:
                    XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
                    XMLNodeproperty = XMLNodethermo.xmlInitNode('properties')
                    Variables(self.case).setNewVariableProperty("property", "constant", XMLNodeproperty, "none", "surface_tension", "Surf_tens")

        XMLNodeNonCondens = self.case.xmlGetNode('non_condensable_list')
        if XMLNodeNonCondens != None:
            nodeList = XMLNodeNonCondens.xmlGetNodeList('variable')
            for i in range(len(nodeList)):
                oldName = "non_condensable" + str(i+1)
                newName = "MassFractionNonCondensableGas_" + str(i+1)
                for node in self.case.xmlGetNodeList('variable', name = oldName):
                    if node != None:
                        node['name'] = newName

        # update name for node property
        dicoO2N = {"drho_dh"                             : "d_rho_d_h",
                   "drho_dP"                             : "d_rho_d_P",
                   "Hsat1"                               : "SaturationEnthalpyLiquid",
                   "Hsat2"                               : "SaturationEnthalpyGas",
                   "TsatK"                               : "SaturationTemperature",
                   "DTSDPDerivative"                     : "d_Tsat_d_P",
                   "DHSDPDerivativeLiquid"               : "d_Hsat_d_P_Liquid",
                   "DHSDPDerivativeGas"                  : "d_Hsat_d_P_Gas"}

        for node in self.case.xmlGetNodeList('property'):
            if node['name'] in dicoO2N.keys():
                old_name = node['name']
                node['name'] = dicoO2N[old_name]


    def __backwardCompatibilityFrom_2_2(self):
        """
        Change XML in order to ensure backward compatibility from 2.2 to 2.4-alpha.
        """
        for node in self.case.xmlGetNodeList('hresolution'):
            if node['model'] == None:
                if node['status'] == 'off':
                    node['model'] = 'off'
                else:
                    node['model'] = 'total_enthalpy'


    def __backwardCompatibilityFrom_4_2(self):
        """
        Change XML in order to ensure backward compatibility from versions prior to 4.3
        Reason: Renaming of wall_temperature as boundary_temperature
        """
        # For versions prior to 5.0,renaming of wall_temperature as boundary_temperature
        for node in self.case.xmlGetNodeList('property'):
            if node['name'] == 'wall_temperature':
                node['name']  = 'boundary_temperature'

    def __backwardCompatibilityFrom_5_0(self):

        # For versions prior to 5.0,renaming of wall_temperature as boundary_temperature
        for node in self.case.xmlGetNodeList('property'):
            if node['name'] == 'wall_temperature':
                node['name']  = 'boundary_temperature'

            if node['name'] == 'wall_friction_velocity':
                self.case.xmlRemoveChild('property',
                                         name='wall_friction_velocity',
                                         field_id='none')


        # Add the choice between SGDH and GGDH turbulent thermal flux models
        cnode = self.case.xmlGetNode('closure_modeling')
        tnode = cnode.xmlGetNode('turbulence')
        tvn   = tnode.xmlGetNode('variables')

        thermo_node = self.case.xmlGetNode('thermophysical_models')
        fnode = thermo_node.xmlGetNode('fields')

        for node in tnode.xmlGetNodeList('field'):
            if node['turb_flux'] == None:
                node['turb_flux'] = 'sgdh'

            # Renaming of Rij tensor
            for node in fnode.xmlGetNodeList('field'):
                fieldId = node['field_id']

                rn = tvn.xmlGetNode('variable',
                                    name="ReynoldsStressXX",
                                    field_id=fieldId)
                if rn != None:
                    rn['name']  = "ReynoldsStress"
                    rn['label'] = "ReynoldsStress"+str(fieldId)
                    rn['dim']   = 6

                    for comp in ["XY", "XZ", "YY", "YZ", "ZZ"]:
                        tvn.xmlRemoveChild('variable',
                                           name="ReynoldsStress"+comp,
                                           field_id=fieldId)

        # Renaming k and espilon
        turb_dico = {'TurbDissip':'epsilon',
                     'TurbKineEner_k':'k'}
        for node in tvn.xmlGetNodeList("variable"):
            fieldId = node['field_id']
            for tv in turb_dico.keys():
                if tv in node['name']:
                    node['name']  = turb_dico[tv]
                    node['label'] = turb_dico[tv]+str(fieldId)



        # Modify the rad transfer xml node name for particles to allow a correct
        # workflow with the RTE SOLVER
        tpnode = self.case.xmlGetNode('thermophysical_models')
        if fnode != None:
            for node in fnode.xmlGetNodeList('field'):
                rn = node.xmlGetNode('radiative_transfer')
                if rn != None:
                    st = rn['status']
                    node.xmlRemoveChild('radiative_transfer')
                    node.xmlInitChildNode('particles_radiative_transfer', status=st)

        # Renaming of Pressure
        vnode = thermo_node.xmlGetNode('variables')
        rdico = {'Enthalpy':'enthalpy',
                 'Pressure':'pressure',
                 'Velocity':'velocity',
                 'VolumeFraction':'volume_fraction'}

        if vnode != None:
            for node in vnode.xmlGetNodeList('variable'):
                vname = node['name']
                if vname in rdico.keys():
                    node['name'] = rdico[vname]
                    for nzi in node.xmlGetNodeList('initial_value'):
                        nf = nzi.xmlGetNode('formula')
                        f  = nzi.xmlGetString('formula')
                        nf.xmlSetTextNode(f.replace(vname, rdico[vname]))

            bcnode = self.case.xmlGetNode('boundary_conditions')
            bc_list = ['inlet', 'wall', 'outlet']
            for bc_type in bc_list:
                for nb in bcnode.xmlGetNodeList(bc_type):
                    for nv in nb.xmlGetNodeList('variable'):
                        if nv['name'] in rdico.keys():
                            nv['name'] = rdico[nv['name']]


    def __backwardCompatibilityCurrentVersion(self):
        """
        Change XML in order to ensure backward compatibility.
        """

        # Retrocompatibility: the use of probes in neptune is now the same as for saturn
        for variableType in ('variable', 'property', 'scalar', 'time_average') :
            for parent in self.case.xmlGetNodeList(variableType):
                parent.xmlRemoveChild('probe_recording')
                if parent.xmlGetNode('no_probe') != None :
                    parent.xmlRemoveChild('no_probe')
                    probes_recording = parent.xmlInitNode('probes_recording', 'status')
                    probes_recording['status'] = 'off'

        timelst = []
        timetup = {}
        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                time_node = node.xmlGetNode("time_start")
                if not time_node:
                    node.xmlSetData('time_start', -1.)
                # construc map for profiles
                if node['name']:
                    timetup[node['name']] = node['label']
                    timelst.append(node['name'])
                    # now name = label
                    node['name'] = node['label']

        for node in self.case.xmlGetNodeList('profiles'):
            if node:
                for nodevar in node.xmlGetNodeList('var_prop'):
                    if nodevar['name'] in timelst:
                        nodevar['name'] = timetup[nodevar['name']]


        # suppress gradient and flux reconstruction if needed
        for node in self.case.xmlGetNodeList('variable'):
            n = node.xmlGetNode('flux_reconstruction')
            if n:
                node.xmlRemoveChild('flux_reconstruction')
            n = node.xmlGetNode('gradient_reconstruction')
            if n:
                node.xmlRemoveChild('gradient_reconstruction')

        # update for cdudn
        self.XMLNodethermo   = self.case.xmlGetNode('thermophysical_models')
        self.__XMLNodefields = self.XMLNodethermo.xmlInitNode('fields')
        for node in self.__XMLNodefields.xmlGetNodeList('field'):
            fieldId = node['field_id']
            XMLWallNode = node.xmlGetNode('wall_model')

            if XMLWallNode:
                mdl = XMLWallNode['model']
                node.xmlRemoveChild('wall_model')

                for n in self.case.xmlGetNodeList('wall'):
                    n.xmlInitChildNode('wall_model', field_id = fieldId, model = mdl)

        for node in self.case.xmlGetNodeList('time_average'):
            if node:
                if node['field_id'] == "None":
                    node['field_id'] = "none"
                if not node['name']:
                    node['name'] = node['label']
                node.xmlInitChildNode('listing_printing')
                node.xmlInitChildNode('postprocessing_recording')

        # ------------------------------------------------------------
        # FIXME: TO REMOVE ONCE NCFD 5.0 is out!
        # For versions prior to 5.0,renaming of wall_temperature as boundary_temperature
        for node in self.case.xmlGetNodeList('property'):
            if node['name'] == 'wall_temperature':
                node['name']  = 'boundary_temperature'

            if node['name'] == 'wall_friction_velocity':
                self.case.xmlRemoveChild('property',
                                         name='wall_friction_velocity',
                                         field_id='none')


        # Add the choice between SGDH and GGDH turbulent thermal flux models
        cnode = self.case.xmlGetNode('closure_modeling')
        tnode = cnode.xmlGetNode('turbulence')

        if tnode != None:
            tvn = tnode.xmlGetNode('variables')

            for node in tnode.xmlGetNodeList('field'):
                if node['turb_flux'] == None:
                    node['turb_flux'] = 'sgdh'

            # Renaming of Rij tensor
            for node in self.__XMLNodefields.xmlGetNodeList('field'):
                fieldId = node['field_id']

                rn = tvn.xmlGetNode("variable",
                                    name="ReynoldsStressXX",
                                    field_id=fieldId)
                if rn != None:
                    rn['name']  = "ReynoldsStress"
                    rn['label'] = "ReynoldsStress"+str(fieldId)
                    rn['dim']   = 6

                    for comp in ["XY", "XZ", "YY", "YZ", "ZZ"]:
                        tvn.xmlRemoveChild("variable",
                                           name="ReynoldsStress"+comp,
                                           field_id=fieldId)

            # Renaming k and espilon
            turb_dico = {'TurbDissip':'epsilon',
                         'TurbKineEner_k':'k'}
            for node in tvn.xmlGetNodeList("variable"):
                fieldId = node['field_id']
                for tv in turb_dico.keys():
                    if tv in node['name']:
                        node['name']  = turb_dico[tv]
                        node['label'] = turb_dico[tv]+str(fieldId)

        # Modify the rad transfer xml node name for particles to allow a correct
        # workflow with the RTE SOLVER
        tpnode = self.case.xmlGetNode('thermophysical_models')
        fnode  = tpnode.xmlGetNode('fields')
        if fnode != None:
            for node in fnode.xmlGetNodeList('field'):
                rn = node.xmlGetNode('radiative_transfer')
                if rn != None:
                    st = rn['status']
                    node.xmlRemoveChild('radiative_transfer')
                    node.xmlInitChildNode('particles_radiative_transfer', status=st)

        # Renaming of Pressure
        vnode = tpnode.xmlGetNode('variables')
        rdico = {'Enthalpy':'enthalpy',
                 'Pressure':'pressure',
                 'Velocity':'velocity',
                 'VolumeFraction':'volume_fraction'}

        if vnode != None:
            for node in vnode.xmlGetNodeList('variable'):
                vname = node['name']
                if vname in rdico.keys():
                    node['name'] = rdico[vname]
                    for nzi in node.xmlGetNodeList('initial_value'):
                        nf = nzi.xmlGetNode('formula')
                        f  = nzi.xmlGetString('formula')
                        nf.xmlSetTextNode(f.replace(vname, rdico[vname]))

            bcnode = self.case.xmlGetNode('boundary_conditions')
            bc_list = ['inlet', 'wall', 'outlet']
            for bc_type in bc_list:
                for nb in bcnode.xmlGetNodeList(bc_type):
                    for nv in nb.xmlGetNodeList('variable'):
                        if nv['name'] in rdico.keys():
                            nv['name'] = rdico[nv['name']]


        # ------------------------------------------------------------

#-------------------------------------------------------------------------------
# XMLinit test case
#-------------------------------------------------------------------------------


class XMLinitTestCaseNeptune(unittest.TestCase):
    """
    """
    def setUp(self):
        """
        This method is executed before all "check" methods.
        """
        import XMLengine
        Toolbox.GuiParam.lang = 'en'
        self.doc = XMLengine.XMLDocument("")
        self.case = XMLengine.Case(None)


    def tearDown(self):
        """
        This method is executed after all "check" methods.
        """
        del self.case
        del self.doc


    def xmlNodeFromString(self, string):
        """Private method to return a xml node from string"""
        return self.doc.parseString(string).root()


    def checkXMLinitInstantiation(self):
        """
        Check whether the Case class could be instantiated
        """
        xmldoc = None
        xmldoc = XMLinit(self.case)
        assert xmldoc != None, 'Could not instantiate XMLinit'


    def checkInitHeading(self):
        """
        Check whether the headings markups could be initialized
        """
        #doc = \
        #'<NeptuneCFD case="" study="" version="1.0">'\
        #'<solution_domain/>'\
        #'</NEPTUNE_CFD_GUI>'

        XMLinit(self.case)

        assert self.case.root() == self.xmlNodeFromString(doc), \
               'Could not use the constructor of the XMLinit class'


def suite():
    testSuite = unittest.makeSuite(XMLinitTestCase, "check")
    return testSuite


def runTest():
    print("XMLinitTestCase to be completed...")
    runner = unittest.TextTestRunner()
    runner.run(suite())


#-------------------------------------------------------------------------------
# End of XMLinit
#-------------------------------------------------------------------------------
