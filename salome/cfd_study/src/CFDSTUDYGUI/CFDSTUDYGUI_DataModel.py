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
Data Model
==========
Definitions of the function that allow to represent the CFC studies in a
tree representation.

SALOME data structure
---------------------
SALOMEDS (SALOME data structure) is a library that provides support for a multi-component
document of SALOME platform. Components can use SALOMEDS to publish their data inside a SALOMEDS
document (Study object). Publishing the data in a common document gives the following advantages
for a custom component:

 - the data becomes available for other components (for processing, visualization, etc.),
   it can accessed using SALOMEDS tools and services;

 - the data becomes automatically persistent (can be saved and restored), as persistence is
   already implemented in SALOMEDS library.

SALOMEDS also provides the mechanism of data persistence for components that do not publish
their data in a common SALOMEDS data structure. This mechanism is described in Implementing
persistence section of the tutorial. Briefly, SALOMEDS provides the following: a component
saves its data in arbitiary format to an external file and returns the name of this file
to SALOMEDS. SALOMEDS serializes this file into a binary stream and includes it into the common
Study file on save operation. When the data must be restored, exactly the same file is created
by SALOMEDS for the component, and the component itself is responsible for loading it.

SALOME Study
------------

A SALOME platform document that contains data of multiple components. The data is organized in
a tree-like structure within the Study. SALOMEDS library supports persistence of Study.
Every branch of the tree is represented by an SObject.

WARNING: a SALOME Study should not be confused with a CFD study.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------
from __future__ import print_function
import os
import re
import string
import logging
import subprocess

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtWidgets import *
from code_saturne import cs_runcase

from omniORB import CORBA

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

from SALOME_NamingServicePy import *
from LifeCycleCORBA import *
import SALOMEDS
import SALOMEDS_Attributes_idl

import SMESH
import salome

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

from CFDSTUDYGUI_Commons import CFD_Code, BinCode, Trace, sgPyQt, sg
from CFDSTUDYGUI_Commons import CaseInProcessStart, CaseInProcessEnd
from CFDSTUDYGUI_Commons import CFD_Saturne, CFD_Neptune
import CFDSTUDYGUI_SolverGUI
import CFDSTUDYGUI_Commons
from CFDSTUDYGUI_CommandMgr import runCommand
from CFDSTUDYGUI_Message import cfdstudyMess
from cs_exec_environment import separate_args

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_DataModel")
log.setLevel(logging.NOTSET)


#-------------------------------------------------------------------------------
# Module name. Attribut "AttributeName" for the related SObject.
#-------------------------------------------------------------------------------

__MODULE_NAME__ = "CFDSTUDY"
__MODULE_ID__   = 10000
__OBJECT_ID__   = 10010

#-------------------------------------------------------------------------------
# Definition of the type of objects for representation in the Object Browser.
# Attribut "AttributeLocalID" for the related SObject.
#-------------------------------------------------------------------------------

dict_object = {}

dict_object["OtherFile"]     = 100000
dict_object["OtherFolder"]   = 100001
dict_object["Study"]         = 100002
dict_object["Case"]          = 100003
dict_object["CaseInProcess"] = 100004

dict_object["DATAFolder"]             = 100010
dict_object["REFERENCEDATAFolder"]    = 100011
dict_object["REFERENCEDATAFile"]      = 100012
dict_object["DATAFile"]               = 100013
dict_object["DRAFTFolder"]            = 100014
dict_object["DATADRAFTFile"]          = 100015
dict_object["DATAPyFile"]             = 100016


dict_object["DATALaunch"]             = 100018
dict_object["DATAfileXML"]            = 100019

dict_object["SRCFolder"]          = 100020
dict_object["SRCFile"]            = 100021
dict_object["SRCDRAFTFile"]       = 100022
dict_object["LOGSRCFile"]         = 100023
dict_object["USERSFolder"]        = 100024
dict_object["USRSRCFile"]         = 100025

dict_object["RESUSubErrFolder"]  = 100029
dict_object["RESUFolder"]     = 100030
dict_object["RESUFile"]       = 100031
dict_object["RESUSubFolder"]  = 100032
dict_object["RESSRCFolder"]   = 100033
dict_object["RESSRCFile"]     = 100034
dict_object["HISTFolder"]     = 100035
dict_object["HISTFile"]       = 100036
dict_object["PRETFolder"]     = 100037
dict_object["SUITEFolder"]    = 100038
dict_object["RESMEDFile"]     = 100039
dict_object["RESXMLFile"]     = 100040
dict_object["POSTPROFolder"]  = 100041
dict_object["RESENSIGHTFile"] = 100042
dict_object["RESUPNGFile"]    = 100043

dict_object["SCRPTFolder"]    = 100050
dict_object["SCRPTLanceFile"] = 100051
dict_object["SCRPTScriptFile"]= 100052
dict_object["SCRPTFile"]      = 100053

dict_object["MESHFolder"]     = 100070
dict_object["MEDFile"]        = 100071
dict_object["DESFile"]        = 100072
dict_object["MESHFile"]       = 100073
dict_object["DATFile"]        = 100074
dict_object["CGNSFile"]       = 100075
dict_object["GeomFile"]       = 100076
dict_object["CaseFile"]       = 100077
dict_object["NeuFile"]        = 100078
dict_object["MSHFile"]        = 100079
dict_object["HexFile"]        = 100080
dict_object["UnvFile"]        = 100081
dict_object["SYRMESHFile"]    = 100082

dict_object["POSTFolder"]     = 100090
dict_object["POSTFile"]       = 100091

#Model objects for COUPLING with SYRTHES CODE
dict_object["CouplingFilePy"]           = 100100
dict_object["RESU_COUPLINGFolder"]      = 100101
dict_object["SYRCaseFolder"]            = 100102
dict_object["SyrthesFile"]              = 100103
dict_object["SyrthesSydFile"]           = 100104
dict_object["CouplingRuncase"]          = 100105
dict_object["RESU_COUPLINGSubFolder"]   = 100106
dict_object["RESUSubFolderSYR"]         = 100107
dict_object["SRCSYRFolder"]             = 100108
dict_object["USRSRCSYRFile"]            = 100109
dict_object["CouplingStudy"]            = 100110
dict_object["OpenSyrthesCaseFile"]      = 100111
#-------------------------------------------------------------------------------
# Definition of the icon of objects to represent in the Object Browser.
# Attribut "AttributePixMap" for the related SObject.
#-------------------------------------------------------------------------------

icon_collection = {}

icon_collection[dict_object["OtherFile"]]      = "CFDSTUDY_UNKNOWN_OBJ_ICON"
icon_collection[dict_object["OtherFolder"]]    = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["Study"]]          = "CFDSTUDY_STUDY_OBJ_ICON"
icon_collection[dict_object["Case"]]           = "CFDSTUDY_CASE_OBJ_ICON"
icon_collection[dict_object["CaseInProcess"]]  = "CFDSTUDY_CASE_IN_PROC_OBJ_ICON"

icon_collection[dict_object["DATAFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["DATAFile"]]       = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["DATAPyFile"]]       = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["DRAFTFolder"]]    = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["REFERENCEDATAFolder"]] = "CFDSTUDY_FOLDER_OBJ_ICON"

icon_collection[dict_object["DATADRAFTFile"]]  = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"

icon_collection[dict_object["REFERENCEDATAFile"]] = "CFDSTUDY_DOCUMENT_OBJ_ICON"

icon_collection[dict_object["DATALaunch"]]     = "CFDSTUDY_EXECUTABLE_OBJ_ICON"
icon_collection[dict_object["DATAfileXML"]]    = "CFDSTUDY_DATA_XML_FILE_OBJ_ICON"

icon_collection[dict_object["SRCFolder"]]      = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["SRCFile"]]        = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["SRCDRAFTFile"]]   = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["LOGSRCFile"]]     = "CFDSTUDY_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["USERSFolder"]]    = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["USRSRCFile"]]     = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"

icon_collection[dict_object["RESUFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESUFile"]]       = "CFDSTUDY_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["RESUSubFolder"]]  = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESUSubErrFolder"]]  = "CFDSTUDY_FOLDER_RED_OBJ_ICON"
icon_collection[dict_object["RESSRCFolder"]]   = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESSRCFile"]]     = "CFDSTUDY_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["HISTFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["HISTFile"]]       = "POST_FILE_ICON"
icon_collection[dict_object["PRETFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["SUITEFolder"]]    = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESMEDFile"]]     = "VISU_OBJ_ICON"
icon_collection[dict_object["RESXMLFile"]]     = "CFDSTUDY_EXECUTABLE_OBJ_ICON"
icon_collection[dict_object["POSTPROFolder"]]  = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESENSIGHTFile"]] = "VISU_OBJ_ICON"
icon_collection[dict_object["RESUPNGFile"]]    = "VIEW_ACTION_ICON"

icon_collection[dict_object["SCRPTFolder"]]    = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["SCRPTLanceFile"]] = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["SCRPTScriptFile"]]= "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["SCRPTFile"]]      = "CFDSTUDY_DOCUMENT_OBJ_ICON"

icon_collection[dict_object["MESHFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["MEDFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["MESHFile"]]       = "CFDSTUDY_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["DESFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["DATFile"]]        = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["CGNSFile"]]       = "MESH_OBJ_ICON"
icon_collection[dict_object["GeomFile"]]       = "MESH_OBJ_ICON"
icon_collection[dict_object["CaseFile"]]       = "MESH_OBJ_ICON"
icon_collection[dict_object["NeuFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["MSHFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["HexFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["UnvFile"]]        = "MESH_OBJ_ICON"
icon_collection[dict_object["SYRMESHFile"]]    = "MESH_OBJ_ICON"

icon_collection[dict_object["POSTFolder"]]     = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["POSTFile"]]       = "CFDSTUDY_DOCUMENT_OBJ_ICON"

#Icons for coupling with SYRTHES CODE
icon_collection[dict_object["SYRCaseFolder"]]         = "SYRTHES_CASE_ICON"
icon_collection[dict_object["SyrthesFile"]]           = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["SyrthesSydFile"]]        = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["CouplingFilePy"]]        = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["CouplingRuncase"]]       = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["RESU_COUPLINGFolder"]]   = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESU_COUPLINGSubFolder"]]= "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["RESUSubFolderSYR"]]      = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["SRCSYRFolder"]]          = "CFDSTUDY_FOLDER_OBJ_ICON"
icon_collection[dict_object["USRSRCSYRFile"]]         = "CFDSTUDY_EDIT_DOCUMENT_OBJ_ICON"
icon_collection[dict_object["CouplingStudy"]]         = "CFDSTUDY_STUDY_OBJ_ICON"

#-------------------------------------------------------------------------------
# ObjectTR is a convenient object for traduction purpose
#-------------------------------------------------------------------------------

ObjectTR = QObject()

#-------------------------------------------------------------------------------
# CORBA
#-------------------------------------------------------------------------------

# init ORB
orb = CORBA.ORB_init([''], CORBA.ORB_ID)

# create naming service instance
naming_service = SALOME_NamingServicePy_i(orb)

# create life cycle CORBA instance
lcc = LifeCycleCORBA(orb)

#-------------------------------------------------------------------------------
# Get Study manager
#-------------------------------------------------------------------------------
# Interface to manipulate SALOME Studies. You will find in this interface the
# methods to create, open, close, and save a SALOME Study. Since a SALOME
# session is multi-document, you will also find methods allowing to navigate
# through a collection of SALOME Studies currently present in a session.

obj = naming_service.Resolve('/myStudyManager')
studyManager = obj._narrow(SALOMEDS.StudyManager)

#-------------------------------------------------------------------------------
# Internal methods
#-------------------------------------------------------------------------------

def _getEngine():
    """
    Gets component engine.

    @return: engine part of the module CFDSTUDY.
    @rtype: C{Engine}
    """
    import CFDSTUDY_ORB
    engine = lcc.FindOrLoadComponent("FactoryServerPy", __MODULE_NAME__)
    return engine


def _getStudy():
    """
    Gets active SALOME Study. C{Study} is a warehouse of data. It can be understood as a document,
    the data storage of the upper level. SALOME Study contains data of multiple components, it's a
    single document for all components. Most of operations on a SALOME Study object are handled by
    C{StudyManager} and C{StudyBuilder} interfaces.

    @return: active Study.
    @rtype: C{Study}
    """
    # getStudyId() -> get study associated to component instance
    # return -1: not initialised (Internal Error)
    #         0: multistudy component instance
    #        >0: study id associated to this instance

    studyId = sgPyQt.getStudyId()
    study = studyManager.GetStudyByID(studyId)
    return study


def _getStudy_Id(studyName) :
    s = studyManager.GetStudyByName(studyName)
    return s._get_StudyId()


def _getNewBuilder():
    study   = _getStudy()
    builder = study.NewBuilder()
    return builder


def _getComponent():
    """
    Returns the component object if CFDSTUDY is active.

    @return: component object if CFDSTUDY is active.
    @rtype: C{Component} or C{None}
    """
    study = _getStudy()
    return study.FindComponent(__MODULE_NAME__)


def _hasChildren(sobj):
    """
    Returns 1 if object has children.

    @type sobj: C{SObject}
    @param sobj: branch of the tree
    @return: 1 if I{sobj} has children, 0 if not.
    @rtype: C{int}
    """
    if sobj:
        study = _getStudy()
        iter  = study.NewChildIterator(sobj)
        while iter.More():
            name = iter.Value().GetName()
            if name:
                return 1
            iter.Next()
    return 0


def _findOrCreateComponent():
    """
    Finds or creates component object, i.e. root of the tree.

    @return: the root C{SObject} for the Object browser representation.
    @rtype: C{SObject}
    """
    study = _getStudy()
    father = study.FindComponent(__MODULE_NAME__)
    if father is None:
        builder = study.NewBuilder()
        father = builder.NewComponent(__MODULE_NAME__)
        attr = builder.FindOrCreateAttribute(father, "AttributeName")
        attr.SetValue(__MODULE_NAME__)
        attr = builder.FindOrCreateAttribute(father, "AttributeLocalID")
        attr.SetValue(__MODULE_ID__)
        attr = builder.FindOrCreateAttribute(father, "AttributePixMap")
        attr.SetPixMap("CFDSTUDY.png")

        try:
            builder.DefineComponentInstance(father, _getEngine())
        except:
            pass

    return father

def _SetCaseLocation(theCasePath):
    log.debug("_SetCaseLocation")
    study         = _getStudy()
    builder       = study.NewBuilder()
    father        = _findOrCreateComponent()
    theStudyPath  = os.path.dirname(theCasePath)

    aCaseName      = os.path.basename(theCasePath)
    studyObject   = FindStudyByPath(theStudyPath)
    if studyObject == None:
        if CFDSTUDYGUI_Commons.isaCFDStudy(theStudyPath):
            studyObject  = builder.NewObject(father)
            attr = builder.FindOrCreateAttribute(studyObject, "AttributeLocalID")
            attr.SetValue(dict_object["Study"])
            attr = builder.FindOrCreateAttribute(studyObject, "AttributePixMap")
            attr.SetPixMap(str(ObjectTR.tr("CFDSTUDY_STUDY_OBJ_ICON")))
            attr = builder.FindOrCreateAttribute(studyObject, "AttributeName")
            attr.SetValue(os.path.basename(theStudyPath))
            attr = builder.FindOrCreateAttribute(studyObject, "AttributeComment")
            attr.SetValue(os.path.dirname(theStudyPath))
    _CreateItem(studyObject,aCaseName)
    caseObject = getSObject(studyObject,aCaseName)
    UpdateSubTree(caseObject)

def _SetStudyLocation(theStudyPath, theCaseNames,theCreateOpt,
                      theCopyOpt, theNameRef = "", theSyrthesOpt =False, theSyrthesCase = "",theNprocs=""):
    """
    Constructs the tree representation of a CFD study (with the
    associated cases) for the Object Browser. All branch of the tree is
    an C{SObject} object.

    @type theStudyPath: C{String}
    @param theStudyPath: unix path of the CFD study.
    @type theCaseNames: C{String}
    @param theCaseNames: unix pathes of the new CFD cases to be build.
    """
    log.debug("_SetStudyLocation")

    iok = True
    if theCopyOpt:
        if not os.path.exists(theNameRef):
            raise ValueError("reference case is not a repository")
    if os.path.exists(theStudyPath) :

        if theCreateOpt:
            mess = cfdstudyMess.trMessage(ObjectTR.tr("STUDY_DIRECTORY_ALREADY_EXISTS"),[""])
            cfdstudyMess.criticalMessage(mess)
            return False
    if theCreateOpt or (not theCreateOpt and theCaseNames!=""):
        iok = _CallCreateScript(theStudyPath, theCreateOpt, theCaseNames,
                                theCopyOpt, theNameRef, theSyrthesOpt, theSyrthesCase)

    study   = _getStudy()
    builder = study.NewBuilder()
    father  = _findOrCreateComponent()
    studyObject = FindStudyByPath(theStudyPath)

    if studyObject == None:
        #obtain name and dir for new study
        lst = os.path.split(theStudyPath)
        aStudyDir = lst[0]
        aStudyName = lst[1]
        if aStudyName == "":
            raise ValueError("Empty Study Name!")
        if aStudyDir == "":
            raise ValueError("Empty Study Directory!")

        studyObject  = builder.NewObject(father)
        attr = builder.FindOrCreateAttribute(studyObject, "AttributeLocalID")
        if CFDSTUDYGUI_Commons.isaSaturneSyrthesCouplingStudy(theStudyPath):
            attr.SetValue(dict_object["CouplingStudy"])
        elif CFDSTUDYGUI_Commons.isaCFDStudy(theStudyPath):
            attr.SetValue(dict_object["Study"])
        else:
            return False
        attr = builder.FindOrCreateAttribute(studyObject, "AttributePixMap")
        attr.SetPixMap(str(ObjectTR.tr("CFDSTUDY_STUDY_OBJ_ICON")))
        attr = builder.FindOrCreateAttribute(studyObject, "AttributeName")
        attr.SetValue(aStudyName)
        attr = builder.FindOrCreateAttribute(studyObject, "AttributeComment")
        attr.SetValue(aStudyDir)

    if iok:
        UpdateSubTree(studyObject)
        if "coupling_parameters.py" in os.listdir(theStudyPath):
            # Update SYRTHES PATH into coupling_parameters.py file
            replaceOrInsertCouplingPath(os.path.join(theStudyPath,
                                                     "coupling_parameters.py"))
        if "runcase" in os.listdir(theStudyPath) and theCreateOpt:
            # Update the coupling study PATH into runcase file
            runcase = cs_runcase.runcase(os.path.join(theStudyPath, "runcase"))
            if theNprocs != "":
                runcase.set_nprocs(theNprocs)
            runcase.save()

        if CFDSTUDYGUI_Commons.isaCFDStudy(theStudyPath):
            theStudy = FindStudyByPath(theStudyPath)
            theStudyCaseNameList = GetCaseNameList(theStudy)
            if theStudyCaseNameList != [] or theStudyCaseNameList != None :
                for theStudyCase in theStudyCaseNameList :
                    runcasePath = os.path.join(theStudyPath, theStudyCase,
                                               "SCRIPTS/runcase")
                    if os.path.exists(runcasePath):
                        runcase = cs_runcase.runcase(runcasePath)
                        runcase.save()
    return iok


def replaceOrInsertCouplingPath(The_coupling_parameters_path):
    """
    coupling_parameters.py file is updated according to syrthes path needed
    """
    f = open(The_coupling_parameters_path,"r")
    l = f.readlines()
    f.close()
    boo = False
    if not os.getenv("SYRTHES4_HOME"):
        mess = cfdstudyMess.trMessage(ObjectTR.tr("COUPLING_SYRTHES4_HOME_MISSING"),[])
        cfdstudyMess.criticalMessage(mess)
        return
    lineToInsert = "sys.path.insert(1,"+"'"+os.path.join(os.path.join(os.environ["SYRTHES4_HOME"],"share"),"syrthes")+os.sep+"')\n"
    for i in l :
        if  i.startswith("sys.path.insert"):
            ind = l.index(i)
            l[ind] = lineToInsert
            boo = True
    if not boo:
        l.insert(1,lineToInsert)
    f = open(The_coupling_parameters_path,"w")
    f.writelines(l)
    f.close()

def _CallCreateScript(theStudyPath, isCreateStudy, theCaseNames,
                      theCopyOpt, theNameRef, theSyrthesOpt, theSyrthesCase):
    """
    Builds new CFD study, and/or new cases on the file system.

    @type theStudyPath: C{String}
    @param theStudyPath: unix path of the CFD study.
    @type isCreateStudy: C{True} or C{False}
    @param isCreateStudy: if C{True} build the new CFD study, if C{False}, only cases have to be build.
    @type theCaseNames: C{String}
    @param theCaseNames: unix pathes of the new CFD cases to be build.
    """
    log.debug("_CallCreateScript")
    mess = ""
    scrpt, c ,mess = BinCode()
    if mess == "" :
        curd = os.getcwd()

        start_dir = ""
        if isCreateStudy:
            fatherdir,etude = os.path.split(theStudyPath)
            start_dir = fatherdir
        else:
            start_dir = theStudyPath

        args = [scrpt]
        args.append('create')

        if isCreateStudy:
            args.append("--study")
            args.append(theStudyPath)
        if theCaseNames != "":
            for i in theCaseNames.split(' '):
                args.append("--case")
                args.append(os.path.join(theStudyPath,i))
        if theCopyOpt:
            args.append("--copy-from")
            args.append(theNameRef)

        if theSyrthesOpt:
            args.append("--syrthes")
            args.append(os.path.join(theStudyPath,theSyrthesCase))

        runCommand(args, start_dir, "")

    else:
        cfdstudyMess.criticalMessage(mess)
    return mess == ""

def updateCasePath(theCasePath):

    mess = ""
    scrpt, c ,mess = BinCode()
    if mess == "" :
        curd = os.getcwd()

        start_dir = ""
        args = [scrpt]
        args.append('create')
        args.append(theCasePath)
        args.append('--import-only')
        runCommand(args, start_dir, "")
    else:
        cfdstudyMess.criticalMessage(mess)
    return mess == ""


def _UpdateStudy():
    """
    Updates CFD study tree of data from the root.
    """
    study   = _getStudy()
    component = study.FindComponent(__MODULE_NAME__)
    if component == None:
        return

    iter  = study.NewChildIterator(component)
    while iter.More():
        _RebuildTreeRecursively(iter.Value())
        iter.Next()


def UpdateSubTree(theObject=None):
    """
    Updates CFD study sub-tree from the argument object.

    @type theObject: C{SObject}
    @param theObject: branch of a tree of data to update.
    """
    log.debug("UpdateSubTree")
    if theObject != None:
        log.debug("_RebuildTreeRecursively -> path: %s" % _GetPath(theObject))
        _RebuildTreeRecursively(theObject)
    else:
        _UpdateStudy()
    # --- update object browser from a thread different of the main thread is not safe !
    studyId = sgPyQt.getStudyId()
    sgPyQt.updateObjBrowser(studyId,1)

def closeCFDStudyTree(theObject):
    """
    Close a CFD Study from the Object browser
    """
    log.debug("closeCFDStudyTree")
    if theObject == None:
        return
    study   = _getStudy()
    builder = study.NewBuilder()
    builder.RemoveObjectWithChildren(theObject)
    return

def _RebuildTreeRecursively(theObject):
    """
    Builds or rebuilds a branch of the tree of data for the Object Browser.

    @type theObject: C{SObject}
    @param theObject: branch of a tree of data
    """
    # SObject is the main constituent of SALOMEDS-based data structure.
    # If you are familiar with CAF (Cascade Application Framework) - the
    # analogy of SObject would be TDF_Label class. It can be understood as
    # a branch of a tree of data, or as a record in a database table. Usually
    # it does not store the data itself, it uses child Attributes - successors
    # of SALOMEDS::GenericAttribute - for storing specific data, properties
    # of the object.
    #
    # type(SObject) -> SALOMEDS._objref_SObject instance
    #
    if theObject == None:
        return
    log.debug("_RebuildTreeRecursively -> %s childs: %s" % (theObject.GetName(), ScanChildNames(theObject,  ".*")))
    theObjectPath = _GetPath(theObject)

    if theObjectPath == None:
        return

    study   = _getStudy()
    builder = study.NewBuilder()
    attr = builder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    # clean the SObject, if the corresponding file or directory
    # does not exist any more in the file system

    if os.path.isfile(theObjectPath) and attr.Value() == dict_object["MEDFile"]:
        return

    if not os.path.isdir(theObjectPath) and not os.path.isfile(theObjectPath):
        builder.RemoveObjectWithChildren(theObject)
        return

    # build a list of file from the file system
    dirList = _GetDirList(theObject)
    # build a list and a dictionary of childs SObject from the current SObject
    objList = []
    objMap  = {}

    iter  = study.NewChildIterator(theObject)
    while iter.More():
        v = iter.Value()
        n = v.GetName()
        objList.append(n)
        objMap[n] = v
        iter.Next()

    objList.sort()
    objIndex = 0
    dirIndex = 0
    # Case with empty list of existing SObject: every SObject must be build
    if len(objList) == 0:
        while dirIndex < len(dirList):
            #append new objects
            if Trace() : print("Whole append new Item: ", dirList[dirIndex])
            _CreateObject(theObject, builder, dirList[dirIndex])
            dirIndex+=1

    # Case with empty list of file: every SObject must be clean
    elif len(dirList) == 0:
        builder.RemoveObjectWithChildren(theObject)
        log.debug("_RebuildTreeRecursively 3: %s childs: %s" % (theObject.GetName(), ScanChildNames(theObject,  ".*")))

    else:
        objEnd = False
        dirEnd = False
        while True:
            objName = objList[objIndex]
            dirName = dirList[dirIndex]

            if dirName < objName:
                if not dirEnd:
                    #append new object
                    if Trace(): print("1 Append new Item: ", dirName)
                    log.debug("_RebuildTreeRecursively 4: dirName = %s objName = %s" %(dirName,objName))
                    _CreateObject(theObject, builder, dirName)
                    dirIndex+=1

                    if objEnd and dirIndex == len(dirList):
                        break
                else:
                    if Trace(): print("1 Remove Item from tree: ", objName)
                    builder.RemoveObjectWithChildren(objMap[objName])
                    objIndex+=1

                    if objIndex == len(objList):
                        break

            elif dirName > objName:
                #remove object if no end
                if not objEnd:
                    if Trace(): print("2 Remove Item from tree: ", objName)
                    builder.RemoveObjectWithChildren(objMap[objName])
                    objIndex+=1

                    if dirEnd and objIndex == len(objList):
                        break

                else:
                    #append new item at the end
                    if Trace(): print("2 Append new Item: ", dirName)
                    log.debug("_RebuildTreeRecursively 5: dirName = %s objName = %s" %(dirName,objName))
                    _CreateObject(theObject, builder, dirName)
                    dirIndex+=1

                    if dirIndex == len(dirList):
                        break

            else:
                # no changes
                _FillObject(objMap[objName], theObject, builder)
                dirIndex+=1
                objIndex+=1

            if dirIndex != len(dirList) or objIndex != len(objList):
                if dirIndex == len(dirList):
                    dirEnd = True
                    dirIndex-=1

                if objIndex == len(objList):
                    objEnd = True
                    objIndex-=1

            if dirIndex == len(dirList) and objIndex == len(objList):
                break

    # recursively calling
    iter  = study.NewChildIterator(theObject)
    while iter.More():
        if iter.Value().GetName():
            _RebuildTreeRecursively(iter.Value())
        iter.Next()
    log.debug("_RebuildTreeRecursively -> %s END" % (theObject.GetName()))


def _CreateObject(theFather, theBuilder, theName):
    """
    Creates a child branch in the tree from the father branch I{theFather}.
    Sets the AttributeName value of this new child with the name theName of the child object
    Calls _FillObject which sets the AttributeLocalID of the child
    _FillObject calls _setIcon which sets the AttributePixMap and AttributeComment of the child object

    Result : an object entry in the Object Browser with AttributeName, AttributeLocalID, AttributePixMap, AttributeComment
    @type theFather: C{SObject}
    @param theFather: branch of the tree to add a child.
    @type theBuilder: C{SUIT_Study}
    @param theBuilder: C{SObject} constructor.
    @type theName: C{String}
    @param theName: AttributeName of the new child branch.
    """
    log.debug("_CreateObject: %s" % theName)
    newChild = theBuilder.NewObject(theFather)
    attr = theBuilder.FindOrCreateAttribute(newChild, "AttributeName")
    attr.SetValue(theName)
    _FillObject(newChild, theFather, theBuilder)


def _CreateItem(theFather,theNewName) :
    """
    Creates a child with name theNewName under theFather root into Object Browser
    @type theFather: C{SObject}
    @type theNewName : C{String}
    """
    log.debug("_CreateItem: NewItem = %s with Parent = %s" % (theNewName,theFather.GetName()))
    if theNewName not in ScanChildNames(theFather,  ".*") :
        theBuilder = _getNewBuilder()
        _CreateObject(theFather, theBuilder, theNewName)


def getNameCodeFromXmlCasePath(XMLCasePath) :
    """
    """
    code = ""
    if os.path.isfile(XMLCasePath):
        fd = os.open(XMLCasePath,os.O_RDONLY)
        f = os.fdopen(fd)
        l1 = f.readline()
        if l1.startswith('''<?xml version="1.0" encoding="utf-8"?><Code_Saturne_GUI''') or l1.startswith('''<?xml version="1.0" encoding="utf-8"?><NEPTUNE_CFD_GUI'''):
            if "Code_Saturne" in l1 :
                code = "Code_Saturne"
            elif "NEPTUNE_CFD" in l1 :
                code = "NEPTUNE_CFD"
        elif l1.startswith('''<?xml version="1.0" encoding="utf-8"?>''') :
            l2 = f.readline()
            if l2.startswith('''<Code_Saturne_GUI''') or l2.startswith('''<NEPTUNE_CFD_GUI'''):
                if "Code_Saturne" in l2 :
                    code = "Code_Saturne"
                elif "NEPTUNE_CFD" in l2 :
                    code = "NEPTUNE_CFD"
            else :
                mess = cfdstudyMess.trMessage(ObjectTR.tr("XML_DATA_FILE"),[XMLCasePath])
                cfdstudyMess.warningMessage(mess)
        else :
            mess = cfdstudyMess.trMessage(ObjectTR.tr("XML_DATA_FILE"),[XMLCasePath])
            cfdstudyMess.warningMessage(mess)
        f.close()
    return code


def _FillObject(theObject, theParent, theBuilder):
    """
    Creates the attribute "AttributeLocalID" for the branch I{theObject}.
    This attribute keeps the type of the I{theObject}.

    @type theObject: C{SObject}
    @param theObject: branch of the tree to add an attribut.
    @type theParent: C{SObject}
    @param theParent: parent of the branch I{theObject}.
    @type theBuilder: C{SUIT_Study}
    @param theBuilder: C{SObject} constructor for create an attribut.
    """
    attr = theBuilder.FindOrCreateAttribute(theParent, "AttributeLocalID")
    parentId = attr.Value()
    name = theObject.GetName()
    objectId = dict_object["OtherFile"]
    path = os.path.join(_GetPath(theParent), name)
    #log.debug("_FillObject: %s" % name)

    # Parent is study
    if parentId == dict_object["Study"] or parentId == dict_object["CouplingStudy"]:
        if Trace(): print("_FillObject : parent is Study ", theParent.GetName())
        #check for case
        if os.path.isdir(path):
            if CFDSTUDYGUI_Commons.isaCFDCase(path):
                objectId = dict_object["Case"]
            else:
                boo = False
                dirList = os.listdir(path)
                for i in dirList:
                    if re.match(".*\.syd$", i) or re.match(".*\.syd_example$", i): boo = True
                if boo :
                    objectId = dict_object["SYRCaseFolder"]
                else:
                    if name == "MESH":
                        objectId = dict_object["MESHFolder"]
                    elif name == "POST":
                        objectId = dict_object["POSTFolder"]
                    else:
                        objectId = dict_object["OtherFolder"]

        if name == "coupling_parameters.py":
            objectId = dict_object["CouplingFilePy"]
        elif name == "runcase":
            objectId = dict_object["CouplingRuncase"]
        elif name == "RESU_COUPLING":
            objectId = dict_object["RESU_COUPLINGFolder"]
    #parent is Syrthes Case
    elif parentId == dict_object["SYRCaseFolder"]:
        if os.path.isdir(path):
            if name == "usr_examples":
                objectId = dict_object["SRCSYRFolder"]
        if name in ["Makefile","syrthes.py","user_cond.c"]:
            objectId = dict_object["SyrthesFile"]
        if re.match(".*\.syd$", name) or re.match(".*\.syd_example$", name) :
            objectId = dict_object["SyrthesSydFile"]
    #parent is Syrthes user examples
    elif parentId == dict_object["SRCSYRFolder"]:
        if re.match(".*\.c$", name):
            objectId = dict_object["USRSRCSYRFile"]

    #parent is Case
    elif parentId == dict_object["Case"]:
        if os.path.isdir(path):
            if name == "DATA":
                objectId = dict_object["DATAFolder"]
            elif name == "SRC":
                objectId = dict_object["SRCFolder"]
            elif name == "RESU":
                objectId = dict_object["RESUFolder"]
            elif name == "SCRIPTS":
                objectId = dict_object["SCRPTFolder"]
            else:
                objectId = dict_object["OtherFolder"]

    # parent is DATA folder
    elif parentId == dict_object["DATAFolder"]:
        if os.path.isdir(path):
            if name == "REFERENCE":
                objectId = dict_object["REFERENCEDATAFolder"]
            if name == "DRAFT":
                objectId = dict_object["DRAFTFolder"]
        else:
            if name[0:10] == "SaturneGUI" or name[0:10] == "NeptuneGUI":
                objectId = dict_object["DATALaunch"]
            elif re.match("^dp_", name) or re.match("^meteo",name) or re.match("^cs_", name):
                objectId = dict_object["DATAFile"]
            elif re.match(".*\.py$", name):
                objectId = dict_object["DATAPyFile"]
            else:
                if os.path.isfile(path):
                    fd = os.open(path , os.O_RDONLY)
                    f = os.fdopen(fd)
                    l1 = f.readline()
                    if l1.startswith('''<?xml version="1.0" encoding="utf-8"?><Code_Saturne_GUI''') or l1.startswith('''<?xml version="1.0" encoding="utf-8"?><NEPTUNE_CFD_GUI'''):
                        objectId = dict_object["DATAfileXML"]
                    elif l1.startswith('''<?xml version="1.0" encoding="utf-8"?>''') :
                        l2 = f.readline()
                        if l2.startswith('''<Code_Saturne_GUI''') or l2.startswith('''<NEPTUNE_CFD_GUI'''):
                            objectId = dict_object["DATAfileXML"]
                    f.close()

    # parent is DRAFT folder
    elif parentId == dict_object["DRAFTFolder"]:
        draftParentFolder = os.path.basename(_GetPath(theParent.GetFather()))
        if os.path.isfile(path):
            if draftParentFolder == "DATA":
                if re.match("^dp_", name) or re.match("^meteo",name) or re.match("^cs_", name):
                    objectId = dict_object["DATADRAFTFile"]
            elif draftParentFolder == "SRC":
                if re.match(".*\.[fF]$", name) or \
                    re.match(".*\.[fF]90$", name) or \
                    re.match(".*\.for$", name) or \
                    re.match(".*\.FOR$", name):
                    objectId = dict_object["SRCDRAFTFile"]
                elif re.match(".*\.c$", name):
                    objectId = dict_object["SRCDRAFTFile"]
                elif re.match(".*\.cxx$", name) or \
                     re.match(".*\.cpp$", name):
                    objectId = dict_object["SRCDRAFTFile"]
                elif re.match(".*\.h$", name) or \
                     re.match(".*\.hxx$", name) or \
                     re.match(".*\.hpp$", name):
                    objectId = dict_object["SRCDRAFTFile"]
        elif os.path.isdir(path):
            objectId = dict_object["OtherFolder"]

    # parent is REFERENCE folder into DATA folder
    elif parentId == dict_object["REFERENCEDATAFolder"]:
        if os.path.isfile(path):
            if re.match("^dp_", name) or re.match("^meteo",name) or re.match("^cs_", name):
                objectId = dict_object["REFERENCEDATAFile"]
        elif os.path.isdir(path):
            objectId = dict_object["OtherFolder"]

    # parent is SCRIPTS folder
    elif parentId == dict_object["SCRPTFolder"]:
        if os.path.isdir(path):
            objectId = dict_object["OtherFolder"]
        else:
            if name[0:7] == "runcase":
                objectId = dict_object["SCRPTLanceFile"]
            else:
                if os.path.isfile(path):
                    fd = os.open(path , os.O_RDONLY)
                    f = os.fdopen(fd)
                    l = f.readline()
                    f.close()
                    if l[0:2] == "#!":
                        objectId = dict_object["SCRPTScriptFile"]
                    else:
                        objectId = dict_object["SCRPTFile"]
                else:
                    objectId = dict_object["OtherFile"]

    # parent is MESH folder
    elif parentId == dict_object["MESHFolder"]:
        if os.path.isdir(path):
            objectId = dict_object["OtherFolder"]
        else:
            if re.match(".*\.des$", name):
                objectId = dict_object["DESFile"]
            elif re.match(".*\.med$", name):
                objectId = dict_object["MEDFile"]
            elif re.match(".*\.dat$", name):
                objectId = dict_object["DATFile"]
            elif re.match(".*\.cgns$", name):
                objectId = dict_object["CGNSFile"]
            elif re.match(".*\.ngeom$", name):
                objectId = dict_object["GeomFile"]
            elif re.match(".*\.case$", name):
                objectId = dict_object["CaseFile"]
            elif re.match(".*\.neu$", name):
                objectId = dict_object["NeuFile"]
            elif re.match(".*\.msh$", name):
                objectId = dict_object["MSHFile"]
            elif re.match(".*\.hex$", name):
                objectId = dict_object["HexFile"]
            elif re.match(".*\.unv$", name):
                objectId = dict_object["UnvFile"]
            elif re.match(".*\.syr$", name):
                objectId = dict_object["SYRMESHFile"]
            else:
                objectId = dict_object["MESHFile"]

    # parent is POST folder
    elif parentId == dict_object["POSTFolder"]:
        if os.path.isdir(path):
            objectId = dict_object["OtherFolder"]
        else:
            objectId = dict_object["POSTFile"]

    # parent is SRC folder
    elif parentId == dict_object["SRCFolder"]:
        if os.path.isfile(path):
            if re.match(".*\.[fF]$", name) or re.match(".*\.[fF]90$", name) \
              or re.match(".*\.for$", name) or re.match(".*\.FOR$", name):
                objectId = dict_object["SRCFile"]
            elif re.match(".*\.c$", name):
                objectId = dict_object["SRCFile"]
            elif re.match(".*\.cpp$", name) or re.match(".*\.cxx$", name):
                objectId = dict_object["SRCFile"]
            elif re.match(".*\.h$", name) or re.match(".*\.hpp$", name) or re.match(".*\.hxx$", name):
                objectId = dict_object["SRCFile"]
            elif re.match(".*\.log$", name):
                objectId = dict_object["LOGSRCFile"]
        elif os.path.isdir(path):
            if name == "REFERENCE" or name == "EXAMPLES" :
                objectId = dict_object["USERSFolder"]
            elif name == "DRAFT":
                objectId = dict_object["DRAFTFolder"]
            else:
                objectId = dict_object["OtherFolder"]

    # parent REFERENCE/base... folder
    elif parentId == dict_object["USERSFolder"]:
        if os.path.isfile(path):
            if re.match(".*\.[fF]$", name) or re.match(".*\.[fF]90$", name) \
              or re.match(".*\.for$", name) or re.match(".*\.FOR$", name):
                objectId = dict_object["USRSRCFile"]
            elif re.match(".*\.c$", name):
                objectId = dict_object["USRSRCFile"]
            elif re.match(".*\.cpp$", name) or re.match(".*\.cxx$", name):
                objectId = dict_object["USRSRCFile"]
            elif re.match(".*\.h$", name) or re.match(".*\.hpp$", name) or re.match(".*\.hxx$", name):
                objectId = dict_object["USRSRCFile"]
            elif re.match(".*\.log$", name):
                objectId = dict_object["LOGSRCFile"]
        elif os.path.isdir(path):
            if name in ("atmo", "base", "cplv", "cfbl", "cogz", \
                        "ctwr", "elec", "fuel", "lagr", "pprt", "rayt"):
                objectId = dict_object["USERSFolder"]
            else:
                objectId = dict_object["OtherFolder"]

    # parent is RESU folder
    elif parentId == dict_object["RESUFolder"]:
        if os.path.isdir(path):
            if "error" in os.listdir(path):
                objectId = dict_object["RESUSubErrFolder"]
            else:
                objectId = dict_object["RESUSubFolder"]

    # parent is RESULT SRC folder
    elif parentId == dict_object["RESSRCFolder"]:
        if os.path.isfile(path):
            if re.match(".*\.[fF]$", name) or re.match(".*\.[fF]90$", name) \
              or re.match(".*\.for$", name) or re.match(".*\.FOR$", name):
                objectId = dict_object["RESSRCFile"]
            elif re.match(".*\.c$", name):
                objectId = dict_object["RESSRCFile"]
            elif re.match(".*\.cpp$", name) or re.match(".*\.cxx$", name):
                objectId = dict_object["RESSRCFile"]
            elif re.match(".*\.h$", name) or re.match(".*\.hpp$", name) or re.match(".*\.hxx$", name):
                objectId = dict_object["RESSRCFile"]

    # parent is RESULT sub folder
    elif parentId == dict_object["RESUSubFolder"] or parentId == dict_object["RESUSubErrFolder"]:
        if os.path.isdir(path):
            if name == "src_neptune" or name == "src_saturne":
                objectId = dict_object["RESSRCFolder"]
            elif name == "monitoring":
                objectId = dict_object["HISTFolder"]
            elif name == "checkpoint":
                objectId = dict_object["SUITEFolder"]
            elif name == "mesh_input":
                objectId = dict_object["PRETFolder"]
            elif name == "partition_output":
                objectId = dict_object["PRETFolder"]
            elif name == "postprocessing":
                objectId = dict_object["POSTPROFolder"]
        else:
            if re.match(".*\.dat$", name) or re.match(".*\.csv$", name):
                objectId = dict_object["HISTFile"]
            elif re.match(".*\.xml$", name):
                objectId = dict_object["RESXMLFile"]
            elif re.match(".*\.log$", name):
                objectId = dict_object["RESUFile"]
            elif re.match("listing$", name):
                objectId = dict_object["RESUFile"]
            elif re.match("error$", name):
                objectId = dict_object["RESUFile"]
            elif re.match(".*\.png$", name):
                objectId = dict_object["RESUPNGFile"]

    elif parentId == dict_object["POSTPROFolder"] :
        if os.path.isfile(path):
            if re.match(".*\.med$", name):
                objectId = dict_object["RESMEDFile"]
            if re.match(".*\.case$", name):
                objectId = dict_object["RESENSIGHTFile"]

    # parent is HIST folder
    elif parentId == dict_object["HISTFolder"]:
        if os.path.isfile(path):
            if re.match(".*\.dat$", name) or re.match(".*\.csv$", name):
                objectId = dict_object["HISTFile"]

    # parent is RESU_COUPLING folder
    elif parentId == dict_object["RESU_COUPLINGFolder"]:
        if os.path.isdir(path):
            objectId = dict_object["RESU_COUPLINGSubFolder"]

    # parent is RESU_COUPLING sub folder
    elif parentId == dict_object["RESU_COUPLINGSubFolder"]:
        if os.path.isdir(path):
            if os.path.isfile(os.path.join(path,"syrthes")):
                objectId = dict_object["RESUSubFolderSYR"]
            else:
                # test if folder is a result cfd folder?
                objectId = dict_object["RESUSubFolder"]

        else:
            if name == "coupling_parameters.py":
                objectId = dict_object["CouplingFilePy"]

    elif parentId == dict_object["RESUSubFolderSYR"]:
        if re.match(".*\.log$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.dat$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.rdt$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.res$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.syr$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.data$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.add$", name):
            objectId = dict_object["RESUFile"]
        if re.match(".*\.c$", name):
            objectId = dict_object["RESUFile"]
        elif re.match("listing$", name):
            objectId = dict_object["RESUFile"]

    if objectId == dict_object["OtherFile"]:
        if re.match(".*\.[fF]$", name) or \
           re.match(".*\.[fF]90$", name) or \
           re.match(".*\.for$", name) or \
           re.match(".*\.FOR$", name):
            if _DetectUSERSObject(theObject) == True:
                if Trace(): print("******************************", path)
                objectId = _DetectSRCObject(theParent)
        elif re.match(".*\.c$", name):
            if _DetectUSERSObject(theObject) == True:
                if Trace(): print("******************************", path)
                objectId = _DetectSRCObject(theParent)
        elif re.match(".*\.cpp$", name) or \
           re.match(".*\.cxx$", name):
            if _DetectUSERSObject(theObject) == True:
                if Trace(): print("******************************", path)
                objectId = _DetectSRCObject(theParent)
        elif re.match(".*\.h$", name) or \
           re.match(".*\.hxx$", name) or \
           re.match(".*\.hpp$", name):
            if _DetectUSERSObject(theObject) == True:
                if Trace(): print("******************************", path)
                objectId = _DetectSRCObject(theParent)

    if objectId == dict_object["OtherFile"]:
        if os.path.isdir(path):
            objectId = dict_object["OtherFolder"]

    log.debug("_FillObject: %s %s" % \
        (name, [k for k, v in dict_object.items() if v == objectId][0]))

    if objectId in (dict_object["OtherFile"],
                    dict_object["OtherFolder"],
                    dict_object["SCRPTFile"],
                    dict_object["MESHFile"],
                    dict_object["DATFile"]):
        study   = _getStudy()
        builder = study.NewBuilder()
        builder.RemoveObjectWithChildren(theObject)
        return

    attr = theBuilder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    attr.SetValue(objectId)

    _SetIcon(theObject, theBuilder)


def _SetIcon(theObject, theBuilder):
    """
    Creates the attribute "AttributePixMap" and "AttributeComment" for the branch I{theObject}.

    @type theObject: C{SObject}
    @param theObject: branch of the tree to add an icon.
    @type theBuilder: C{SUIT_Study}
    @param theBuilder: C{SObject} constructor for create an attribut.
    """
    attr = theBuilder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    id = int(attr.Value())
    if icon_collection[id] == "":
        return
    attr = theBuilder.FindOrCreateAttribute(theObject, "AttributePixMap")
    attr.SetPixMap(str(ObjectTR.tr(icon_collection[id])))
    #check path for link and create new attribute
    if id != dict_object["Case"]:
        path = _GetPath(theObject)
        if os.path.islink(path):
            attr = theBuilder.FindOrCreateAttribute(theObject, "AttributeComment")
            attr.SetValue("->" + os.path.realpath(path))


def _GetPath(theObject):
    """
    Returns the unix path of the branch I{theObject}.

    @type theObject: C{SObject}
    @param theObject: branch of the tree to add an icon.
    @return: unix path of the branch I{theObject}
    @rtype: C{String}
    """
    # check for null object
    # check object from others component
    # check if CFDSTUDY component object

    if _getComponent() == None:
        return ""

    if not theObject or \
           theObject.GetFatherComponent().GetID() != _getComponent().GetID() or \
           theObject.GetID() == _getComponent().GetID():
        return ""

    study   = _getStudy()
    builder = study.NewBuilder()
    path = str(theObject.GetName())
    attr = builder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    if attr.Value() == dict_object["Study"] or attr.Value() == dict_object["CouplingStudy"]:
        dir = builder.FindOrCreateAttribute(theObject, "AttributeComment")
        return os.path.join(dir.Value(), path)

    father = theObject.GetFather()
    attr = builder.FindOrCreateAttribute(father, "AttributeLocalID")
    path = os.path.join(_GetPath(father), path)

    return path


def _GetDirList(theObject):
    """
    Returns the unix pathes of the directories which are child of the branch I{theObject}.

    @type theObject: C{SObject}
    @param theObject: branch of the tree.
    @return: list of unix pathes of directory.
    @rtype: C{List} of C{String}
    """
    study   = _getStudy()
    builder = study.NewBuilder()

    path = _GetPath(theObject)
    attr = builder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    lst = []
    if os.path.isdir(path):
        lst = os.listdir(path)
    lst.sort()
    return lst


def _DetectUSERSObject(theObject):
    """
    Search if the branch I{theObject} represents the USERS folder.

    @type theObject: C{SObject}
    @param theObject: branch of the tree.
    @return: C{True} if the I{theObject} represents the USERS folder
    @rtype: C{True} or C{False}
    """
    study   = _getStudy()
    builder = study.NewBuilder()
    cur = theObject.GetFather()
    attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")

    while True:
        if attr.Value() == dict_object["USERSFolder"]:
            return True
        elif attr.Value() == dict_object["Study"]:
            return False

        cur = cur.GetFather()
        attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")

    return False


def _DetectSRCObject(theObject):
    """
    Returns the type of the branch I{theObject} which represents
    the files in the SRC folder.

    @type theObject: C{SObject}
    @param theObject: branch of the tree.
    @return: type of the I{theObject} which represents files in the SRC folder.
    @rtype: C{int}
    """
    study   = _getStudy()
    builder = study.NewBuilder()
    cur = theObject.GetFather()
    attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")

    while True:
        if attr.Value() == dict_object["SRCFolder"]:
            return dict_object["USRSRCFile"]
        if attr.Value() == dict_object["SRCFolder"]:
            return dict_object["USRSRCFile"]

        cur = cur.GetFather()
        attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")

    return dict_object["USRSRCFile"]


def GetCase(theObject):
    """
    Returns the case to which belongs the I{theObject}.

    @type theObject: C{SObject}
    @param theObject: file or folder we want to know the case.
    @return: case to which belongs the I{theObject}.
    @rtype: C{SObject}
    """
    if theObject == None:
        return None

    study   = _getStudy()
    builder = study.NewBuilder()
    cur = theObject

    while cur:
        attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")
        if Trace():
            print("attr:",attr)
            print("Value", attr.Value())
        value = attr.Value()
        if value == dict_object["Case"]:
            return cur
        elif value == dict_object["Study"] or value == dict_object["CouplingStudy"] \
             or value == __MODULE_ID__ \
             or value == 0:
            return None

        cur = cur.GetFather()

    return None


def GetFirstStudy():
    """
    Returns the first CFD study loaded in the Object Browser.

    @return: first study of the tree.
    @rtype: C{SObject}
    """
    study = _getStudy()

    component = _getComponent()
    if component == None:
        return None

    iter  = study.NewChildIterator(component)
    return iter.Value()


def GetStudyByObj(theObject):
    """
    Returns the CFD study to which belongs the I{theObject}.

    @type theObject: C{SObject}
    @param theObject: file or folder we want to know the father's CFD study.
    @return: study to which belongs the I{theObject}.
    @rtype: C{SObject}
    """
    if theObject == None:
        return None

    study   = _getStudy()
    builder = study.NewBuilder()
    cur = theObject

    while cur:
        attr = builder.FindOrCreateAttribute(cur, "AttributeLocalID")
        value = attr.Value()
        if value == dict_object["Study"] or value == dict_object["CouplingStudy"]:
            return cur
        elif value == __MODULE_ID__ or value == 0:
            return None

        cur = cur.GetFather()

    return None


def FindStudyByPath(theStudyPath):
    """
    Returns a CFD study described by the unix path I{theStudyPath}.

    @type theStudyPath: C{String}
    @param theStudyPath: unix path of the CFD study.
    @return: the CFD study.
    @rtype: C{SObject} or C{None}
    """
    component = _getComponent()
    if component == None:
        return None

    study = _getStudy()
    builder = study.NewBuilder()

    iter  = study.NewChildIterator(component)
    while iter.More():
        attr = builder.FindOrCreateAttribute(iter.Value(), "AttributeLocalID")
        if attr.Value() == dict_object["Study"] or attr.Value() == dict_object["CouplingStudy"] :
            #compare study path
            aCurStudyPath = _GetPath(iter.Value())
            if aCurStudyPath == theStudyPath:
                return iter.Value()
        iter.Next()

    return None

def FindCaseByPath(theCasePath):
    """
    Returns a CFD study described by the unix path I{theCasePath}.

    @type theCasePath: C{String}
    @param theCasePath: unix path of the CFD study.
    @return: the CFD study.
    @rtype: C{SObject} or C{None}
    """
    component = _getComponent()
    if component == None:
        return None

    study = _getStudy()
    builder = study.NewBuilder()
    studyCfdObject = FindStudyByPath(os.path.dirname(theCasePath))
    iter  = study.NewChildIterator(studyCfdObject)
    while iter.More():
        attr = builder.FindOrCreateAttribute(iter.Value(), "AttributeLocalID")
        if attr.Value() == dict_object["Case"]:
            #compare case path
            aCurCasePath = _GetPath(iter.Value())
            if aCurCasePath == theCasePath:
                return iter.Value()
        iter.Next()

    return None

def GetCaseNameList(theStudy):
    """
    Returns the list of the existing cases from a CFD study in the Object Browser.
    Used into slotAddCase to verify the existing cases
    @type theStudy: C{SObject}
    @param theStudy: CFD study data in the Object Browser.
    @return: list of names of the loaded CFD studies.
    @rtype: C{List} or C{String}
    """
    CaseList = []

    study   = _getStudy()
    builder = study.NewBuilder()

    attr = builder.FindOrCreateAttribute(theStudy, "AttributeLocalID")
    if attr.Value() != dict_object["Study"] :
        if attr.Value() != dict_object["CouplingStudy"]:
            return CaseList

    iter  = study.NewChildIterator(theStudy)
    while iter.More():
        anObject = iter.Value()
        attr = builder.FindOrCreateAttribute(anObject, "AttributeLocalID")
        if attr.Value() == dict_object["Case"]:
            CaseList.append(anObject.GetName())
        iter.Next()

    return CaseList


def GetCaseList(theStudy):
    """
    Returns a list of data which are cases folder in the Object Browser.

    @type theStudy: C{SObject}
    @param theStudy: CFD study data in the Object Browser.
    @return: list of branch which are CFD cases.
    @rtype: C{list} of C{SObject}
    """
    CaseList = []

    study   = _getStudy()
    builder = study.NewBuilder()

    attr = builder.FindOrCreateAttribute(theStudy, "AttributeLocalID")
    if attr.Value() != dict_object["Study"] :
        if attr.Value() != dict_object["CouplingStudy"]:
            return CaseList

    iter  = study.NewChildIterator(theStudy)

    while iter.More():
        anObject = iter.Value()
        attr = builder.FindOrCreateAttribute(anObject, "AttributeLocalID")
        if attr.Value() == dict_object["Case"]:
            CaseList.append(anObject)
        iter.Next()

    return CaseList

def getXmlCaseNameList(theCase):
    """
    Returns a list of xml file names from case aCase
    """
    XmlCaseNameList = []
    if not checkType(theCase, dict_object["Case"]):
        return XmlCaseNameList
    aChildList = []

    aChildList = ScanChildren(theCase, "^DATA$")
    if len(aChildList) != 1:
        # no DATA folder
        print("There are no data folder in selected by user case")
        return

    aDataObj =  aChildList[0]
    aDataPath = _GetPath(aDataObj)
    study   = _getStudy()
    builder = study.NewBuilder()

    iter  = study.NewChildIterator(aDataObj)

    while iter.More():
        aName = iter.Value().GetName()
        if aName != "" :
            if "XML" in subprocess.check_output(["file",_GetPath(iter.Value())]):
                XmlCaseNameList.append(iter.Value().GetName())
        iter.Next()
    return XmlCaseNameList


def ScanChildren(theObject, theRegExp):
    """
    Returns a list of children data from a parent branch data.
    The list of the children is filtered whith a regular expression.

    @type theObject: C{SObject}
    @param theObject: parent data.
    @type theRegExp: C{String}
    @param theRegExp: regular expression to filter children data.
    @return: list of branch of children data.
    @rtype: C{list} of C{SObject}
    """
    ChildList = []

    study   = _getStudy()
    builder = study.NewBuilder()

    iter  = study.NewChildIterator(theObject)

    while iter.More():
        aName = iter.Value().GetName()
        if not aName == "" and re.match(theRegExp, aName):
            ChildList.append(iter.Value())
        iter.Next()

    return ChildList


def ScanChildNames(theObject, theRegExp):
    """
    Returns a list of children data names from a parent branch data.
    The list of the children is filtered whith a regular expression.

    @type theObject: C{SObject}
    @param theObject: parent data.
    @type theRegExp: C{String}
    @param theRegExp: regular expression to filter children data.
    @return: list name of branch of children data.
    @rtype: C{list} of C{String}
    """
    NameList = []

    study   = _getStudy()
    builder = study.NewBuilder()

    iter  = study.NewChildIterator(theObject)

    while iter.More():
        aName = iter.Value().GetName()
        if not aName == "" and re.match(theRegExp, aName):
            NameList.append(aName)
        iter.Next()

    #log.debug("ScanChildNames: %s -> %s" % (theObject.GetName(), NameList))
    return NameList


def getType(theObject):
    if theObject == None:
        return None
    study   = _getStudy()
    builder = study.NewBuilder()
    attr = builder.FindOrCreateAttribute(theObject, "AttributeLocalID")
    return attr.Value()

def hasTheSameType(ListObject):

    if ListObject == []:
        return False
    typListBool = True
    typListBoolRESUSub = None
    typList = []
    typ     = getType(ListObject[0])
    typList.append(typ)
    if len(ListObject)> 1:
        for SObject in ListObject[1:]:
            typListBool = typListBool and getType(SObject) == typ
            typList.append(getType(SObject))
        if not typListBool:
            typListBoolRESUSub = True
            for ty in typList:
                typListBoolRESUSub = typListBoolRESUSub and (ty == dict_object["RESUSubFolder"] or ty == dict_object["RESUSubErrFolder"] or ty == dict_object["RESU_COUPLINGSubFolder"])
    if typListBoolRESUSub != None:
        if typListBoolRESUSub == True:
            return typListBoolRESUSub
    else :
        return typListBool

def isACFDSTUDYListObject(ListObject):
    if ListObject == []:
        return False
    typListBool = True
    for sobj in ListObject :
        typListBool = typListBool and sobj.GetFatherComponent().GetName() == "CFDSTUDY"
    return typListBool

def isASmeshListObject(ListObject):
    if ListObject == []:
        return False
    typListBool = True
    for sobj in ListObject :
        if sobj.GetFatherComponent().GetName() == "Mesh":
            if getMeshFromMesh(sobj) == None:
                meshGroupObject,group = getMeshFromGroup(sobj)
                typListBool = typListBool and meshGroupObject != None
            else:
                typListBool = True
        else:
           return False
    return typListBool

def checkType(theObject, theType):
    """
    Checks if I{theObject} has the type ("AttributeLocalID") I{theType}.

    @type theObject: C{SObject}
    @param theObject: object from the Object Browser.
    @type theType: C{String}
    @param theType: type of the object in the Object Browser.
    @rtype: C{True} or C{False}
    @return: C{True} if C{theObject} has the type I{theType}.
    """
    if theObject == None or theType == None:
        return False
    if theObject != None and theType!= None :
        return getType(theObject) == theType

def checkPreMEDType(theObject):
    """
    Checks if I{theObject} is a mesh file, that can be converted to med format.

    @type theObject: C{SObject}
    @param theObject: object from the Object Browser.
    @rtype: C{True} or C{False}
    @return: C{True} if C{theObject} is a mesh file, that can be converted to med.
    """
    return checkType(theObject, dict_object["DESFile"]) or \
           checkType(theObject, dict_object["CGNSFile"]) or \
           checkType(theObject, dict_object["GeomFile"]) or \
           checkType(theObject, dict_object["CaseFile"]) or \
           checkType(theObject, dict_object["NeuFile"]) or \
           checkType(theObject, dict_object["MSHFile"]) or \
           checkType(theObject, dict_object["HexFile"]) or \
           checkType(theObject, dict_object["UnvFile"])


def checkCaseLaunchGUI(theCase):
    """
    Checks if I{theCase} has the script to start GUI in the DATA folder.

    @type theCase: C{SObject}
    @param theCase: object from the Object Browser.
    @rtype: C{True} or C{False}
    @return: C{True} if C{theCase} has the script to start GUI in the DATA folder.
    """
    if not checkType(theCase, dict_object["Case"]):
        return False

    aChildList = ScanChildren(theCase, "^DATA$")
    if not len(aChildList) == 1:
        # no DATA folder
        print("There are no data folder in selected by user case")
        return False

    aDataObj =  aChildList[0]
    aDataPath = _GetPath(aDataObj)

    import sys
    if CFD_Code() == "Code_Saturne":
        if sys.platform.startswith("win"):
            aChildList = ScanChildren(aDataObj, "^SaturneGUI.bat$")
        else:
            aChildList = ScanChildren(aDataObj, "^SaturneGUI$")
    elif CFD_Code() == "NEPTUNE_CFD":
        if sys.platform.startswith("win"):
            aChildList = ScanChildren(aDataObj, "^NeptuneGUI.bat$")
        else:
            aChildList = ScanChildren(aDataObj, "^NeptuneGUI$")
    if not len(aChildList) == 1:
        if Trace(): print("There are no SaturneGUI or NeptuneGUI in selected by user case")
        return False

    return True


def checkCode(theCase):
    """
    Checks if I{theCase} is code_saturne or neptune_cfd.

    @type theCase: C{SObject}
    @param theCase: object from the Object Browser.
    @rtype: C{CFD_Saturne} or C{CFD_Neptune}
    @return: C{True} if C{theCase} has the script to start GUI in the DATA folder.
    """
    if not checkType(theCase, dict_object["Case"]):
        return CFD_Code()

    aChildList = ScanChildren(theCase, "^DATA$")
    if not len(aChildList) == 1:
        # no DATA folder
        print("There is no data folder in selected case")
        return CFD_Code()

    aDataObj =  aChildList[0]
    aDataPath = _GetPath(aDataObj)

    import sys
    if sys.platform.startswith("win"):
        aChildList = ScanChildren(aDataObj, "^SaturneGUI.bat$")
        if len(aChildList) == 1:
            return CFD_Saturne
        aChildList = ScanChildren(aDataObj, "^NeptuneGUI.bat$")
        if len(aChildList) == 1:
            return CFD_Neptune
    else:
        aChildList = ScanChildren(aDataObj, "^SaturneGUI$")
        if len(aChildList) == 1:
            return CFD_Saturne
        aChildList = ScanChildren(aDataObj, "^NeptuneGUI$")
        if len(aChildList) == 1:
            return CFD_Neptune
    ### try in SCRIPTS
    aChildList = ScanChildren(theCase, "^SCRIPTS$")
    if not len(aChildList) == 1:
        # no SCRIPTS folder
        print("There is no scripts folder in selected case")
        return CFD_Code()

    aDataObj =  aChildList[0]
    aDataPath = _GetPath(aDataObj)
    aChildList = ScanChildren(aDataObj, "^runcase$")
    if len(aChildList) == 1:
        path = os.path.join(aDataPath, "runcase")
        try:
            f = open(path, mode = 'r')
            lines = f.readlines()
            f.close()

            for i in range(len(lines) - 1, -1, -1):
                line = lines[i]
                # Skip comment and empty lines
                if len(line) == 0:
                    continue
                if line[0] == '#' or line[0:4] == 'rem ':
                    continue
                j = line.find('#')
                if j > -1:
                    line = line[0:j]
                args = separate_args(line)
                if args.count('run') == 1:
                    if args.index('run') == 1: # "<package_name> run"
                        for name in ('code_saturne', 'neptune_cfd'):
                            if not sys.platform.startswith('win'):
                                test_name = '\\' + name
                            else:
                                test_name = name
                            if args[0].find(test_name) == 0:
                                if name == 'code_saturne':
                                    return CFD_Saturne
                                else:
                                    return CFD_Neptune
        except IOError:
            print(("Error: can not open or read %s\n" % path))
            return CFD_Code()


def isLinkPathObject(theObject):
    """
    Checks if I{theObject} represents a unix symbolic link.

    @type theObject: C{SObject}
    @param theObject: object from the Object Browser.
    @rtype: C{True} or C{False}
    @return: C{True} if C{SObject} represents a unix symbolic link.
    """
    if theObject == None:
        return False

    study   = _getStudy()
    builder = study.NewBuilder()
    attr = builder.FindOrCreateAttribute(theObject, "AttributeComment")
    return re.match("^->", attr.Value())


def setCaseInProcess(theCasePath, isInProcess):
    """
    Udpates the case icon with I{Case} or I{CaseInProcess} in the Object Browser.

    @type theCasePath: C{String}
    @param theCasePath: absolute path of the case.
    @type isInProcess: C{True} or C{False}
    @param isInProcess: if C{True}, shows the I{CaseInProcess} icon.
    """
    aStudyPath, aCaseName = os.path.split(theCasePath)
    aStudyObj = FindStudyByPath(aStudyPath)
    if not aStudyPath:
        if Trace():
            print("Study by case path not found")
        return

    #get case object
    lst = ScanChildren(aStudyObj, aCaseName)
    if len(lst) != 1:
        if Trace():
            print("Invalid number of cases under study")
        return

    aCaseObj = lst[0]

    study   = _getStudy()
    builder = study.NewBuilder()

    attr = builder.FindOrCreateAttribute(aCaseObj, "AttributePixMap")
    if isInProcess:
        attr.SetPixMap(str(ObjectTR.tr(icon_collection[dict_object["CaseInProcess"]])))
    else:
        attr.SetPixMap(str(ObjectTR.tr(icon_collection[dict_object["Case"]])))


def getSObject(theParent,Name) :
    Sobjlist = ScanChildren(theParent,  ".*")
    SObj = None
    for i in Sobjlist :
        if i.GetName() == Name :
            SObj = i
    return SObj


#def publishInStudySalome(SO_father, objName, idElem):
    #"""
    #Publish objName into Object Browser under SO_father with the AttributeLocalID idElem
    #listPublishedId is used into PublishedIntoObjectBrowser method and caracterize entries
    #PublishedIntoObjectBrowser method adds entries into Salome Object Browser.
    #These entries do not provide from an Unix cfd study directory, and are idendified into the object browser
    #by a localId Attribute from the python list listPublishedId
    #"""
    #study = _getStudy()
    #builder = study.NewBuilder()
    #studyObject = builder.NewObject(SO_father)
    #attr = builder.FindOrCreateAttribute(studyObject, "AttributeName")
    #attr.SetValue(objName)
    #attr = builder.FindOrCreateAttribute(studyObject, "AttributeLocalID")
    #attr.SetValue(idElem)
    #_SetIcon(studyObject, builder)
    #log.debug("publishInStudySalome: %s" % ScanChildNames(SO_father,  ".*"))
    #return studyObject


def getOrLoadObject(item):
    """
    Get the CORBA object associated with the SObject `item`, eventually by
    first loading it with the corresponding engine.
    """
    object = item.GetObject()
    if object is None: # the engine has not been loaded yet
        sComponent = item.GetFatherComponent()
        #self.loadComponentEngine(sComponent)
        study   = _getStudy()
        builder = study.NewBuilder()
        engine = _getEngine()
        if engine is None:
            print("Cannot load component ", __MODULE_NAME__)
        object = item.GetObject()
    return object

def getMeshFromMesh(meshSobjItem) :
    """
    return: The SALOMEDS._objref_SObject instance of the mesh, if the meshSobjItem is a sobj of a mesh, None if not
    """
    meshItem = None
    obj = getOrLoadObject(meshSobjItem)
    if obj != None:
        mesh = obj._narrow(SMESH.SMESH_Mesh)
        if mesh != None:
            meshItem = salome.ObjectToSObject(mesh)
    return meshItem

def SetAutoColor(meshSobjItem) :
    obj = getOrLoadObject(meshSobjItem)
    if obj is not None:
        mesh = obj._narrow(SMESH.SMESH_Mesh)
        if mesh is not None:
            mesh.SetAutoColor(1)

def getMeshFromGroup(meshGroupItem):
    """
    Get the mesh item owning the mesh group `meshGroupItem`.

    :type   meshGroupItem: SObject
    :param  meshGroupItem: Mesh group belonging to the searched mesh.

    :return: The SALOMEDS._objref_SObject instance corresponding to the mesh group or None
             and SMESH._objref_SMESH_Group instance  or None if it was not
             found.
    """
    group = None
    meshItem = None
    obj = getOrLoadObject(meshGroupItem)
    #obj = self.editor.getOrLoadObject(meshGroupItem) # version 515

    if obj is not None:
        group = obj._narrow(SMESH.SMESH_GroupBase)
        if group != None: # The type of the object is ok
            meshObj = group.GetMesh()
            meshItem = salome.ObjectToSObject(meshObj)
    return meshItem, group
