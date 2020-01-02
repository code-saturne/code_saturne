# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
Common
======

"""
from code_saturne.Base.QtCore    import *

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, re, subprocess
import logging

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------
from CFDSTUDYGUI_Message import cfdstudyMess
# Get SALOME PyQt interface
import SalomePyQt
sgPyQt = SalomePyQt.SalomePyQt()

# Get SALOME Swig interface
import libSALOME_Swig
sg = libSALOME_Swig.SALOMEGUI_Swig()

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Global variables
#-------------------------------------------------------------------------------

CFD_Saturne = "Code_Saturne"
CFD_Neptune = "NEPTUNE_CFD"

# ObjectTR is a convenient object for traduction purpose
ObjectTR = QObject()

# Main variable for solver
_CFD_Code = None #By default

# True or false for log tracing
_Trace = False  #

# If True all stdout redirected to MassageWindow
_LogModeOn = True

#---Enumerations---
#Event type for indicate of case in process
CaseInProcessStart  = -1000
CaseInProcessEnd    = -1001
UpdateScriptFolder  = -1002

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_Commons")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Functions definitions
#-------------------------------------------------------------------------------

def CFD_Code():
    global _CFD_Code
    return _CFD_Code


def _SetCFDCode(var):
    log.debug("_SetCFDCode : var = %s" % var)
    global _CFD_Code
    _CFD_Code = var


def Trace():
    global _Trace
    return _Trace


def LogModeOn():
    global _LogModeOn
    _LogModeOn = True


def LogModeOff():
    global _LogModeOn
    _LogModeOn = False


# check for avalable type of solver

def CheckCFD_CodeEnv(code):
    """
    This method try to found the config file of the CFD I{code}.

    @param code: name of the searching code (CFD_Saturne or CFD_Neptune).
    @type theType: C{String}
    @rtype: C{True} or C{False}
    @return: C{True} if the searching code is found.
    """
    mess = ""
    prefix = ""
    bindir = ""

    if code not in [CFD_Saturne, CFD_Neptune]:
        mess = cfdstudyMess.trMessage(ObjectTR.tr("CFDSTUDY_INVALID_SOLVER_NAME"),[code,CFD_Saturne,CFD_Neptune])
        iok= False
        return iok,mess

    if code == CFD_Saturne:
        try:
            from code_saturne.cs_package import package
            iok = True
        except ImportError as e:
            mess = cfdstudyMess.trMessage(ObjectTR.tr("INFO_DLG_INVALID_ENV"),[code]) + e.__str__()
            if "cs_package" in e.__str__():
                mess = mess + cfdstudyMess.trMessage(ObjectTR.tr("CHECK_CODE_PACKAGE"),["cs_package",code])
            elif "code_saturne" in e.__str__():
                mess = mess + cfdstudyMess.trMessage(ObjectTR.tr("CHECK_PYTHON_PATH"),[])
            iok = False
    elif code == CFD_Neptune:
        try:
            from neptune_cfd.nc_package import package
            iok = True
        except ImportError as e:
            mess = cfdstudyMess.trMessage(ObjectTR.tr("INFO_DLG_INVALID_ENV"),[code]) + e.__str__()
            if "nc_package" in e.__str__():
                mess = mess + cfdstudyMess.trMessage(ObjectTR.tr("CHECK_CODE_PACKAGE"),["nc_package",code])
            elif "neptune_cfd" in e.__str__():
                mess = mess + cfdstudyMess.trMessage(ObjectTR.tr("CHECK_PYTHON_PATH"),[])
            iok = False
    else:
        raise ValueError("Invalid name of solver!")

    if iok:
        pkg = package()
        prefix = pkg.get_dir('prefix')
        log.debug("CheckCFD_CodeEnv -> prefix = %s" % (prefix))

        bindir = pkg.get_dir('bindir')
        log.debug("CheckCFD_CodeEnv -> prefix = %s" % (bindir))

        if not os.path.exists(prefix):
            mess = cfdstudyMess.trMessage(ObjectTR.tr("ENV_DLG_INVALID_DIRECTORY"),[prefix])
            iok = False
        else:
            if not os.path.exists(bindir):
                mess = cfdstudyMess.trMessage(ObjectTR.tr("ENV_DLG_INVALID_DIRECTORY"),[bindir])
                iok = False

    log.debug("CheckCFD_CodeEnv -> %s = %s" % (code, iok))
    log.debug("CheckCFD_CodeEnv -> %s: %s" % (code, mess))
    return iok, mess


def BinCode():
    b = ""
    c = ""
    mess = ""
    # default package is code_saturne (for convert...)
    from code_saturne.cs_package import package
    pkg = package()

    if CFD_Code() == CFD_Saturne:
        bindir = pkg.get_dir('bindir')
        b = os.path.join(bindir, "code_saturne")
    elif CFD_Code() == CFD_Neptune:
        from neptune_cfd.nc_package import package
        pkg = package()
        bindir = pkg.get_dir('bindir')
        b = os.path.join(bindir, "neptune_cfd")
    c = pkg.get_preprocessor()
    log.debug("BinCode -> \n    %s\n    %s" % (b, c))
    return b, c, mess

def isaCFDCase(theCasePath):
    log.debug("isaCFDCase")
    iok = True
    dirList = []
    if os.path.isdir(theCasePath):
        try:#python3
            dirList = os.walk(theCasePath).__next__()[1]
        except: #python27
            dirList = os.walk(theCasePath).next()[1]
        if (dirList.count("DATA") and \
           dirList.count("SRC")  and \
           dirList.count("SCRIPTS")):
            if not (dirList.count("RESU")):
                subprocess.call(["mkdir","-p",os.path.join(theCasePath,"RESU")])
        else:
            iok = False
    return iok

def isaCFDStudy(theStudyPath):
    log.debug("isaCFDStudy")
    dirList = []
    if os.path.isdir(theStudyPath):
        try:#python3
            dirList = os.walk(theStudyPath).__next__()[1]
        except: #python27
            dirList = os.walk(theStudyPath).next()[1]
        for i in dirList:
            if i not in ["MESH"] :
                if isaCFDCase(os.path.join(theStudyPath,i)) :
                    return True
    return False

def isSyrthesCase(theCasePath):
    log.debug("isSyrthesCase")
#a minima
    iok = True
    if os.path.isdir(theCasePath):
        dirList = os.listdir(theCasePath)
        if not dirList.count("Makefile") and not dirList.count("syrthes.py"):
            iok = False
    return iok

def isaSaturneSyrthesCouplingStudy(theStudyPath):
    log.debug("isaSaturneSyrthesCouplingStudy")
    iok = False
    hasCFDCase     = False
    hasSyrthesCase = False
    if not os.path.isdir(theStudyPath):
        mess = cfdstudyMess.trMessage(ObjectTR.tr("MUST_BE_A_DIRECTORY"),[theStudyPath])
        cfdstudyMess.criticalMessage(mess)
        return False
    dirList = os.listdir(theStudyPath)
    if not (dirList.count("RESU_COUPLING") and dirList.count("coupling_parameters.py") and dirList.count("runcase")):
        return False
    for i in dirList:
        ipath = os.path.join(theStudyPath,i)
        if os.path.isdir(ipath):
            if i not in ["MESH", "RESU_COUPLING"]:
                if isaCFDCase(ipath):
                    hasCFDCase = True
                if isSyrthesCase(ipath):
                    hasSyrthesCase = True
    if hasCFDCase and hasSyrthesCase:
        iok = True
    return iok

#-------------------------------------------------------------------------------
# Classes definitions
#-------------------------------------------------------------------------------

class LoggingAgent:
    def __init__(self, stream ):
        self.stream = stream


    def write( self, Text ):
        global _LogModeOn
        global sgPyQt

        #self.stream.write( Text )

        if len(Text) == 0:
            return

        lst = re.split( "\n", Text )
        for s in lst:
            if not len(s) == 0:
                sgPyQt.message( re.sub('<','&lt;',re.sub( '>', '&gt;', s)), False )


    def close(self):
        return self.stream


class LoggingMgr:
    def __init__(self ):
        pass


    def start( self, sys_obj):
        self.AgentOut = LoggingAgent( sys_obj.stdout )
        sys_obj.stdout = self.AgentOut

        self.AgentErr = LoggingAgent( sys_obj.stderr )
        sys_obj.stderr = self.AgentErr


    def finish( self, sys_obj):

        if self.AgentOut != None:
            sys_obj.stdout = self.AgentOut.close()

        if self.AgentErr != None:
            sys_obj.stderr = self.AgentErr.close()
