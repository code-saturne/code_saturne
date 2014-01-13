# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
from PyQt4.QtCore import QObject
from PyQt4.QtGui  import QMessageBox

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, re
import logging

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

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
_Trace = False  #True

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
        mess = ObjectTR.tr("CFDSTUDY_INVALID_SOLVER_NAME").arg(code).arg(CFD_Saturne).arg(CFD_Neptune)
        iok= False

    if code == CFD_Saturne:
        try:
            from cs_package import package
            iok = True
        except ImportError,e:
            mess = ObjectTR.tr("INFO_DLG_INVALID_ENV").arg(code) + e.__str__()
            if "cs_package" in e.__str__():
                mess = mess + " ; Check for cs_package file in Code_Saturne python package"
            elif "code_saturne" in e.__str__():
                mess = mess + " ; Check PYTHONPATH then your installation "
            iok = False
    elif code == CFD_Neptune:
        try:
            from nc_package import package
            iok = True
        except ImportError,e:
            mess = ObjectTR.tr("INFO_DLG_INVALID_ENV").arg(code) + e.__str__()
            if "nc_package" in e.__str__():
                mess = mess + " ; Check for nc_package file in NEPTUNE_CFD python package"
            elif "neptune_cfd" in e.__str__():
                mess = mess + " ; Check PYTHONPATH then your installation "
            iok = False
    else:
        raise ApplicationError, "Invalid name of solver!"

    if iok:
        pkg = package()
        prefix = pkg.get_dir('prefix')
        log.debug("CheckCFD_CodeEnv -> prefix = %s" % (prefix))

        bindir = pkg.get_dir('bindir')
        log.debug("CheckCFD_CodeEnv -> prefix = %s" % (bindir))

        if not os.path.exists(prefix):
            mess1 = ObjectTR.tr("ENV_DLG_INVALID_DIRECTORY")
            mess = mess + mess1.arg(prefix)
            iok = False
        else:
            if not os.path.exists(bindir):
                mess2 =  ObjectTR.tr("ENV_DLG_INVALID_DIRECTORY")
                mess = mess + mess2.arg(bindir)
                iok = False

    log.debug("CheckCFD_CodeEnv -> %s = %s" % (code, iok))
    log.debug("CheckCFD_CodeEnv -> %s: %s" % (code, mess))
    return iok, mess


def BinCode():
    b = ""
    c = ""
    mess = ""
    if CFD_Code() == CFD_Saturne:
        from cs_package import package
        pkg = package()
        bindir = pkg.get_dir('bindir')
        b = os.path.join(bindir, "code_saturne")
    elif CFD_Code() == CFD_Neptune:
        from nc_package import package
        pkg = package()
        bindir = pkg.get_dir('bindir')
        b = os.path.join(bindir, "neptune_cfd")

    c = pkg.get_preprocessor()
    log.debug("BinCode -> \n    %s\n    %s" % (b, c))
    return b, c, mess

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
