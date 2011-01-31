# -*- coding: utf-8 -*-
#============================================================================
#
#     This file is part of CFDSTUDY the plug-in for Salome
#     of Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     CFDSTUDY is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     CFDSTUDY is distributed in the hope that it will be
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
#============================================================================

"""
Common
======

"""

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
log = logging.getLogger("CFDSTUDYGUI")
#log.setLevel(logging.DEBUG)
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
    if code == CFD_Saturne:
        try:
            from code_saturne import cs_package
            pkg = cs_package.package()
            prefix = pkg.prefix
            bindir = pkg.bindir
            iok = True
        except:
            iok = False
    elif code == CFD_Neptune:
        try:
            from ncs import nc_config
            prefix = nc_config.dirs.prefix
            bindir = nc_config.dirs.bindir
            iok = True
        except:
            iok = False
    else:
        raise ApplicationError, "Invalid name of solver!"

    if iok:
        if not os.path.exists(prefix):
            iok = False
        if not os.path.exists(bindir):
            iok = False

    if iok:
        if not os.path.isfile(os.path.join(bindir, "code_saturne")) and \
           not os.path.isfile(os.path.join(bindir, "cs")) and \
           not os.path.isfile(os.path.join(bindir, "nc")):
            iok = False

    log.debug("CheckCFD_CodeEnv -> %s = %s" % (code, iok))
    return iok


def BinCode():
    if CFD_Code() == CFD_Saturne:
        from code_saturne import cs_package
        pkg = cs_package.package()
        bindir = pkg.bindir
        if os.path.isfile(os.path.join(bindir, "code_saturne")):
            b = os.path.join(bindir, "code_saturne")
        elif os.path.isfile(os.path.join(bindir, "cs")):
            b = os.path.join(bindir, "cs")
    elif CFD_Code() == CFD_Neptune:
        from ncs import nc_config
        bindir = nc_config.dirs.bindir
        b = os.path.join(bindir, "nc")

    c = os.path.join(bindir, "cs_preprocess")
    log.debug("BinCode -> \n    %s\n    %s" % (b, c))
    return b, c


#def runCommand(cmd, start_directory, prefix, * args):
#    """
#    Run command cmd and asynchronize put it to LogWindow
#    Each string logged with prefix
#    """
#    import subprocess
#    import os
#
#    if start_directory != None and start_directory != "":
#        os.chdir(start_directory)
#
#    try:
#        pipe = subprocess.Popen(cmd, bufsize = 0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#
#        while True:
#            text = pipe.stdout.readline()
#            if not text:
#                break
#
#            sgPyQt.message( prefix + text, False )
#
#    except OSError, e:
#        sgPyQt.message( prefix + "Exception had occured during script execution " + e.__str__(), True )

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
