#============================================================================
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2008 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne Kernel is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne Kernel is distributed in the hope that it will be
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
#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, os, shutil

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

import Autovalidation.Common as Common
import Autovalidation.Case as Case
import Autovalidation.Parser as Parser

class Study:
    
    def __init__(self, parser, studyLabel):
        
        self.cases = []
        #
        # if study doesn't exist, create it with "cree_sat"
        #
        localDir = Common.localDirectory
        studyPath = localDir+"/"+studyLabel.upper()

        reportListing = studyLabel + "_listing.report"
        fa = open(reportListing,'w')
        fa.write('-------------------------\n')
        fa.write('Listing files comparison \n')
        fa.write('-------------------------\n')
        fa.close()

        reportChrono = studyLabel + "_chrono.report"
        fa = open(reportChrono,'w')
        fa.write('-------------------------\n')
        fa.write('Chrono files comparison  \n')
        fa.write('-------------------------\n')
        fa.close()

        if not os.path.isdir(studyLabel.upper()) :
            proc = os.popen("cs_create -nogui -study "+studyLabel)
            proc.close()
            shutil.rmtree(studyPath+"/CASE1")
        #
        # lien ou copie des maillages et des scripts du cas de reference
        #
        refMeshPath = Common.referencePath+"/"+studyLabel.upper()+"/MESH"
        try:
            meshesList = os.listdir(refMeshPath)
            for mesh in meshesList :
                os.symlink(refMeshPath+"/"+mesh,studyPath+"/MESH/"+mesh)
        except:
            pass

        refPostPath = Common.referencePath+"/"+studyLabel.upper()+"/POST"
        try:
            scriptsList = os.listdir(refPostPath)
            for script in scriptsList :
                shutil.copyfile(refPostPath+"/"+script,studyPath+"/POST/"+script)
        except:
            pass
        #
        # On recupere les variables a comparer
        variables = parser.getVariablesDefinition(studyLabel)
        #
        # create cases list
        casesLabels = parser.getCasesLabels(studyLabel)
        if casesLabels != [] and casesLabels != None:
            for caseLabel  in casesLabels:
                case = Case.Case(parser, studyLabel, caseLabel, variables, reportListing, reportChrono )
                self.cases.append(case)
        else :
            print "Error : no case in "+studyLabel+" study"
            sys.exit(1)
        
    def getCases(self):
        return self.cases
