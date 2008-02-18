#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
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

import os, sys
import re
import string

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

class Listing :
    """
    Describe a listing
    """

    def __init__(self, casePath):
        """
        constructor
        """
        #
        # Verifier qu'il n'y a qu'un listing dans RESU
        listingPath = casePath+"/RESU"

        filesNames = os.listdir(listingPath)
        for fileName in filesNames :
            if fileName.find('listing_n') >= 0 :
                os.remove(listingPath+"/"+fileName)    

        counter = 0
        filesNames = os.listdir(listingPath)
        for fileName in filesNames :
            if fileName.find('listing') >= 0 :
                listingName = fileName 
                counter += 1
        if counter > 1 :
            print ("Error : Several listing files in :" + listingPath)
            print ("        Delete old files and run again")  
            print sys.exit(1)
        if counter < 1 :
            print ("Error : No listing file in :" + listingPath)
            print sys.exit(1)

        try:
            self.listingFile = file(listingPath+'/'+listingName, mode='r')
        except IOError:
            print "Error : opening "+listingPath+'/'+listingName
            print sys.exit(1)            


    def getMinMaxVariables(self, variables):

        clipMin = {}
        clipMax = {}
        varMin = {}
        varMax = {}
        variablesListingNames = {}

        for variable in variables :
            #
            # in listing, variable name len is < 12
            if len(variable) > 12 :
                variableListingName = variable[:12]
            else : 
                variableListingName = variable

            variablesListingNames[variable] = variableListingName
            clipMin[variable] = None
            clipMax[variable] = None
            varMin[variable] = None
            varMax[variable] = None
            
        while 1:
            line = self.listingFile.readline()
            if (line == ""):
                break
            
            for variable in variables :
                
                keyword = "v  " + variablesListingNames[variable]
                kw = re.compile(keyword)

                lineVar = kw.match(line)
                if lineVar :
                    tmp = line.split()
                    
                    if len(tmp) == 6 :
                        varMin[variable] = float(tmp[2])
                        varMax[variable] = float(tmp[3])

                        if tmp[4] != '--' :
                            clipMin[variable] = int(tmp[4])

                        if tmp[5] != '--' :
                            clipMax[variable] = int(tmp[5])

                    else :                        
                        print " Warning : Variable found but bad listing structure"

        return  varMin, varMax, clipMin, clipMax
    

