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

class Chrono :
    """
    Describe a chrono
    """

    def __init__(self, casePath, variables):
        """
        constructor
        """
        self.variables = variables
        self.casePath = casePath
        #
        # Verifier qu'il n'y a qu'un chr dans RESU
        self.chrPath = casePath+"/RESU"

        filesNames = os.listdir(self.chrPath)

        counter = 0
        for fileName in filesNames :
            if fileName.find('CHR.ENSIGHT') >= 0 :
                self.chrName = fileName
                counter += 1
        if counter > 1 :
            print ("Error : Several chrono files in :" + self.chrPath)
            print ("        Delete old files and run again")  
            print sys.exit(1)
        if counter < 1 :
            print ("Error : No chrono file in :" + self.chrPath)
            print sys.exit(1)

        try:
            chrFile = file(self.chrPath+'/'+self.chrName+'/CHR.case', mode='r')
        except IOError:
            print "Error : opening "+self.chrPath+'/'+self.chrName+'/CHR.case'
            print sys.exit(1)            
        
        variablesChrNames={}
        for variable in self.variables :
            #
            # in CHR.case, variable name len is < 19
            if len(variable) > 19 :
                variableChrName = variable[:19]
            else : 
                variableChrName = variable
            #
            # in CHR.case, replace special characters
            characters=["(",")","[","]","+","-","@"," ","!",
                        "#","*","^","$","/"]
            
            for character in characters :
                variableChrName.replace(character,"_")
                
            variablesChrNames[variable] = variableChrName

        line = ""
        while line.find("VARIABLE") < 0 :
            line = chrFile.readline()
            if (line == ""):
                print "Error : reading "+self.chrPath+'/'+self.chrName+'/CHR.case'
                sys.exit(1)
            
        self.variablesChrFiles={}
        self.typ = {}

        while line.find("TIME") < 0:
            line = chrFile.readline()
            if (line == ""):
                print "Error : reading "+self.chrPath+'/'+self.chrName+'/CHR.case'
                sys.exit(1)

            for variable in self.variables :
                numCar = line.find(variablesChrNames[variable])
                if numCar > -1 :
                    self.typ[variable] = line.split(":")[0]
                    nameTab = line[numCar:].split()
                    self.variablesChrFiles[variable]=nameTab[1][:nameTab[1].find("*")]
                      
        while line.find("number of steps:") < 0:
            line = chrFile.readline()
            if (line == ""):
                print "Error : reading "+self.chrPath+'/'+self.chrName+'/CHR.case'
                sys.exit(1)
        stepNb = int(line.split(":")[1])
        filenumber = '%4.4i' % stepNb

        for variable in self.variables :
            self.variablesChrFiles[variable]+=str(filenumber)

        while line.find("time values:") < 0:
            line = chrFile.readline()
            if (line == ""):
                print "Error : reading "+self.chrPath+'/'+self.chrName+'/CHR.case'
                sys.exit(1)

#        self.time = float(line.split(":")[1])
        
        while 1 :
            line = chrFile.readline()
            if (line == ""):
                break
            try:
                self.time = float(line.split()[0])
            except:
                pass
            

    def getTime(self):
        return self.time


    def getTyp(self,variable):
        return self.typ[variable]
            

    def getValues(self, variable):
        chrFicName = self.variablesChrFiles[variable]
        path = self.chrPath+'/'+self.chrName+'/'+self.variablesChrFiles[variable]
        try:
            chrVar = file(path, mode='r')
        except IOError:
            print "Error : opening "+path
            print sys.exit(1)

        tabl = []
        counter = 0

        while 1:
            text = chrVar.readline()
            counter += 1		
            if (text == ""):
                break
                            
            if counter > 4:
                tmp = text.split()
                try:
                    tabl.append(float(tmp[0]))
                except:
                    pass
                            
        chrVar.close()                        
        return tabl


    def getStructure(self,variable):
        chrFicName = self.variablesChrFiles[variable]
        path = self.chrPath+'/'+self.chrName+'/'+self.variablesChrFiles[variable]
        print path
        try:
            chrVar = file(path, mode='r')
        except IOError:
            print "Error : opening "+path
            print sys.exit(1)
    
        structure = []

        line = chrVar.readline()
        while line.find("part") < 0:
            line = chrVar.readline()
            if (line == ""):
                break

        while 1:
            line = chrVar.readline()		
            if (line == ""):
                break

            valpart=[]
            counter = 0
            typ = ""
            while line.find("part") < 0:
                line = chrVar.readline()
                if (line == ""):
                    break

                tmp=line.split()
                try:
                    float(tmp[0])
                    counter += 1
                except:
                    if counter > 0 and typ != "":
                        valpart.append([typ,counter])
                    typ = tmp[0]
                    counter = 0
                    
            if typ != "part" and typ !="" :
                valpart.append([typ,counter])
            structure.append(valpart)
                
            if (line == ""):
                break

        return structure


    def addVariable(self, values, typ, name, varModel):
        #
        # on ajoute les variables supplementaires dans CHR.case
        #
        print self.chrPath+'/'+self.chrName+'/CHR.case'

        try:
            chrFile = file(self.chrPath+'/'+self.chrName+'/CHR.case', mode='r')
        except IOError:
            print "Error : opening "+self.chrPath+'/'+self.chrName+'/CHR.case'
            print sys.exit(1)

        chrFileTmp = file(self.chrPath+'/'+self.chrName+'/CHR.case.tmp', mode='w')

        line = ""
        while line.find("VARIABLE") < 0 :
            line = chrFile.readline()
            if (line == ""):
                print "Error : reading "+self.chrPath+'/'+self.chrName+'/CHR.case'
                sys.exit(1)
            chrFileTmp.write(line)

        if len(name) > 19 :
            line = typ+":"+" 1 "+name[:19]+" "+"chr."+name[:19]+"\n"
        else:
            line = typ+":"+" 1 "+name+" "+" "+"chr."+name+"\n"

        chrFileTmp.write(line)

        while 1 :
            line = chrFile.readline()
            if (line == ""):
                break
            chrFileTmp.write(line)

        chrFile.close()
        chrFileTmp.close()
        os.rename(self.chrPath+'/'+self.chrName+'/CHR.case.tmp',self.chrPath+'/'+self.chrName+'/CHR.case')
            
        #
        # creation du chr.delta
        #
        structure = self.getStructure(varModel)
        
        if len(name) > 19 :
            chrVar = file(self.chrPath+'/'+self.chrName+'/chr.'+name[:19], mode='w')
        else:
            chrVar = file(self.chrPath+'/'+self.chrName+'/chr.'+name, mode='w')

        chrVar.write('Variable : '+name+'\n')
        
        numPart = 1
        debut = 0
        for part in structure:
            chrVar.write("part\n")
            chrVar.write("%10i\n" % numPart)
            for elt in part:
                chrVar.write(elt[0]+"\n")
                for i in range(debut, debut+elt[1]) :
                    chrVar.write("%12.5e\n" % values[i])
                debut = elt[1]
            numPart +=1
