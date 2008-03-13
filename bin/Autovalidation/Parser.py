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
import sys
from xml.dom import minidom
#from neptune import case

class Parser:
    """ Parser -- class to parse XML file."""

    def __init__ (self, XMLFileName):
        "Constructor --"
        try:
            self.doc =  minidom.parse (XMLFileName)
        except:
            print "No file or syntax error"
            sys.exit(1)
            
        self.root = self.doc.firstChild

        name = self.root.nodeName

        if name != "autovalid":
            print XMLFileName + " :Wrong XML file"
            print sys.exit(1)


    def getReferenceVersion(self):
        return self.getDataFromNode(self.root,"referenceversion")


    def getReferencePath(self):
        return self.getDataFromNode(self.root,"referencepath")
    

    def getStudiesLabels(self):
        studiesLabels=[]

        studyNodes = self.root.getElementsByTagName("study")

        for node in studyNodes :
            label = str(node.attributes["label"].value)
            status = str(node.attributes["status"].value)
            if status == 'on':
                studiesLabels.append(label)

        return studiesLabels


    def getVariablesDefinition(self, studyLabel):
        variables = {}
        
        studyNode = self.getStudyNode(studyLabel)

        variableNodes = studyNode.getElementsByTagName("variable")

        for node in variableNodes :
            label = str(node.attributes["label"].value)
            status = str(node.attributes["status"].value)
            if status == 'on' :
                variables[label] = float(self.getDataFromNode(node,"tolerance"))

        return variables   
        

    def getCasesLabels(self, studyLabel):
        casesLabels=[]

        studyNode = self.getStudyNode(studyLabel)
    
        caseNodes = studyNode.getElementsByTagName("case")

        for node in caseNodes :
            label = str(node.attributes["label"].value)
            status = str(node.attributes["status"].value)
            compute = str(node.attributes["compute"].value)
            if status == 'on':
                casesLabels.append(label)
       
        return casesLabels


    def getIfComputeCase(self, studyLabel, caseLabel):
        ifComputeCase = False
        
        studyNode = self.getStudyNode(studyLabel)
        caseNode = self.getCaseNode(studyNode,caseLabel)

        compute = str(caseNode.attributes["compute"].value)
        if compute == 'on':
            ifComputeCase = True

        return ifComputeCase
            

    def getScriptPostName(self, studyLabel, caseLabel):
        studyNode = self.getStudyNode(studyLabel)
        caseNode = self.getCaseNode(studyNode,caseLabel)

        postNodes = caseNode.getElementsByTagName("post")
        returnLabel = None

        if postNodes != None :
            label = str(postNodes[0].attributes["label"].value)
            status = str(postNodes[0].attributes["status"].value)
            if status == 'on':
                returnLabel = label

        return returnLabel


    def getProcNumber(self, studyLabel, caseLabel):
        studyNode = self.getStudyNode(studyLabel)
        caseNode = self.getCaseNode(studyNode,caseLabel)

        try:
            procNumber = int(self.getDataFromNode(caseNode,"nproc"))
        except:
            procNumber = 1

        return procNumber

    
    
    def getDataFromNode(self, node, childName):
        data = None
        list = node.getElementsByTagName(childName)
              
        if (list.length == 1):
            current = list.item(0)
            data = current.firstChild.data
                    
        return data


    def getStudyNode(self, studyLabel):        
        studyNodes = self.root.getElementsByTagName("study")

        for node in studyNodes :
            label = str(node.attributes["label"].value)
            if label == studyLabel :
                studyNode = node

        return studyNode

    
    def getCaseNode(self, studyNode, caseLabel):        
        caseNodes = studyNode.getElementsByTagName("case")

        for node in caseNodes :
            label = str(node.attributes["label"].value)
            if label == caseLabel :
                caseNode = node

        return caseNode

       
        
