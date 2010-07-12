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
This module contains the following classes and function:
- preprocessorFile
- Informations
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, re, os, string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
import Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PreprocessingInformationsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Function for select the Preprocessor listing
#-------------------------------------------------------------------------------

def preprocessorFile(parent, initdir):
    """
    Verify if the choses of file is correct
    """
    file_name = ""
    title = tr("Select a Code_Saturne Preprocessor listing")
    filetypes = "Preprocessor listing (listpre.*);;All Files (*)"
    filt = "All files (*)"
    file_name = QFileDialog.getOpenFileName(parent, title, initdir, filetypes, filt)
    file_name = str(file_name)

    if file_name:
        f = open(file_name, 'r')
        lines = f.readlines()
        f.close()

        j=0
        for i in range(len(lines)):
            index = string.rfind(lines[i], "ECS   version")
            if index != -1:
                j = i
                break

        if j == 0:
            title = tr("Informations")
            msg = tr("Warning : the selected file is not a correct file.\n\n"\
                     "Verify your selection")
            QMessageBox.information(parent, title, msg)
            preprocessorFile(initdir) # ???

    return file_name



def tr(text):
    """
    Translation
    """
    return text

#-------------------------------------------------------------------------------
# Informations class
#-------------------------------------------------------------------------------

class Informations:
    def __init__(self, file, chain):

        self.chain = chain
        if self.chain not in ('faces', 'cells'):
            raise ValueError("Informations class is called with a wrong parameter 'chain'")

        lines = self.readFile(file)
        if not lines:
            raise ValueError("Code_Saturne Preprocessor listing language unknown.")

        refList, groupList = self.getListes(lines)

        self.refList = refList
        self.groupList = groupList

        self.updateListes()


    def readFile(self, file):
        if not file:
            return []
        else:
            f = open(file, 'r')
            lines = f.readlines()
            f.close()

            lang = ""
            for i in range(len(lines)):
                #index = re.search("processeur", lines[i])
                index = string.rfind(lines[i], "processeur")
                if index > 0:
                    lang = 'fr'
                    break

            if lang == 'fr':
                self.str1 = "finition des familles de faces et cellules"
                if self.chain == 'faces':
                   self.str2 = "Nombre de faces de bord"
                elif self.chain == 'cells':
                   self.str2 = "Nombre de cellules"
                self.str3 = 'Famille'
                self.str4 = 'Couleur'
                self.str5 = 'Groupe'
                self.str6 = 'Définition des couleurs et groupes en fonction des familles'
            else:
                self.str1 = "Definition of face and cell families"
                if self.chain == 'faces':
                    self.str2 = "Number of boundary faces"
                elif self.chain == 'cells':
                    self.str2 = "Number of cells"
                self.str3 = 'Family'
                self.str4 = 'Color'
                self.str5 = 'Group'
                self.str6 = 'Definition of colors and groups beyond families'

            return lines


    def getListes(self, lines):
        refList = []
        groupList = []
        for i in range(len(lines)):
            index = re.search(self.str1, lines[i])
            if index != None:
                j = i
                break

        for n in range(j,len(lines),1):
            if re.search(self.str6, lines[n]): break
            index = re.search(self.str3, lines[n])
            if index != None:
                familyList =[]
                fam = re.split(self.str3, lines[n])
                for f in fam[1:]: numfam = re.split('\n', f)[0]
                for m in range(n+1,len(lines),1):
                    if re.search(self.str2, lines[m]) != None:
                        p = m
                        for p in range(p-1,p-m,-1):
                            if re.search(self.str4, lines[p]) != None:
                                coul = re.split(self.str4 + ' ', lines[p])
                                for cl in coul[1:]:
                                    ref = re.split('\n', cl)[0]
                                    if ref:
                                        if ref not in refList: refList.append(ref)
                            if re.search(self.str5, lines[p]) != None:
                                gr = re.split(self.str5 + ' ', lines[p])
                                for g in gr[1:]:
                                    group = re.split('\n', g)[0]
                                    group = re.split('"',group)[1]
                                    if group:
                                        if group not in groupList: groupList.append(group)
                            if re.search(self.str3, lines[p]) != None:
                                n = m-1
                                break

        return refList, groupList


    # TODO: delete this method
    def updateListes(self):
        return self.refList, self.groupList


    def getLocalizations(self):
        list = self.refList + self.groupList
        return map(str, list)


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
