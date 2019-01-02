# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module contains the following classes and function:
- preprocessorFile
- Informations
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, re, os, logging
import os.path

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore import *
from code_saturne.Base.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import getopenfilename

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PreprocessingInformationsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Function to select the Preprocessor log
#-------------------------------------------------------------------------------

def preprocessorFile(parent, initdir):
    """
    Verify if the choses of file is correct
    """
    file_name = ""
    title = tr("Select a preprocessor log")
    filetypes = "Preprocessor log (*.log);;All Files (*)"
    filt = "All files (*)"
    initdir = os.path.join(initdir, 'check_mesh')
    file_name, _selfilter = getopenfilename(parent, title, initdir, filetypes)
    file_name = str(file_name)

    if file_name:
        f = open(file_name, 'r')
        lines = f.readlines()
        f.close()

        j=-1
        for i in range(len(lines)):
            index = lines[i].rfind("Code_Saturne")
            if index != -1:
                j = i
                break

        if j == -1:
            title = tr("Informations")
            msg = tr("Warning : the selected file is not a correct file.\n\n"\
                     "Verify your selection")
            QMessageBox.information(parent, title, msg)

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
            raise ValueError("Preprocessor log language unknown.")

        refList, groupList = self.getLists(lines)

        self.refList = refList
        self.groupList = groupList


    def readFile(self, file):
        if not file:
            return []
        else:
            f = open(file, 'r')
            lines = f.readlines()
            f.close()

            lang = ""
            for i in range(len(lines)):
                index = lines[i].find("familles de faces et cellules")
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
                self.str4 = 'Groupe'
                self.str5 = 'Configuration locale du cas'
            else:
                self.str1 = "Definition of face and cell families"
                if self.chain == 'faces':
                    self.str2 = "Number of boundary faces"
                elif self.chain == 'cells':
                    self.str2 = "Number of cells"
                self.str3 = 'Family'
                self.str4 = 'Group'
                self.str5 = 'Local case configuration'

            return lines


    def getLists(self, lines):
        refList = []
        groupList = []
        j = len(lines)
        for i in range(len(lines)):
            index = re.search(self.str1, lines[i])
            if index != None:
                j = i
                break

        for n in range(j,len(lines),1):
            if re.search(self.str5, lines[n]): break
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
                                gr = re.split(self.str4 + ' ', lines[p])
                                for g in gr[1:]:
                                    group = re.split('\n', g)[0]
                                    group = re.split('"',group)[1]
                                    if group:
                                        if group not in groupList: groupList.append(group)
                            if re.search(self.str3, lines[p]) != None:
                                n = m-1
                                break

        return refList, groupList


    def getLocalizations(self):
        lst = self.refList + self.groupList
        return list(map(str, lst))


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
