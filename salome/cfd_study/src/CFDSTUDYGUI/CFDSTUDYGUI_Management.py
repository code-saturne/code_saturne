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
Management of windows in the SALOME GUI
==========
called by CFDSTUDYGUI_SolverGUI.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, logging


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_Management")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Class definitions
#-------------------------------------------------------------------------------

class Mapper:
    def __init__(self, d1, d2 = {} ):
        self.d1 = d1
        self.d2 = d2
    def __getitem__(self,expr):
        try:
            return eval(expr, self.d1, self.d2)
        except SyntaxError as m:
            print("Syntax Error '%s' in the mapper" % expr)
        raise

#-------------------------------------------------------------------------------

class CFDGUI_Management:
    """
    Dock windows are managed by CFDGUI_Management class
    CFDGUI_Management.d_CfdCases.append([dock,mw,aStudyCFD,aCaseCFD,xmlFileName,sobjXML])
    """
    def __init__(self):
      self.dockPosInListe                  = 0
      self.mwCFDPosInListe                 = 1
      self.studyCFDPosInListe              = 2
      self.caseCFDPosInListe               = 3
      self.xmlCFDFileNamePosInListe        = 4
      self.sobjXmlPosInListe               = 5

      self.nbelem                          = 6

      self.dock                            = None
      self.aMwCFD                          = None
      self.aStudyCFD                       = None
      self.aCaseCFD                        = None
      self.aXmlCFDFile                     = None
      self.sobjXml                         = None

      self.d_CfdCases                      = []


    def set_d_CfdCases(self,
                       dock, mwCFD,
                       aStudyCFD, aCaseCFD,
                       axmlCFDFile, sobjXml):
      """
      Add a new Solver GUI in the SALOME desktop.
      """

      self.d_CfdCases.append([dock, mwCFD,
                                       aStudyCFD, aCaseCFD,
                                       axmlCFDFile, sobjXml])

      self.dock          = dock
      self.aMwCFD        = mwCFD
      self.aStudyCFD     = aStudyCFD
      self.aCaseCFD      = aCaseCFD
      self.aXmlCFDFile   = axmlCFDFile
      self.sobjXml       = sobjXml

      log.debug("set_d_CfdCases \n\tdock = %s\n\tmwCFD = %s\n\taStudyCFD = %s\n\taCaseCFD = %s\n\taxmlCFDFile = %s" % \
                 (dock, mwCFD, aStudyCFD, aCaseCFD, axmlCFDFile))


    def checkDockWindowsLists(self):
      if self.d_CfdCases!= []:
          return True
      else:
          return False


    def getDockListe(self):
      """
      return a liste which contains all the CFD DockWidget instances opened in GUI
      """
      listeDock = []
      if self.checkDockWindowsLists():
          for liste in self.d_CfdCases:
              listeDock.append(liste[self.dockPosInListe])
      return listeDock


    def getElem(self, elempos):

        d = {}
        if elempos not in list(range(self.nbelem)):
            return d
        if self.checkDockWindowsLists():
            for liste in self.d_CfdCases:
                d[liste[elempos]] = self.d_CfdCases.index(liste)
        return d


    def getDocks(self):
        """
        return a liste l - called in CFDSTUDYGUI_SolverGUI
        """
        return self.getElem(self.dockPosInListe)


    def getDockListeWithCFDStudyAndCaseNames(self, studyCFDName, caseName):
        """
        a CFD case can have more xml files
        """
        l = []
        if not self.checkDockWindowsLists():
            return l
        for liste in self.d_CfdCases:
            if liste[self.studyCFDPosInListe].GetName() == studyCFDName \
              and liste[self.caseCFDPosInListe].GetName() == caseName:
                l.append(liste)
        return l


    def getListeWithCFDNames(self, studyCFDName, caseName, xmlName):
        l = []
        for liste in self.d_CfdCases:
            if liste[self.studyCFDPosInListe].GetName() == studyCFDName \
              and liste[self.caseCFDPosInListe].GetName() == caseName \
              and liste[self.xmlCFDFileNamePosInListe] == xmlName:
                return liste
        return l


    def getDockWithCFDNames(self, studyCFDName, caseName, xmlName):
        l = self.getListeWithCFDNames(studyCFDName, caseName, xmlName)
        if l != []:
            return l[self.dockPosInListe]
        else:
            return None


    def getStudyCaseXmlNames(self, mw):
        log.debug("getStudyCaseXmlNames mw = %s" % mw)
        if self.checkDockWindowsLists():
            for l in self.d_CfdCases:
                if l[self.mwCFDPosInListe] == mw:
                    return l[self.studyCFDPosInListe].GetName(), \
                            l[self.caseCFDPosInListe].GetName(), \
                            l[self.xmlCFDFileNamePosInListe]
        return None, None, None


    def getCase(self, mw):
        if self.checkDockWindowsLists():
            for l in self.d_CfdCases:
                if l[self.mwCFDPosInListe] == mw:
                    return l[self.caseCFDPosInListe]
        return None


    def hideDocks(self):
        if not self.checkDockWindowsLists():
            return
        for liste in self.d_CfdCases:
            if liste[self.dockPosInListe] != None:
                liste[self.dockPosInListe].hide()
                liste[self.dockPosInListe].toggleViewAction().setVisible(False)


    def showDocks(self):
        if not self.checkDockWindowsLists():
            return
        for liste in self.d_CfdCases:
            if liste[self.dockPosInListe] != None:
                liste[self.dockPosInListe].show()
                liste[self.dockPosInListe].setVisible(True)
                liste[self.dockPosInListe].toggleViewAction().setVisible(True)


    def findElem(self, xmlName, caseName, studyCFDName):
        boo = False
        if self.checkDockWindowsLists():
            for l in self.d_CfdCases:
                if l[self.xmlCFDFileNamePosInListe] == xmlName:
                    if l[self.caseCFDPosInListe].GetName() == caseName:
                        if l[self.studyCFDPosInListe].GetName() == studyCFDName:
                            l[self.dockPosInListe].show()
                            l[self.dockPosInListe].raise_()
                            l[self.dockPosInListe].setVisible(True)
                            l[self.dockPosInListe].toggleViewAction().setVisible(True)
                            boo = True
        return boo

    def findDock(self, xmlName, caseName, studyCFDName):
        boo = False
        if self.checkDockWindowsLists():
            for l in self.d_CfdCases:
                if l[self.xmlCFDFileNamePosInListe] == xmlName:
                    if l[self.caseCFDPosInListe].GetName() == caseName:
                        if l[self.studyCFDPosInListe].GetName() == studyCFDName:
                            boo = True
        return boo


    def delDockfromStudyAndCaseNames(self, dsk, studyCFDName, caseName):
        """
        Delete all the opened dock windows from a study name and a case name
        """
        liste = self.getDockListeWithCFDStudyAndCaseNames(studyCFDName, caseName)
        if liste == []:
            return
        for l in liste:
            dockcfd = l[self.dockPosInListe]
            if dockcfd != None:
                dsk.removeDockWidget(dockcfd)
                dockcfd.setParent(None)
                dockcfd.close()
        # remove the liste which contains the removed docks in the main list self.d_CfdCases
            self.d_CfdCases.remove(l)


    def delDock(self, dsk, studyCFDName, caseName, xmlName):
        """
        Delete the opened dock window from a study name, a case name, a xml file name
        """
        liste = self.getListeWithCFDNames(studyCFDName, caseName, xmlName)
        if liste == []:
            return
        dockcfd = liste[self.dockPosInListe]
        if dockcfd != None:
            dsk.removeDockWidget(dockcfd)
            dockcfd.setParent(None)
            dockcfd.close()
        # remove the liste which contains the removed docks in the dictionary
        self.d_CfdCases.remove(liste)


    def cleanAllDock(self, dsk):
        """
        clean all dock windows of cfd cases and clean attached main liste;
        called when closing salome study and remaining into the desktop
        """
        log.debug("cleanAllDock")
        if self.d_CfdCases == [] : return
        for liste_object in self.d_CfdCases :
            dockcfd = liste_object[self.dockPosInListe]
            if dockcfd != None:
                dsk.removeDockWidget(dockcfd)
                dockcfd.setParent(None)
                dockcfd.close()
        # clean the associated dictionary
        self.d_CfdCases.clear()

