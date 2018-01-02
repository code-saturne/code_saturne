# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

import os, sys, string, logging


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
        except SyntaxError, m:
            print "Syntax Error '%s' in the mapper" % expr
        raise

#-------------------------------------------------------------------------------

class CFDGUI_Management:
    """
    Dock windows are managed by CFDGUI_Management class
    CFDGUI_Management.d_CfdCases[studyId].append([dock,mw.dockWidgetBrowser,mw,aStudyCFD,aCaseCFD,xmlFileName,sobjXML])
    """
    def __init__(self):
      self.dockPosInListe                  = 0
      self.dockWBPosInListe                = 1
      self.mwCFDPosInListe                 = 2
      self.studyCFDPosInListe              = 3
      self.caseCFDPosInListe               = 4
      self.xmlCFDFileNamePosInListe        = 5
      self.sobjXmlPosInListe               = 6

      self.nbelem                          = 7

      self.dock                            = None
      self.dockWB                          = None
      self.aMwCFD                          = None
      self.aStudyCFD                       = None
      self.aCaseCFD                        = None
      self.aXmlCFDFile                     = None
      self.sobjXml                         = None

      self.d_CfdCases                      = {}


    def set_d_CfdCases(self, studyId,
                       dock, dockWB, mwCFD,
                       aStudyCFD, aCaseCFD,
                       axmlCFDFile, sobjXml):
      """
      Add a new Solver GUI in the SALOME desktop.
      """
      if studyId not in self.d_CfdCases.keys():
          self.d_CfdCases[studyId] = []

      self.d_CfdCases[studyId].append([dock, dockWB, mwCFD,
                                       aStudyCFD, aCaseCFD,
                                       axmlCFDFile, sobjXml])

      self.dock          = dock
      self.dockWB        = dockWB
      self.aMwCFD        = mwCFD
      self.aStudyCFD     = aStudyCFD
      self.aCaseCFD      = aCaseCFD
      self.aXmlCFDFile   = axmlCFDFile
      self.sobjXml       = sobjXml

      log.debug("set_d_CfdCases \n\tdock = %s\n\tdockWB = %s\n\tmwCFD = %s\n\taStudyCFD = %s\n\taCaseCFD = %s\n\taxmlCFDFile = %s" % \
                 (dock, dockWB, mwCFD, aStudyCFD, aCaseCFD, axmlCFDFile))


    def checkDockWindowsLists(self, studyId):
      if studyId in self.d_CfdCases.keys():
          return True
      else:
          return False


    def getdockWB(self, studyId, dock):
      dockWB = None
      if self.checkDockWindowsLists(studyId):
          d = self.getDocks(studyId)
          if dock in d.keys():
              ind = d[dock]
              dockWB = self.d_CfdCases[studyId][ind][self.dockWBPosInListe]
      return dockWB


    def getdock(self, studyId, dockWB):
      dock = None
      if self.checkDockWindowsLists(studyId):
          d = self.getDocksWB(studyId)
          if dockWB in d.keys():
              ind = d[dockWB]
              dock = self.d_CfdCases[studyId][ind][self.dockPosInListe]
      return dock


    def getDockListes(self, studyId):
      dockListe = []
      dockListeWB = []
      if self.checkDockWindowsLists(studyId):
          for liste in self.d_CfdCases[studyId]:
              dockListe.append(liste[self.dockPosInListe])
              dockListeWB.append(liste[self.dockWBPosInListe])
      return dockListe, dockListeWB


    def getElem(self, studyId, elempos):
        d = {}
        if elempos not in range(self.nbelem):
            return d
        if studyId in self.d_CfdCases.keys():
            for liste in self.d_CfdCases[studyId]:
                d[liste[elempos]] = self.d_CfdCases[studyId].index(liste)
        return d


    def getDocks(self, studyId):
        """
        return a dictionary d
        """
        return self.getElem(studyId, self.dockPosInListe)


    def getDocksWB(self, studyId):
        """
        """
        return self.getElem(studyId, self.dockWBPosInListe)

    def getDockWithCFDStudyAndCaseNames(self, studyId, studyCFDName, caseName):
        l = []
        if self.d_CfdCases == {} :
            return l
        for liste in self.d_CfdCases[studyId]:
            if liste[self.studyCFDPosInListe].GetName() == studyCFDName \
              and liste[self.caseCFDPosInListe].GetName() == caseName:
                l.append(liste)
        return l

    def getDockWithCFDNames(self, studyId, studyCFDName, caseName, xmlName):
        l = []
        for liste in self.d_CfdCases[studyId]:
            if liste[self.studyCFDPosInListe].GetName() == studyCFDName \
              and liste[self.caseCFDPosInListe].GetName() == caseName \
              and liste[self.xmlCFDFileNamePosInListe] == xmlName:
                l = liste
        return l


    def getStudyCaseXmlNames(self, studyId, mw):
        log.debug("getStudyCaseXmlNames mw = %s" % mw)
        if studyId in self.d_CfdCases.keys():
            for l in self.d_CfdCases[studyId]:
                if l[self.mwCFDPosInListe] == mw:
                    return l[self.studyCFDPosInListe].GetName(), \
                            l[self.caseCFDPosInListe].GetName(), \
                            l[self.xmlCFDFileNamePosInListe]
        return None, None, None


    def getCase(self, studyId, mw):
        if studyId in self.d_CfdCases.keys():
            for l in self.d_CfdCases[studyId]:
                if l[self.mwCFDPosInListe] == mw:
                    return l[self.caseCFDPosInListe]
        return None


    def hideDocks(self,studyId):
        if not self.checkDockWindowsLists(studyId):
            return
        for liste in self.d_CfdCases[studyId]:
            for pos in [self.dockPosInListe, self.dockWBPosInListe]:
                if liste[pos] != None:
                    liste[pos].hide()
                    liste[pos].toggleViewAction().setVisible(False)


    def showDocks(self, studyId):
        if not self.checkDockWindowsLists(studyId):
            return
        for liste in self.d_CfdCases[studyId]:
            for pos in [self.dockPosInListe, self.dockWBPosInListe]:
                if liste[pos] != None:
                    liste[pos].show()
                    liste[pos].setVisible(True)
                    liste[pos].toggleViewAction().setVisible(True)


    def findElem(self, xmlName, caseName, studyCFDName):
        boo = False
        for studyId in self.d_CfdCases.keys():
            for l in self.d_CfdCases[studyId]:
                if l[self.xmlCFDFileNamePosInListe] == xmlName:
                    if l[self.caseCFDPosInListe].GetName() == caseName:
                        if l[self.studyCFDPosInListe].GetName() == studyCFDName:
                            for pos in [self.dockPosInListe,self.dockWBPosInListe]:
                                l[pos].show()
                                l[pos].raise_()
                                l[pos].setVisible(True)
                                l[pos].toggleViewAction().setVisible(True)
                                boo = True
        return boo

    def findDock(self, xmlName, caseName, studyCFDName):
        boo = False
        for studyId in self.d_CfdCases.keys():
            for l in self.d_CfdCases[studyId]:
                if l[self.xmlCFDFileNamePosInListe] == xmlName:
                    if l[self.caseCFDPosInListe].GetName() == caseName:
                        if l[self.studyCFDPosInListe].GetName() == studyCFDName:
                            boo = True
        return boo

    def showDockWindows(self, studyId, xmlName, caseName, studyCFDName):
        for l in self.d_CfdCases[studyId]:
            if l[self.xmlCFDFileNamePosInListe] == xmlName:
                if l[self.caseCFDPosInListe].GetName() == caseName:
                    if l[self.studyCFDPosInListe].GetName() == studyCFDName:
                        for pos in [self.dockPosInListe, self.dockWBPosInListe]:
                            l[pos].show()
                            l[pos].raise_()
                            l[pos].setVisible(True)
                            l[pos].toggleViewAction().setVisible(True)


    def getMW(self, studyId, dock):
        """
        return mW CFD window attached to dock in the liste d_CfdCases[StudyId]
        """
        d = self.getDocks(studyId)
        if d != {}:
            if dock in d.keys():
                return self.d_CfdCases[studyId][d[dock]][self.mwCFDPosInListe]
        else:
            return None


    def delDockfromStudyAndCaseNames(self, dsk, studyId, studyCFDName, caseName):
        """
        Delete all the opened dock windows from a study name and a case name
        """
        liste = self.getDockWithCFDStudyAndCaseNames(studyId, studyCFDName, caseName)
        if liste == []:
            return
        for ll in liste :
            dockcfd, docwb = ll[self.dockPosInListe], ll[self.dockWBPosInListe]
            for dock in [dockcfd, docwb]:
                if dock != None:
                    dsk.removeDockWidget(dock)
                    dock.setParent(None)
                    dock.close()
            # remove the liste which contains the removed docks in the dictionary
            self.d_CfdCases[studyId].remove(ll)

    def delDock(self, dsk, studyId, studyCFDName, caseName, xmlName):
        """
        Delete the opened dock window from a study name, a case name, a xml file name
        """
        liste = self.getDockWithCFDNames(studyId, studyCFDName, caseName, xmlName)
        if liste == []:
            return
        dockcfd, docwb = liste[self.dockPosInListe], liste[self.dockWBPosInListe]
        for dock in [dockcfd, docwb]:
            if dock != None:
                dsk.removeDockWidget(dock)
                dock.setParent(None)
                dock.close()
        # remove the liste which contains the removed docks in the dictionary
        self.d_CfdCases[studyId].remove(liste)


    def cleanAllDock(self, dsk):
        """
        clean all dock windows of cfd cases and clean attached dictionary;
        called when closing salome study and remaining into the desktop
        """
        if self.d_CfdCases == {} : return
        for liste in self.d_CfdCases.values() :
            for liste_object in liste :
                dockcfd, docwb = liste_object[self.dockPosInListe], liste_object[self.dockWBPosInListe]
                for dock in [dockcfd, docwb]:
                    if dock != None:
                        dsk.removeDockWidget(dock)
                        dock.setParent(None)
                        dock.close()
        # clean the associated dictionary
        self.d_CfdCases.clear()


    def tabifyDockWindows(self,dsk,studyId):
        """
        tabify all opened CFD windows and window CFD Browser
        force le regroupement en onglets des fenetres d'etudes CFD
        """
        docListe, docListeWB = self.getDockListes(studyId)

        if len(docListe) > 1:
            for i in range(1,len(docListe)):
                dsk.tabifyDockWidget(docListe[0], docListe[i])

        if len(docListeWB) > 1:
            for i in range(1,len(docListeWB)):
                dsk.tabifyDockWidget(docListeWB[0], docListeWB[i])

#-------------------------------------------------------------------------------
