# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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
log.setLevel(logging.DEBUG)
#log.setLevel(logging.NOTSET)

#MP1 2011/02/14 nouvelle classe CFDGUI_Management permettant a terme de remplacer _d_DockWindows, _d_DockWindowsBrowser, #_d_DockWindowsRuncase
# CFDGUI_Management.d_CfdCases = {SalomeStudyId : [[DockWindow, DockWindowBrowser, MwCFD, AStudyCFD, ACaseCFD, AxmlFileName], ...]}
#CFDGUI_Management.d_CfdCases[studyId].append([dock,mw.dockWidgetBrowser,mw,aStudyCFD,aCaseCFD,xmlFileName,sobjXML,adockWRuncase,])
#MP1 Fin

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

_d_DockWindowsRuncase = {}

#-------------------------------------------------------------------------------
# Class definitions
#-------------------------------------------------------------------------------
class Mapper :
    def __init__(self, d1, d2 = {} ):
        self.d1 = d1
        self.d2 = d2
    def __getitem__(self,expr):
        try:
            return eval(expr, self.d1, self.d2)
        except SyntaxError, m :
            print "Syntax Error '%s' in the mapper" % expr
        raise
#-------------------------------------------------------------------------------

class CFDGUI_Management:
    model_windowCFD = """\
List of CFD STUDY CASES in SALOME StudyId : %(StudyID)d
Dock window Name             : %(dockName)s
Dock window Browser Name     : %(dockWBName)s
MainView CFD Name            : %(mwCFDName)s
Study CFD Name               : %(aStudyCFDName)s
Case CFD Name                : %(aCaseCFDName)s
XML CFD File Name            : %(axmlCFDFileName)s
Run text edit window         : %(adockWRuncaseName)s
"""

#################################################################################
    def __init__(self):
      """
      """
      self.dockPosInListe                  = 0
      self.dockWBPosInListe                = 1
      self.mwCFDPosInListe                 = 2
      self.studyCFDPosInListe              = 3
      self.caseCFDPosInListe               = 4
      self.xmlCFDFileNamePosInListe        = 5
      self.sobjXmlPosInListe               = 6
      self.dockWRuncasePosInListe          = 7

      self.nbelem                          = 8

      self.dockName                        = ""
      self.dockWBName                      = ""
      self.aMwCFDName                      = ""
      self.aStudyCFDName                   = ""
      self.aCaseCFDName                    = ""
      self.aXmlCFDFileName                 = ""
      self.aDockWRuncaseName               = ""

      self.dock                            = None
      self.dockWB                          = None
      self.aMwCFD                          = None
      self.aStudyCFD                       = None
      self.aCaseCFD                        = None
      self.aXmlCFDFile                     = None
      self.sobjXml                         = None
      self.aDockWRuncase                   = None

      self.d_CfdCases = {}
      self.studyId                         = None

#################################################################################

    def format(self,studyId):
      StudyID                              = studyId
      dockName                             = self.dockName
      dockWBName                           = self.dockWBName
      mwCFDName                            = self.aMwCFDName
      aStudyCFDName                        = self.aStudyCFDName
      aCaseCFDName                         = self.aCaseCFDName
      axmlCFDFileName                      = self.aXmlCFDFileName
      adockWRuncaseName                    = self.aDockWRuncaseName
      return self.model_windowCFD % Mapper(locals())

#################################################################################

    def set_d_CfdCases(self,studyId,dock,dockWB,mwCFD,aStudyCFD,aCaseCFD,axmlCFDFile,sobjXml,adockWRuncase) :
      if studyId not in self.d_CfdCases.keys() :
          self.d_CfdCases[studyId] = []
      self.studyId = studyId
      self.d_CfdCases[studyId].append([dock,dockWB,mwCFD,aStudyCFD,aCaseCFD,axmlCFDFile,sobjXml,adockWRuncase])
      self.dock                            = dock
      self.dockWB                          = dockWB
      self.aMwCFD                          = mwCFD
      self.aStudyCFD                       = aStudyCFD
      self.aCaseCFD                        = aCaseCFD
      self.aXmlCFDFile                     = axmlCFDFile
      self.sobjXml                         = sobjXml
      self.aDockWRuncase                   = adockWRuncase

#################################################################################

    def getdockWB(self,studyId,dock):
      dockWB = None
      if self.checkDockWindowsLists(studyId) :
        d = self.getDocks(studyId)
        if dock in d.keys() :
          ind = d[dock]
          dockWB = self.d_CfdCases[studyId][ind][self.dockWBPosInListe]
      return dockWB

#################################################################################

    def getdock(self,studyId,dockWB):
      dock = None
      if self.checkDockWindowsLists(studyId) :
        d = self.getDocksWB(studyId)
        if dockWB in d.keys() :
          ind = d[dockWB]
          dock = self.d_CfdCases[studyId][ind][self.dockPosInListe]
      return dock

#################################################################################

    def checkDockWindowsLists(self,studyId):
      """
      """
      if studyId in self.d_CfdCases.keys() :
          return True
      else:
          return False

#################################################################################

    def printDockListe(self,dockListe) :
      for dock in dockListe :
        print "dockListe = ",dock.windowTitle()

#################################################################################

    def getDockListes(self,studyId) :
      """
      """
      dockListe = []
      dockListeWB = []
      if self.checkDockWindowsLists(studyId):
        for liste in self.d_CfdCases[studyId] :
          dockListe.append(liste[self.dockPosInListe])
          dockListeWB.append(liste[self.dockWBPosInListe])
      return dockListe,dockListeWB

#################################################################################

    def print_d_CfdCases(self) :

      print "self.d_CfdCases = ",self.d_CfdCases
      for studySalome in self.d_CfdCases.keys() :
        for l_winValue in self.d_CfdCases[studySalome] :
          if l_winValue[self.dockPosInListe] != None :
            self.dockName           = l_winValue[self.dockPosInListe].windowTitle()
          if l_winValue[self.dockWBPosInListe] != None :
            self.dockWBName         = l_winValue[self.dockWBPosInListe].windowTitle()
          if l_winValue[self.mwCFDPosInListe] != None :
            self.aMwCFDName         = l_winValue[self.mwCFDPosInListe].windowTitle()
          if l_winValue[self.studyCFDPosInListe] != None :
            self.aStudyCFDName      = l_winValue[self.studyCFDPosInListe].GetName()
          if l_winValue[self.caseCFDPosInListe] != None :
            self.aCaseCFDName       = l_winValue[self.caseCFDPosInListe].GetName()
          if l_winValue[self.xmlCFDFileNamePosInListe] != None :
            self.aXmlCFDFileName    = l_winValue[self.xmlCFDFileNamePosInListe]
          if l_winValue[self.dockWRuncasePosInListe] != None :
            self.aDockWRuncaseName  = l_winValue[self.dockWRuncasePosInListe].windowTitle()
          print self.format(studySalome)
          print "  "
          

#################################################################################

    def getElem(self,studyId,elempos) :
      """
      """
      d = {}
      if elempos not in range(self.nbelem) : return
      if self.d_CfdCases != {} :
        if studyId in self.d_CfdCases.keys() :
          for liste in self.d_CfdCases[studyId] :
            d[liste[elempos]] = self.d_CfdCases[studyId].index(liste)
      return d

#################################################################################

    def getDocks(self,studyId) :
      """
      return a dictionary d 
      """
      d = self.getElem(studyId,self.dockPosInListe)
      return d

#################################################################################

    def getDocksWB(self,studyId) :
      """
      """
      d = self.getElem(studyId,self.dockWBPosInListe)
      return d

#################################################################################

    def getDockWithCFDNames(self,studyId,studyCFDName, caseName, xmlName) :
      l = []
      for liste in self.d_CfdCases[studyId] :
        if liste[self.studyCFDPosInListe].GetName() == studyCFDName and liste[self.caseCFDPosInListe].GetName() == caseName and liste[self.xmlCFDFileNamePosInListe] == xmlName :
          l = liste
      return l

#################################################################################

    def getListElem(self,studyId,elempos) :
      """
      """
      l = []
      if elempos not in range(self.nbelem) : return
      if self.d_CfdCases != {} :
        if studyId in self.d_CfdCases.keys() :
          for liste in self.d_CfdCases[studyId] :
            l.append(liste[elempos])
      return l

#################################################################################

    def getListSobj(self,studyId) :
      """
      return a list of Sobj corresponding to the opened dock window for xml CFD file  
      """
      liste = self.getListElem(studyId,self.sobjXmlPosInListe)
      return liste

#################################################################################

    def getStudyCaseXmlNames(self,studyId,mw) :
      if self.d_CfdCases.keys() != [] :
        if studyId in self.d_CfdCases.keys() :
          for l in self.d_CfdCases[studyId] :
            if l[self.mwCFDPosInListe] == mw :
              return l[self.studyCFDPosInListe].GetName(),l[self.caseCFDPosInListe].GetName(),l[self.xmlCFDFileNamePosInListe]
      return None,None,None

#################################################################################

    def getCase(self,studyId,mw) :
      if self.d_CfdCases.keys() != [] :
        if studyId in self.d_CfdCases.keys() :
          for l in self.d_CfdCases[studyId] :
            if l[self.mwCFDPosInListe] == mw :
              return l[self.caseCFDPosInListe]
      return None

#################################################################################

    def hideDocks(self,studyId) :
      if not self.checkDockWindowsLists(studyId) : return
      for liste in self.d_CfdCases[studyId] :
        for pos in [self.dockPosInListe,self.dockWBPosInListe,self.dockWRuncasePosInListe] :
          if liste[pos] != None :
            liste[pos].hide()
            liste[pos].toggleViewAction().setVisible(False)

#################################################################################

    def showDocks(self,studyId) :
      if not self.checkDockWindowsLists(studyId) : return
      for liste in self.d_CfdCases[studyId] :
        for pos in [self.dockPosInListe,self.dockWBPosInListe,self.dockWRuncasePosInListe] :
          if liste[pos] != None :
            
            liste[pos].show()
            liste[pos].setVisible(True)
            liste[pos].toggleViewAction().setVisible(True)
           
#################################################################################

    def findElem(self,xmlName, caseName, studyCFDName):
        
      boo = False
      if self.d_CfdCases.keys() != [] :
        for studyId in self.d_CfdCases.keys() :
          for l in self.d_CfdCases[studyId] :
            if l[self.xmlCFDFileNamePosInListe] == xmlName :
              if l[self.caseCFDPosInListe].GetName() == caseName :
                if l[self.studyCFDPosInListe].GetName() == studyCFDName :
                  for pos in [self.dockPosInListe,self.dockWBPosInListe]:
                    l[pos].show()
                    l[pos].raise_()
                    l[pos].setVisible(True)
                    l[pos].toggleViewAction().setVisible(True)
                    boo = True
      return boo

#################################################################################

    def showDockWindows(self,studyId,xmlName, caseName, studyCFDName) :

      for l in self.d_CfdCases[studyId] : 
        if l[self.xmlCFDFileNamePosInListe] == xmlName :
          if l[self.caseCFDPosInListe].GetName() == caseName :
            if l[self.studyCFDPosInListe].GetName() == studyCFDName :
              for pos in [self.dockPosInListe,self.dockWBPosInListe]:
                l[pos].show()
                l[pos].raise_()
                l[pos].setVisible(True)
                l[pos].toggleViewAction().setVisible(True)
          return

#################################################################################

    def getDockId(self,studyId,dock) :
      """
      return position Id of the list attached to dock in the liste d_CfdCases[StudyId]
      """
      d = self.getDocks(studyId)
      if d != {} :
        if dock in d.keys() :
          return d[dock]
      else :
        return None

#################################################################################

    def getMW(self,studyId,dock) :
      """
      return mW CFD window attached to dock in the liste d_CfdCases[StudyId]
      """
      d = self.getDocks(studyId)
      if d != {} :
        if dock in d.keys() :
          return self.d_CfdCases[studyId][d[dock]][self.mwCFDPosInListe]
      else :
        return None

#################################################################################

    def getSobjXml(self,studyId,dock) :
      """
      return mW CFD window attached to dock in the liste d_CfdCases[StudyId]
      """
      d = self.getDocks(studyId)
      if d != {} :
        if dock in d.keys() :
          return self.d_CfdCases[studyId][d[dock]][self.sobjXmlPosInListe]
      else :
        return None

#################################################################################

    def getCaseCFD(self,studyId,dock) :
        """
        return mW CFD window attached to dock in the liste d_CfdCases[StudyId]
        """
        d = self.getDocks(studyId)
        if d != {} :
            if dock in d.keys() :
                return self.d_CfdCases[studyId][d[dock]][self.caseCFDPosInListe]
        else :
            return None

#################################################################################

    def delDock(self,dsk,studyId,studyCFDName, caseName, xmlName) :

      liste = self.getDockWithCFDNames(studyId,studyCFDName, caseName, xmlName)
      if liste == [] : return
      dockcfd,docwb = liste[self.dockPosInListe],liste[self.dockWBPosInListe]
      for dock in [dockcfd,docwb] :
        if dock != None :
          dsk.removeDockWidget(dock)
          dock.setParent(None)
          dock.close()
      # remove the liste which contains the removed docks in the dictionary
      self.d_CfdCases[studyId].remove(liste)

#################################################################################

    def tabifyDockWindows(self,dsk,studyId):
      """
      tabify all opened CFD windows and window CFD Browser
      force le regroupement en onglets des fenetres d'etudes CFD
      """
      
      docListe,docListeWB = self.getDockListes(studyId)

      if _d_DockWindowsRuncase.has_key(studyId) :
        docListe = docListe+_d_DockWindowsRuncase[studyId]
      if len(docListe) > 1:
        for i in range(1,len(docListe)):
          dsk.tabifyDockWidget(docListe[0], docListe[i])

      if len(docListeWB) > 1:
        for i in range(1,len(docListeWB)):
          dsk.tabifyDockWidget(docListeWB[0], docListeWB[i])

#################################################################################

