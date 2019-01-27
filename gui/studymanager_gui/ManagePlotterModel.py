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
This module manages plots :

This module defines the following classes:
- ManagePlotterModel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLmodel import ModelTest

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class ManagePlotterModel(Model):
    """
    Class for Variables and Scalar model initialization.
    """
    def __init__(self, case):
        """
        Constructor
        """
        self.case = case
        self.repo = self.case.xmlGetNode("studymanager").xmlGetString('repository')


    def defaultInitialValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['subplot_name']      = ""
        default['subplot_legstatus'] = "on"
        default['figure_name']       = "Figure"
        default['figure_title']      = ""
        default['figure_fmt']        = "pdf"
        default['color']             = "bedf"

        return default


    def getStudyList(self):
        """
        Return list of studies in xml file
        """
        list_study = self.case.xmlGetNodeList("study")
        lst = []
        for node in list_study:
            lst.append(node['label'])
        return lst


    def getSubplotList(self, study):
        """
        Return list of subplots for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = []
        for node in study_node.xmlGetNodeList("subplot"):
            lst.append(node['id'])
        return lst


    def getFigureList(self, study):
        """
        Return list of figures for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = []
        for node in study_node.xmlGetNodeList("figure"):
            lst.append(node['id'])
        return lst


    def getMeasurementList(self, study):
        """
        Return list of measurements for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = []
        for node in study_node.xmlGetNodeList("measurement"):
            lst.append(node['file'])
        return lst


    def getCaseList(self, study):
        """
        Return list of case for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = []
        for nn in study_node.xmlGetNodeList("case"):
            lst.append(int(nn['id']))
        return lst


    def addSubplot(self, study):
        """
        add a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        idx = len(study_node.xmlGetNodeList("subplot"))
        study_node.xmlInitChildNode("subplot", id = idx)
        return idx


    def delSubplot(self, study, idx):
        """
        suppress a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        node = study_node.xmlGetNode("subplot", id = idx)
        node.xmlRemoveNode()
        for nn in study_node.xmlGetNodeList("subplot"):
            try:
                if int(nn['id']) > idx:
                    nn['id'] = str(int(nn['id']) - 1)
            except:
                pass

        for nn in study_node.xmlGetNodeList("figure"):
            lst = nn['idlist']
            new_lst = ''
            if lst:
                for idl in lst:
                    if idl != " ":
                        new_idl = idl
                        if int(idl) > idx:
                            new_idl = str(int(idl) - 1)
                        if new_lst != "":
                            new_lst = new_lst + " " + str(new_idl)
                        else:
                            new_lst = str(new_idl)
                    nn['idlist'] = new_lst

        for nn in study_node.xmlGetNodeList("plot"):
            lst = nn['fig']
            new_lst = ''
            if lst:
                for idl in lst:
                    if idl != " ":
                        new_idl = idl
                        if int(idl) > idx:
                            new_idl = str(int(idl) - 1)
                        if new_lst != "":
                            new_lst = new_lst + " " + str(new_idl)
                        else:
                            new_lst = str(new_idl)
                    nn['fig'] = new_lst


    def getSubplotTitle(self, study, idx):
        """
        Return title for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        name = subplot['title']
        if not name:
            name = self.defaultInitialValues()['subplot_name']
        return name


    def setSubplotTitle(self, study, idx, title):
        """
        Set title for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        if title != "":
            subplot['title'] = title


    def getSubplotXLabel(self, study, idx):
        """
        Return X label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        name = subplot['xlabel']
        if not name:
            name = ""
        return name


    def setSubplotXLabel(self, study, idx, title):
        """
        Set X label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        subplot['xlabel'] = title


    def getSubplotYLabel(self, study, idx):
        """
        Return Y label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        name = subplot['ylabel']
        if not name:
            name = ""
        return name


    def setSubplotYLabel(self, study, idx, title):
        """
        Set Y label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        subplot['ylabel'] = title


    def getSubplotLegPos(self, study, idx):
        """
        Return legend position for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        name = subplot['legpos']
        if not name:
            name = ""
        return name


    def setSubplotLegPos(self, study, idx, title):
        """
        Set legend position for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        subplot['legpos'] = title


    def getSubplotLegStatus(self, study, idx):
        """
        Return legend status for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        status = subplot['legstatus']
        if not status:
            status = self.defaultInitialValues()['subplot_legstatus']
            self.setSubplotLegStatus(study, idx, status)
        return status


    def setSubplotLegStatus(self, study, idx, status):
        """
        Set legend status for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        subplot['legstatus'] = status


    def getSubplotXLim(self, study, idx):
        """
        Return X limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        status = subplot['xlim']
        if not status:
            status = ""
        return status


    def setSubplotXLim(self, study, idx, pos):
        """
        Set X limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        if pos != "":
            subplot['xlim'] = pos


    def getSubplotYLim(self, study, idx):
        """
        Return Y limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        status = subplot['ylim']
        if not status:
            status = ""
        return status


    def setSubplotYLim(self, study, idx, pos):
        """
        Set Y limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("subplot", id = idx)
        if pos != "":
            subplot['ylim'] = pos


#-------------------------------------------------------------------------------
# Figure
#-------------------------------------------------------------------------------
    def addFigure(self, study):
        """
        add a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        idx = len(study_node.xmlGetNodeList("figure"))
        study_node.xmlInitChildNode("figure", id = idx)
        return idx


    def delFigure(self, study, idx):
        """
        suppress a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        node = study_node.xmlGetNode("figure", id = idx)
        node.xmlRemoveNode()
        for nn in study_node.xmlGetNodeList("figure"):
            try:
                if int(nn['id']) > idx:
                    nn['id'] = str(int(nn['id']) - 1)
            except:
                pass


    def getFigureName(self, study, idx):
        """
        Return name for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        name = subplot['name']
        if not name:
            name = self.defaultInitialValues()['figure_name']
            self.setFigureName(study, idx, name)
        return name


    def setFigureName(self, study, idx, name):
        """
        Set name for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        subplot['name'] = name


    def getFigureTitle(self, study, idx):
        """
        Return title for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        name = subplot['title']
        if not name:
            name = self.defaultInitialValues()['figure_title']
            self.setFigureTitle(study, idx, name)
        return name


    def setFigureTitle(self, study, idx, name):
        """
        Set title for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        if name != "":
            subplot['title'] = name
        elif subplot['title']:
            subplot['title'] = name


    def getFigureRow(self, study, idx):
        """
        Return number of row for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        row = subplot['nbrow']
        if not row:
            row = ""
        return row


    def setFigureRow(self, study, idx, row):
        """
        Set number of row for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        if row != "":
            subplot['nbrow'] = row


    def getFigureColumn(self, study, idx):
        """
        Return number of column for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        col = subplot['nbcol']
        if not col:
            col = ""
        return col


    def setFigureColumn(self, study, idx, col):
        """
        Set number of column for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        if col != "":
            subplot['nbcol'] = col


    def getFigureFormat(self, study, idx):
        """
        Return format for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        name = subplot['format']
        if not name:
            name = self.defaultInitialValues()['figure_fmt']
            self.setFigureFormat(study, idx, name)
        return name


    def setFigureFormat(self, study, idx, fmt):
        """
        Set format for a figure
        """
        self.isInList(fmt, ("png", "pdf"))
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        subplot['format'] = fmt


    def getFigureIdList(self, study, idx):
        """
        Return subplot id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        lst = subplot['idlist']
        if not lst:
            lst = ""
        return lst


    def setFigureIdList(self, study, idx, idlist):
        """
        Set subplot id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        subplot = study_node.xmlGetNode("figure", id = idx)
        subplot['idlist'] = idlist


#-------------------------------------------------------------------------------
# Measurement
#-------------------------------------------------------------------------------
    def addMeasurementFile(self, study, fle):
        """
        add a measurement file for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        study_node.xmlInitChildNode("measurement", file=fle)


    def delMeasurementFile(self, study, idx):
        """
        suppress a measurement file for a study
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[idx]
        node.xmlRemoveNode()


    def addMeasurementPlot(self, study, measurement_idx):
        """
        add a plot for a measurement file
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        lst2 = study_node.xmlGetNodeList("plot")
        idx = len(lst2)
        node.xmlInitChildNode("plot", id = idx)
        return idx


    def deleteMeasurementPlot(self, study, measurement_idx, idx):
        """
        suppress a plot for a measurement file
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        nn = node.xmlGetNode("plot", id = idx)
        nn.xmlRemoveNode()

        idx = 0
        for n in study_node.xmlGetNodeList("plot"):
            n['id'] = idx
            idx = idx + 1


    def getMeasurementPlotList(self, study, measurement_name):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        measurement_node = study_node.xmlGetNode("measurement", file = measurement_name)
        lst = []
        for node in measurement_node.xmlGetNodeList("plot"):
            lst.append(node['id'])
        return lst


    def getMeasurementColor(self, study, measurement_idx, idx):
        """
        Return plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        color = plot['color']
        if not color:
            color = self.defaultInitialValues()['color']
            self.setMeasurementColor(study, measurement_idx, idx, color)
        return color


    def setMeasurementColor(self, study, measurement_idx, idx, color):
        """
        Set plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['color'] = color


    def getMeasurementFormat(self, study, measurement_idx, idx):
        """
        Return plot format
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        fmt = plot['fmt']
        if not fmt:
            fmt = ""
        return fmt


    def setMeasurementFormat(self, study, measurement_idx, idx, fmt):
        """
        Set plot format
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['fmt'] = fmt


    def getMeasurementLegend(self, study, measurement_idx, idx):
        """
        Return plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        legend = plot['legend']
        if not legend:
            legend = ""
        return legend


    def setMeasurementLegend(self, study, measurement_idx, idx, legend):
        """
        Set plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['legend'] = legend


    def getMeasurementXcol(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['xcol']
        if not val:
            val = ""
        return val


    def setMeasurementXcol(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['xcol'] = val


    def getMeasurementYcol(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['ycol']
        if not val:
            val = ""
        return val


    def setMeasurementYcol(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['ycol'] = val


    def getMeasurementWidth(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['linewidth']
        if not val:
            val = ""
        return val


    def setMeasurementWidth(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['linewidth'] = val


    def getMeasurementMarker(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['markersize']
        if not val:
            val = ""
        return val


    def setMeasurementMarker(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['markersize'] = val


    def getMeasurementXerr(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['xerr']
        if not val:
            val = ""
        return val


    def setMeasurementXerr(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['xerr'] = val


    def getMeasurementXerrp(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['xerrp']
        if not val:
            val = ""
        return val


    def setMeasurementXerrp(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['xerrp'] = val


    def getMeasurementYerr(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['yerr']
        if not val:
            val = ""
        return val


    def setMeasurementYerr(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['yerr'] = val


    def getMeasurementYerrp(self, study, measurement_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        val = plot['yerrp']
        if not val:
            val = ""
        return val


    def setMeasurementYerrp(self, study, measurement_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['yerrp'] = val


    def getMeasurementIdList(self, study, measurement_idx, idx):
        """
        Return subplot id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        lst = plot['fig']
        if not lst:
            lst = ""
        return lst


    def setMeasurementIdList(self, study, measurement_idx, idx, idlist):
        """
        Set subplot id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("measurement")
        node = lst[measurement_idx]
        plot = node.xmlGetNode("plot", id = idx)
        plot['fig'] = idlist


#-------------------------------------------------------------------------------
# Case
#-------------------------------------------------------------------------------
    def getCaseName(self, study, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = idx)
        return nn['label']


    def getCaseDataList(self, study, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = idx)
        lst = []
        for n in nn.xmlGetNodeList("data"):
            lst.append(n['file'])
        return lst


    def addCaseDataFile(self, study, case_idx, name):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = case_idx)
        nn.xmlInitChildNode("data", file=name)


    def delCaseDataFile(self, study, case_idx, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = case_idx)
        lst = nn.xmlGetNodeList("data")
        node = lst[idx]
        node.xmlRemoveNode()


    def addAssociatedCasePlot(self, study, case_idx, data_idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = case_idx)
        node = nn.xmlGetNodeList("data")[data_idx]
        idx = len(study_node.xmlGetNodeList("plot"))
        node.xmlInitChildNode("plot", id = idx)
        return idx


    def delAssociatedCasePlot(self, study, case_idx, data_idx, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNode("case", id = case_idx)
        node = nn.xmlGetNodeList("data")[data_idx]
        lst = node.xmlGetNodeList("plot")
        n = lst[idx]
        n.xmlRemoveNode()

        idx = 0
        for n in study_node.xmlGetNodeList("plot"):
            n['id'] = idx
            idx = idx + 1


    def getCasePlotList(self, study, idx, data_idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        case = study_node.xmlGetNode("case", id = idx)
        node = case.xmlGetNodeList("data")[data_idx]
        lst = []
        for n in node.xmlGetNodeList("plot"):
            lst.append(n['id'])
        return lst


    def getAssociatedCaseColor(self, study, case_idx, data_idx, idx):
        """
        Return plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        color = plot['color']
        if not color:
            color = self.defaultInitialValues()['color']
            self.setAssociatedCaseColor(study, case_idx, data_idx, idx, color)
        return color


    def setAssociatedCaseColor(self, study, case_idx, data_idx, idx, color):
        """
        Set plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['color'] = color


    def getAssociatedCaseFormat(self, study, case_idx, data_idx, idx):
        """
        Return plot format
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        fmt = plot['fmt']
        if not fmt:
            fmt = ""
        return fmt


    def setAssociatedCaseFormat(self, study, case_idx, data_idx, idx, fmt):
        """
        Set plot format
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['fmt'] = fmt


    def getAssociatedCaseLegend(self, study, case_idx, data_idx, idx):
        """
        Return plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        legend = plot['legend']
        if not legend:
            legend = ""
        return legend


    def setAssociatedCaseLegend(self, study, case_idx, data_idx, idx, legend):
        """
        Set plot color
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['legend'] = legend


    def getAssociatedCaseXcol(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['xcol']
        if not val:
            val = ""
        return val


    def setAssociatedCaseXcol(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['xcol'] = val


    def getAssociatedCaseYcol(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['ycol']
        if not val:
            val = ""
        return val


    def setAssociatedCaseYcol(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['ycol'] = val


    def getAssociatedCaseWidth(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['linewidth']
        if not val:
            val = ""
        return val


    def setAssociatedCaseWidth(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['linewidth'] = val


    def getAssociatedCaseMarker(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['markersize']
        if not val:
            val = ""
        return val


    def setAssociatedCaseMarker(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['markersize'] = val


    def getAssociatedCaseXerr(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['xerr']
        if not val:
            val = ""
        return val


    def setAssociatedCaseXerr(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['xerr'] = val


    def getAssociatedCaseXerrp(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['xerrp']
        if not val:
            val = ""
        return val


    def setAssociatedCaseXerrp(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['xerrp'] = val


    def getAssociatedCaseYerr(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['yerr']
        if not val:
            val = ""
        return val


    def setAssociatedCaseYerr(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['yerr'] = val


    def getAssociatedCaseYerrp(self, study, case_idx, data_idx, idx):
        """
        Return plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        val = plot['yerrp']
        if not val:
            val = ""
        return val


    def setAssociatedCaseYerrp(self, study, case_idx, data_idx, idx, val):
        """
        Set plot val
        """
        study_node = self.case.xmlGetNode("study", label = study)
        lst = study_node.xmlGetNodeList("case")
        node = lst[case_idx]
        n = node.xmlGetNodeList("data")[data_idx]
        plot = n.xmlGetNode("plot", id = idx)
        plot['yerrp'] = val


    def getCasePlotId(self, study, case_idx, data_idx, idx):
        """
        Return subplot id list for a case
        """
        study_node = self.case.xmlGetNode("study", label = study)
        case = study_node.xmlGetNode("case", id = case_idx)
        node = case.xmlGetNodeList("data")[data_idx]
        n = node.xmlGetNodeList("plot")[idx]
        return n['id']


    def getCaseIdList(self, study, case_idx, data_idx, plot_id):
        """
        Return subplot id list for a case
        """
        study_node = self.case.xmlGetNode("study", label = study)
        case = study_node.xmlGetNode("case", id = case_idx)
        node = case.xmlGetNodeList("data")[data_idx]
        n = node.xmlGetNode("plot", id = plot_id)
        lst = n['fig']
        if not lst:
            lst = ""
        return lst


    def setCaseIdList(self, study, case_idx, data_idx, plot_id, idlist):
        """
        Set subplot id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        case = study_node.xmlGetNode("case", id = case_idx)
        node = case.xmlGetNodeList("data")[data_idx]
        plot = node.xmlGetNode("plot", id = plot_id)
        plot['fig'] = idlist

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
