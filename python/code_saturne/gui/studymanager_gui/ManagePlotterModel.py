# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

import sys, os.path, unittest

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

        self.plot_keys = ('color', 'format', 'legend',
                          'xcol', 'ycol', 'xscale', 'yscale', 'xplus', 'yplus',
                          'xerr', 'xerrp', 'yerr', 'yerrp')


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
        default['ycol']              = "0"

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
        study_node = self.case.xmlGetNode("study", label=study)
        lst = []
        for node in study_node.xmlGetNodeList("subplot"):
            lst.append(node['id'])
        return lst


    def getFigureList(self, study):
        """
        Return list of figures for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        lst = []
        for idx, node in enumerate(study_node.xmlGetNodeList("figure")):
            lst.append(idx)
        return lst


    def getMeasurementList(self, study):
        """
        Return list of measurements for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        lst = []
        for node in study_node.xmlGetNodeList("measurement"):
            lst.append((node['file'], node['path']))
        return lst


    def getCaseList(self, study):
        """
        Return list of cases for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        lst = []
        for idx, nn in enumerate(study_node.xmlGetNodeList("case")):
            lst.append(idx)
        return lst


    def addSubplot(self, study):
        """
        add a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        id_max = -1
        for nn in study_node.xmlGetNodeList("subplot"):
            idx = int(nn['id'])
            if id_max < idx:
                id_max = idx
        id_max = id_max+1
        study_node.xmlInitChildNode("subplot", id=id_max)
        return id_max


    def delSubplot(self, study, idx):
        """
        suppress a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        node = study_node.xmlGetNode("subplot", id=idx)
        node.xmlRemoveNode()

        for nn in study_node.xmlGetNodeList("figure"):
            lst = nn['idlist']
            new_lst = ''
            if lst:
                for idl in lst.split():
                    if idl != idx:
                        new_lst += idl + " "
            nn['idlist'] = new_lst.strip()

        for nn in study_node.xmlGetNodeList("plot"):
            lst = nn['spids']
            new_lst = ''
            new_lst = ''
            if lst:
                for idl in lst.split():
                    if idl != idx:
                        new_lst += idl + " "
            nn['spids'] = new_lst.strip()


    def getSubplotIdByIdx(self, study, idx):
        """
        Return id for a subplot, based on tis position in list
        Returns id after set (which may be the same as before if
        the requested id is already used).
        """

        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNodeByIdx("subplot", idx)
        return subplot['id']


    def setSubplotId(self, study, old_idx, new_idx):
        """
        Set id for a subplot
        Returns id after set (which may be the same as before if
        the requested id is already used).
        """

        if new_idx == old_idx:
            return False

        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=old_idx)
        if new_idx != "" and new_idx not in self.getSubplotList(study):
            subplot['id'] = new_idx
        else:
            return False

        # Renumber references

        for nn in study_node.xmlGetNodeList("figure"):
            lst = nn['idlist']
            new_lst = ''
            if lst:
                for idl in lst.split():
                    if idl != old_idx:
                        new_lst += idl + " "
                    else:
                        new_lst += new_idx + " "
            nn['idlist'] = new_lst.strip()

        for nn in study_node.xmlGetNodeList("plot"):
            lst = nn['spids']
            new_lst = ''
            new_lst = ''
            if lst:
                for idl in lst.split():
                    if idl != old_idx:
                        new_lst += idl + " "
                    else:
                        new_lst += new_idx + " "
            nn['spids'] = new_lst

        return True


    def getSubplotTitle(self, study, idx):
        """
        Return title for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        name = subplot['title']
        if not name:
            name = self.defaultInitialValues()['subplot_name']
        return name


    def setSubplotTitle(self, study, idx, title):
        """
        Set title for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        if title != "":
            subplot['title'] = title


    def getSubplotXLabel(self, study, idx):
        """
        Return X label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        name = subplot['xlabel']
        if not name:
            name = ""
        return name


    def setSubplotXLabel(self, study, idx, title):
        """
        Set X label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        subplot['xlabel'] = title


    def getSubplotYLabel(self, study, idx):
        """
        Return Y label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        name = subplot['ylabel']
        if not name:
            name = ""
        return name


    def setSubplotYLabel(self, study, idx, title):
        """
        Set Y label for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        subplot['ylabel'] = title


    def getSubplotLegPos(self, study, idx):
        """
        Return legend position for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        name = subplot['legpos']
        if not name:
            name = ""
        return name


    def setSubplotLegPos(self, study, idx, title):
        """
        Set legend position for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        subplot['legpos'] = title


    def getSubplotLegStatus(self, study, idx):
        """
        Return legend status for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        status = subplot['legstatus']
        if not status:
            status = self.defaultInitialValues()['subplot_legstatus']
            self.setSubplotLegStatus(study, idx, status)
        return status


    def setSubplotLegStatus(self, study, idx, status):
        """
        Set legend status for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        subplot['legstatus'] = status


    def getSubplotXLim(self, study, idx):
        """
        Return X limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        status = subplot['xlim']
        if not status:
            status = ""
        return status


    def setSubplotXLim(self, study, idx, pos):
        """
        Set X limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        if pos != "":
            subplot['xlim'] = pos


    def getSubplotYLim(self, study, idx):
        """
        Return Y limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        status = subplot['ylim']
        if not status:
            status = ""
        return status


    def setSubplotYLim(self, study, idx, pos):
        """
        Set Y limit for a subplot
        """
        study_node = self.case.xmlGetNode("study", label=study)
        subplot = study_node.xmlGetNode("subplot", id=idx)
        if pos != "":
            subplot['ylim'] = pos


#-------------------------------------------------------------------------------
# Figure
#-------------------------------------------------------------------------------

    def addFigure(self, study):
        """
        add a fiure for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        idx = len(study_node.xmlGetNodeList("figure"))
        nn = study_node.xmlInitChildNode("figure", idx=-1)
        del(nn['idx'])
        for idx, n in enumerate(study_node.xmlGetNodeList("figure")):
            if nn == n:
                return idx

        return -1  # Should not reach here


    def delFigure(self, study, idx):
        """
        suppress a subplot for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        node = study_node.xmlGetNodeByIdx("figure", idx)
        node.xmlRemoveNode()


    def getFigureName(self, study, idx):
        """
        Return name for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        name = figure['name']
        if not name:
            name = self.defaultInitialValues()['figure_name']
            self.setFigureName(study, idx, name)
        return name


    def setFigureName(self, study, idx, name):
        """
        Set name for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        figure['name'] = name


    def getFigureTitle(self, study, idx):
        """
        Return title for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        name = figure['title']
        if not name:
            name = self.defaultInitialValues()['figure_title']
            self.setFigureTitle(study, idx, name)
        return name


    def setFigureTitle(self, study, idx, name):
        """
        Set title for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        if name != "":
            figure['title'] = name
        elif figure['title']:
            figure['title'] = name


    def getFigureRows(self, study, idx):
        """
        Return number of row for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        row = figure['nbrow']
        if not row:
            row = ""
        return row


    def setFigureRows(self, study, idx, row):
        """
        Set number of row for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        if row != "":
            figure['nbrow'] = row


    def getFigureColumns(self, study, idx):
        """
        Return number of column for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        col = figure['nbcol']
        if not col:
            col = ""
        return col


    def setFigureColumns(self, study, idx, col):
        """
        Set number of column for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        if col != "":
            figure['nbcol'] = col


    def getFigureSize(self, study, idx):
        """
        Return size for a figure
        """
        study_node = self.case.xmlGetNode("study", label = study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        col = figure['figsize']
        if not col:
            col = ""
        return col


    def setFigureSize(self, study, idx, col):
        """
        Set number of column for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        if col != "":
            figure['figsize'] = col


    def getFigureFormat(self, study, idx):
        """
        Return format for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        name = figure['format']
        if not name:
            name = self.defaultInitialValues()['figure_fmt']
            self.setFigureFormat(study, idx, name)
        return name


    def setFigureFormat(self, study, idx, fmt):
        """
        Set format for a figure
        """
        self.isInList(fmt, ("png", "pdf"))
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        figure['format'] = fmt


    def getFigureIdList(self, study, idx):
        """
        Return figure id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        lst = figure['idlist']
        if not lst:
            lst = ""
        return lst


    def setFigureIdList(self, study, idx, idlist):
        """
        Set figure id list for a figure
        """
        study_node = self.case.xmlGetNode("study", label=study)
        figure = study_node.xmlGetNodeByIdx("figure", idx)
        figure['idlist'] = idlist


#-------------------------------------------------------------------------------
# Measurement or Case Data
#-------------------------------------------------------------------------------

    def dataKeys(self):
        """
        Return list of keys used for plot definitions.
        """
        return ('color', 'format', 'legend', 'xcol', 'ycol', 'xscale', 'yscale',
                'xplus', 'yplus', 'xerr', 'xerrp', 'yerr', 'yerrp')


    def dataDictDefaults(self):
        """
        Return dictionnary default with key/value definitions.
        """

        defaults = self.defaultInitialValues()
        d = {}
        for k in self.plot_keys:
            d[k] = ''
            if k in defaults:
                d[k] = defaults[k]

        return d


    def addDataPlot(self, data):
        """
        add a plot for a data file
        """
        nn = data.xmlInitChildNode("plot", idx=-1)
        del(nn['idx'])

        lst = data.xmlGetNodeList("plot")
        for idx, n in enumerate(lst):
            if n == nn:
                return idx

        return -1 # should not arrive here


    def deleteDataPlot(self, data, idx):
        """
        suppress a plot for a measurement file
        """
        nn = data.xmlGetNodeByIdx("plot", idx)
        if nn:
            nn.xmlRemoveNode()


    def getDataPlotList(self, data):
        """
        """
        lst = []
        for idx, n in enumerate(data.xmlGetNodeList("plot")):
            lst.append(idx)
        return lst


    def getDataDict(self, data, idx):
        """
        Return dictionnary with key/value definitions.
        """
        plot = data.xmlGetNodeByIdx("plot", idx)

        d = {}
        for k in self.plot_keys:
            d[k] = plot[k]

        return d


    def setDataDict(self, data, idx, values):
        """
        Return dictionnary with key/value definitions.
        """
        plot = data.xmlGetNodeByIdx("plot", idx)

        d = {}
        for k in self.plot_keys:
            if k in values:
                plot[k] = values[k]
                if plot[k] == '':
                    del(plot[k])


    def getDataIdList(self, data, idx):
        """
        Return subplot id list for a plot of a data
        """
        plot = data.xmlGetNodeByIdx("plot", idx)
        lst = plot['spids']
        if not lst:
            lst = ""
        return lst


    def setDataIdList(self, data, idx, idlist):
        """
        Set subplot id list for a plot of a data
        """
        plot = data.xmlGetNodeByIdx("plot", idx)
        plot['spids'] = idlist


#-------------------------------------------------------------------------------
# Measurement
#-------------------------------------------------------------------------------

    def addMeasurementFile(self, study, fle):
        """
        add a measurement file for a study
        """
        study_node = self.case.xmlGetNode("study", label=study)
        study_node.xmlInitChildNode("measurement", file=fle)


    def delMeasurementFile(self, measurement):
        """
        suppress a measurement file for a study
        """
        measurement.xmlRemoveNode()


    def getMeasurementNode(self, study, name):
        """
        """
        study_node = self.case.xmlGetNode("study", label=study)
        mp, mf = os.path.split(name)
        m_nodes = study_node.xmlGetNodeList("measurement", file=mf, path=mp)
        l = len(m_nodes)
        if l > 0:
            measurement_node = m_nodes[0]
            for i in range(1, l):
                measurement_node.xmlMergeNode(m_nodes[i])
        else:
            measurement_node = study_node.xmlInitNode("measurement",
                                                      file=mf, path=mp)
        return measurement_node


    def getMeasurementPlotList(self, measurement):
        """
        """
        lst = []
        for idx, node in enumerate(measurement.xmlGetNodeList("plot")):
            lst.append(idx)
        return lst


#-------------------------------------------------------------------------------
# Case
#-------------------------------------------------------------------------------

    def getCaseName(self, study, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNodeByIdx("case", idx)
        name = nn['label']
        run_id = nn['run_id']
        if run_id is not None and run_id != '':
            name += '/' + run_id
        return name


    def getCaseNode(self, study, case):
        """
        """
        study_node = self.case.xmlGetNode("study", label=study)

        case_name, run_id = os.path.split(case)

        case_node = None
        if case_name != '':
            case_node = study_node.xmlGetNode("case",
                                              label=case_name,
                                              run_id=run_id)
        else:
            case_name = run_id
            case_node = study_node.xmlGetNode("case",
                                              label=case_name,
                                              run_id='')
            if case_node is None:
                case_node = study_node.xmlGetNode("case",
                                                  label=case_name)

        return case_node


    def getCaseDataNode(self, study, case, name):
        """
        """
        case_node = self.getCaseNode(study, case)

        d_nodes = case_node.xmlGetNodeList("data", file=name)
        l = len(d_nodes)
        if l > 0:
            data_node = d_nodes[0]
            for i in range(1, l):
                data_node.xmlMergeNode(d_nodes[i])
        else:
            data_node = case_node.xmlInitNode("data",
                                              file=name)
        return data_node


    def getCaseDataList(self, study, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNodeByIdx("case", idx)
        lst = []
        for n in nn.xmlGetNodeList("data"):
            lst.append(n['file'])
        return lst


    def addCaseDataFile(self, case_node, name):
        """
        """
        nn = case_node.xmlInitChildNode("data", dest="", file=name, idx=-1)
        del(nn['idx'])


    def delCaseDataFile(self, study, case_idx, idx):
        """
        """
        study_node = self.case.xmlGetNode("study", label = study)
        nn = study_node.xmlGetNodeByIdx("case", case_idx)
        lst = nn.xmlGetNodeList("data")
        node = lst[idx]
        node.xmlRemoveNode()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
