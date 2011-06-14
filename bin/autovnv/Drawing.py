#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne Scripts, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
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

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, string, logging
from string import *

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

import matplotlib

#-------------------------------------------------------------------------------
# matplotlib config
#-------------------------------------------------------------------------------

# Backend: Agg creates PNG output using the high quality Anti-Grain Geometry
#          library
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['legend.fontsize'] = 'medium'
rcParams['axes.labelsize']  = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['savefig.dpi']     = 200
rcParams['figure.figsize'] = (4,4)
rcParams['font.family']     = 'sans-serif'
rcParams['text.usetex']     = True

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)
#log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
# Curve
#-------------------------------------------------------------------------------

class Plot(object):
    """
    Storage of data for a single curve.
    """
    def __init__ (self, node, parser, file):
        """
        Constructor of a curve.
        @type node: C{DOM Element} instance
        @param node: xml node <plot>
        @type file: C{String}
        @param file: label of the file which contains data of the curve
        """
        self.id       = 0
        self.xspan    = []
        self.yspan    = []
        self.ismesure = False
        self.errorbar = []

        # Open file of data

        self.f = open(file, 'r')

        # Read mandatory attributes

        self.id     = int(node.attributes["fig"].value)
        self.legend = node.attributes["legend"].value
        ycol        = int(node.attributes["ycol"].value)

        # Read optional attributes
        try:
            self.fmt = node.attributes["fmt"].value
        except:
            self.fmt = ""

        try:
            xcol = int(node.attributes["xcol"].value)
        except:
            xcol = None
        try:
            xplus = float(node.attributes["xplus"].value)
        except:
            xplus = 0
        try:
            xfois = float(node.attributes["xfois"].value)
        except:
            xfois = 1
        try:
            yplus = float(node.attributes["yplus"].value)
        except:
            yplus = 0
        try:
            yfois = float(node.attributes["yfois"].value)
        except:
            yfois = 1

        self.uploadData(xcol, ycol, xplus, xfois, yplus, yfois)

        try:
            errorbar = [int(s) for s in node.attributes["errorbar"].value.split()]
        except:
            errorbar = None

        self.uploadErrorBar(errorbar)

        self.f.close()

        self.cmd = parser.getPltCommands(node)


    def setMeasurement(self, bool):
        self.ismesure = bool
        if self.ismesure and self.fmt == "":
            sef.fmt = 'o'


    def measurement(self):
        return self.ismesure


    def uploadData(self, xcol, ycol, xplus, xfois, yplus, yfois):
        """
        Upload and parse listing
        """
        j = 0

        for line in self.f.readlines():
            if line[0] != '#' and line != "\n":
                j += 1
                if xcol:
                    self.xspan.append(float(split(line)[xcol-1])*xfois + xplus)
                else:
                    self.xspan.append(j)

                self.yspan.append(float(split(line)[ycol-1])*yfois + yplus)


    def uploadErrorBar(self, errorbar):
        """
        load and parse listing for Measurement incertitude
        """
        if errorbar == None:
            self.errorbar = None

        elif len(errorbar) == 2:
            error = [ [], [] ]
            for line in self.f.readlines():
                if line[0] != '#' and line != "\n":
                    error[0].append(float(split(line)[errorbar[0]-1]))
                    error[1].append(float(split(line)[errorbar[1]-1]))
            self.errorbar = error

        elif len(errorbar) == 1:
            error = []
            for line in self.f.readlines():
                if line[0] != '#' and line != "\n":
                    error.append(float(split(line)[errorbar[0]-1]))
            self.errorbar = error

#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

class Subplot(object):
    """
    Management of a single subplot (frame with several curves).
    """
    def __init__ (self, node, parser, curves):
        """
        Constructor of a plot.
        """
        self.title     = ""
        self.xlabel    = ""
        self.ylabel    = ""
        self.id        = 0
        self.legstatus = "No"
        self.legpos    = ""
        self.xrange    = ""
        self.yrange    = ""
        self.curves    = []

        self.id = int(node.attributes["fig"].value)

        try:
            self.xlabel = node.attributes["xlabel"].value
        except:
            self.xlabel = ""

        try:
            self.ylabel = node.attributes["ylabel"].value
        except:
            self.ylabel = ""

        try:
            self.title = node.attributes["title"].value
        except:
            self.title = ""

        try:
            self.legstatus = node.attributes["legstatus"].value
        except:
            self.legstatus = 'Yes'

        try:
            self.legpos = [float(s) for s in node.attributes["legpos"].value.split()]
        except:
            self.legpos = None

        try:
            self.xrange = [float(s) for s in node.attributes["xrange"].value.split()]
        except:
            self.xrange = None

        try:
            self.yrange = [float(s) for s in node.attributes["yrange"].value.split()]
        except:
            self.yrange = None

        for curve in curves:
            if curve.id == self.id:
                self.curves.append(curve)

        self.cmd = parser.getPltCommands(node)

        # build of legends
        self.legends = []

        for curve in self.curves:
            if curve.id == self.id:
                if curve.measurement():
                    if curve.errorbar != None:
                        self.legends.append(curve.legend + " errorbar")
                        self.legends.append(curve.legend + " errorbar")
                        self.legends.append(curve.legend)
                    else:
                        self.legends.append(curve.legend)
                else:
                    self.legends.append(curve.legend)

#-------------------------------------------------------------------------------
# Figure
#-------------------------------------------------------------------------------

class Figure(object):
    """
    Management of a single figure (layout of several plots).
    """
    def __init__ (self, node, parser, curves, n_plots, plots):
        """
        Constructor of a figure.
        """
        self.file_name = node.attributes["name"].value
        self.l_plots = [int(s) for s in node.attributes["fig"].value.split()]

        self.tags = ("title", "figsize", "dpi", "facecolor", "edgecolor", "linewidth")
        for tag in self.tags:
            try:
                self.__dict__[tag] = node.attributes[tag].value
            except:
                self.__dict__[tag] = None

        self.cmd = parser.getPltCommands(node)

        self.o_plots = []
        for i in self.l_plots:
            for p in plots:
                if p.id == i:
                    self.o_plots.append(p)

        # Figure options
        if self.figsize:
            t = tuple(int(s) for s in self.figsize[1:-1].split(','))
            plt.figure(figsize=t)
        if self.dpi:
            plt.figure(dpi=int(self.dpi))
        if self.linewidth:
            plt.figure(linewidth=float(self.linewidth))
        if self.facecolor:
            plt.figure(facecolor=self.facecolor)
        if self.edgecolor:
            plt.figure(edgecolor=self.edgecolor)
        if self.title:
            plt.suptitle(self.title, fontsize=10)


    def layout(self, n_plots):
        """
        Automatic layout, based on the number of subplots
        Parameters with their defaults:
        left = 0.125 the left side of the subplots of the ﬁgure
        right = 0.9 the right side of the subplots of the ﬁgure
        bottom = 0.1 the bottom of the subplots of the ﬁgure
        top = 0.9 the top of the subplots of the ﬁgure
        wspace = 0.2 the amount of width reserved for blank space between subplots
        hspace = 0.2 the amount of height reserved for white space between subplots
        """
        if n_plots < 3:
            nbrow = 1
            nbcol = n_plots
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 2
            ri = 0.8
            le = 0.15
            hs = 0.35
            ws = 0.5
            plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)
        elif n_plots > 24:
            rcParams['font.size'] = 4
            rcParams['lines.markersize']= 1
            nbrow = 5
            nbcol = (n_plots-1)/nbrow+1
            hs=0.5
            if nbcol > 7:
                hs = 1.0
            elif nbcol > 6:
                hs = 0.85
            elif nbcol > 5:
                hs = 0.65
            ri = 0.9
            le = 0.1
            ws = 1.25
            plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)
        elif n_plots > 12:
            nbrow = 4
            rcParams['font.size'] = 5
            rcParams['lines.markersize'] = 2
            nbcol = (n_plots-1) / nbrow+1
            hs=0.45
            if nbcol > 5:
                hs = 0.8
            elif nbcol > 4:
                hs = 0.6
            ri = 0.9
            le = 0.1
            ws = 0.99
            plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)
        elif n_plots > 6:
            nbrow = 3
            nbcol = (n_plots-1) / nbrow+1
            rcParams['lines.markersize'] = 3
            rcParams['font.size'] = 6
            hs = 0.4
            if nbcol > 3:
                hs = 0.5
            ri = 0.9
            le = 0.1
            ws = 0.8
            plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)
        elif n_plots > 2:
            nbrow = 2
            nbcol = (n_plots-1) / nbrow+1
            hs = 0.3
            if nbcol > 2:
              hs = 0.45
            ri = 0.85
            le = 0.1
            ws = 0.6
            plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 4

        return nbrow, nbcol, hs, ri, le, ws


    def save(self, f):
        """method used to save the figure in a png format"""
        plt.savefig(f)
        plt.close()

#-------------------------------------------------------------------------------
# Plotter
#-------------------------------------------------------------------------------

class Plotter(object):
    """
    Manager of the matplotlib commands.
    """
    def __init__ (self, parser):
        """
        Constructor of the plotter.
        @type parser: C{Parser} instance
        @param parser: parser of the xml file
        """
        self.parser = parser


    def plot_study(self, study_label, study_object):
        """
        Method used to plot all plots from a I{study_label} (all cases).
        @type study_label: C{String}
        @param study_label: label of a study
        """
        # initialisation for each Study
        self.n_plots = 0
        self.curves  = []
        self.plots   = []
        self.figures = []

        # Read the parser for the Measurments Files
        nodes_list, files = self.parser.getMeasurement(study_label)
        for i in range(len(files)):
            nodes = nodes_list[i]
            file  = files[i]
            for node in nodes:
                curve = Plot(node, self.parser, file)
                curve.setMeasurement(True)
                self.curves.append(curve)
                self.n_plots = max(self.n_plots, curve.id)

        # Read the files of results
        for case in study_object.Cases:
            if case.plot == "on" and case.is_run != "KO":
                for node in self.parser.getChilds(case.node, "data"):
                    plots, file, dest, repo = self.parser.getResult(node)
                    f = os.path.join(self.parser.getDestination(),
                                     study_label,
                                     case.label, "RESU",
                                     dest, file)
                    if not os.path.isfile(f):
                        raise ValueError, "This file does not exist: %s" % f
                    for node in plots:
                        curve = Plot(node, self.parser, f)
                        curve.setMeasurement(False)
                        self.curves.append(curve)
                        self.n_plots = max(self.n_plots, curve.id)

        for node in self.parser.getPlots(study_label):
            self.plots.append(Subplot(node, self.parser, self.curves))

        for p in self.plots:
            self.n_plots = max(self.n_plots, p.id)

        for node in self.parser.getFigures(study_label):
            self.figures.append(Figure(node, self.parser, self.curves, self.n_plots, self.plots))

        for figure in self.figures:
            self.plot_figure(figure)
            f = os.path.join(self.parser.getDestination(),
                             study_label,
                             "POST",
                             figure.file_name)

            # additional matplotlib raw commands for figure
            for cmd in figure.cmd:
                f = open("./tmp.py", "w")
                f.write(cmd)
                f.close()
                execfile("./tmp.py")
                os.remove("./tmp.py")

            figure.save(f)
            study_object.matplotlib_figures.append(f)


    def __draw_curve(self, curve, n_fig, j):
        xspan = curve.xspan
        yspan = curve.yspan

        if curve.measurement():
            if curve.errorbar:
                lines = plt.errorbar(xspan, yspan, yerr=curve.errorbar,
                                     hold=True, fmt=curve.fmt)
            else:
                lines = plt.plot(xspan, yspan, curve.fmt)
        else:
            if curve.fmt:
                lines = plt.plot(xspan, yspan, curve.fmt)
            else:
                if n_fig[curve.id -1] < 8:
                    lines = plt.plot(xspan, yspan)
                elif n_fig[curve.id -1] < 15:
                    lines = plt.plot(xspan, yspan, '--')
                else:
                    lines = plt.plot(xspan, yspan, ':')

            # additional matplotlib raw commands for line2D
            line = lines[0]
            for cmd in curve.cmd:
                f = open("./tmp.py", "w")
                f.write(cmd)
                f.close()
                execfile("./tmp.py")
                os.remove("./tmp.py")
        plt.hold(True)


    def __draw_legend(self, p, bool, hs, ws, ri, le):
        if p.xlabel:
            plt.xlabel(p.xlabel)
        if p.ylabel:
            plt.ylabel(p.ylabel)
        if p.title:
            plt.title(p.title)
        if p.xrange:
            plt.xlim((p.xrange[0], p.xrange[1]))
        if p.yrange:
            plt.ylim((p.yrange[0], p.yrange[1]))

        if p.legends:
            if p.legstatus == 'Yes':
                if p.legpos == None:
                    plt.legend(p.legends, bbox_to_anchor=(1.02, 1), \
                                loc=2, borderaxespad=0., markerscale=0.5)
                else:
                    id_loc = 2
                    if   p.legpos[0] > 0.9 and p.legpos[1] > 0.9:
                        id_loc = 1
                    elif p.legpos[0] > 0.9 and p.legpos[1] < 0.1:
                        id_loc = 4
                    elif p.legpos[0] < 0.1 and p.legpos[1] < 0.1:
                        id_loc = 3
                    elif p.legpos[0] < 0.1 and p.legpos[1] > 0.9:
                        id_loc = 2
                    elif p.legpos[1] > 0.9:
                        id_loc = 2
                    elif p.legpos[1] < 0.1:
                        id_loc = 3 
                    elif p.legpos[0] > 0.9:
                        id_loc = 4
                    elif p.legpos[0] < 0.1:
                        id_loc = 3
                    plt.legend(p.legends,
                               bbox_to_anchor=(p.legpos[0], p.legpos[1]),
                               loc=id_loc, borderaxespad=0., markerscale=0.5)
                    if bool:
                        ri2 = min(0.9, ri*1.1)
                        le2 = max(0.125, le/1.1)
                        ws2 = min(0.3, ws/1.8)
                        plt.subplots_adjust(hspace=hs, wspace=ws2,
                                            right=ri2, left=le2)


    def plot_figure(self, figure):
        """
        Plotter of a single figure with several subplots.
        """
        # Layout
        nbrow, nbcol, hs, ri, le, ws = figure.layout(self.n_plots)
        log.debug("plot_figure --> layout: n_plots: %s nbcol: %s nbrow: %s l_plots: %s" % \
                  (self.n_plots, nbcol, nbrow, figure.l_plots))

        # list of numbers of curves per plot
        n_fig = [0] * self.n_plots

        # single plot in the figure
        if len(figure.l_plots) == 1:
            log.debug("plot_figure --> single plot draw")

            # build of the plots
            for curve in self.curves:
                plt.subplot(nbcol, nbrow, curve.id)
                n_fig[curve.id-1] += 1
                self.__draw_curve(curve, n_fig, curve.id-1)

        # several plots in the figure
        else:
            # build of the plots
            for j in range(len(figure.l_plots)):
                log.debug("plot_figure --> several plot draw: id = %s" % j)
                plt.subplot(nbcol, nbrow, j+1)

                for curve in self.curves:
                    if figure.l_plots[j] == curve.id:
                        n_fig[j] += 1
                        self.__draw_curve(curve, n_fig, j)

        plt.hold(False)

        # single plot in the figure: legends layout
        if len(figure.l_plots) == 1:
            for i in range(self.n_plots):
                p = figure.o_plots[i]
                plt.subplot(nbcol, nbrow, i+1)
                log.debug("plot_figure --> single plot legends")
                self.__draw_legend(p, False, hs, ws, ri, le)

        else:
            # several plots in the figure: legends layout
            for i in range(len(figure.l_plots)):
                j = figure.l_plots[i]
                p = figure.o_plots[i]
                plt.subplot(nbcol, nbrow, j)
                log.debug("plot_figure --> several plot legends: id = %s" % j)
                self.__draw_legend(p, True, hs, ws, ri, le)

        # additional matplotlib raw commands for subplot
        for i in range(len(figure.l_plots)):
            j = figure.l_plots[i]
            #p = figure.o_plots[i]
            ax = plt.subplot(nbcol, nbrow, j)
            for cmd in p.cmd:
                f = open("./tmp.py", "w")
                f.write(cmd)
                f.close()
                execfile("./tmp.py")
                os.remove("./tmp.py")

#-------------------------------------------------------------------------------
