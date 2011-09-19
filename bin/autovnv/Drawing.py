# -*- coding: utf-8 -*-
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
rcParams['figure.dpi']      = 200
rcParams['figure.figsize']  = (4,4)
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
# Curve or plot
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
        self.subplot  = 0
        self.xspan    = []
        self.yspan    = []
        self.ismesure = False
        self.cmd      = []

        # Open file of data

        self.f = open(file, 'r')

        # Read mandatory attributes
        self.subplot = int(parser.getAttribute(node, "fig"))
        ycol         = int(parser.getAttribute(node, "ycol"))

        # Read optional attributes
        try:
            self.legend = parser.getAttribute(node, "legend")
        except:
            self.legend = ""

        try:
            self.fmt = parser.getAttribute(node, "fmt")
        except:
            self.fmt = ""

        try:
            xcol = int(parser.getAttribute(node, "xcol"))
        except:
            xcol = None

        try:
            xplus = float(parser.getAttribute(node, "xplus"))
        except:
            xplus = 0

        try:
            xfois = float(parser.getAttribute(node, "xfois"))
        except:
            xfois = 1

        try:
            yplus = float(parser.getAttribute(node, "yplus"))
        except:
            yplus = 0

        try:
            yfois = float(parser.getAttribute(node, "yfois"))
        except:
            yfois = 1

        self.uploadData(xcol, ycol, xplus, xfois, yplus, yfois)

        try:
            xerr = [int(s) for s in parser.getAttribute(node, "xerr").split()]
        except:
            xerr = None

        try:
            yerr = [int(s) for s in parser.getAttribute(node, "yerr").split()]
        except:
            yerr = None

        self.f.close()

        # Error Bar
        self.f = open(file, 'r')
        self.xerr = self.uploadErrorBar(xerr)
        self.yerr = self.uploadErrorBar(yerr)
        self.f.close()

        # List of additional matplotlib commands
        for k, v in parser.getAttributes(node).items():
            if k not in ('fig', 'fmt', 'legend', 'xcol', 'ycol', \
                         'xplus', 'yplus', 'xfois', 'yfois', 'xerr', 'yerr'):
                self.cmd.append("plt.setp(lines, " + k + "=" + v + ")")

        self.cmd += parser.getPltCommands(node)


    def setMeasurement(self, bool):
        self.ismesure = bool
        if self.ismesure and self.fmt == "":
            self.fmt = 'ko'


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
            return None

        elif len(errorbar) == 2:
            error = [ [], [] ]
            for line in self.f.readlines():
                if line[0] != '#' and line != "\n":
                    error[0].append(float(split(line)[errorbar[0]-1]))
                    error[1].append(float(split(line)[errorbar[1]-1]))
            return error

        elif len(errorbar) == 1:
            error = []
            for line in self.f.readlines():
                if line[0] != '#' and line != "\n":
                    error.append(float(split(line)[errorbar[0]-1]))
            return error

#-------------------------------------------------------------------------------
# SubPlot
#-------------------------------------------------------------------------------

class Subplot(object):
    """
    Management of a single subplot (frame with several curves).
    """
    def __init__ (self, node, parser, curves):
        """
        Constructor of a plot.
        """
        # Read mandatory attribute
        self.id = int(parser.getAttribute(node, "id"))

        # Read optional attributes
        self.xlabel    = parser.getAttribute(node, "xlabel", False)
        self.ylabel    = parser.getAttribute(node, "ylabel", False)
        self.title     = parser.getAttribute(node, "title", False)
        self.legstatus = parser.getAttribute(node, "legstatus", False)

        try:
            self.legpos = [float(s) for s in parser.getAttribute(node, "legpos").split()]
        except:
            self.legpos = None

        try:
            self.xlim = [float(s) for s in parser.getAttribute(node, "xlim").split()]
        except:
            self.xlim = None

        try:
            self.ylim = [float(s) for s in parser.getAttribute(node, "ylim").split()]
        except:
            self.ylim = None

        # Store the list of curves to be plotted in this subplot
        self.curves = []
        for curve in curves:
            if curve.subplot == self.id:
                self.curves.append(curve)

        # List of additional matplotlib commands
        self.cmd = []
        for k, v in parser.getAttributes(node).items():
            if k not in ('id', 'xlabel', 'ylabel', 'title', 'legstatus', \
                         'legpos', 'xlim', 'ylim'):
                self.cmd.append("plt.subplot(" + k + "=" + v + ")")

        self.cmd += parser.getPltCommands(node)

#-------------------------------------------------------------------------------
# Figure
#-------------------------------------------------------------------------------

class Figure(object):
    """
    Management of a single figure (layout of several subplots).
    """
    def __init__ (self, node, parser, curves, n_plots, subplots):
        """
        Constructor of a figure.
        """
        self.file_name = node.attributes["name"].value
        self.l_subplots = [int(s) for s in node.attributes["idlist"].value.split()]

        self.tags = ("title", "nbrow", "nbcol")
        for tag in self.tags:
            self.__dict__[tag] = parser.getAttribute(node, tag, False)

        self.cmd = []
        for k, v in parser.getAttributes(node).items():
            if k not in ("name", "idlist", "title", "nbrow", "nbcol"):
                self.cmd.append("plt.figure(" + k + "=" + v + ")")

        self.cmd += parser.getPltCommands(node)

        self.o_subplots = []
        for i in self.l_subplots:
            for p in subplots:
                if p.id == i:
                    self.o_subplots.append(p)


    def options(self):
        """
        Additional figure options.
        """
        if self.title:
            plt.suptitle(self.title, fontsize=10)


    def layout(self):
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
        nbr = len(self.l_subplots)

        if nbr < 3:
            nbrow = 1
            nbcol = nbr
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 2
            ri = 0.8
            le = 0.15
            hs = 0.35
            ws = 0.5
        elif nbr > 24:
            rcParams['font.size'] = 4
            rcParams['lines.markersize']= 1
            nbrow = 5
            nbcol = (nbr-1)/nbrow+1
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
        elif nbr > 12:
            nbrow = 4
            rcParams['font.size'] = 5
            rcParams['lines.markersize'] = 2
            nbcol = (nbr-1) / nbrow+1
            hs=0.45
            if nbcol > 5:
                hs = 0.8
            elif nbcol > 4:
                hs = 0.6
            ri = 0.9
            le = 0.1
            ws = 0.99
        elif nbr > 6:
            nbrow = 3
            nbcol = (nbr-1) / nbrow+1
            rcParams['lines.markersize'] = 3
            rcParams['font.size'] = 6
            hs = 0.4
            if nbcol > 3:
                hs = 0.5
            ri = 0.9
            le = 0.1
            ws = 0.8
        elif nbr > 2:
            nbrow = 2
            nbcol = (nbr-1) / nbrow+1
            hs = 0.3
            if nbcol > 2:
              hs = 0.45
            ri = 0.85
            le = 0.1
            ws = 0.6
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 4

        plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le)

        if self.nbrow:
            nbrow = int(self.nbrow)
        if self.nbcol:
            nbcol = int(self.nbcol)

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
                self.n_plots = max(self.n_plots, curve.subplot)

        # Read the files of results
        for case in study_object.Cases:
            if case.plot == "on" and case.is_run != "KO":
                for node in self.parser.getChilds(case.node, "data"):
                    plots, file, dest, repo = self.parser.getResult(node)

                    if dest:
                        d = dest
                    elif repo:
                        d = repo

                    f = os.path.join(self.parser.getDestination(),
                                     study_label,
                                     case.label, "RESU",
                                     d, file)

                    if not os.path.isfile(f):
                        f = os.path.join(self.parser.getDestination(),
                                         study_label,
                                         case.label, "RESU",
                                         d, "monitoring", file)
                        if not os.path.isfile(f):
                            raise ValueError, "This file does not exist: %s" % f

                    for node in plots:
                        curve = Plot(node, self.parser, f)
                        curve.setMeasurement(False)
                        self.curves.append(curve)
                        self.n_plots = max(self.n_plots, curve.subplot)

        subplots = []
        for node in self.parser.getSubplots(study_label):
            subplots.append(Subplot(node, self.parser, self.curves))

        for p in subplots:
            self.n_plots = max(self.n_plots, p.id)

        for node in self.parser.getFigures(study_label):
            self.figures.append(Figure(node, self.parser, self.curves, self.n_plots, subplots))

        for figure in self.figures:
            figure.options()

            # additional matplotlib raw commands for figure
            for cmd in figure.cmd:
                c = open("./tmp.py", "w")
                c.write(cmd)
                c.close()
                try:
                    execfile("./tmp.py")
                except:
                    print "Error with the matplotlib command: %s" % cmd
                os.remove("./tmp.py")

            # Plot curve
            self.plot_figure(figure)

            # save the figure
            f = os.path.join(self.parser.getDestination(),
                             study_label,
                             "POST",
                             figure.file_name)

            figure.save(f)

            # store the name of the figure for the build of
            # the detailed report
            study_object.matplotlib_figures.append(f)


    def __draw_curve(self, ax, curve, n_fig):
        """
        Draw a single curve.
        """
        xspan = curve.xspan
        yspan = curve.yspan

        if curve.measurement():
            if curve.xerr:
                lines = ax.errorbar(xspan, yspan,
                                    xerr=curve.xerr,
                                    fmt=curve.fmt,
                                    label=curve.legend)
            if curve.yerr:
                lines = ax.errorbar(xspan, yspan,
                                    yerr=curve.yerr,
                                    fmt=curve.fmt,
                                    label=curve.legend)
            else:
                lines = ax.plot(xspan, yspan, curve.fmt, label=curve.legend)
        else:
            if curve.fmt:
                lines = ax.plot(xspan, yspan, curve.fmt, label=curve.legend)
            else:
                if n_fig[curve.subplot -1] < 8:
                    lines = ax.plot(xspan, yspan, label=curve.legend)
                elif n_fig[curve.subplot -1] < 15:
                    lines = ax.plot(xspan, yspan, '--', label=curve.legend)
                else:
                    lines = ax.plot(xspan, yspan, ':', label=curve.legend)

        # additional matplotlib raw commands for line2D
        line = lines[0]
        for cmd in curve.cmd:
            c = open("./tmp.py", "w")
            c.write(cmd)
            c.close()
            try:
                execfile("./tmp.py")
            except:
                print "Error with the matplotlib command: %s" % cmd
            os.remove("./tmp.py")

        plt.hold(True)


    def __draw_axis(self, ax, p):
        if p.xlabel:
            ax.set_xlabel(p.xlabel)
        if p.ylabel:
            ax.set_ylabel(p.ylabel)
        if p.xlim:
            ax.set_xlim((p.xlim[0], p.xlim[1]))
        if p.ylim:
            ax.set_ylim((p.ylim[0], p.ylim[1]))


    def __draw_legend(self, ax, p, bool, hs, ws, ri, le):
        if p.title:
            ax.set_title(p.title)

        if p.legstatus == "on":
            handles, labels = ax.get_legend_handles_labels()

            if p.legpos == None:
                ax.legend(handles,
                          labels,
                          bbox_to_anchor=(1.02, 1),
                          loc=2,
                          borderaxespad=0.,
                          markerscale=0.5)
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
                ax.legend(handles,
                          labels,
                          bbox_to_anchor=(p.legpos[0], p.legpos[1]),
                          loc=id_loc,
                          borderaxespad=0.,
                          markerscale=0.5)

                if bool:
                    ri2 = min(0.9, ri*1.1)
                    le2 = max(0.125, le/1.1)
                    ws2 = min(0.3, ws/1.8)
                    plt.subplots_adjust(hspace=hs, wspace=ws2, right=ri2, left=le2)


    def plot_figure(self, figure):
        """
        Plotter of a single figure with several subplots.
        """
        # Layout
        nbrow, nbcol, hs, ri, le, ws = figure.layout()
        log.debug("plot_figure --> layout: n_plots: %s nbcol: %s nbrow: %s l_subplots: %s" % \
                  (self.n_plots, nbcol, nbrow, figure.l_subplots))

        # list of numbers of curves per subplot
        n_fig = [0] * self.n_plots

        # draw curves in the right subplot
        for j in figure.l_subplots:
            idx = figure.l_subplots.index(j)
            log.debug("plot_figure --> plot draw: id = %s" % j)
            ax = plt.subplot(nbrow, nbcol, idx + 1)

            for curve in self.curves:
                if j == curve.subplot:
                    n_fig[j-1] += 1
                    self.__draw_curve(ax, curve, n_fig)

        plt.hold(False)

        # axis and legend
        bool = len(figure.l_subplots) > 1

        for i in range(len(figure.l_subplots)):
            j = figure.l_subplots[i]
            p = figure.o_subplots[i]
            ax = plt.subplot(nbrow, nbcol, i + 1)
            self.__draw_axis(ax, p)
            self.__draw_legend(ax, p, bool, hs, ws, ri, le)

        # additional matplotlib raw commands for subplot
        for i in range(len(figure.l_subplots)):
            j = figure.l_subplots[i]
            p = figure.o_subplots[i]
            ax = plt.subplot(nbrow, nbcol, i + 1)
            for cmd in p.cmd:
                c = open("./tmp.py", "w")
                c.write(cmd)
                c.close()
                try:
                    execfile("./tmp.py")
                except:
                    print "Error with the matplotlib command: %s" % cmd
                os.remove("./tmp.py")

#-------------------------------------------------------------------------------
