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

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, string, logging
from string import *

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

try:
    import matplotlib
except ImportError:
    print("Warning: import matplotlib failed.")
    pass

matplotlib.use("Agg")

#-------------------------------------------------------------------------------
# matplotlib config
#-------------------------------------------------------------------------------

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
# additional colors
#-------------------------------------------------------------------------------

bedf      = (0.    , 0.3569, 0.7333)
dbedf     = (0.0353, 0.2078, 0.4784)
oedf      = (1.    , 0.6275, 0.1843)
doedf     = (0.9961, 0.3451, 0.0824)
yedf      = (0.7686, 0.8392, 0.    )
gedf      = (0.3137, 0.6196, 0.1843)
dgedf     = (0.    , 0.3922, 0.    )
violet    = (0.5803, 0.    , 0.8275)
turquoise = (0.    , 0.8078, 0.8196)
brown     = (0.5451, 0.2706, 0.0745)
red       = (1.    , 0.    , 0.    )

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)
#log.setLevel(logging.DEBUG)

#===============================================================================
# Plot class
#===============================================================================

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
        self.subplots  = []
        self.xspan    = []
        self.yspan    = []
        self.ismesure = False
        self.cmd      = []

        # Open file of data

        self.f = open(file, 'r')

        # Read mandatory attributes
        self.subplots = [int(s) for s in parser.getAttribute(node,"spids").split()]
        ycol         = int(parser.getAttribute(node, "ycol"))

        # Read optional attributes
        self.legend = parser.getAttribute(node, "legend", "")

        self.fmt = parser.getAttribute(node, "fmt", "")
        xcol  = int(parser.getAttribute(node,   "xcol",  0))
        xplus = float(parser.getAttribute(node, "xplus", 0))
        xscale = float(parser.getAttribute(node, "xscale", 1))
        yplus = float(parser.getAttribute(node, "yplus", 0))
        yscale = float(parser.getAttribute(node, "yscale", 1))

        self.uploadData(xcol, ycol, xplus, xscale, yplus, yscale)

        try:
            xerr = [int(s) for s in parser.getAttribute(node, "xerr").split()]
        except:
            xerr = None

        try:
            yerr = [int(s) for s in parser.getAttribute(node, "yerr").split()]
        except:
            yerr = None

        try:
            xerrp = float(parser.getAttribute(node, "xerrp"))
        except:
            xerrp = None

        try:
            yerrp = float(parser.getAttribute(node, "yerrp"))
        except:
            yerrp = None

        self.f.close()

        # Error Bar
        self.f = open(file, 'r')
        self.xerr = self.uploadErrorBar(xerr, xerrp, xcol)
        self.f.close()

        self.f = open(file, 'r')
        self.yerr = self.uploadErrorBar(yerr, yerrp, ycol)
        self.f.close()

        # List of additional matplotlib commands
        for k, v in parser.getAttributes(node).items():
            if k not in ('spids', 'fmt', 'legend', 'xcol', 'ycol', \
                         'xplus', 'yplus', 'xscale', 'yscale', \
                         'xerr', 'yerr', 'xerrp', 'yerrp', 'id'):
                self.cmd.append("plt.setp(lines, " + k + "=" + v + ")")

        self.cmd += parser.getPltCommands(node)

    #---------------------------------------------------------------------------

    def setMeasurement(self, bool):
        self.ismesure = bool
        if self.ismesure and self.fmt == "":
            self.fmt = 'ko'

    #---------------------------------------------------------------------------

    def measurement(self):
        return self.ismesure

    #---------------------------------------------------------------------------

    def uploadData(self, xcol, ycol, xplus, xscale, yplus, yscale):
        """
        Upload and parse data
        """
        j = 0

        for line in self.f.readlines():
            line = line.lstrip()
            if line and line[0] != '#':
                j += 1
                line = line.replace(", ", " ") # compatibility with CSV
                line = line.lstrip()

                # for CSV files, try to detect a header to skip it
                if j == 1:
                    try:
                        val = float(line.split()[0])
                    # if it can not be converted to float
                    except ValueError as not_float:
                        continue

                if xcol:
                    self.xspan.append(float(line.split()[xcol-1])*xscale + xplus)
                else:
                    self.xspan.append(j)

                self.yspan.append(float(line.split()[ycol-1])*yscale + yplus)

    #---------------------------------------------------------------------------

    def uploadErrorBar(self, errorbar, errorp, col):
        """
        load and parse data for Measurement uncertainty
        """
        if errorbar == None and errorp == None:
            return None

        elif errorbar:
            if errorp:
                print("Warning: ambiguous definitions of error bars for "
                      "one data set, percentage and set of error values.\n"
                      "The error definition by percentage will be ignored.")

            if len(errorbar) == 2:
                error = [ [], [] ]
                j = 0
                for line in self.f.readlines():
                    line = line.lstrip()
                    if line and line[0] != '#':
                        j += 1
                        line = line.replace(", ", " ") # compatibility with CSV
                        line = line.lstrip()
                        # for CSV files, try to detect a header to skip it
                        if j == 1:
                            try:
                                val = float(line.split()[0])
                            except ValueError as not_float:
                                continue
                        error[0].append(float(line.split()[errorbar[0]-1]))
                        error[1].append(float(line.split()[errorbar[1]-1]))
                return error

            elif len(errorbar) == 1:
                error = []
                j = 0
                for line in self.f.readlines():
                    line = line.lstrip()
                    if line and line[0] != '#':
                        j += 1
                        line = line.replace(", ", " ") # compatibility with CSV
                        line = line.lstrip()
                        # for CSV files, try to detect a header to skip it
                        if j == 1:
                            try:
                                val = float(line.split()[0])
                            except ValueError as not_float:
                                continue
                        error.append(float(line.split()[errorbar[0]-1]))
                return error
        elif errorp:
            if col == 0:
                print("Error: can not compute errors by percentage of an "
                      "unspecified data set (column number missing).\n")
                sys.exit(1)
            else:
                error = []
                j = 0
                for line in self.f.readlines():
                    line = line.lstrip()
                    if line and line[0] != '#':
                        j += 1
                        line = line.replace(", ", " ") # compatibility with CSV
                        line = line.lstrip()
                        # for CSV files, try to detect a header to skip it
                        if j == 1:
                            try:
                                val = float(line.split()[0])
                            except ValueError as not_float:
                                continue
                        error.append(errorp/100.*float(line.split()[col-1]))
                return error

#===============================================================================
# Probes class
# Curve or plot for monitoring point (probe)
#===============================================================================

class Probes(object):
    """
    Curve from a probe.
    """
    def __init__ (self, file_name, fig, ycol):
        """
        Constructor of a curve.
        """
        self.subplots = [int(fig)]
        self.xspan    = []
        self.yspan    = []
        self.cmd      = []
        self.legend   = "Probe " + str(ycol - 1)
        self.fmt      = ""
        self.ismesure = False
        self.xerr     = None
        self.yerr     = None

        xcol = 1

        f = open(file_name, 'r')

        for line in f.readlines():
            line = line.lstrip()
            if line and line[0] != '#':
                self.xspan.append(float(line.split()[xcol - 1]))
                self.yspan.append(float(line.split()[ycol - 1]))

        f.close()

    #---------------------------------------------------------------------------

    def measurement(self):
        return self.ismesure

#===============================================================================
# Subplot class
#===============================================================================

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
        self.xlabel    = parser.getAttribute(node, "xlabel", "")
        self.ylabel    = parser.getAttribute(node, "ylabel", "")
        self.title     = parser.getAttribute(node, "title",  "")
        self.legstatus = parser.getAttribute(node, "legstatus", "off")

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
            for subplot in curve.subplots:
                if subplot == self.id:
                    self.curves.append(curve)

        # List of additional matplotlib commands
        self.cmd = []
        for k, v in parser.getAttributes(node).items():
            if k not in ('id', 'xlabel', 'ylabel', 'title', 'legstatus', \
                         'legpos', 'xlim', 'ylim'):
                self.cmd.append("plt.subplot(" + k + "=" + v + ")")

        self.cmd += parser.getPltCommands(node)

#===============================================================================
# Figure class
#===============================================================================

class Figure(object):
    """
    Management of a single figure (layout of several subplots).
    """
    def __init__ (self, node, parser, curves, subplots, default_fmt):
        """
        Constructor of a figure.
        """
        self.file_name = node.attributes["name"].value
        self.fmt = parser.getAttribute(node, "format", default_fmt)

        for tag in ("title", "nbrow", "nbcol"):
            self.__dict__[tag] = parser.getAttribute(node, tag, False)

        self.figsize = None
        self.dpi = None

        for k, v in parser.getAttributes(node).items():
            if k == "figsize":
                v_spl = v.strip("() ").split(",")
                self.figsize = tuple([float(co) for co in v_spl])
            elif k == "dpi":
                self.dpi = int(v)

        # Store te list of subplot objects associated to the current figure
        self.subplots = []
        for id in [int(s) for s in node.attributes["idlist"].value.split()]:
            for p in subplots:
                if p.id == id:
                    self.subplots.append(p)

    #---------------------------------------------------------------------------

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
        nbr = len(self.subplots)

        if nbr < 3:
            nbrow = 1
            nbcol = nbr
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 2
            ri = 0.9
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
            le = 0.15
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
            le = 0.15
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
            le = 0.15
            ws = 0.8
        elif nbr > 2:
            nbrow = 2
            nbcol = (nbr-1) / nbrow+1
            hs = 0.3
            if nbcol > 2:
              hs = 0.45
            ri = 0.9
            le = 0.15
            ws = 0.6
            rcParams['font.size'] = 6
            rcParams['lines.markersize']= 4

        if self.nbrow:
            nbrow = int(self.nbrow)
        if self.nbcol:
            nbcol = int(self.nbcol)

        return nbrow, nbcol, hs, ri, le, ws

#===============================================================================
# Plotter class
#===============================================================================

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

        # initialisation for each Study
        self.curves  = []
        self.figures = []

    #---------------------------------------------------------------------------

    def __number_of_column(self, file_name):
        """
        Compute the number of column of the data file.
        """
        nbr = 0
        f = open(file_name, 'r')
        for line in f.readlines():
            line = line.lstrip()
            if line and line[0] != '#':
                nbr = len(line.split())
                break
        f.close()
        return nbr

    #---------------------------------------------------------------------------

    def plot_study(self, study_label, study_object, disable_tex, default_fmt):
        """
        Method used to plot all plots from a I{study_label} (all cases).
        @type study_label: C{String}
        @param study_label: label of a study

        """
        # disable tex in Matplotlib (use Mathtext instead)
        rcParams['text.usetex'] = not disable_tex

        # Read the parser for the Measurements Files
        nodes_list, files = self.parser.getMeasurement(study_label)
        for i in range(len(files)):
            nodes = nodes_list[i]
            file  = files[i]
            for node in nodes:
                curve = Plot(node, self.parser, file)
                curve.setMeasurement(True)
                self.curves.append(curve)

        # Read the files for probes
        for case in study_object.cases:
            if case.plot == "on" and case.is_run != "KO":
                for node in self.parser.getChildren(case.node, "probes"):
                    file_name, dest, fig = self.parser.getProbes(node)

                    f = os.path.join(self.parser.getDestination(),
                                     study_label,
                                     case.label, case.resu,
                                     dest, "monitoring", file_name)
                    if not os.path.isfile(f):
                        raise ValueError("\n\nThis file does not exist: %s\n (call with path: %s)\n" % (file_name, f))

                    for ycol in range(2, self.__number_of_column(f) + 1):
                        curve = Probes(f, fig, ycol)
                        self.curves.append(curve)

        # Read the files of results of cases
        for case in study_object.cases:
            if case.plot == "on" and case.is_run != "KO":
                for node in self.parser.getChildren(case.node, "data"):
                    plots, file_name, dest, repo = self.parser.getResult(node)

                    if dest:
                        d = dest
                        dd = self.parser.getDestination()
                    elif repo:
                        d = repo
                        dd = self.parser.getRepository()

                    if case.subdomains:
                        dom_list = case.subdomains
                    else:
                        dom_list = [""]

                    iok = False
                    for ddd in dom_list:

                        for sd in (".", "monitoring", "profiles"):

                            f = os.path.join(dd,
                                             study_label,
                                             case.label, case.resu,
                                             d, ddd, sd, file_name)

                            if os.path.isfile(f):
                                iok = True
                                break

                    if not iok:
                        raise ValueError("\n\nThis file does not exist: %s\n (last call with path: %s)\n" % (file_name, f))

                    for nn in plots:
                        curve = Plot(nn, self.parser, f)
                        curve.setMeasurement(False)
                        self.curves.append(curve)

        # Read the files of results of postpro
        script, label, nodes, args = self.parser.getPostPro(study_label)
        for i in range(len(label)):
            if script[i]:
                for node in self.parser.getChildren(nodes[i], "data"):
                    plots, file_name, dest, repo = self.parser.getResult(node)

                    dd = self.parser.getDestination()
                    f = os.path.join(dd, study_label, "POST", file_name)

                    if not os.path.isfile(f):
                        raise ValueError("\n\nThis file does not exist: %s\n (call with path: %s)\n" % (file_name, f))

                    for nn in plots:
                        curve = Plot(nn, self.parser, f)
                        curve.setMeasurement(False)
                        self.curves.append(curve)

        subplots = []
        for node in self.parser.getSubplots(study_label):
            subplots.append(Subplot(node, self.parser, self.curves))

        # Build the list of figures to handle
        for node in self.parser.getFigures(study_label):
            self.figures.append(Figure(node,
                                       self.parser, self.curves,
                                       subplots,
                                       default_fmt))

        # create one figure and use it for all figures
        # this figure is the current figure
        fig = plt.figure()

        for figure in self.figures:
            # Plot curve
            self.plot_figure(figure)

            f = os.path.join(self.parser.getDestination(),
                             study_label,
                             "POST",
                             figure.file_name)

            # store the name of the figure for the build of
            # the detailed report without the png or pdf extension.
            study_object.matplotlib_figures.append(f)

            # save the figure
            self.__save(f, figure)

        # close current figure
        plt.close()

    #---------------------------------------------------------------------------

    def __draw_curve(self, ax, curve, p):
        """
        Draw a single curve.
        """

        # data sets x and y
        xspan = curve.xspan
        yspan = curve.yspan

        # data sets errors on x and y
        xerr = curve.xerr
        yerr = curve.yerr

        # draw curve with error bars
        if xerr or yerr:
            lines = ax.errorbar(xspan, yspan,
                                xerr=xerr,
                                yerr=yerr,
                                fmt=curve.fmt,
                                label=curve.legend)
        # draw curve only
        else:
            lines = ax.plot(xspan, yspan, curve.fmt, label=curve.legend)


        # additional matplotlib raw commands for line2D
        line = lines[0]
        for cmd in curve.cmd:
            try:
                exec(cmd)
            except:
                print("Error with the matplotlib command: %s" % cmd)

    #---------------------------------------------------------------------------

    def __draw_axis(self, ax, p):
        if p.xlabel:
            ax.set_xlabel(p.xlabel)
        if p.ylabel:
            ax.set_ylabel(p.ylabel)
        if p.xlim:
            ax.set_xlim((p.xlim[0], p.xlim[1]))
        if p.ylim:
            ax.set_ylim((p.ylim[0], p.ylim[1]))

    #---------------------------------------------------------------------------

    def __draw_legend(self, ax, p, bool, hs, ws, ri, le):
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
                    plt.subplots_adjust(hspace=hs,
                                        wspace=ws2,
                                        right=ri2,
                                        left=le2)

    #---------------------------------------------------------------------------

    def plot_figure(self, figure):
        """
        Plotter of a single figure with several subplots.
        """

        # clears the entire current figure
        plt.clf()

        # get current figure and set size
        if figure.figsize:
            fig = plt.gcf()
            fig.set_size_inches(figure.figsize)

        # Layout
        nbrow, nbcol, hs, ri, le, ws = figure.layout()
        log.debug("plot_figure --> layout: nbcol: %s nbrow: %s " % (nbcol, nbrow))
        plt.subplots_adjust(hspace=hs, wspace=ws, right=ri, left=le, bottom=0.2)

        # draw curves in the right subplot
        for p in figure.subplots:
            idx = figure.subplots.index(p)
            log.debug("plot_figure --> plot draw: id = %s" % p.id)
            ax = plt.subplot(nbrow, nbcol, idx + 1)

            for curve in p.curves:
                self.__draw_curve(ax, curve, p)

        # title of subplot, axis and legend
        bool = len(figure.subplots) > 1

        for p in figure.subplots:
            idx = figure.subplots.index(p)
            ax = plt.subplot(nbrow, nbcol, idx + 1)

            if p.title:
                ax.set_title(p.title)

            self.__draw_axis(ax, p)

            if len(p.curves) > 0:
                self.__draw_legend(ax, p, bool, hs, ws, ri, le)

        # title of the figure
        if figure.title:
           plt.suptitle(figure.title, fontsize=12)

        # additional matplotlib raw commands for subplot
        for p in figure.subplots:
            idx = figure.subplots.index(p)
            ax = plt.subplot(nbrow, nbcol, idx + 1)
            for cmd in p.cmd:
                try:
                    exec(cmd)
                except:
                    print("Error with the matplotlib command: %s" % cmd)

    #---------------------------------------------------------------------------

    def __save(self, f, figure):
        """method used to save the figure"""
        fmt = figure.fmt
        f = f + "." + fmt
        if fmt == "png":
            dpi = 800
            if figure.dpi:
                dpi = figure.dpi
            plt.savefig(f, format=fmt, dpi=dpi)
        else:
            plt.savefig(f, format=fmt)

#-------------------------------------------------------------------------------
