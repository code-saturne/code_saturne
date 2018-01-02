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

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, string, logging, string
from math import sqrt

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

import vtk
from matplotlib import colors as mpl_colors

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)
#log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
# Scalar
#-------------------------------------------------------------------------------

class Scalar(object):
    """
    Storage of data for a single image
    """
    def __init__ (self, node, parser, fn):
        """
        Constructor of a image.
        @type node: C{DOM Element} instance
        @param node: xml node <plot>
        @type file: C{String}
        @param file: label of the file which contains data of the curve
        """
        # Open file of data
        self.f_name = fn

        f = open(self.f_name, 'r')

        # scalar
        #-------

        # Read mandatory attributes
        self.name     = parser.getAttribute(node, "name")
        self.variable = parser.getAttribute(node, "variable")

        # Read optional attributes

        self.normal = parser.getAttributeTuple(node, "normal", (0, 0, 1))
        self.center = parser.getAttributeTuple(node, "center", ())
        self.stretch = parser.getAttributeTuple(node, "stretch", ())
        self.time_step = float(parser.getAttribute(node, "time-step", -1.))
        self.size = map(int, parser.getAttributeTuple(node, "size", (500, 400)))
        self.zoom = float(parser.getAttribute(node, "zoom", 1.))
        self.wireframe = parser.getAttribute(node, "wireframe", "off")
        self.smult = float(parser.getAttribute(node, "smult", 1.))

        # color scale and color legend
        #-----------------------------

        n_scale = parser.getChild(node, "scale")

        if n_scale:
            self.scale = True

            self.color_map = parser.getAttribute(n_scale, "color", "")
            if self.color_map not in ("hsv", "gray", "hot", "flag", "jet",  \
                                      "blue_to_yellow", "spring", "summer", \
                                      "winter", "autumn"):
                self.color_map = ""

            self.color_srange = parser.getAttributeTuple(n_scale, "range", ())
            self.color_coord = parser.getAttributeTuple(n_scale, "coord", ())
            self.color_levels = int(parser.getAttribute(n_scale, "levels", 10))
            self.color_height = float(parser.getAttribute(n_scale, "height", 0))
            self.color_width = float(parser.getAttribute(n_scale, "width", 0))
            self.color_position = parser.getAttribute(n_scale, "position", "East")
            self.legend = parser.getAttribute(n_scale, "legend", "")
            self.legend_fontsize = int(parser.getAttribute(n_scale, "fontsize", 20))
            self.legend_format = parser.getAttribute(n_scale, "format", "%4.4f")
        else:
            self.scale           = False
            self.color_map       = ""
            self.color_levels    = 10
            self.color_srange    = ()
            self.color_position  = "East"
            self.color_height    = 0
            self.color_width     = 0
            self.color_coord     = ()
            self.legend          = ""
            self.legend_fontsize = 20
            self.legend_format   = "%4.4f"


        # title
        #------

        n_title = parser.getChild(node, "title")

        if n_title:
            self.title = parser.getAttribute(n_title, "label", "")
            self.title_fontsize = int(parser.getAttribute(n_title, "fontsize", 20))
            self.title_coord = parser.getAttributeTuple(n_title, "coord", ())
        else:
            self.title = ""
            self.title_fontsize = 20
            self.title_coord = ()

        # axes
        #-----

        n_axes = parser.getChild(node, "axes")

        if n_axes:
            self.axes = True
            self.axes_fontsize = int(parser.getAttribute(n_axes, "fontsize", 20))
            self.axes_format = parser.getAttribute(n_axes, "format", "%6.3g")
            self.axes_levels = parser.getAttribute(n_axes, "levels", 0)
            self.axes_xlabel = parser.getAttribute(n_axes, "xlabel", "")
            self.axes_ylabel = parser.getAttribute(n_axes, "ylabel", "")
            self.axes_zlabel = parser.getAttribute(n_axes, "zlabel", "")
        else:
            self.axes = False
            self.axes_fontsize = 20
            self.axes_format = "%6.3g"
            self.axes_levels = 0
            self.axes_xlabel = ""
            self.axes_ylabel = ""
            self.axes_zlabel = ""

        # contours
        #------

        n_cont = parser.getChild(node, "contours")

        if n_cont:
            self.cont = parser.getAttribute(n_cont, "status", "off")
            self.cont_nval = int(parser.getAttribute(n_cont, "nval", 10))
            self.cont_range = parser.getAttributeTuple(n_cont, "range", ())
            self.cont_color = parser.getAttributeTuple(n_cont, "color", "k")
        else:
            self.cont = "off"
            self.cont_nval = 10
            self.cont_range = ()
            self.cont_color = "k"

        # raw commands
        #-------------

        self.cmd = parser.getVTKCommands(node)

        f.close()

#-------------------------------------------------------------------------------
# Plotter
#-------------------------------------------------------------------------------

class PlotVTK(object):
    """
    Manager of the vtk commands.
    """
    def __init__ (self, parser):
        """
        Constructor of the plotter.
        @param parser: parser of the xml file
        """
        self.parser = parser


    def plot_study(self, study_label, study_object):
        """
        Method used to plot all plots from a I{study_label} (all cases).
        @param study_label: label of a study
        """
        # initialisation for each Study
        self.figures = []

        # Read the files of results
        for case in study_object.Cases:
            if case.plot == "on" and case.is_run != "KO":
                for node in self.parser.getChildren(case.node, "resu"):
                    plots, file_name, dest, repo = self.parser.getResult(node)

                    if dest:
                        d = dest
                        dd = self.parser.getDestination()
                    elif repo:
                        d = repo
                        dd = self.parser.getRepository()

                    f = os.path.join(dd,
                                     study_label,
                                     case.label, "RESU",
                                     d, "postprocessing", file_name)

                    if not os.path.isfile(f):
                        raise ValueError("\n\nThis file does not exist: %s\n\n" % f)

                    for nn in plots:
                        self.figures.append(Scalar(nn, self.parser, f))

        for options in self.figures:
            # Plot
            figure = Builder(options)

            # save the figure
            f = os.path.join(self.parser.getDestination(),
                             study_label,
                             "POST",
                             options.name)
            figure.screenshot(f)

            # store the name of the figure for the build of
            # the detailed report
            study_object.vtk_figures.append(f)

#-------------------------------------------------------------------------------
# Draw the scene
#-------------------------------------------------------------------------------

class Builder(object):
    def __init__(self, options):
        """VTK window manager"""
        self.opt  = options

        ds = self.readData(self.opt.time_step)
        cellData = ds.GetCellData()

        if not cellData.HasArray(self.opt.variable):
            print("Error: variable %s not found." % self.opt.variable)
            sys.exit(1)

        # get dimension of data field
        dim = cellData.GetArray(self.opt.variable).GetNumberOfComponents()

        convert = self.convertCell2Point(ds)

        # Multiply given variable by a real number using vtk array calculator
        calc = vtk.vtkArrayCalculator()
        calc.SetInputConnection(convert.GetOutputPort())
        if dim == 1:
            calc.AddScalarVariable("var", self.opt.variable, 0)
            attribute = "SCALARS"
        else:
            calc.AddVectorVariable("var", self.opt.variable, 0 , 1, 2)
            attribute = "VECTORS"

        calc.SetFunction("var*%f"%self.opt.smult)
        resultName = self.opt.variable+"transf"
        calc.SetResultArrayName(resultName)

        aa = vtk.vtkAssignAttribute()
        aa.SetInputConnection(calc.GetOutputPort())
        aa.Assign(resultName, attribute, "POINT_DATA")
        aa.Update()

        convert = aa
        self.opt.variable = resultName

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(1., 1., 1.)
        self.cam = self.ren.GetActiveCamera()
        self.ren.ResetCamera()

        # mapper for the plane cut
        self.cutMapper = vtk.vtkPolyDataMapper()

        # define the cut plane and the cut input
        cut = self.cutPlane(convert)

        if self.opt.stretch:
            t = vtk.vtkTransform()
            m = vtk.vtkTransformPolyDataFilter()

            t.Scale(self.opt.stretch[0],
                    self.opt.stretch[1],
                    self.opt.stretch[2])

            m.SetInputConnection(cut.GetOutputPort())
            m.SetTransform(t)

            self.cutMapper.SetInputConnection(m.GetOutputPort())
        else:
            self.cutMapper.SetInputConnection(cut.GetOutputPort())

        if not self.opt.color_map:
            self.cutMapper.CreateDefaultLookupTable()
            lut = self.cutMapper.GetLookupTable()
        else:
            lut = eval("self." + self.opt.color_map + "()")

        # scalar map of the magnitude of vector fields
        lut.SetVectorModeToMagnitude()

        self.colorDataSetByArray(convert.GetOutput(), lut, self.cutMapper)
        grid = vtk.vtkLODActor()

        grid.SetMapper(self.cutMapper)

        if self.opt.wireframe == "on":
            grid.GetProperty().SetRepresentationToWireframe()

        self.ren.AddActor(grid)

        if self.opt.cont == 'on':
            # mapper for the contours
            self.contoursMapper = vtk.vtkPolyDataMapper()
            # define contours on cut plane
            contours = self.contoursFilterOnCutPlane(convert, cut, dim)

            # apply same transformation as to the cut
            if self.opt.stretch:
                t = vtk.vtkTransform()
                cm = vtk.vtkTransformPolyDataFilter()

                t.Scale(self.opt.stretch[0],
                        self.opt.stretch[1],
                        self.opt.stretch[2])

                cm.SetInputConnection(contours.GetOutputPort())
                cm.SetTransform(t)

                self.contoursMapper.SetInputConnection(cm.GetOutputPort())
            else:
                self.contoursMapper.SetInputConnection(contours.GetOutputPort())

            self.colorDataSetByArray(convert.GetOutput(), lut, self.contoursMapper)

            # create actor and add it to the renderer
            cont_grid = vtk.vtkLODActor()
            cont_grid.SetMapper(self.contoursMapper)
            rgb = mpl_colors.ColorConverter().to_rgb(self.opt.cont_color)
            cont_grid.GetProperty().SetColor(rgb)
            self.ren.AddActor(cont_grid)

        if self.opt.axes:
            if self.opt.stretch:
                axes = self.addAxes(m)
            else:
                axes = self.addAxes(cut)
            axes.SetCamera(self.cam)
            self.ren.AddActor(axes)

        if self.opt.title:
            title = self.setTitle()
            self.ren.AddActor(title)

        if self.opt.scale:
            legend = self.addLegend(lut)
            self.ren.AddActor(legend)

        # Camera parameters

        if self.opt.axes and self.opt.zoom == 1.:
            self.cam.Zoom(0.9)
        elif self.opt.zoom != 1.:
            self.cam.Zoom(self.opt.zoom)

        if self.opt.center:
            self.cam.SetFocalPoint(self.opt.center[0],
                                   self.opt.center[1],
                                   self.opt.center[2])
            self.cam.SetPosition(self.opt.center[0] - self.opt.normal[0],
                                 self.opt.center[1] - self.opt.normal[1],
                                 self.opt.center[2] - self.opt.normal[2])
        else:
            self.cam.SetFocalPoint(0,0,0)
            self.cam.SetPosition(-1 * self.opt.normal[0],
                                 -1 * self.opt.normal[1],
                                 -1 * self.opt.normal[2])

        if self.opt.axes:
            axes.SetXAxisVisibility(1)
            axes.SetYAxisVisibility(1)
            axes.SetZAxisVisibility(1)
            if self.opt.normal[1] == 0 and self.opt.normal[2] == 0:
                axes.SetYAxisVisibility(0)
                self.cam.SetViewUp(0, 0, 1)
            elif self.opt.normal[0] == 0 and self.opt.normal[2] == 0:
                axes.SetZAxisVisibility(0)
                self.cam.SetViewUp(0, 0, 1)
            elif self.opt.normal[0] == 0 and self.opt.normal[1] == 0:
                axes.SetYAxisVisibility(0)
                self.cam.SetViewUp(0, 1, 0)

        cam = self.cam

        for cmd in self.opt.cmd:
            c = open("./tmp.py", "w")
            c.write(cmd)
            c.close()
            try:
                execfile("./tmp.py")
            except:
                print("Error with the vtk command: %s" % cmd)
            os.remove("./tmp.py")

        self.win = vtk.vtkRenderWindow()
        self.win.SetSize(self.opt.size[0], self.opt.size[1])
        self.win.AddRenderer(self.ren)
        self.win.Render()


    def readData(self, ntime = -1.):
        """
        Return the mesh vtkDataSet with the associated scalars for a single time step.
        """
        p = os.path.abspath(self.opt.f_name)
        path, case = os.path.dirname(p), os.path.basename(p)
        ens = vtk.vtkGenericEnSightReader()
        ens.SetByteOrderToLittleEndian()
        ens.SetFilePath(path)
        ens.SetCaseFileName(case)
        ens.ReadAllVariablesOn()
        ens.Update()

        times = ens.GetTimeSets().GetItem(0)

        # last result by default
        if ntime == -1.:
            ntime = times.GetSize() - 1

        ens.SetTimeValue(times.GetTuple1(ntime))
        ens.Update()
        return ens.GetOutput().GetBlock(0) # Now multiblock objet with 0 for volume, 1 for bord !


    def convertCell2Point(self, cellDataSet):
        """Apply filter vtkCellDataToPointData"""
        convert = vtk.vtkCellDataToPointData()
        try:
            convert.SetInputData(cellDataSet)
        except AttributeError:
            convert.SetInput(cellDataSet)
        convert.Update()
        return convert


    def cutPlane(self, ptDataSet):
        """Extract a slice from the 3D computational domain."""
        # The (implicit) plane is used to do the cutting
        plane = vtk.vtkPlane()
        if self.opt.center:
            plane.SetOrigin(self.opt.center[0],
                            self.opt.center[1],
                            self.opt.center[2])
        else:
            plane.SetOrigin(ptDataSet.GetOutput().GetCenter())

        plane.SetNormal(self.opt.normal[0],
                        self.opt.normal[1],
                        self.opt.normal[2])

        # The cutter is set up to process each contour value over all cells
        # (SetSortByToSortByCell). This results in an ordered output of polygons
        # which is key to the compositing.
        cut = vtk.vtkCutter()
        cut.SetInputConnection(ptDataSet.GetOutputPort())
        cut.SetCutFunction(plane)
        #cut.GenerateCutScalarsOn()
        #cut.SetSortByToSortByCell()
        return cut

    def contoursFilterOnCutPlane(self, convert, cut, dim):
        """Create contours on a slice from the 3D computational domain."""
        # create the contour filter
        contours = vtk.vtkContourFilter()

        # scalar field
        if dim == 1:
            # set the plane as input
            contours.SetInputConnection(cut.GetOutputPort())
            # select the variable
            contours.SetInputArrayToProcess(0,0,0,0,self.opt.variable)
        # contours of the magnitude of vector fields
        elif dim > 1:
            vectorNorm = vtk.vtkVectorNorm()
            vectorNorm.SetInputConnection(cut.GetOutputPort())
            vectorNorm.SetInputArrayToProcess(0,0,0,0,self.opt.variable)
            contours.SetInputConnection(vectorNorm.GetOutputPort())

        # define the iso values
        self.isoValuesByArray(convert.GetOutput(), contours)

        contours.Update()

        return contours


    def textProperty(self, fontsize = 20):
        """Return properties for a text."""
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0., 0., 0.)
        tprop.SetFontSize(fontsize)
        tprop.SetFontFamilyToArial()
        tprop.BoldOff()
        return tprop


    def setTitle(self):
        """Put a text at the top of the figure."""
        tprop = self.textProperty(fontsize = self.opt.title_fontsize)
        tprop.SetVerticalJustificationToTop()
        tprop.SetJustificationToCentered()
        mapper = vtk.vtkTextMapper()
        mapper.SetInput(self.opt.title)
        mapper.SetTextProperty(tprop)
        actor = vtk.vtkActor2D()
        actor.SetMapper(mapper)
        actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        if self.opt.title_coord:
            actor.GetPositionCoordinate().SetValue(self.opt.title_coord[0],
                                                   self.opt.title_coord[1])
        else:
            actor.GetPositionCoordinate().SetValue(0.5, 0.99)

        return actor


    def addLegend(self, lut):
        """Add a color bar and a color legend"""
        tprop = self.textProperty(fontsize = self.opt.legend_fontsize)
        legend = vtk.vtkScalarBarActor()
        legend.SetLookupTable(lut)
        legend.SetLabelTextProperty(tprop)
        legend.GetLabelTextProperty().BoldOff()
        legend.SetLabelFormat(self.opt.legend_format)
        legend.SetNumberOfLabels(self.opt.color_levels)
        legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()

        if self.opt.legend:
            legend.SetTitle(self.opt.legend)
            legend.SetTitleTextProperty(tprop)
            legend.GetTitleTextProperty().BoldOff()

        if self.opt.color_position in ('South', 'North'):
            legend.SetOrientationToHorizontal()
        else:
            legend.SetOrientationToVertical()

        __locations = {'North': (0.1,  0.9 , 0.8, 0.1),
                       'South': (0.1,  0.01, 0.8, 0.1),
                       'West':  (0.01, 0.1,  0.1, 0.8),
                       'East':  (0.9,  0.09, 0.1, 0.8)}

        cbloc = __locations[self.opt.color_position]

        if not self.opt.color_width:
            legend.SetWidth(cbloc[2])
        else:
            legend.SetWidth(self.opt.color_width)

        if not self.opt.color_height:
            legend.SetHeight(cbloc[3])
        else:
            legend.SetHeight(self.opt.color_height)

        if not self.opt.color_coord:
            legend.GetPositionCoordinate().SetValue(cbloc[0], cbloc[1])
        else:
            legend.GetPositionCoordinate().SetValue(self.opt.color_coord[0],
                                                    self.opt.color_coord[1])

        return legend


    def addAxes(self, obj):
        """Add axes"""
        tprop = self.textProperty(fontsize = self.opt.axes_fontsize)

        axes = vtk.vtkCubeAxesActor2D()
        try:
            axes.SetInputData(obj.GetOutput())
        except AttributeError:
            axes.SetInput(obj.GetOutput())
        axes.SetLabelFormat(self.opt.axes_format)
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(1.5)
        axes.SetCornerOffset(0.0)

        axes.GetProperty().SetColor(tprop.GetColor())
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)

        if self.opt.axes_levels:
            axes.SetNumberOfLabels()
        if self.opt.axes_xlabel:
            axes.SetXLabel(self.opt.axes_xlabel)
        if self.opt.axes_ylabel:
            axes.SetYLabel(self.opt.axes_ylabel)
        if self.opt.axes_zlabel:
            axes.SetZLabel(self.opt.axes_zlabel)

        return axes


    def colorDataSetByArray(self, data, lut, mapper):
        """Build the color map."""
        array = None
        flag = True
        mapper.Modified()

        if data.GetCellData() and data.GetCellData().HasArray(self.opt.variable):
            array = data.GetCellData().GetArray(self.opt.variable)
            if array:
                mapper.SetScalarModeToUseCellFieldData()

        if not array and data.GetPointData() and data.GetPointData().HasArray(self.opt.variable):
            array = data.GetPointData().GetArray(self.opt.variable)
            if array:
                mapper.SetScalarModeToUsePointFieldData()

        if not array:
            mapper.SetScalarModeToDefault()
            mapper.SetInterpolateScalarsBeforeMapping(1)
            print("Variable name not found in data set.")
            flag = False

        if flag:
            mapper.SetLookupTable(lut)
            if self.opt.color_srange:
                mapper.SetScalarRange(self.opt.color_srange[0],
                                      self.opt.color_srange[1])
            else:
                comp = -1
                mapper.SetScalarRange(array.GetRange(comp))
            mapper.SetInterpolateScalarsBeforeMapping(1)
            mapper.SelectColorArray(array.GetName())


    def isoValuesByArray(self, data, contours):
        """Generate the isovalues from data range."""
        array = None

        if data.GetCellData() and data.GetCellData().HasArray(self.opt.variable):
            array = data.GetCellData().GetArray(self.opt.variable)

        if not array and data.GetPointData() and data.GetPointData().HasArray(self.opt.variable):
            array = data.GetPointData().GetArray(self.opt.variable)

        if not array:
            print("Error: variable %s not found." % self.opt.variable)
            sys.exit(1)

        if self.opt.cont_range:
            contours.GenerateValues(self.opt.cont_nval, self.opt.cont_range)
        else:
            comp = -1
            contours.GenerateValues(self.opt.cont_nval, array.GetRange(comp))


    def screenshot(self, f, mode = "PNG"):
        """Save the image to png format."""
        if mode not in ["PostScript", "JPEG", "PNG", "TIFF"]:
            raise ValueError("format unknown: %s" % mode)

        # capture what is being shown on self.renderer
        scr = vtk.vtkWindowToImageFilter()
        scr.SetInput(self.win)
        scr.Update()

        # write data to file
        writer = eval("vtk.vtk" + mode + "Writer()")
        writer.SetInputConnection(scr.GetOutputPort())
        writer.SetFileName(f + ".png")
        writer.Write()

#-------------------------------------------------------------------------------
# Colormaps from Scitools
# http://code.google.com/p/scitools/
#
# Copyright (c) 2007-2009, Hans Petter Langtangen <hpl@simula.no> and
# Simula Resarch Laboratory.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#
#    * Neither the name of Simula Research Laboratory nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#-------------------------------------------------------------------------------

    def hsv(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.0, 1.0)
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.SetNumberOfColors(m)
        lut.Build()
        return lut


    def gray(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.0, 0.0)
        lut.SetSaturationRange(0.0, 0.0)
        lut.SetValueRange(0.0, 1.0)
        lut.SetNumberOfColors(m)
        lut.Build()
        return lut


    def hot(self, m=256):
        lut = vtk.vtkLookupTable()
        inc = 0.01175
        lut.SetNumberOfColors(256)
        i = 0
        r = 0.0; g = 0.0; b = 0.0
        while r <= 1.:
            lut.SetTableValue(i, r, g, b, 1)
            r += inc;  i += 1
        r = 1.
        while g <= 1.:
            lut.SetTableValue(i, r, g, b, 1)
            g += inc;  i += 1
        g = 1.
        while b <= 1:
            if i == 256: break
            lut.SetTableValue(i, r, g, b, 1)
            b += inc;  i += 1
        lut.Build()
        return lut


    def flag(self, m=64):
        """Alternating red, white, blue, and black color map.

        - flag(m)
          'm' must be a multiple of 4
        """
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        # the last parameter alpha is set to 1 by default
        # in method declaration
        for i in range(0,m,4):
            lut.SetTableValue(i,1,0,0,1)   # red
            lut.SetTableValue(1+i,1,1,1,1) # white
            lut.SetTableValue(2+i,0,0,1,1) # blue
            lut.SetTableValue(3+i,0,0,0,1) # black
        lut.Build()
        return lut


    def jet(self, m=256):
        # blue, cyan, green, yellow, red, black
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        lut.SetHueRange(0.667,0.0)
        lut.Build()
        return lut


    def blue_to_yellow(self, m=200):
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        for i in range(m):
            frac = i / float(m / 2.0 - 1.0)
            if (frac <= 1):
                r = frac
                g = r
                b = 1
            else:
                r = 1
                g = r
                b = 2 - frac
            # SetTableValue(indx, red, green, blue, alpha)
            lut.SetTableValue(i, r, g, b, 1)
        lut.Build()
        return lut


    def spring(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        lut.SetHueRange(0.0, 0.17)
        lut.SetSaturationRange(0.5, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.Build()
        return lut


    def summer(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        lut.SetHueRange(0.47, 0.17)
        lut.SetSaturationRange(1.0, 0.6)
        lut.SetValueRange(0.5, 1.0)
        lut.Build()
        return lut


    def winter(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        lut.SetHueRange(0.8, 0.42)
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(0.6, 1.0)
        lut.Build()
        return lut


    def autumn(self, m=256):
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(m)
        lut.SetHueRange(0.0, 0.15)
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.Build()
        return lut

#-------------------------------------------------------------------------------
