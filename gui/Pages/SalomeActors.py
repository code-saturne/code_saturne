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
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore     import *
from code_saturne.Base.QtGui      import *
from code_saturne.Base.QtWidgets  import *

import vtk

import libSalomePy
import salome
import SMESH
import CFDSTUDYGUI_DataModel

try:
    from libvtkRenderingCorePython import *
    from libvtkFiltersSourcesPython import *
except:
    try:
        from vtkRenderingCorePython import *
        from vtkFiltersSourcesPython import *
    except:
        # for compatibility with salome 6.6
        try:
            from libvtkRenderingPython import *
        except:
            from vtkRenderingPython import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SalomeActors")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ActorBase(object):
    def __init__(self, render, label, color):
        self.render = render

        self.colors = {
               'Marine Blue': [0.1, 0.2, 0.4],
               'Deep Blue':   [0.2, 0.2, 1.0],
               'Red':         [1.0, 0.0, 0.0],
               'Yellow':      [1.0, 1.0, 0.0],
               'Green':       [0.0, 1.0, 0.0],
               'Orange':      [1.0, 0.5, 0.0],
               'Magenta':     [1.0, 0.0, 1.0],
               'Grey':        [0.5, 0.5, 0.5],
               'Black':       [0.0, 0.0, 0.0],
               'Cyan':        [0.0, 1.0, 1.0],
               'Turquoise':   [0.0, 0.5, 0.5],
               'White':       [1.0, 1.0, 1.0],
               'Light green': [0.5, 1.0, 0.5],
                }

        if color:
            self.color = color
        else:
            self.color = "Yellow"

        self.__isSelected = False

        self.points = vtk.vtkPoints()
        self.points.SetDataTypeToFloat()

        self.property = vtk.vtkProperty()
        self.actor = vtk.vtkActor()
        self.actor.SetProperty(self.property)

        self.label = label
        self.labels = vtk.vtkVectorText()
        self.labels.SetText(label)
        self.labelProp = vtk.vtkProperty()
        self.labelActor = vtk.vtkFollower()


    def setOpacity(self, opacity):
        if opacity >= 0.0 and opacity <= 1.0:
            self.property.SetOpacity(opacity)
            self.labelProp.SetOpacity(opacity)


    def getOpacity(self):
        return self.property.GetOpacity()


    def setColor(self, color):
        if color in self.colors.keys():
            self.color = color
            self.property.SetColor(self.colors[color])
            self.labelProp.SetColor(self.colors[color])


    def getColor(self):
        return self.color


    def getColors(self,):
        return self.colors.keys()


    def getActor(self):
        return self.actor


    def getLabelActor(self):
        return self.labelActor


    def getRadius(self):
        return self.radius


    def getLabel(self):
        return self.label


    def setLabel(self, label):
        self.label = label
        self.labels.SetText(label)


    def select(self):
        self.actor.GetProperty().SetAmbient(1)
        self.labelActor.GetProperty().SetAmbient(1)
        #self.labelProp.SetAmbient(1)
        self.__isSelected = True


    def unSelect(self):
        self.actor.GetProperty().SetAmbient(0)
        self.labelActor.GetProperty().SetAmbient(0)
        #self.labelProp.SetAmbient(0)
        self.__isSelected = False


    def selected(self):
        return self.__isSelected

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbeActor(ActorBase):
    def __init__(self, render, label, point, radius, color = ""):
        """
        """
        ActorBase.__init__(self, render, label, color)

        self.points.InsertNextPoint(point[0], point[1], point[2])

        self.radius = radius
        self.__sphere = vtk.vtkSphereSource()
        self.__sphere.SetThetaResolution(10)
        self.__sphere.SetPhiResolution(10)
        self.__sphere.SetRadius(self.radius)

        self.polyData = vtk.vtkPolyData()
        self.polyData.SetPoints(self.points)

        self.__glyph = vtk.vtkGlyph3D()
        try:
            self.__glyph.SetInputData(self.polyData)
        except:
            # for compatibility with salome 6.6
            self.__glyph.SetInput(self.polyData)

        try:
            self.__glyph.SetSourceConnection(self.__sphere.GetOutputPort())
        except:
            # for compatibility with salome 6.6
            self.__glyph.SetSource(self.__sphere.GetOutput())

        self.__mapper = vtk.vtkPolyDataMapper()
        self.__mapper.SetInputConnection(self.__glyph.GetOutputPort(0))

        self.actor.SetMapper(self.__mapper)
        self.render.AddActor(self.actor)

        self.__labelMapper = vtk.vtkPolyDataMapper()
        self.__labelMapper.SetInputConnection(self.labels.GetOutputPort(0))
        self.labelActor.SetMapper(self.__labelMapper)
        self.labelActor.SetProperty(self.labelProp)
        scale = 1.5 * self.radius
        self.labelActor.SetScale(scale, scale, scale)
        factor = 1.2
        self.labelActor.AddPosition(point[0] + factor * self.radius,
                                    point[1] + factor * self.radius,
                                    point[2] + factor * self.radius)

        self.render.AddActor2D(self.labelActor)
        self.labelActor.SetCamera(self.render.GetActiveCamera())

        self.setColor(self.color)
        self.setOpacity(self.getOpacity())


    def setRadius(self, radius):
        if radius > 0:
            self.radius = radius
            self.__sphere.SetRadius(self.radius)
            self.__sphere.Update()

            scale = 1.5 * self.radius
            self.labelActor.SetScale(scale, scale, scale)


    def updateLocation(self, point):
        idx = 0
        factor = 1.2
        self.points.InsertPoint(idx, point[0], point[1], point[2])
        self.points.Modified()
        self.labelActor.SetPosition(point[0] + factor * self.radius,
                                    point[1] + factor * self.radius,
                                    point[2] + factor * self.radius)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProfilActor(ActorBase):
    def __init__(self, render, label, p1, p2, radius, color = ""):
        """
        """
        ActorBase.__init__(self, render, label, color)

        pt1 = self.points.InsertNextPoint(p1[0], p1[1], p1[2])
        pt2 = self.points.InsertNextPoint(p2[0], p2[1], p2[2])

        self.radius = radius
        r = vtk.vtkFloatArray()
        r.InsertNextValue(self.radius)

        self.__array = vtk.vtkCellArray()
        self.__array.InsertNextCell(2)
        self.__array.InsertCellPoint(pt1)
        self.__array.InsertCellPoint(pt2)

        self.polyData = vtk.vtkPolyData()
        self.polyData.SetPoints(self.points)
        self.polyData.SetLines(self.__array)
        self.polyData.GetCellData().SetScalars(self.radius)

        # TODO: to finish
        #edgeColors = vtk.vtkIdTypeArray()
        #edgeColors.InsertNextValue(2)
        #self.polyData.GetCellData().AddArray(edgeColors)

        self.__tubes = vtk.vtkTubeFilter()
        try:
            self.__tubes.SetInputData(self.polyData)
        except:
            # for compatibility with salome 6.6
            self.__tubes.SetInput(self.polyData)
        self.__tubes.SetVaryRadius(0)
        self.__tubes.SetRadius(self.radius)
        self.__tubes.CappingOff()
        self.__tubes.SetNumberOfSides(12)
        self.__tubes.SetGenerateTCoordsToOff()

        self.__mapper = vtk.vtkPolyDataMapper()
        self.__mapper.SetInputConnection(self.__tubes.GetOutputPort(0))

        self.actor.SetMapper(self.__mapper)
        self.render.AddActor(self.actor)

        self.__labelMapper = vtk.vtkPolyDataMapper()
        self.__labelMapper.SetInputConnection(self.labels.GetOutputPort(0))
        self.labelActor.SetMapper(self.__labelMapper)
        self.labelActor.SetProperty(self.labelProp)
        scale = 7.5 * self.radius
        self.labelActor.SetScale(scale, scale, scale)
        factor = 1.2
        c0 = 0.5 * (p1[0] + p2[0])
        c1 = 0.5 * (p1[1] + p2[1])
        c2 = 0.5 * (p1[2] + p2[2])
        self.labelActor.AddPosition(c0 + factor * self.radius,
                                    c1 + factor * self.radius,
                                    c2 + factor * self.radius)

        self.render.AddActor2D(self.labelActor)
        self.labelActor.SetCamera(self.render.GetActiveCamera())

        self.setColor(self.color)
        self.setOpacity(self.getOpacity())


    def setRadius(self, radius):
        if radius > 0:
            self.radius = radius
            self.polyData.GetCellData().SetScalars(self.radius)
            self.__tubes.SetRadius(self.radius)
            self.__tubes.Update()

            scale = 1.5 * self.radius
            self.labelActor.SetScale(scale, scale, scale)


    # TODO: to finish
    #def updateLocation(self, point):
        #idx = 0
        #factor = 1.2
        #self.points.InsertPoint(idx, point[0], point[1], point[2])
        #self.points.Modified()
        #self.labelActor.SetPosition(point[0] + factor * self.radius,
                                      #point[1] + factor * self.radius,
                                      #point[2] + factor * self.radius)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class Actors(object):
    """
    Manage the collection of Generic Actors.
    """
    def __init__(self):
        """
        http://www.google.com/codesearch/p?hl=fr#a-DtlzUst6U/Examples/Annotation/Python/textOrigin.py&q=annotatePick&d=3
        """
        #super(Actors, self).__init__()

        self.render = libSalomePy.getRenderer()
        self.__interactor = libSalomePy.getRenderWindowInteractor()
        self.picker = self.__interactor.GetPicker()
        self.__interactor.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)

        self.radius = 0.025
        self.p = []


    def setTableView(self, view):
        self.__view = view


    def refreshView(self):
        salome.sg.UpdateView()
        #libSalomePy.fitAll()


    def add(self, p):
        """
        Add a new single actor.

        @type p: C{ProbeActor}
        @param p: new ProbeActor to store in the collection
        """
        self.p.append(p)
        self.refreshView()


    def actors(self):
        return self.p


    def labels(self):
        l = []
        for p in self.p:
            l.append(p.getLabel())
        return l


    def remove(self, label):
        log.debug("remove -> %s in %s" % (label, self.labels()))
        idx = self.labels().index(label)
        self.p[idx].setLabel("")
        self.render.RemoveActor(self.p[idx].getLabelActor())
        self.render.RemoveActor(self.p[idx].getActor())
        t = self.p.pop(idx)
        del t
        # Renumber the probes or profils
        for i in range(len(self.p)):
            self.p[i].setLabel(str(i+1))
        log.debug("remove renumbering -> %s" % (self.labels(),))
        self.refreshView()


    def removeActors(self):
        for p in self.p:
            self.render.RemoveActor(p.getLabelActor())
            self.render.RemoveActor(p.getActor())
        self.p = []


    def setRadius(self, r):
        self.radius = r
        for a in self.actors():
            a.setRadius(r)
        self.refreshView()


    def getRadius(self):
        return self.radius


    def setVisibility(self, b):
        for a in self.actors():
            a.getActor().SetVisibility(b)
            a.getLabelActor().SetVisibility(b)
        self.refreshView()


    def getVisibility(self):
        try:
            r = self.actors()[0].getActor().GetVisibility()
        except:
            r = 1

        return r


    def select(self, label):
        log.debug("select -> %s in %s" % (label, self.labels()))
        idx = self.labels().index(label)
        self.p[idx].select()
        self.refreshView()


    def unSelect(self, label):
        idx = self.labels().index(label)
        self.p[idx].unSelect()
        self.refreshView()


    def unSelectAll(self):
        for p in self.p:
            p.unSelect()
        self.refreshView()


    def leftButtonPressEvent(self, obj=None, event=""):
        x, y = obj.GetEventPosition()
        if self.picker.PickProp(x, y, self.render):
            actor = self.picker.GetActor()
            if obj.GetShiftKey() > 0:
                self.__select(actor)
            else:
                self.__unSelectAll()
                self.__select(actor)
        else:
            self.__unSelectAll()


    def __actors(self):
        l = []
        for p in self.p:
            l.append(p.getActor())
        return l


    def __select(self, actor):
        idx = self.__actors().index(actor)
        self.p[idx].select()

        if self.__view:
            self.__view.setSelectionMode(QAbstractItemView.MultiSelection)
            self.__view.clearSelection()
            for idx in range(len(self.p)):
                if self.p[idx].selected():
                    if idx not in self.__view.selectionModel().selectedRows():
                        self.__view.selectRow(idx)
            self.__view.setSelectionMode(QAbstractItemView.ExtendedSelection)


    def __unSelectAll(self):
        self.unSelectAll()

        if self.__view:
            self.__view.clearSelection()


    def __del__(self):
        self.__unSelectAll()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProbeActors(Actors):
    """
    Manage the collection of Probes Actors.
    """
    def addProbe(self, label, point):
        """
        Create a new Actor.
        """
        log.debug("addProbe -> %s: %s" % (label, point))

        if label not in self.labels():
            self.add(ProbeActor(self.render, label, point, self.radius))


    def updateLocation(self, label, point):
        log.debug("updateLocation -> %s: %s" % (label, point))
        idx = self.labels().index(label)
        self.p[idx].updateLocation(point)
        self.refreshView()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class ProfilActors(Actors):
    """
    Manage the collection of Profiles Actors.
    """
    def addProfil(self, label, point1, point2):
        log.debug("addProfil -> %s: %s %s" % (label, point1, point2))
        self.add(ProfilActor(self.render, label, point1, point2, self.radius))


    def updateLocation(self, label, point1, point2):
        log.debug("updateLocation -> %s: %s %s" % (label, point1, point2))
        idx = self.labels().index(label)
        self.p[idx].updateLocation(point1, point2)
        self.refreshView()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


#def getMesh(label):
    #for aStudy in CFDSTUDYGUI_DataModel.GetStudyList():
        #print("study = " +  aStudy + " " + aStudy.GetName())
        #aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aStudy, label)

    #print("tototototo = " + " ".join(aChildList))
    #aDataObj = aChildList[0]
    #print(aDataObj.GetName())


#def BoundingBox(my_mesh):
    #xyz_max = [-1e20, -1e20, -1e20]
    #xyz_min = [1e20, 1e20, 1e20]
    #for id in my_mesh.GetNodesId():
        #xyz_id = my_mesh.GetNodeXYZ(id)
        #for i in (0, 1, 2):
            #xyz_max[i] = max(xyz_id[i], xyz_max[i])
            #xyz_min[i] = min(xyz_id[i], xyz_min[i])

    #return xyz_min, xyz_max

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------



