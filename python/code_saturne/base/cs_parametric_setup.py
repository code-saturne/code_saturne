#!/usr/bin/env python3

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
Update XML file with common parametric arguments.
"""

import os, sys
from argparse import ArgumentParser

#-------------------------------------------------------------------------------
# Process the command line
#-------------------------------------------------------------------------------

def arg_parser(argv):
    """
    Build argument parser for command line arguments.
    """

    parser = ArgumentParser(description="Filter a case setup model.")

    parser.add_argument("--notebook", "--nb", dest="notebook", type=str,
                        action="append",
                        metavar="<var>:<val>",
                        help="Notebook parameters.")

    parser.add_argument("-m", "--mesh", dest="mesh", type=str,
                        action="append",
                        help="Name of the new mesh")

    parser.add_argument("--mi", "--mesh_input", dest="mesh_input", type=str,
                        help="Name of the new mesh_input")

    parser.add_argument("-a", "--perio-angle", dest="rotationAngle", type=float,
                        help="Periodicity angle")

    parser.add_argument("-r", "--restart", dest="RestartRun", type=str,
                        help="Run to restart from (in same case)")

    parser.add_argument("--different-restart-mesh", dest="DiffRestartMesh",
                        action="store_true",
                        help="Restart from run on a different mesh. " \
                        "Will check for mesh in checkpoint folder.")

    parser.add_argument("--different-restart-mesh-path", dest="DiffRestartMeshPath", type=str,
                        help="Restart from run on a different mesh. " \
                        "Provide original mesh with this argument.")

    parser.add_argument("-n", "--iter-num", dest="iterationsNumber", type=int,
                        help="New iteration number")

    parser.add_argument("--tmax", dest="tmax", type=float,
                        help="Time based stop criterion")

    parser.add_argument("--iter-dt", dest="TimeStep", type=float,
                        help="New time step")

    parser.add_argument("-u", "--unsteady", dest="TimeModel", action="store_true",
                        help="Unsteady time model (steady by default)")

    parser.add_argument("--imrgra", dest="imrgra", type=int,
                        help="Gradient reconstruction")

    parser.add_argument("--blencv", dest="blencv", type=str, nargs='*',
                        metavar="<var>:<val>",
                        help="Blencv")

    parser.add_argument("--update-bc-criteria", dest="bc_criteria",
                        type=str, action="append",
                        metavar="<bc_label>:<selection_criteria>",
                        help="Update boundary zone selection criteria.")

    return parser

#-------------------------------------------------------------------------------
# Class which provides a set of predefined methods to modify a
# code_saturne xml file
#-------------------------------------------------------------------------------

class case_setup_filter(object):
    """
    Class to modify a code_saturne xml file
    """

    #---------------------------------------------------------------------------

    def __init__(self, case=None, pkg=None):
        """
        Constructor.
        @param case: XML case object
        @param pkg: package
        """

        # package
        self.pkg = None
        if pkg:
            self.pkg = pkg
        else:
            from code_saturne.base.cs_package import package
            self.pkg = package()

        from code_saturne.model.XMLengine import Case

        self.case = case

        self.outputModel   = None
        self.bcModel       = None
        self.bndModel      = None
        self.restartModel  = None
        self.timeStepModel = None
        self.meshModel     = None
        self.notebookModel = None
        self.numParamModel = None
        self.numParamEquationModel = None

    #---------------------------------------------------------------------------

    def initMeshModel(self):
        """
        Initialize the SolutionDomainModel
        """

        if self.meshModel is None:
            from code_saturne.model.SolutionDomainModel import SolutionDomainModel
            self.meshModel = SolutionDomainModel(self.case)

    #---------------------------------------------------------------------------

    def initOutputModel(self):
        """
        Method to initialize the outputModel data structure which handles
        output parameters including monitoring points
        """
        if self.outputModel is None:
            from code_saturne.model.OutputControlModel import OutputControlModel
            self.outputModel = OutputControlModel(self.case)

    #---------------------------------------------------------------------------

    def initBcModel(self):
        """
        Initialize the BC model
        """

        if self.bcModel is None:
            from code_saturne.model.LocalizationModel import LocalizationModel
            self.bcModel = LocalizationModel("BoundaryZone", self.case)

    #---------------------------------------------------------------------------

    def initBndModel(self):
        """
        Initialize the boundary model.
        """

        if self.bndModel is None:
            from code_saturne.model.Boundary import Boundary
            self.bndModel = Boundary

    #---------------------------------------------------------------------------

    def initRestartModel(self):
        """
        Initialize the restart model
        """

        if self.restartModel is None:
            from code_saturne.model.StartRestartModel import StartRestartModel
            self.restartModel = StartRestartModel(self.case)

    #---------------------------------------------------------------------------

    def initTimeStepModel(self):
        """
        Initialize the time step model.
        """

        if self.timeStepModel is None:
            if self.case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
                from code_saturne.model.TimeStepModelNeptune import TimeStepModel
            else:
                from code_saturne.model.TimeStepModel import TimeStepModel
            self.timeStepModel = TimeStepModel(self.case)

    #---------------------------------------------------------------------------

    def initNotebookModel(self):
        """
        Initialize the notebook model.
        """

        if self.notebookModel is None:
            from code_saturne.model.NotebookModel import NotebookModel
            self.notebookModel = NotebookModel(self.case)

    #---------------------------------------------------------------------------

    def initNumParamModel(self):
        """
        Initialize the numerical parameters model
        """
        if self.numParamModel is None:
            from code_saturne.model.NumericalParamGlobalModel import NumericalParamGlobalModel
            self.numParamModel = NumericalParamGlobalModel(self.case)

    #---------------------------------------------------------------------------

    def initNumParamEquationModel(self):
        """
        Initialize the numerical parameters equation model
        """
        if self.numParamEquationModel is None:
            from code_saturne.model.NumericalParamEquationModel import NumericalParamEquationModel
            self.numParamEquationModel = NumericalParamEquationModel(self.case)

    #---------------------------------------------------------------------------

    def cleanMeshList(self):
        """
        Empty mesh list.
        """
        self.initMeshModel()

        for m in self.meshModel.getMeshList():
            self.meshModel.delMesh(m)

    #---------------------------------------------------------------------------

    def addMesh(self, mesh, path=None):
        """
        Add a mesh to the xml file.
        @param mesh : name of the mesh file
        @param path : path to the mesh file if file is not in MESH folder.
        """
        self.initMeshModel()

        self.meshModel.addMesh((mesh, path))

    #---------------------------------------------------------------------------

    def setMeshInput(self, mesh_input):
        """
        Define a mesh input as mesh to use.
        @param mesh_input : name/path of mesh_input file
        """
        self.initMeshModel()

        self.cleanMeshList()
        self.meshModel.setMeshInput(mesh_input)

    #---------------------------------------------------------------------------

    def rotateMesh(self, rotationAngle):
        """
        Define a mesh rotation operation
        @param rotationAngle : angle of rotation.
        """
        self.initMeshModel()

        self.meshModel.setRotationAngle(0, rotationAngle)

    #---------------------------------------------------------------------------

    def defineProbesFromCSV(self, probes_file):
        """
        Method to add probes based on a csv file
        @param probes_file: file containing the coordinates of the probes
        to add.
        """
        self.initOutputModel()

        self.outputModel.ImportProbesFromCSV(probes_file)

    #---------------------------------------------------------------------------

    def defineProbeFromCoordinates(self, x, y, z):
        """
        Add a probe using coordinates.
        @param x: first coordinate
        @param y: second coordinate
        @param z: third coordinate
        """
        self.initOutputModel()

        self.outputModel.addMonitoringPoint(x, y, z)

    #---------------------------------------------------------------------------

    def removeExistingProbes(self):
        """
        Remove all existing probes in an xml file
        """
        self.initOutputModel()

        nProbes = self.outputModel.getNumberOfMonitoringPoints()

        for i in range(nProbes):
            self.outputModel.deleteMonitoringPoint(str(nProbes-i))

    #---------------------------------------------------------------------------

    def setProbesFrequency(self, time_freq=None, iter_freq=None):
        """
        Modify the output frequency of monitoring points.
        @param time_freq: set physical time frequency in [s]
        @param iter_freq: set iteration frequency in [iter]
        """
        if time_freq and iter_freq:
            raise Exception("You cannot define both time and iteration frequency")

        elif time_freq:
            self.outputModel.setMonitoringPointType('Frequency_h_x')
            self.outputModel.setMonitoringPointFrequencyTime(time_freq)

        elif iter_freq:
            self.outputModel.setMonitoringPointType('Frequency_h')
            self.outputModel.setMonitoringPointFrequency(iter_freq)

        else:
            raise Exception("You provided no frequency")

    #---------------------------------------------------------------------------

    def setBcType(self, bc_label, new_type):
        """
        Set boundary condition type.
        @param bc_label: name of the boundary condition
        @param new_type: boundary condition type.
                         Examples: 'inlet', 'wall'
        """

        self.initBcModel()

        self.bcModel.setNature(bc_label, new_type)

    #---------------------------------------------------------------------------

    def setBcLocalization(self, bc_label, localization):
        """
        Set selection criteria for the boundary definition.
        @param bc_label: name of the boundary condition.
        @param localization : selection criteria.
        """

        self.initBcModel()

        self.bcModel.setLocalization(bc_label, localization)

    #---------------------------------------------------------------------------

    def getBoundary(self, bc_label, bc_type):
        """
        Get boundary object based on name and type.
        @param bc_label: name of the boundary condition
        @param bc_type : boundary condition type (inlet, wall)
        """

        self.initBndModel()

        bnd = self.bndModel(bc_type, bc_label, self.case)

        return bnd

   #---------------------------------------------------------------------------

    def setInletVelocity(self, bc_label, vel):
        """
        Set velocity value on an inlet.
        @param bc_label: name of the boundary condition
        @param vel     : value to set
        """

        bnd = self.getBoundary(bc_label, "inlet")

        bnd.setVelocityChoice("norm")
        bnd.setVelocity(vel)
        bnd.setDirectionChoice("normal")

        pass

    #---------------------------------------------------------------------------

    def setInletScalar(self, bc_label, scalar_name, val):
        """
        Set scalar injection value on an inlet
        @param bc_label   : name of the boundary condition
        @param scalar_name: scalar name
        @param val        : value to set
        """

        bnd = self.getBoundary(bc_label, "inlet")

        bnd.setScalarChoice(scalar_name, "dirichlet")
        bnd.setScalarValue(scalar_name, "dirichlet", val)

        pass

    #---------------------------------------------------------------------------

    def setRestartPath(self, checkpoint_path, diff_restart_mesh=False):
        """
        Set restart path in xml file.
        @param path: path to checkpoint folder
        """

        self.initRestartModel()

        # If restart is done using a new mesh which is in the checkpoint folder
        if diff_restart_mesh:
            _mp = os.path.join(checkpoint_path, "mesh_input.csm")
            self.setDifferentRestartMesh(_mp)

        self.restartModel.setRestartPath(checkpoint_path)
        pass


    #---------------------------------------------------------------------------

    def setDifferentRestartMesh(self, restart_mesh_path):
        """
        Set a different mesh for restart data than the one used for current
        computation.
        @param restart_mesh_path: path to original mesh
        """

        self.initRestartModel()

        self.restartModel.setRestartMeshPath(restart_mesh_path)

        pass
    #---------------------------------------------------------------------------

    def removeRestartPath(self):
        """
        Remove restart option from xml file
        """

        self.initRestartModel()

        self.restartModel.setRestartPath(None)
        self.restartModel.setRestartMeshPath(None)
        pass

    #---------------------------------------------------------------------------

    def setTimeIterationsNumber(self, n_iter):
        """
        Set simulation stop criterion to a given number of iterations.
        @param n_iter: number of iterations to compute.
        """
        self.initTimeStepModel()

        self.timeStepModel.setStopCriterion('iterations_add', n_iter)

    #---------------------------------------------------------------------------

    def setMaxTime(self, t_max):
        """
        Set simulation stop criterion to a given physical time.
        @param t_max: physical time to achieve.
        """
        self.initTimeStepModel()

        self.timeStepModel.setStopCriterion('maximum_time', t_max)

    #---------------------------------------------------------------------------

    def setTimeStep(self, dt):
        """
        Set the computation time step value.
        @param dt: time step value
        """
        self.initTimeStepModel()

        self.timeStepModel.setTimeStep(dt)

    #---------------------------------------------------------------------------

    def setUnsteadySimulation(self):
        """
        Define unsteady simulation.
        """
        self.initTimeStepModel()

        self.timeStepModel.setTimePassing(0)

    #---------------------------------------------------------------------------

    def changeNotebookParameter(self, name, val):
        """
        Change the value of a notebook parameter.
        @param name : name of the parameter
        @param val  : new value to assign
        """
        self.initNotebookModel()

        if name not in self.notebookModel.getVarNameList():
            raise Exception("%s is not in the notebook!" % name)

        self.notebookModel.setVariableValue(val=val, var=name)

    #---------------------------------------------------------------------------

    def setGradientReconstruction(self, imrgra):
        """
        Set gradient reconstruction method.
        @param imrgra : Method indicator (int)
        """
        self.initNumParamModel()

        self.numParamModel.setGradientReconstruction(imrgra)

    #---------------------------------------------------------------------------

    def setBlendingFactor(self, varName, blencv):
        """
        Set the blending factor for a given variable.
        @param varName : name of the variable
        @param blencv  : Blending factor value
        """

        self.initNumParamEquationModel()

        self.numParamEquationModel.setBlendingFactor(varName, blencv)

    #---------------------------------------------------------------------------

    def find_nearest_checkpoint(self, search_path, restart_time):
        """
        Find inside a checkpoint folder the restart file which time
        is closest (and smaller) to a given value.
        @param search_path : path to checkpoint folder
        @param restart_time: restart time which is sought for
        """

        from code_saturne.model.StartRestartModel import getRestartInfo
        from glob import glob

        # Check where we are before launching
        origin = os.getcwd()

        # Go to case folder
        os.chdir(search_path)

        # Get list of available checkpoints
        spaths = glob("previous_dump_*")
        spaths.sort()
        spaths.append(".")

        # initialize delta and closest checkpoint
        delta = 100000000000.0
        checkpoint = None

        # Loop on all checkpoint dumps looking for closest dump done with a time
        # less or equal to the desired restart time
        for p in spaths:
            ret = getRestartInfo(package=self.pkg, restart_path=p)

            time_p = ret[2]
            if time_p <= restart_time:
                dt = restart_time - time_p
                if dt < delta:
                    checkpoint = p
                    delta      = dt

        # Return path to the found checkpoint. None if none found
        restart_checkpoint = None
        if checkpoint:
            restart_checkpoint = os.path.join(search_path, checkpoint)

        os.chdir(origin)

        return restart_checkpoint

#-------------------------------------------------------------------------------
# Function which will update case in memory based on the case_setup_filter class
#-------------------------------------------------------------------------------

def update_case_model(case, options, pkg):
    """
    Update a Case model in memory with given options
    """

    xml_controller = case_setup_filter(case, pkg=pkg)

    # Mesh parameters
    # ---------------

    if options.mesh:
        xml_controller.cleanMeshList()
        for m in options.mesh:
            xml_controller.addMesh(m, None)

    if options.mesh_input:
        xml_controller.setMeshInput(options.mesh_input)

    if options.rotationAngle:
        xml_controller.rotateMesh(options.rotationAngle)

    # Boundary conditions parameters
    # ------------------------------

    if options.bc_criteria:
        for bc in options.bc_criteria:
            bc_label, bc_criteria = bc.split(':')
            xml_controller.setBcLocalization(bc_label, bc_criteria)

    # Time parameters
    # ---------------

    if options.iterationsNumber:
        xml_controller.setTimeIterationsNumber(options.iterationsNumber)

    if options.tmax:
        xml_controller.setMaxTime(options.tmax)

    if options.TimeStep:
        xml_controller.setTimeStep(options.TimeStep)

    if options.TimeModel:
        xml_controller.setUnsteadySimulation()

    # Numerical parameters
    # --------------------

    if options.imrgra:
        xml_controller.setGradientReconstruction(options.imrgra)

    if options.blencv:
        for elt in options.blencv:
            fname, factor_str = elt.split(':')
            factor = float(factor_str)
            xml_controller.setBlendingFactor(fname, factor)

    # Notebook
    # --------

    if options.notebook:
        for var in options.notebook:
            vname, value = var.split(':')
            xml_controller.changeNotebookParameter(vname, value)

    # Restart parameters
    # ------------------

    if options.RestartRun:
        restart_path = os.path.join('RESU', options.RestartRun, 'checkpoint')
        xml_controller.setRestartPath(restart_path, options.DiffRestartMesh)

    if options.DiffRestartMeshPath:
        xml_controller.setDifferentRestartMesh(options.DiffRestartMeshPath)

#-------------------------------------------------------------------------------
# Main function which modifies the case setup based on given argument list
#-------------------------------------------------------------------------------

def update_case_setup(case, args, pkg):

    # Read input arguments
    # --------------------

    if args == []:
        args.append("--help")

    parser = arg_parser(args)
    (options, rargs) = parser.parse_known_args(args)

    if rargs:
        print("cs_parametric_setup: unknown arguments", rargs)
        print("-------------------")

    if case:
        xml_controller = case_setup_filter(case, pkg=pkg)
        update_case_model(xml_controller.case, options, pkg)

    else:
        print("parsed options:", options)

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    retval = update_case_setup(argv=sys.argv[1:])

    sys.exit(retval)

#-------------------------------------------------------------------------------
