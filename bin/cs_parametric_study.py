#!/usr/bin/env python

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2021 EDF S.A.
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
This module
"""

# ------------------------------------------------------------------------------
# Import required Python modules
# ------------------------------------------------------------------------------

import os, sys
from optparse import OptionParser

# ------------------------------------------------------------------------------
# Function used to process the command line
# ------------------------------------------------------------------------------
def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """

    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("--update", dest="updateXml", default = False,
                      action= 'store_true',
                      help="Update the xml file.")

    parser.add_option("--create-from-ref", dest="createFromRef", type="string",
                      help="Reference case to create a copy of")

    parser.add_option("-c", "--case", dest="case", type="string",
                      help="Directory of the current case")

    parser.add_option("-p", "--param", dest="param", type="string",
                      help="Name f the file of parameters")

    parser.add_option("--notebook", "--nb", dest="notebook", type="string",
                      action="append",
                      help="Notebook parameters. Format is <var>:<val>")

    parser.add_option("-m", "--mesh", dest="mesh", type="string",
                      action="append",
                      help="Name of the new mesh")

    parser.add_option("--mi", "--mesh_input", dest="mesh_input", type="string",
                      help="Name of the new mesh_input")

    parser.add_option("-a", "--perio-angle", dest="rotationAngle", type="float",
                      help="Periodicity angle")

    parser.add_option("-r", "--restart", dest="RestartRun", type="string",
                      help="Run to restart from (in same case)")

    parser.add_option("-n", "--iter-num", dest="iterationsNumber", type="int",
                      help="New iteration number")

    parser.add_option("--tmax", dest="tmax", type="float",
                      help="Time based stop criterion")

    parser.add_option("--iter-dt", dest="TimeStep", type="float",
                      help="New time step")

    parser.add_option("-u", "--unsteady", dest="TimeModel", action="store_true",
                      help="Unsteady time model (steady by default)")

    parser.add_option("--imrgra", dest="imrgra", type="int",
                      help="Gradient reconstruction")

    parser.add_option("--blencv", dest="blencv", type="float",
                      help="Blencv")

    (options, args) = parser.parse_args(argv)

    if args != []:
        print(" -------------------------- ")
        print(" Ignored arguments : ")
        print(" ",args)
        print(" -------------------------- ")

    return options

# ------------------------------------------------------------------------------
# Class which provides a set of predefined methods to modify a
# code_saturne xml file
# ------------------------------------------------------------------------------
class cs_modify_xml(object):
    """
    Toy class to modify a Code_Saturne xml file
    """

    # ---------------------------
    def __init__(self, xml_file, pkg=None):
        """
        Init method.
        @param xml_file: path to xml file
        """

        # package
        self.pkg = None
        if pkg:
            self.pkg = pkg
        else:
            from code_saturne.cs_package import package
            self.pkg = package()

        from code_saturne.model.XMLengine import Case

        try:
            case = Case(package = pkg, file_name = filepath)
        except:
            print("Error while reading parameters files.")
            print("This file is not in accordance with XML specifications.")
            sys.exit(1)

        if self.pkg.name == 'code_saturne':
            from code_saturne.model.XMLinitialize import XMLinit
        else:
            from model.XMLinitializeNeptune import XMLinitNeptune as XMLinit

        self.xml = xml_file

        # internal functions of code_saturne to properly initialize the
        # xml structure
        self.case = Case(file_name=xml_file)
        self.case['xmlfile'] = xml_file
        self.case.xmlCleanAllBlank(self.case.xmlRootNode())
        XMLinit(self.case).initialize()

        self.outputModel   = None
        self.bcModel       = None
        self.bndModel      = None
        self.restartModel  = None
        self.timeStepModel = None
        self.meshModel     = None
        self.notebookModel = None
        self.numParamModel = None
    # ---------------------------

    # ---------------------------
    def saveXml(self):
        """
        Save the xml file
        """

        self.case.xmlSaveDocument()
    # ---------------------------

    # ---------------------------
    def initMeshModel(self):
        """
        Initialize the SolutionDomainModel
        """

        if self.meshModel == None:
            from code_saturne.model.SolutionDomainModel import SolutionDomainModel
            self.meshModel = SolutionDomainModel(self.case)
    # ---------------------------

    # ---------------------------
    def initOutputModel(self):
        """
        Method to initialize the outputModel data structure which handles
        output parameters including monitoring points
        """
        if self.outputModel == None:
            from code_saturne.model.OutputControlModel import OutputControlModel
            self.outputModel = OutputControlModel(self.case)
    # ---------------------------

    # ---------------------------
    def initBcModel(self):
        """
        Initialize the BC model
        """

        if self.bcModel == None:
            from code_saturne.model.LocalizationModel import LocalizationModel
            self.bcModel = LocalizationModel("BoundaryZone", self.case)
    # ---------------------------

    # ---------------------------
    def initBndModel(self):
        """
        Initialize the boundary model.
        """

        if self.bndModel == None:
            from code_saturne.model.Boundary import Boundary
            self.bndModel = Boundary
    # ---------------------------

    # ---------------------------
    def initRestartModel(self):
        """
        Initialize the restart model
        """

        if self.restartModel == None:
            from code_saturne.model.StartRestartModel import StartRestartModel
            self.restartModel = StartRestartModel(self.case)
    # ---------------------------

    # ---------------------------
    def initTimeStepModel(self):
        """
        Initialize the time step model.
        """

        if self.timeStepModel == None:
            from code_saturne.model.TimeStepModel import TimeStepModel
            self.timeStepModel = TimeStepModel(self.case)
    # ---------------------------

    # ---------------------------
    def initNotebookModel(self):
        """
        Initialize the notebook model.
        """

        if self.notebookModel == None:
            from code_saturne.model.NotebookModel import NotebookModel
            self.notebookModel = NotebookModel(self.case)
    # ---------------------------

    # ---------------------------
    def initNumParamModel(self):
        """
        Initialize the numerical parameters model
        """
        if self.numParamModel == None:
            from code_saturne.model.NumericalParamGlobalModel import NumericalParamGlobalModel
            self.numParamModel = NumericalParamGlobalModel(self.case)
    # ---------------------------

    # ---------------------------
    def cleanMeshList(self):
        """
        Empty mesh list.
        """
        self.initMeshModel()

        for m in self.meshModel.getMeshList():
            self.meshModel.delMesh(m)
    # ---------------------------

    # ---------------------------
    def addMesh(self, mesh, path=None):
        """
        Add a mesh to the xml file.
        @param mesh : name of the mesh file
        @param path : path to the mesh file if file is not in MESH folder.
        """
        self.initMeshModel()

        self.meshModel.addMesh((mesh, path))
    # ---------------------------

    # ---------------------------
    def setMeshInput(self, mesh_input):
        """
        Define a mesh input as mesh to use.
        @param mesh_input : name/path of mesh_input file
        """
        self.initMeshModel()

        self.cleanMeshList()
        self.meshModel.setMeshInput(mesh_input)
    # ---------------------------

    # ---------------------------
    def rotateMesh(self, rotationAngle):
        """
        Define a mesh rotation operation
        @param rotationAngle : angle of rotation.
        """
        self.initMeshModel()

        self.meshModel.setRotationAngle(0, rotationAngle)
    # ---------------------------

    # ---------------------------
    def defineProbesFromCSV(self, probes_file):
        """
        Method to add probes based on a csv file
        @param probes_file: file containing the coordinates of the probes
        to add.
        """
        self.initOutputModel()

        self.outputModel.ImportProbesFromCSV(probes_file)
    # ---------------------------

    # ---------------------------
    def defineProbeFromCoordinates(self, x, y, z):
        """
        Add a probe using coordinates.
        @param x: first coordinate
        @param y: second coordinate
        @param z: third coordinate
        """
        self.initOutputModel()

        self.outputModel.addMonitoringPoint(x, y, z)
    # ---------------------------

    # ---------------------------
    def removeExistingProbes(self):
        """
        Remove all existing probes in an xml file
        """
        self.initOutputModel()

        nProbes = self.outputModel.getNumberOfMonitoringPoints()

        for i in range(nProbes):
            self.outputModel.deleteMonitoringPoint(str(nProbes-i))
    # ---------------------------

    # ---------------------------
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
    # ---------------------------

    # ---------------------------
    def setBcType(self, bc_label, new_type):
        """
        Set boundary condition type.
        @param bc_label: name of the boundary condition
        @param new_type: boundary condition type.
                         Examples: 'inlet', 'wall'
        """

        self.initBcModel()

        self.bcModel.setNature(bc_label, new_type)
    # ---------------------------

    # ---------------------------
    def setBcLocalization(self, bc_label, localization):
        """
        Set selection criteria for the boundary definition.
        @param bc_label: name of the boundary condition.
        @param localization : selection criteria.
        """

        self.initBcModel()

        self.bcModel.setLocalization(bc_label, localization)
    # ---------------------------

    # ---------------------------
    def getBoundary(self, bc_label, bc_type):
        """
        Get boundary object based on name and type.
        @param bc_label: name of the boundary condition
        @param bc_type : boundary condition type (inlet, wall)
        """

        self.initBndModel()

        bnd = self.bndModel(bc_type, bc_label, self.case)

        return bnd

    # ---------------------------

    # ---------------------------
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
    # ---------------------------

    # ---------------------------
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
    # ---------------------------

    # ---------------------------
    def setRestartPath(self, checkpoint_path):
        """
        Set restart path in xml file.
        @param path: path to checkpoint folder
        """

        self.initRestartModel()

        self.restartModel.setRestartPath(checkpoint_path)
        _mp = os.path.join(checkpoint_path, 'mesh_input.csm')
        self.restartModel.setRestartMeshPath(_mp)
        pass
    # ---------------------------

    # ---------------------------
    def removeRestartPath(self):
        """
        Remove restart option from xml file
        """

        self.initRestartModel()

        self.restartModel.setRestartPath(None)
        self.restartModel.setRestartMeshPath(None)
        pass
    # ---------------------------

    # ---------------------------
    def setTimeIterationsNumber(self, n_iter):
        """
        Set simulation stop criterion to a given number of iterations.
        @param n_iter: number of iterations to compute.
        """
        self.initTimeStepModel()

        self.timeStepModel.setStopCriterion('iterations_add', n_iter)
    # ---------------------------

    # ---------------------------
    def setMaxTime(self, t_max):
        """
        Set simulation stop criterion to a given physical time.
        @param t_max: physical time to achieve.
        """
        self.initTimeStepModel()

        self.timeStepModel.setStopCriterion('maximum_time', t_max)
    # ---------------------------

    # ---------------------------
    def setTimeStep(self, dt):
        """
        Set the computation time step value.
        @param dt: time step value
        """
        self.initTimeStepModel()

        self.timeStepModel.setTimeStep(dt)
    # ---------------------------

    # ---------------------------
    def setUnsteadySimulation(self):
        """
        Define unsteady simulation.
        """
        self.initTimeStepModel()

        self.timeStepModel.setTimePassing(0)
    # ---------------------------

    # ---------------------------
    def changeNotebookParameter(self, name, val):
        """
        Change the value of a notebook parameter.
        @param name : name of the parameter
        @param val  : new value to assign
        """
        self.initNotebookModel()

        if name not in self.notebookModel.getVarNamesList():
            raise Exception("%s is not in the notebook!" % name)

        self.notebookModel.setVariableValue(val=val, var=name)
    # ---------------------------

    # ---------------------------
    def setGradientReconstruction(self, imrgra):
        """
        Set gradient reconstruction method.
        @param imrgra : Method indicator (int)
        """
        self.initNumParamModel()

        self.numParamModel.setGradientReconstruction(imrgra)
    # ---------------------------

    # ---------------------------
    def setBlendingFactor(self, varName, blencv):
        """
        Set the blending factor for a given variable.
        @param varName : name of the variable
        @param blencv  : Blending factor value
        """
        self.initNumParamModel()

        self.numParamModel.setBlendingFactor(varName, blencv)
    # ---------------------------

    # ---------------------------
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
        # less or equal to the wished restart time
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
    # ---------------------------


# ------------------------------------------------------------------------------
# Class which allows the execution of code_saturne script commands
# ------------------------------------------------------------------------------

class cs_exec(object):

    def __init__(self, pkg=None):

        # package
        self.pkg = None
        if pkg:
            self.pkg = pkg
        else:
            from code_saturne.cs_package import package
            self.pkg = package()

    # ---------------------------
    def run_cs_command(self, args):
        """
        Execute code_saturne with a given command line.
        @param args: list of arguments composing the command line
        """

        from code_saturne.cs_script import master_script
        script = master_script(argv=args, package=self.pkg)

        script.execute()
        pass
    # ---------------------------

    # ---------------------------
    def run_cs_case(self, case_dir, xmlfile="setup.xml", nprocs=1):
        """
        Run a code_saturne case.
        @param case_dir: path to the case folder
        @param xmlfile : name of the xml file (default is setup.xml)
        @param nprocs  : number of cores to uses (default is 1)
        """

        # Check where we are before launching
        origin = os.getcwd()

        # Go to case folder
        os.chdir(os.path.join(case_dir, "DATA"))

        # Set the arguments list
        args = ["run", "--param", xmlfile, "--nprocs", str(nprocs)]

        self.run_cs_command(args)

        os.chdir(origin)

    # ---------------------------

    # ---------------------------
    def generate_new_case(self, new_case=None, ref_case=None):
        """
        Generate a new case based on a reference case.
        @param new_case: name of new case
        @param ref_case: path to reference case
        """

        args = ["create", "-c"]

        # Set case name
        if new_case:
            args.append(new_case)
        else:
            args.append("CASE1")

        # Set ref case if needed
        if ref_case:
            args.extend(["--copy-from", ref_case])

        self.run_cs_command(args)
    # ---------------------------


# ------------------------------------------------------------------------------
# Function which will update an xml file based on the cs_modify_xml class
# ------------------------------------------------------------------------------

def update_xml_file(pkg, filepath, options):

    # Sanity check
    if not os.path.isfile(filepath):
        print("File %s does not exist!" % filepath)
        sys.exit(1)

    xml_controller = cs_modify_xml(filepath, pkg=pkg)

    # --------------------
    # Mesh parameters

    if options.mesh != []:
        xml_controller.cleanMeshList()
        for m in options.mesh:
            xml_controller.addMesh(m, None)

    if options.mesh_input:
        xml_controller.setMeshInput(options.mesh_input)

    if options.rotationAngle:
        xml_controller.rotateMesh(options.rotationAngle)
    # --------------------

    # --------------------
    # Time parameters
    if options.iterationsNumber:
        xml_controller.setTimeIterationsNumber(options.iterationsNumber)

    if options.tmax:
        xml_controller.setMaxTime(options.tmax)

    if options.TimeStep:
        xml_controller.setTimeStep(options.TimeStep)

    if options.timeModel:
        xml_controller.setUnsteadySimulation()
    # --------------------

    # --------------------
    # Numerical parameters
    if options.imrgra:
        xml_controller.setGradientReconstruction(options.imrgra)

    if (options.blencv):
        for elt in options.blencv:
            fname, factor = elt.split(':')
            xml_controller.setBlendingFactor(fname, factor)
    # --------------------

    # --------------------
    # Notebook
    if len(options.notebook) > 0:
        for var in options.notebook:
            vname, value = var.split(':')
            xml_controller.changeNotebookParameter(vname, value)
    # --------------------

    # --------------------
    # Restart parameters
    if options.RestartRun:
        r = StartRestartModel(case)
        restart_path = os.path.join('RESU', options.RestartRun, 'checkpoint')

        xml_controller.setRestartPath(restart_path)
    # --------------------

    case.xmlSaveDocument()


# ------------------------------------------------------------------------------
# Run main function which modifies the case parameters
# ------------------------------------------------------------------------------

def main(args, pkg):

    # ---------------------------
    # Read input arguments
    # ---------------------------
    if args == []:
        args.append("--help")
    opts = process_cmd_line(args)

    # ---------------------------
    # Check if needs to create a case from a reference one
    # ---------------------------
    if opts.createFromRef:
        cmd_mgr = cs_exec(pkg)
        cmd_mgr.generate_new_case(opts.case, opts.createFromRef)

        # Set default name if needed
        if opts.case == None:
            opts.case = "CASE1"

    # ---------------------------
    # Sanity check of case or xml in input parameters
    # ---------------------------
    if opts.case == None and opts.param == None:
        print("No case name nor parameters file provided")
        sys.exit(1)

    fp = ""
    if opts.case:
        fp = os.path.join(opts.case, "DATA")
        if opts.param:
            fp = os.path.join(fp, opts.param)
        else:
            fp = os.path.join(fp, "setup.xml")

    elif opts.param:
        fp = opts.param

    if fp == "":
        print("No case name nor parameters file provided")
        sys.exit(1)

    # ---------------------------
    # Update xml if needed
    # ---------------------------
    if opts.updateXml:
        update_xml_file(pkg, fp, opts)



if __name__=="__main__":

    pkg = None
    args = None

    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        args = ["--help"]

    main(args=args, pkg=pkg)
# ------------------------------------------------------------------------------
# END
# ------------------------------------------------------------------------------
