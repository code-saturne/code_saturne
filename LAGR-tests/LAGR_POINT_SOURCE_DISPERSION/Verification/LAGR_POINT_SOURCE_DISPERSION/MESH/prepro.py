#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
import string
import re
from optparse import OptionParser

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """Processes the passed command line arguments."""
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-c", "--case", dest="case", type="string",
                      help="Directory of the current case")

    parser.add_option("-p", "--param", dest="param", type="string",
                      help="Name of the file of parameters")

    parser.add_option("-m", "--mesh", dest="mesh", type="string",
                      help="Name of the new mesh")

    parser.add_option("--mi", "--mesh_input", dest="mesh_input", type="string",
                      help="Name of the new mesh_input")

    parser.add_option("-n", "--iter-num", dest="iterationsNumber", type="int",
                      help="New iteration number")

    parser.add_option("-t", "--iter-dt", dest="TimeStep", type="float",
                      help="New time step")

    parser.add_option("-u", "--unsteady", dest="TimeModel", action="store_true",
                      help="Unsteady time model (steady by default)")

    parser.add_option("-a", "--perio-angle", dest="rotationAngle", type="float",
                      help="Periodicity angle")

    parser.add_option("-i", "--imrgra", dest="imrgra", type="int",
                      help="Gradient reconstruction")

    parser.add_option("-b", "--blencv", dest="blencv", type="float",
                      help="Blencv")

    parser.add_option("--nb", "--notebook", dest="notebook", type="string", action="append",
                      default=[],
                      help="Modify a notebook variable with '--nb <var_name>:<value>'")

    (options, args) = parser.parse_args(argv)

    return options

#-------------------------------------------------------------------------------

def main(options):
    from cs_package import package
    from model.XMLengine import Case
    from model.XMLinitialize import XMLinit
    from model.SolutionDomainModel import SolutionDomainModel
    from model.TimeStepModel import TimeStepModel
    from model.NumericalParamGlobalModel import NumericalParamGlobalModel
    from model.NumericalParamEquationModel import NumericalParamEquationModel
    from model.NotebookModel import NotebookModel
    from studymanager.cs_studymanager_parser import Parser

    fp = os.path.join(options.case, "DATA", options.param)
    if os.path.isfile(fp):
        try:
            case = Case(package = package(), file_name = fp)
        except:
            print("Parameters file reading error.\n")
            print("This file is not in accordance with XML specifications.")
            sys.exit(1)

        case['xmlfile'] = fp
        case.xmlCleanAllBlank(case.xmlRootNode())
        XMLinit(case).initialize()

        s = SolutionDomainModel(case)

        if options.mesh:
            l = s.getMeshList()
            s.delMesh(l[0])
            s.addMesh((options.mesh, None))

        if (options.TimeModel):
            t = TimeStepModel(case)
            t.setTimePassing(0)
        if options.mesh_input:
            print('Changing mesh_input with:', options.mesh_input)
            s.setMeshInput(options.mesh_input)

        if options.rotationAngle:
            s.setRotationAngle(0, options.rotationAngle)

        if (options.iterationsNumber):
            t = TimeStepModel(case)
            t.setStopCriterion('iterations', options.iterationsNumber)

        if (options.TimeStep):
            t = TimeStepModel(case)
            t.setTimeStep(options.TimeStep)

        if (options.imrgra):
            p = NumericalParamGlobalModel(case)
            p.setGradientReconstruction(options.imrgra)

        if (options.blencv):
            p = NumericalParamEquationModel(case)
            for fname in ["velocity", "temperature"]:
                p.setBlendingFactor(fname, options.blencv)

        if len(options.notebook) > 0:
            nb = NotebookModel(case)
            nb_idx   = nb.getVarList()
            nb_names = nb.getVarNameList()
            for var in options.notebook:
                name, value = var.split(':')
                if name not in nb_names:
                    raise Exception("%s is not a notebook variable!"%name)
                else:
                    nb.setVariableValue(val=value, var=name)

        case.xmlSaveDocument()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)

#-------------------------------------------------------------------------------
