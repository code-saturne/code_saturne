#-------------------------------------------------------------------------------
#   This file is part of the Code_Saturne CFD tool.
#
#   Copyright (C) 2011-2019 EDF S.A.
#
#   The Code_Saturne CFD tool is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne CFD tool is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Preprocessor; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import os, sys

from module_generator import Generator, Module, Service
from module_generator import ASTERComponent, CPPComponent

# Context definition
# ------------------

asterdir = os.path.expanduser("/home/aster/NEW11")

# SALOME environment be must loaded beforehand

KERNEL_ROOT_DIR = os.getenv("KERNEL_ROOT_DIR")
if KERNEL_ROOT_DIR is None:
    print("KERNEL_ROOT_DIR must be defined.\nPlease load SALOME environment.\n")
    sys.exit(1)

context = {'update':1,
           "prerequisites":os.path.join(os.path.dirname(KERNEL_ROOT_DIR),
                                        'prerequis-V6.3.1.sh'),
           "kernel":KERNEL_ROOT_DIR}

# Definition of all the in/out streams
# ------------------------------------

nbIter = ("NBPDTM", "CALCIUM_integer","I")
nbSubIter = ("NBSSIT", "CALCIUM_integer","I")

indSync =("ISYNCP", "CALCIUM_integer","I")
frequency = ("NTCHRO", "CALCIUM_integer","I")

timeStep = ("PDTREF", "CALCIUM_double","I")
prevTime = ("TTINIT", "CALCIUM_double","I")
epsilon = ("EPSILO", "CALCIUM_double","I")

newTimeStep = ("DTCALC", "CALCIUM_double","I")
nbDyn = ("NB_DYN", "CALCIUM_integer","I")
nbFor = ("NB_FOR", "CALCIUM_integer","I")

nodesColor = ("COLNOD", "CALCIUM_integer","I")
facesColor = ("COLFAC", "CALCIUM_integer","I")
nodesXYZ = ("COONOD", "CALCIUM_double","I")
facesXYZ = ("COOFAC", "CALCIUM_double","I")

astForces = ("FORAST", "CALCIUM_double","I")
astConvergence = ("ICVAST", "CALCIUM_integer","I")

astTimeStep = ("DTAST", "CALCIUM_double","I")
astDisplacement = ("DEPAST", "CALCIUM_double","I")
astVelocity = ("VITAST", "CALCIUM_double","I")

# Code_Aster streams

aster_outstream = [astDisplacement, astVelocity, astTimeStep]

aster_instream = [nbFor, nbDyn,
                  nodesXYZ, facesXYZ,
                  nodesColor, facesColor,
                  astForces,
                  nbIter, nbSubIter,
                  epsilon, astConvergence,
                  indSync, frequency,
                  prevTime, timeStep, newTimeStep]

# Code_Saturne streams

saturne_instream = [astDisplacement, astVelocity, astTimeStep]

saturne_outstream = [nbFor, nbDyn,
                     nodesXYZ, facesXYZ,
                     nodesColor, facesColor,
                     astForces,
                     nbIter, nbSubIter,
                     epsilon, astConvergence,
                     indSync, frequency,
                     prevTime, timeStep, newTimeStep]

# Creation of the different components
# ------------------------------------

# code_aster component

aster_service = Service("op0117",
                        instream = aster_instream,
                        outstream = aster_outstream)

c1 = ASTERComponent("FSI_ASTER",
                    services = [aster_service],
                    aster_dir = asterdir,
                    kind = "exe",
                    exe_path = "./aster_by_yacs.sh")

# Code_Saturne component

saturne_defs = """extern "C" {
  void cs_calcium_set_component(int comp_id, void *comp);
  void cs_calcium_set_verbosity(int n_echo);
  void cs_run(void);
}"""

saturne_body = """cs_calcium_set_component(0, component);
cs_calcium_set_verbosity(verbosity);
cs_run();"""

saturne_service = Service("run",
                          inport = [("app_name", "string"),
                                    ("verbosity",  "long")],
                          outport = [("retval", "long")],
                          instream  = saturne_instream,
                          outstream = saturne_outstream,
                          defs = saturne_defs,
                          body = saturne_body)

c2 = CPPComponent("FSI_SATURNE",
                  services = [saturne_service],
                  kind = "exe",
                  exe_path = "./run_solver")

# Creation of the FSI module
# --------------------------

m = Module("FSI", components=[c1,c2,c3], prefix="fsi_arch")
g = Generator(m, context)

g.generate()

# We don't want to compile and install the module by default
if False:
    g.bootstrap()
    g.configure()
    g.make()
    g.install()
    g.make_appli("fsi_appli", restrict=["KERNEL","GUI","YACS"])
