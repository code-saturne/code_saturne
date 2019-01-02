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

import sys,os
from omniORB import CORBA
from FSI_ASTER_module import FSI_ASTER

if __name__ == '__main__':

  print(sys.argv)
  orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
  poa = orb.resolve_initial_references("RootPOA")
  print("ORB and POA initialized " + str(orb) + ' ' + str(poa))
  sys.stdout.flush()
  sys.stderr.flush()

  container=orb.string_to_object(os.getenv("SALOME_CONTAINER"))
  containerName=os.getenv("SALOME_CONTAINERNAME")
  instanceName=os.getenv("SALOME_INSTANCE")

  compo=FSI_ASTER(orb,poa,container,containerName, instanceName, "FSI_ASTER")
  comp_o = compo._this()
  comp_iors = orb.object_to_string(comp_o)
  print("ior aster " + str(comp_iors))

  sys.stdout.flush()
  sys.stderr.flush()

  #activate the POA
  poaManager = poa._get_the_POAManager()
  poaManager.activate()

  orb.run()
  print("fin du composant aster standalone")

