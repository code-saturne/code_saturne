import sys,os
from omniORB import CORBA
from FSI_ASTER_module import FSI_ASTER

if __name__ == '__main__':

  print sys.argv
  orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
  poa = orb.resolve_initial_references("RootPOA")
  print "ORB and POA initialized",orb,poa
  sys.stdout.flush()
  sys.stderr.flush()

  container=orb.string_to_object(os.getenv("SALOME_CONTAINER"))
  containerName=os.getenv("SALOME_CONTAINERNAME")
  instanceName=os.getenv("SALOME_INSTANCE")

  compo=FSI_ASTER(orb,poa,container,containerName, instanceName, "FSI_ASTER")
  comp_o = compo._this()
  comp_iors = orb.object_to_string(comp_o)
  print "ior aster",comp_iors

  sys.stdout.flush()
  sys.stderr.flush()

  #activate the POA
  poaManager = poa._get_the_POAManager()
  poaManager.activate()

  orb.run()
  print "fin du composant aster standalone"

