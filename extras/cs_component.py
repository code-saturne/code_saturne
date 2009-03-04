import os
from module_generator import Generator,Module,Service,CPPComponent

context={'update':1,
         "prerequisites":"/home/salome/SALOME5/V5_1_0/envSalome-V5_1_0.sh",
         "kernel":"/home/salome/SALOME5/V5_1_0/KERNEL_V5_1_0",
        }

cwd=os.getcwd()

defs="""
void cs_calcium_set_component(int comp_id, void *comp);
void cs_run(void);
"""

body="""
cs_calcium_set_component(0, component);
cs_run();
"""
c1 = CPPComponent("CFD_Single",
                  services = [Service("run",
                                      inport=[("exec_dir", "string"),
                                              ("library", "string"),
                                              ("args", "string")],
                                      outport=[("retval", "long")],
                                      defs=defs, body=body,
                                      ),
                              ],
                  includes="-I/usr/include",
                  kind="exe",
                  exe_path=os.path.join(cwd, "cs14.exe"),
                  )

c2 = CPPComponent("CFD_Aster",
                  services = [Service("run",
                                      inport=[("exec_dir", "string"),
                                              ("library", "string"),
                                              ("args", "string")],
                                      outport=[("retval", "long")],
                                      instream =[("depsat","CALCIUM_double","I"),
                                                 ("epsilo","CALCIUM_double","I"),
                                                 ("dtcalc","CALCIUM_double","I"),
                                                 ("ttinit","CALCIUM_double","I"),
                                                 ("pdtref","CALCIUM_double","I"),
                                                 ("nbpdtm","CALCIUM_integer","I"),
                                                 ("nbssit","CALCIUM_integer","I"),
                                                 ("isyncp","CALCIUM_integer","I"),
                                                 ("ntchro","CALCIUM_integer","I"),
                                                 ("icvext","CALCIUM_integer","I")],
                                      outstream=[("dtsat", "CALCIUM_double","I"),
                                                 ("forsat","CALCIUM_double","I"),
                                                 ("almaxi","CALCIUM_double","I"),
                                                 ("coonod","CALCIUM_double","I"),
                                                 ("coofac","CALCIUM_double","I"),
                                                 ("icv",   "CALCIUM_integer","I"),
                                                 ("dongeo","CALCIUM_integer","I"),
                                                 ("colnod","CALCIUM_integer","I"),
                                                 ("colfac","CALCIUM_integer","I")],
                                      defs=defs, body=body,
                                      ),
                              ],
                  includes="-I/usr/include",
                  kind="exe",
                  exe_path=os.path.join(cwd, "cs14.exe"),
                  )

g=Generator(Module("CFDRUN", components=[c1, c2], prefix="./install"), context)
g.generate()
g.bootstrap()
g.configure()
g.make()
g.install()
g.make_appli("appli", restrict=["KERNEL", "GUI", "YACS"])

