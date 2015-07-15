==============================================================
*Code_Saturne* howto: **run a coupling study**
==============================================================

----------------
Introduction
----------------

This tutorial provides a course for a coupling study with *Code_Saturne* and *SYRTHES*.

Run a coupling study from *SALOME_CFD* in graphical mode is not available yet.

----------------
Open SALOME
----------------

We considere, that we use SALOME_CFD to compute a coupling study. Then, to
open salome shell the command is ``~/salome/appli_x_y_z/salome shell``.
Now, we are in salome environment but in a shell. You can, for example, run
command like ``code_saturne -h`` for general option or
``code_saturne info -g user`` for the *Code_Saturne* user's guide.

-------------------------------------
Create study
-------------------------------------

To create the study :
    - go to the main directory for the study. For example, **my_cfd_study** whith the command : ``cd my_cfd_study``.
    - create the study : ``code_sturne create --study=my_study --case=saturne_case --syrthes=syrthes``.


Now, in the directory **my_cfd_study**, a directory **my_study** has been create with five subdirectory
and two files:

     ======================  =========================================================
     Directory / file        Content
     ======================  =========================================================
     MESH                    directory for meshes files
     POST                    directory for post processing files
     RESU_COUPLING           directory for study results
     saturne_case            directory for *Code_Saturne* case (DATA, SRC)
     syrthes                 directory for *SYRTHES*
     coupling_parameters.py  file for coupling parameters (number of processors, ...)
     runcase                 executable file
     ======================  =========================================================

-------------------------------------
Define calculation
-------------------------------------

Use the function **exit** with the command ``exit`` in the salome shell to log out.
Now, we can launch salome in GUI mode, with the command : ``~/salome/appli_x_y_z/salome``.

We are in salome environment and we can :
    - create CAO and meshes,
    - define CFD calculation parameters,
    - define SYRTHES calculation parameters.

-------------------------------------
Run computation
-------------------------------------

When calculation parameters are defined, close salome GUI and launch again **salome shell**.
Then, go to the directory **my_study** and edit **coupling_parameters.py** file.

Finally, to launch the computation just run **runcase** file with the command ``./runcase``.

