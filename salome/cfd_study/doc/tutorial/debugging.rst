==================================================================
*Code_Saturne* howto: **compile user source in debugging mode**
==================================================================

----------------
Introduction
----------------

We recommend to compile user sources with a debug build of
*code_Saturne* to verify and control array bounds, etc. Two versions of
*Code_Saturne* and *Neptune_CFD* are present in *SALOME_CFD*:

    - a production build in *CFD_MODULE*,
    - a debug build as a tool.

Only optimize version is available with GUI now.

------------------
Use debug version
------------------

To use debugging build, launch a salome shell with the command: ``~/salome/appli_x_y_z/salome shell``.

We define a variable in the environment for a direct access to *Code_Saturne* and *Neptune_CFD*. To use
both codes use:

    - ``$CFDSTUDY_DBG_ROOT_DIR/bin/code_saturne``,
    - ``$CFDSTUDY_DBG_ROOT_DIR/bin/neptune_cfd``.

