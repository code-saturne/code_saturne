==================================================================
*Code_Saturne* howto: **compile user source in debugging mode**
==================================================================

----------------
Introduction
----------------

We recommand to compile user sources with de debugging version of
*code_Saturne* to verify and control array bounds, etc. Two versions of
*Code_Saturne* and *Neptune_CFD* are present in *SALOME_CFD* :

    - an optimize version in *CFD_MODULE*,
    - a debug version as a tool.

Only optimize version is available with GUI now.

------------------
Use debug version
------------------

To use debugging version, launch a salome shell with the command : ``~/salome/appli_x_y_z/salome shell``.

We define a variable in the environment for a direct access to *Code_Saturne* and *Neptune_CFD*. To use
both code use :

    - ``$CFDSTUDY_DBG_ROOT_DIR/bin/code_saturne``,
    - ``$CFDSTUDY_DBG_ROOT_DIR/bin/neptune_cfd``.

