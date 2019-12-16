===========================
Generalities about CFDSTUDY
===========================

----------------
Introduction
----------------

The **CFDSTUDY** is a component for the SALOME platform. The purpose of this
software is to provide an interface between CFD (Computational Fluid Dynamics)
softwares *Code_Saturne* and NEPTUNE_CFD with other modules of the platform.

*Code_Saturne* and NEPTUNE_CFD are CFD softwares from EDF R&D. *Code_Saturne*
could be freely downloaded from ``http://code-saturne.org``.

This document provides a tutorial for the use of CFDSTUDY with *Code_Saturne*.
For a *Code_Saturne* tutorial itself and more information about user issues,
please consult the software documentation.

Note: CFDSTUDY is a pure Python module of SALOME. Only the GUI part of the module
is implemented without its engine counterpart. Therefore, the dump functionality is
not available yet.

------------------------------
Reference functionalities
------------------------------

The main purpose of **CFDSTUDY** is to embed the GUI of *Code_Saturne* inside the
SALOME desktop and to make easier the setup of a case. For that, when the module is
loaded, several hooks are available:

- Menubar

    - **File > CFD code** menu:

        - Save CFD Data file
        - Save as CFD Data file

    - **CFDSTUDY** menu:

        - Set CFD study Location: allow to chose an existing study, or to create a new one.

            .. image:: images/CFDSTUDY_location.png
              :align: center
              :width: 10cm

        - Update Object Browser: refresh the directories list in the Object Browser.

        - Tools: display the file of parameter, open a terminal.

    - **Help** menu:

        - **CFDSTUDY module User's guide** menu: display this document in html format.
        - **CFD module** menu: direct access to Code_Saturne (and NEPTUNE_CFD):
                - theoretical guide,
                - tutorial guide,
                - user guide,
                - doxygen documentation.


- Toolbar (from the left to the right):

        .. image:: images/CFDSTUDY_toolbar.png
          :align: center


    - Set CFD study Location
    - Add a new case in a study
    - Launch the GUI
    - Run a case
    - Save CFD Data file
    - Save as CFD Data file
    - Close a GUI
    - Undo
    - Redo

- GUI: additional functionalities are available:

    - Groups of boundary faces can be selected in the Object Browser or graphically,
    - Groups of cells can be selected in the Object Browser or graphically,
    - Monitoring points can be displayed in the VTK viewver.

- Object Browser: several actions are available through a specific contextual menu (open by *Right click*)

    - Study directory:

    .. image:: images/CFDSTUDY_context_menu_study.png
      :align: center


    - Mesh file:

    .. image:: images/CFDSTUDY_context_menu_mesh.png
      :align: center


    - Case directory:

    .. image:: images/CFDSTUDY_context_menu_case.png
      :align: center

    - *SaturneGUI* file:

    .. image:: images/CFDSTUDY_context_menu_new_gui.png
      :align: center

    - File of parameters:

    .. image:: images/CFDSTUDY_context_menu_xml.png
      :align: center

    - File of functions in the *SRC* directory:

    .. image:: images/CFDSTUDY_context_menu_src.png
      :align: center

    - File of functions in the *REFERENCE* directory:

    .. image:: images/CFDSTUDY_context_menu_ref.png
      :align: center

    - File of functions in the *DRAFT* directory:

    .. image:: images/CFDSTUDY_context_menu_draft.png
      :align: center

    - Script of *runcase* file:

    .. image:: images/CFDSTUDY_context_menu_runcase.png
      :align: center

    - Results directories in the *RESU* directory:

    .. image:: images/CFDSTUDY_context_menu_resu.png
       :align: center

