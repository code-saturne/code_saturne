Tableau de reference
====================

Nature -> Flag C++
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Nature
     - Flag C++
   * - inlet
     - CS_BOUNDARY_INLET
   * - wall
     - CS_BOUNDARY_WALL
   * - outlet
     - CS_BOUNDARY_OUTLET
   * - symmetry
     - CS_BOUNDARY_SYMMETRY
   * - free_inlet_outlet
     - CS_BOUNDARY_FREE_INLET_OUTLET
   * - imposed_p_outlet
     - CS_BOUNDARY_IMPOSED_P_OUTLET
   * - free_surface
     - CS_BOUNDARY_FREE_SURFACE
   * - groundwater
     - CS_BOUNDARY_GROUNDWATER

Defaults
--------

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Parametre
     - Valeur par defaut
   * - bc_nature
     - inlet
   * - velocity_choice
     - norm
   * - velocity_direction
     - normal
   * - turbulence_choice
     - hydraulic_diameter
   * - scalar_bc_type
     - dirichlet
