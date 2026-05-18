Thermal Fields
==============

Champs thermiques auxiliaires crees selon la configuration.
        KineticST si terme source cinetique actif.
        Yv si air humide.
        Temperature si energie interne resolue.

.. contents::
   :local:
   :depth: 1


Champs pour terme source cinetique thermique
--------------------------------------------

**Condition** : ``has_kinetic_st``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``kinetic_energy_thermal_st`` (Kinetic Energy Thermal ST)
     - 1
     - property
     - cells
     - postprocesse
   * - ``rho_k_prev`` (rho_k_prev)
     - 1
     - internal
     - cells
     - -
   * - ``inner_face_velocity`` (Inner Face Velocity)
     - 3
     - internal
     - interior_faces
     - -
   * - ``boundary_face_velocity`` (Boundary Face Velocity)
     - 3
     - internal
     - boundary_faces
     - -


Fraction massique vapeur eau pour air humide
--------------------------------------------

**Condition** : ``ieos_moist_air``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``yv`` (Water Vapor Mass Fraction)
     - 1
     - property
     - cells
     - postprocesse


Champs compressibles gradient pression et Cp
--------------------------------------------

**Condition** : ``ieos_not_none_and_temp_or_energy``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``algo:pressure_gradient`` (Pressure Gradient)
     - 3
     - internal
     - cells
     - -
   * - ``algo:pressure_increment_gradient`` (Pressure Increment Gradient)
     - 3
     - internal
     - cells
     - -
   * - ``isobaric_heat_capacity`` (Isobaric Heat Capacity)
     - 1
     - property
     - cells
     - -


Temperature comme propriete si energie interne resolue
------------------------------------------------------

**Condition** : ``thermal_internal_energy``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``temperature`` (Temperature)
     - 1
     - property
     - cells
     - postprocesse


Nombre CFL pour la variable thermique
-------------------------------------

**Condition** : ``cflt_active``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``cfl_t`` (CFL Thermal)
     - 1
     - postprocess
     - cells
     - postprocesse


Nombre CFL pour la pression
---------------------------

**Condition** : ``cflp_active``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``cfl_p`` (CFL Pressure)
     - 1
     - postprocess
     - cells
     - postprocesse

