Properties Fields
=================

Champs proprietes physiques.
        Density et viscosity toujours presents.
        CFL cache en regime permanent.
        TotalPressure seulement si non compressible.

.. contents::
   :local:
   :depth: 1


Masse volumique cellules et faces de bord
-----------------------------------------

**Condition** : ``always``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``density`` (Density)
     - 1
     - property
     - cells
     - ``CS_ENUMF_(rho)``, postprocesse
   * - ``boundary_density`` (Boundary Density)
     - 1
     - property
     - boundary_faces
     - ``CS_ENUMF_(rho_b)``


Viscosite moleculaire et turbulente
-----------------------------------

**Condition** : ``always``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``molecular_viscosity`` (Laminar Viscosity)
     - 1
     - property
     - cells
     - ``CS_ENUMF_(mu)``
   * - ``turbulent_viscosity`` (Turb Viscosity)
     - 1
     - property
     - cells
     - ``CS_ENUMF_(mu_t)``, cache si CS_TURB_NONE


Nombres adimensionnels CFL et Fourier
-------------------------------------

**Condition** : ``always``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``courant_number`` (CFL)
     - 1
     - property
     - cells
     - cache si permanent
   * - ``fourier_number`` (Fourier Number)
     - 1
     - property
     - cells
     - cache si permanent


Nombre de Courant volumique pour VoF
------------------------------------

**Condition** : ``vof_active``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``volume_courant_number`` (CourantNbVol)
     - 1
     - property
     - cells
     - cache si permanent


Pression totale statique + dynamique
------------------------------------

**Condition** : ``non_compressible``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``total_pressure`` (Total Pressure)
     - 1
     - property
     - cells
     - restart: auxiliary


Constante de Smagorinsky au carre pour LES dynamique
----------------------------------------------------

**Condition** : ``model_eq_CS_TURB_LES_SMAGO_DYN``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``smagorinsky_constant^2`` (Csdyn2)
     - 1
     - property
     - cells
     - -

