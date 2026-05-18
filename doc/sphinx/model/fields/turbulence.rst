Turbulence Fields
=================

Champs turbulence crees selon le modele actif.
        Chaque famille de modele a ses propres variables.

.. contents::
   :local:
   :depth: 1


Modeles k-epsilon et k-epsilon-PL
---------------------------------

**Condition** : ``itytur_eq_2``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``k`` (Turb Kinetic Energy)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(k)``
   * - ``epsilon`` (Turb Dissipation)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(eps)``


Modeles Rij (tenseur de Reynolds complet)
-----------------------------------------

**Condition** : ``second_order``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``rij`` (Rij)
     - 6
     - variable
     - cells
     - ``CS_ENUMF_(rij)``, couple
   * - ``epsilon`` (Turb Dissipation)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(eps)``, exclu si CS_TURB_LES_TAUSGS


EBRSM : champ alpha pour equation elliptique
--------------------------------------------

**Condition** : ``model_eq_CS_TURB_RIJ_EPSILON_EBRSM``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``alpha`` (Alphap)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(alp_bl)``, sans convection, sans terme instationnaire, idircl=0


Modeles v2f : k, epsilon + phi
------------------------------

**Condition** : ``itytur_eq_5``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``k`` (Turb Kinetic Energy)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(k)``
   * - ``epsilon`` (Turb Dissipation)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(eps)``
   * - ``phi`` (Phi)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(phi)``


v2f-phi : f_bar equation elliptique
-----------------------------------

**Condition** : ``model_eq_CS_TURB_V2F_PHI``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``f_bar`` (f_bar)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(f_bar)``, sans convection, sans terme instationnaire, idircl=0


v2f-BL-v2/k : alpha equation elliptique
---------------------------------------

**Condition** : ``model_eq_CS_TURB_V2F_BL_V2K``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``alpha`` (Alpha)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(alp_bl)``, sans convection, sans terme instationnaire, idircl=0


k-omega SST : k + taux de dissipation specifique omega
------------------------------------------------------

**Condition** : ``model_eq_CS_TURB_K_OMEGA``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``k`` (Turb Kinetic Energy)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(k)``
   * - ``omega`` (Omega)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(omg)``


Spalart-Allmaras : modele a une equation
----------------------------------------

**Condition** : ``model_eq_CS_TURB_SPALART_ALLMARAS``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``nu_tilda`` (NuTilda)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(nusa)``

