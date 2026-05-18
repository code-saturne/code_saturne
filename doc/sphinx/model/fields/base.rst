Base Fields
===========

Champs de base crees dans toutes les simulations.
        Velocity et Pressure sont toujours presents.
        VoidFraction uniquement si algorithme VoF actif.

.. contents::
   :local:
   :depth: 1


Champ vitesse toujours present
------------------------------

**Condition** : ``always``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``velocity`` (Velocity)
     - 3
     - variable
     - cells
     - ``CS_ENUMF_(vel)``, couple


Champ pression toujours present
-------------------------------

**Condition** : ``always``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``pressure`` (Pressure)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(p)``, sans convection


Fraction volumique pour algorithme VoF
--------------------------------------

**Condition** : ``vof_active``

.. list-table:::
   :widths: 20 10 15 15 40
   :header-rows: 1

   * - Champ
     - Dim
     - Type
     - Localisation
     - Notes
   * - ``void_fraction`` (Void Fraction)
     - 1
     - variable
     - cells
     - ``CS_ENUMF_(void_f)``, sans diffusion

