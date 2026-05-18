Variables thermiques
====================

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Variable
     - Formule diffusivite
     - Plan thermodynamique
     - Flag is_temperature
   * - ``none``
     - -
     - -
     - -
   * - ``temperature``
     - lambda
     - PT
     - 1
   * - ``enthalpy``
     - lambda/cp
     - PH
     - 0
   * - ``internal_energy``
     - lambda/cp
     - PH
     - 0


Echelles de temperature
-----------------------

* ``none`` : Pas de modele thermique actif
* ``kelvin`` : Temperature en Kelvin (K)
* ``celsius`` : Temperature en Celsius (C)
