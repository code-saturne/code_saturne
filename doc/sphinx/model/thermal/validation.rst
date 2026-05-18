Regles de validation
====================

Parametre ``is_temperature``
-----------------------------

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Min
     - Max
     - Description
   * - 0
     - 1
     - 0 = scalaire passif/enthalpie, 1 = temperature

Parametre ``turbulent_schmidt``
--------------------------------

Doit etre strictement positif (> 0).

Lagrangien ordre 2
------------------

Le schema Lagrangien au second ordre est interdit si :

* Variable thermique = ``temperature`` ET echelle = ``kelvin``
* Variable thermique = ``enthalpy``
