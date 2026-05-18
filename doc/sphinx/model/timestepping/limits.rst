Limites du pas de temps
=======================

.. list-table::
   :widths: 25 15 15 45
   :header-rows: 1

   * - Parametre
     - Defaut
     - Min
     - Description
   * - dt_ref
     - 0.1
     - > 0
     - Pas de temps de reference
   * - dt_min
     - 1e-5
     - >= 0
     - Pas de temps minimum
   * - dt_max
     - 1e4
     - > 0
     - Pas de temps maximum
   * - cfl_max
     - 1.0
     - > 0
     - Nombre de Courant maximum
   * - courant_max
     - 1.0
     - > 0
     - Courant maximum
   * - fourier_max
     - 10.0
     - > 0
     - Nombre de Fourier maximum

Contrainte : ``dt_min < dt_ref < dt_max``
