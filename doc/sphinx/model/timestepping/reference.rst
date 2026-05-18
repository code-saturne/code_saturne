Tableau de reference
====================

Methode d'extrapolation -> theta
---------------------------------

.. list-table::
   :widths: 40 20 40
   :header-rows: 1

   * - Methode
     - Theta
     - Proprietes
   * - none (iviext=0)
     - 0.0
     - viscosity, heat_capacity, scalar_diffusivity
   * - linear (iviext=1)
     - 0.5
     - viscosity, heat_capacity, scalar_diffusivity
   * - second_order (iviext=2)
     - 1.0
     - viscosity, heat_capacity, scalar_diffusivity

Ordre terme source -> theta
----------------------------

.. list-table::
   :widths: 40 20 40
   :header-rows: 1

   * - Ordre
     - Theta
     - Termes sources
   * - explicit (isno2t=0)
     - 0.0
     - NS, turbulence, scalaires
   * - semi_implicit_1 (isno2t=1)
     - 0.5
     - NS, turbulence, scalaires
   * - semi_implicit_2 (isno2t=2)
     - 1.0
     - NS, turbulence, scalaires
