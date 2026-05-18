Schemas Theta
=============

Le schema theta definit le niveau d'extrapolation dans le temps :

* **theta = 0** : explicite
* **theta = 0.5** : extrapolation en n+1/2
* **theta = 1** : extrapolation en n+1

Configurations par propriete
-----------------------------


viscosity
^^^^^^^^^

Viscosity extrapolation scheme

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Methode
     - Theta
     - Description
   * - none
     - 0.0
     - Explicite
   * - linear
     - 0.5
     - Extrapole en n+1/2
   * - second_order
     - 1.0
     - Extrapole en n+1


heat_capacity
^^^^^^^^^^^^^

Heat capacity (Cp) extrapolation scheme

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Methode
     - Theta
     - Description
   * - none
     - 0.0
     - Explicite
   * - linear
     - 0.5
     - Extrapole en n+1/2
   * - second_order
     - 1.0
     - Extrapole en n+1


navier_stokes_source
^^^^^^^^^^^^^^^^^^^^

Navier-Stokes source terms

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Ordre terme source
     - Theta
     - Description
   * - explicit
     - 0.0
     - Termes sources explicites
   * - semi_implicit_1
     - 0.5
     - Semi-implicite theta=0.5
   * - semi_implicit_2
     - 1.0
     - Semi-implicite theta=1


turbulence_source
^^^^^^^^^^^^^^^^^

Turbulent source terms

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Ordre terme source
     - Theta
     - Description
   * - explicit
     - 0.0
     - Termes sources explicites
   * - semi_implicit_1
     - 0.5
     - Semi-implicite theta=0.5
   * - semi_implicit_2
     - 1.0
     - Semi-implicite theta=1


scalar_source
^^^^^^^^^^^^^

Scalar source terms

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Ordre terme source
     - Theta
     - Description
   * - explicit
     - 0.0
     - Termes sources explicites
   * - semi_implicit_1
     - 0.5
     - Semi-implicite theta=0.5
   * - semi_implicit_2
     - 1.0
     - Semi-implicite theta=1


scalar_diffusivity
^^^^^^^^^^^^^^^^^^

Scalar diffusivity extrapolation

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Methode
     - Theta
     - Description
   * - none
     - 0.0
     - Explicite
   * - linear
     - 0.5
     - Extrapole en n+1/2
   * - second_order
     - 1.0
     - Extrapole en n+1

