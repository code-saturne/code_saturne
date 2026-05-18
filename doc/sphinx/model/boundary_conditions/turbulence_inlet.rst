Turbulence a l'entree
======================

Mapping modele -> variables BC
--------------------------------

Ce tableau remplace 200 lignes de if/else if hardcodes dans
``cs_gui_boundary_conditions.cpp``. Il est genere depuis
``BoundaryConditionsRules.xml``.

.. list-table::
   :widths: 25 15 25 35
   :header-rows: 1

   * - Modele
     - N valeurs
     - Fonction MEG
     - Variables
   * - Rij-EBRSM
     - 8
     - turbulence_rij_ebrsm
     - ``rij``, ``epsilon``, ``alpha``
   * - Rij-SSG
     - 7
     - turbulence_rije
     - ``rij``, ``epsilon``
   * - Rij-epsilon
     - 7
     - turbulence_rije
     - ``rij``, ``epsilon``
   * - Spalart-Allmaras
     - 1
     - turbulence_spalart
     - ``nu_tilda``
   * - k-epsilon
     - 2
     - turbulence_ke
     - ``k``, ``epsilon``
   * - k-epsilon-PL
     - 2
     - turbulence_ke
     - ``k``, ``epsilon``
   * - k-omega-SST
     - 2
     - turbulence_kw
     - ``k``, ``omega``
   * - v2f-BL-v2/k
     - 4
     - turbulence_v2f
     - ``k``, ``epsilon``, ``phi``, ``alpha``


Methodes de specification de turbulence
------------------------------------------

* ``hydraulic_diameter`` : Diametre hydraulique (loi de paroi lisse)
* ``turbulent_intensity`` : Intensite turbulente en pourcentage
* ``formula`` : Formule MEG utilisateur
