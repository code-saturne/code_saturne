Tableau de reference
====================

Modeles et leurs parametres
----------------------------

.. list-table::
   :widths: 25 15 25 10 25
   :header-rows: 1

   * - Modele
     - Famille
     - Variables requises
     - WF défaut
     - WF autorisées
   * - LES_Smagorinsky
     - LES
     - -
     - -
     - -
   * - LES_WALE
     - LES
     - -
     - -
     - -
   * - LES_dynamique
     - LES
     - -
     - -
     - -
   * - Rij-EBRSM
     - RANS
     - rij, epsilon, alpha
     - 7
     - 0, 7
   * - Rij-SSG
     - RANS
     - rij, epsilon
     - 3
     - 0, 2, 3, 4
   * - Rij-epsilon
     - RANS
     - rij, epsilon
     - 3
     - 0, 2, 3, 4
   * - Spalart-Allmaras
     - RANS
     - nu_tilda
     - 2
     - 2
   * - k-epsilon
     - RANS
     - k, epsilon
     - 3
     - 0, 2, 3, 4
   * - k-epsilon-PL
     - RANS
     - k, epsilon
     - 3
     - 0, 2, 3, 4
   * - k-omega-SST
     - RANS
     - k, omega
     - 3
     - 0, 2, 3, 4, 7
   * - mixing_length
     - -
     - -
     - 2
     - 2
   * - off
     - -
     - -
     - 0
     - 0
   * - v2f-BL-v2/k
     - RANS
     - k, epsilon, phi, alpha
     - 0
     - 0
   * - v2f-phi
     - RANS
     - k, epsilon, phi
     - 3
     - 0, 3

Wall Functions — Mapping valeur → description
---------------------------------------------

.. list-table::
   :widths: 10 30 60
   :header-rows: 1

   * - Valeur
     - Nom
     - Description
   * - 0
     - Pas de loi de paroi
     - Résolution jusqu'à la paroi (y+ ~ 1). Requis pour v2f et certains LES.
   * - 2
     - Loi logarithmique standard
     - Loi log classique. Nécessite y+ ~ 30-300.
   * - 3
     - All y+ (two-layer)
     - Fonctionne pour tout maillage (y+ faible ou élevé). Recommandé par défaut.
   * - 4
     - Scalable wall function
     - Variante de la loi log, robuste pour y+ variable.
   * - 7
     - Two-scale model (EBRSM)
     - Spécifique au modèle Rij-EBRSM. Modèle à deux échelles de vitesse.

.. note::
   La valeur par défaut est réinitialisée automatiquement quand l'utilisateur
   change de modèle de turbulence dans la GUI (correction du bug wall function).
