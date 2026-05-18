Overview
========

``FieldsRules.xml`` centralise toutes les regles de creation de champs
de Code_Saturne. Il est lu par :

* Le **kernel C++** via ``cs_fields_rules_manager``
* La **GUI Python** via ``FieldsRulesManager`` (ElementTree)

Ce fichier remplace les blocs ``if/else if`` hardcodes dans
``cs_setup.cpp`` et ``cs_parameters.cpp``.

Modules
-------

* **Base** : Champs de base crees dans toutes les simulations.
        Velocity et Pressure sont toujours presents.
        VoidFraction uniquement si algorithme VoF actif.
* **Turbulence** : Champs turbulence crees selon le modele actif.
        Chaque famille de modele a ses propres variables.
* **Properties** : Champs proprietes physiques.
        Density et viscosity toujours presents.
        CFL cache en regime permanent.
        TotalPressure seulement si non compressible.
* **Thermal** : Champs thermiques auxiliaires crees selon la configuration.
        KineticST si terme source cinetique actif.
        Yv si air humide.
        Temperature si energie interne resolue.
