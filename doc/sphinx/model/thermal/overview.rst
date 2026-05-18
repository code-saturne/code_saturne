Overview
========

``ThermalRules.xml`` centralise les regles du modele thermique de Code_Saturne.

Il est lu par :

* Le **kernel C++** via ``cs_thermal_rules_manager``
* La **GUI Python** via ``ThermalRulesManager`` (ElementTree)

Variables thermiques supportees
--------------------------------

* ``none``
* ``temperature``
* ``enthalpy``
* ``internal_energy``


Equations d'etat
----------------

* ``none``
* ``ideal_gas``
* ``moist_air``
* ``real_gas``


Plans thermodynamiques
----------------------

* ``PT``
* ``PH``
* ``PS``
