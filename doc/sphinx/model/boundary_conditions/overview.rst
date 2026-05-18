Overview
========

``BoundaryConditionsRules.xml`` centralise les regles des conditions
aux limites de Code_Saturne. Il est lu par :

* Le **kernel C++** via ``cs_boundary_conditions_rules_manager``
* La **GUI Python** via ``BoundaryConditionsRulesManager`` (ElementTree)

Le principal apport est le mapping modele de turbulence -> variables BC
inlet, qui remplace 200 lignes de if/else if hardcodes dans
``cs_gui_boundary_conditions.cpp``.

Sections principales
--------------------

* ``Definitions`` : types valides (natures, choix vitesse, turbulence, scalaires)
* ``ValidationRules`` : mapping turbulence -> BC inlet, contraintes
* ``Defaults`` : valeurs par defaut
* ``Mappings`` : correspondances vers les flags C++
