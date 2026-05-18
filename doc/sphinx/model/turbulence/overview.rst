Overview
========

``TurbulenceRules.xml`` est la source de verite unique pour les regles
de turbulence dans Code_Saturne. Il est lu par :

* Le **kernel C++** via ``cs_turbulence_rules_manager``
* La **GUI Python** via ``TurbulenceRulesManager`` (ElementTree)

Sections principales
--------------------

* ``Definitions`` : types valides (modeles, variables, proprietes)
* ``TurbulenceModels`` : configuration numerique de chaque modele
* ``ModelGroups`` : classification RANS / LES / Rij / TwoEquations
* ``ValidationRules`` : regles de validation GUI et kernel
* ``Defaults`` : valeurs par defaut
* ``Dependencies`` : compatibilite avec d'autres physiques
