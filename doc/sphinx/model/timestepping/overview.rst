Overview
========

``TimeSteppingRules.xml`` centralise les regles de discretisation temporelle.

Il est lu par :

* Le **kernel C++** via ``cs_time_stepping_rules_manager``
* La **GUI Python** via ``TimeSteppingRulesManager`` (ElementTree)

Schemas temporels supportes
----------------------------

* ``steady`` : Regime permanent
* ``first_order`` : Schema Euler implicite 1er ordre
* ``second_order`` : Schema Crank-Nicolson 2eme ordre


Methodes d'extrapolation
------------------------

* ``none`` : Pas d'extrapolation (theta=0)
* ``linear`` : Extrapolation lineaire (theta=0.5)
* ``second_order`` : Extrapolation second ordre (theta=1)
