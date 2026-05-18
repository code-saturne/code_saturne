Regles de validation
====================

Parametre theta
---------------

* Doit etre dans [0, 1]
* Ne pas modifier manuellement : initialise automatiquement depuis le XML
* Valeur sentinelle : -999 (declanche l'auto-init)

Pas de temps
------------

* ``dt_ref > 0``
* ``dt_min >= 0``
* ``dt_max > 0``
* ``dt_max > dt_min``
