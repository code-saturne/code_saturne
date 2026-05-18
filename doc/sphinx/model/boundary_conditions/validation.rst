Regles de validation
====================

Vitesse
-------

* ``velocity_norm`` : doit etre strictement positive

Diametre hydraulique
--------------------

* ``hydraulic_diameter`` : doit etre strictement positif

Intensite turbulente
---------------------

* ``turbulent_intensity`` : doit etre dans [0, 100] (pourcentage)

Nature de BC
------------

Doit etre une des valeurs definies dans ``Definitions/TypeDef[@name='BCNature']``.
