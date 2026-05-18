Champs crees par variable thermique
====================================

Temperature
-----------

* ``temperature`` (variable, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)

Enthalpie
---------

* ``enthalpy`` (variable, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)
* ``thermal_diffusivity`` (propriete, cellules)

Energie interne
---------------

* ``total_energy`` (variable, cellules)
* ``temperature`` (propriete, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)

Champs conditionnels
--------------------

**Si** ``has_kinetic_st`` :

* ``kinetic_energy_thermal_st`` (propriete)
* ``rho_k_prev`` (interne)
* ``inner_face_velocity`` (faces internes, dim=3)
* ``boundary_face_velocity`` (faces de bord, dim=3)

**Si** ``ieos == moist_air`` :

* ``yv`` (propriete, fraction massique vapeur d'eau)
