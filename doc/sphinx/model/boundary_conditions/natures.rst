Types de conditions aux limites
================================

Nature de BC
------------

* ``inlet`` : Entree de fluide (vitesse imposee)
* ``wall`` : Paroi (vitesse nulle ou glissement)
* ``outlet`` : Sortie de fluide (pression imposee)
* ``symmetry`` : Plan de symetrie
* ``free_inlet_outlet`` : Entree/sortie libre
* ``imposed_p_outlet`` : Sortie a pression imposee
* ``free_surface`` : Surface libre
* ``groundwater`` : Ecoulement en milieu poreux
* ``undefined`` : Type non defini


Choix de vitesse a l'entree
-----------------------------

* ``norm`` : Norme de vitesse imposee
* ``flow1`` : Debit massique impose
* ``flow2`` : Debit volumique impose
* ``norm_formula`` : Norme via formule MEG
* ``flow1_formula`` : Debit massique via formule MEG
* ``flow2_formula`` : Debit volumique via formule MEG


Direction de la vitesse
------------------------

* ``normal`` : Direction normale a la face
* ``coordinates`` : Direction par coordonnees cartesiennes
* ``formula`` : Direction via formule MEG
