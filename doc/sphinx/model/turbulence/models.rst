Turbulence Models
=================

.. contents::
   :local:
   :depth: 2


RANS Models
-----------


Rij-EBRSM
^^^^^^^^^

Variables requises : ``rij``, ``epsilon``, ``alpha``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


Rij-SSG
^^^^^^^

Variables requises : ``rij``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


Rij-epsilon
^^^^^^^^^^^

Variables requises : ``rij``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


Spalart-Allmaras
^^^^^^^^^^^^^^^^

Variables requises : ``nu_tilda``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS


k-epsilon
^^^^^^^^^

Variables requises : ``k``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations


k-epsilon-PL
^^^^^^^^^^^^

Variables requises : ``k``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations


k-omega-SST
^^^^^^^^^^^

Variables requises : ``k``, ``omega``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations


v2f-BL-v2/k
^^^^^^^^^^^

Variables requises : ``k``, ``epsilon``, ``phi``, ``alpha``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS


v2f-phi
^^^^^^^

Variables requises : ``k``, ``epsilon``, ``phi``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS


LES Models
----------


LES_Smagorinsky
^^^^^^^^^^^^^^^

Proprietes requises : ``turbulent_viscosity``

Groupes : LES


LES_WALE
^^^^^^^^

Proprietes requises : ``turbulent_viscosity``

Groupes : LES


LES_dynamique
^^^^^^^^^^^^^

Proprietes requises : ``smagorinsky_constant^2``, ``turbulent_viscosity``

Groupes : LES


RijModels Models
----------------


Rij-EBRSM
^^^^^^^^^

Variables requises : ``rij``, ``epsilon``, ``alpha``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


Rij-SSG
^^^^^^^

Variables requises : ``rij``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


Rij-epsilon
^^^^^^^^^^^

Variables requises : ``rij``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, RijModels


TwoEquations Models
-------------------


k-epsilon
^^^^^^^^^

Variables requises : ``k``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations


k-epsilon-PL
^^^^^^^^^^^^

Variables requises : ``k``, ``epsilon``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations


k-omega-SST
^^^^^^^^^^^

Variables requises : ``k``, ``omega``

Proprietes requises : ``turbulent_viscosity``

Groupes : RANS, TwoEquations

