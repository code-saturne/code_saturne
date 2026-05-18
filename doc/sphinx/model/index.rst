Code_Saturne -- Model Rules Documentation
=========================================

Documentation auto-generee depuis les fichiers XML de regles metier.

Architecture
------------

Les regles metier sont centralisees dans des fichiers XML partages entre
le kernel C++ et la GUI Python (pattern Repository + Layers).

.. code-block:: text

   data/model/
   |-- TurbulenceRules.xml           <- Modeles de turbulence
   |-- ThermalRules.xml              <- Modele thermique
   |-- TimeSteppingRules.xml         <- Schemas temporels
   |-- BoundaryConditionsRules.xml   <- Conditions aux limites
   `-- FieldsRules.xml               <- Creation des champs

Modules documentes
------------------

.. toctree::
   :maxdepth: 1

   turbulence/index
   thermal/index
   timestepping/index
   boundary_conditions/index
   fields/index

* :ref:`genindex`
* :ref:`search`
