Setting Up a Computation
========================

.. toctree::
   :maxdepth: 2

   api/cs_ug_mesh_prepare
   api/cs_ug_mesh_select_c
   api/advanced_ale
   api/advanced_atmospheric
   api/advanced_cavitation
   api/advanced_coal_and_gas_combution
   api/advanced_compressible
   api/advanced_conjugate_heat_transfer
   api/advanced_coupling
   api/advanced_electric_arcs
   api/advanced_particle_tracking
   api/advanced_radiative_thermal
   api/advanced_solidification
   api/advanced_specific_physics
   api/advanced_turbomachinery
   api/cs_ug_cdo_gwf
   api/low_level_boundary_condition_definitions
   api/cs_ug_output

Calculation environment
-----------------------

- Notebook
- Time tables

Mesh
----

- Preprocessing
- Immersed Boundaries
- Volume zones
- Boundary zones

See :doc:`api/cs_ug_mesh_prepare` and :doc:`api/cs_ug_mesh_select_c`.

Calculation features
--------------------

Standard
~~~~~~~~

- Thermal model
- Body forces (gravity, Coriolis)
- Turbulence models
- Species transport
- Fans

Specific physics
~~~~~~~~~~~~~~~~

- Atmospheric flows — :doc:`api/advanced_atmospheric`
- Electric arc / Joule effect — :doc:`api/advanced_electric_arcs`
- Gas and coal combustion — :doc:`api/advanced_coal_and_gas_combution`
- Particle tracking (Lagrangian) — :doc:`api/advanced_particle_tracking`
- Radiative heat transfer — :doc:`api/advanced_radiative_thermal`
- Compressible flows — :doc:`api/advanced_compressible`
- Turbomachinery — :doc:`api/advanced_turbomachinery`
- ALE (deformable mesh) — :doc:`api/advanced_ale`
- Cavitation — :doc:`api/advanced_cavitation`
- Groundwater (CDO) — :doc:`api/cs_ug_cdo_gwf`
- Conjugate heat transfer — :doc:`api/advanced_conjugate_heat_transfer`
- Solidification — :doc:`api/advanced_solidification`

Volume conditions
-----------------

Source terms, porosity, head losses and initialization applied to volume zones.

Boundary conditions
-------------------

Nature and type of each boundary zone.
See :doc:`api/low_level_boundary_condition_definitions`.

Coupling parameters
-------------------

- Internal coupling
- Conjugate heat transfer coupling
- Fluid structure interaction

Time settings
-------------

- Start/Restart
- Time step option
- Velocity-Pressure algorithm
- Stopping criterion

Numerical parameters
--------------------

- Equation parameters
- Gradient reconstruction
- Solver and preconditioning

Postprocessing
--------------

- Calculator
- Time averages
- Volume solution control
- Surface solution control
- Profiles
- Balance by zone

See :doc:`api/cs_ug_output`.

Performance settings
--------------------

- Partitioning
- Input/Output
- MPI algorithms

