Master (not on release branches yet)
------------------------------------

### Numerics:

- Reshape the interface with the MUMPS library (sparse direct solver)
  to make the usage of common MUMPS options easier. Take benefit from
  optimizations available in the latest 5.6.2 version.

### User changes:

- Compressible flows: remove uscfx1 and uscfx2 user-defined functions.
  Standard functions such as cs_user_parameters can be used instead.

Release 8.1.0 (2023-12-13)
--------------------------

### Physical modeling

- Remove obsolete and experimental heavy fuel combustion model
  (`cs_glob_physical_model[CS_COMBUSTION_FUEL]` / `ippmod(icfuel)`).
  That model was not accessible using the GUI, and not covered in the
  validation suite, with no known successful use in at least 10 years.

- Remove old pulverized coal combustion model with Lagrangian reciprocal
  approach (`cs_glob_physical_model[CS_COMBUSTION_PCLC]` / `ippmod(icpl3c)`).
  That model was not accessible using the GUI, and not covered in the
  validation suite, so probably unused since 10 years at least.

- Cooling towers model: Adding rain drift velocity equation solving.
  New keyword solve_rain_velocity has been added to allow the solving of
  rain drops drift velocity equation. The model has been adapted from
  existing one related to coal combustion.

### User changes:

- GUI: add access to various "per variable" equation parameters.
  * Handling of convection schemes is revamped, allowing selection
    of "true SOLU" (with upwind gradient), blended, and TVD/NVD schemes.
  * Most individual gradient options are now available in a new tab.

- Add high-level objects to handle open (inlet/outlet) boundary conditions
  in the general boundary conditions mechanism.
  * This is an evolution the the handling of BC's with the GUI, now
    usable also with user-defined functions.
  * This will allow sharing contexts wich may be used by CDO-type
    functions such as cs_dof_function_t in a shared manner between
    variables.
  * The coal combustion model now uses this mechanism for inlet conditions.

- run.cfg: the user can now always use global sections to define entries,
  with the syntax [key] for a key among "run_prologue", "run_epilogue",
  "compute_prologue" and "compute_epilogue".
  By doing so the defined entries are either used if nothing is defined
  for the computation resource or are prepended to the resource based options.
  This allows defining case specific and resource independent entries.

- GUI: when adding a postprocessing mesh, it is associated to the default
  writer automatically (and this can be modified as usual).

- Postprocessing: add the possibility of viewing interior face data at
  face centers, so as to have the data at correct location but keep it light
  on rendering (using a "Point Gaussian" view under ParaView for example).

- When running mesh quality checks, interior face values are now output
  at those locations, if the user has defined a postprocessing mesh based
  on interior faces or interior face centers.
  * The "per cell" and "per vertex" maxima are removed.

- Modifications in the way to set MUMPS solvers using cs_sles_param_t
  structure. Advanced settings can now be set thanks to the two functions
  cs_param_sles_mumps() and cs_param_sles_mumps_advanced(). More advanced
  settings can still be specified using the user-defined function
  cs_user_sles_mumps_hook

### Numerics:

- Use local fixed-point algorithm for vector and tensor least-squares
  gradient reconstruction at boundary.
  * Improves the performance of this operator, without changing
    the formulation.
  * To revert to the previous algorithm, set the following
     environment variable:
    `CS_GRADIENT_LSQ_BOUNDARY=legacy`.

- Do not activate gradient clipping by default even for least-squares
  gradients; allow its use (as an option) for all gradient types.

### Architectural changes:

- Remove support for freesteam (IAPWS-IF97 steam tables) library,
  as this library has not been maintained since 2013, and its build
  system is based on scons with Python 2 files, so is not usable
  anymore. Use EOS (if available) or CoolProp instead.

### Bug fixes:

- Preprocessor: workaround for CGNS incompatibility with files generated with
  CGNS 3.4.1 and with 3.4 version info (confused with 3.4.0, whose polygon and
  polyhedra representations are inccompatible). That version should be avoided
  but it seems that Star-CCM+ 2022 uses it.

- Preprocessor: fix reading of polygons in MED files, broken in v7.3.

- Fix numerical settings for Navier-Stokes problems with CDO schemes. User
  modifications may not be taken into account and one remains to the default
  settings (e.g. Hodge parameters for the viscous term)

Release 8.0.0 (2023-06-30)
--------------------------

### User changes:

- Add user postprocessing functions to compute mean and/or integrals of
  scalars over selections of boundary faces or over boundary zones.

- Add a user postprocessing function to compute total pressure over a
  selection of boundary faces.

- Postprocessing over a slice (plane, disc, annulus) can now be done
  * Slices are defined using a user function.
  * User can compute integral or mean values using scalar and/or vectorial
    weights of a scalar.
  * This uses medcoupling as an underlying building brick, and thus only
    available if the prerequisite is available.

- Add the functionnality of "time tables".
  * User can now define a time table, based on a data file, either using the
    GUI or an advanced user function.
  * The time table can then be used directly in mathematical formulae in the
    GUI (if the table is defined in the GUI), or in the different user
    functions.

- Add basic immersed boundary method features based on the
  neptune_cfd approach.

- CoolProp support now actually works with tabular properties.

- Cartesian mesher can now handle multiple blocks
  * This requires modification of basic API, specifically by adding
    an additional input argument for some functions (name of id).
  * A cartesian block can now have a name, which will lead to a renaming
    of boundary groups such as instead of "X0" the user will have
    "<name>_X0".

- The API for the definition by array has been modified to handle more
  complex cases. This induces a modification of some user settings,
  for instance when defininig a property, a boundary condition or a
  source term by array. This mainly occurs when using CDO schemes.

- The user can now specify a custom thermodynamic pressure.

### Architectural changes:

- Update and improve build configuration automation, especially
  when using the `--with-salome` option with SALOME 9.9 or 9.10 CAS builds.
  * For ParaView Catalyst, the `--enable-dlopen-rtld-global` option is no
    longer necessary, as the matching option is handled automatically
    when loading the matching plugin.

- Add and deploy several functions in cs_array.[c,h] to centralize
  simple operations on arrays (fill with zero, set a constant value on
  the whole array or partly, copy partly or fully an array into
  another one and so on). These operations take into acount OpenMP if
  available.

Release 7.3.0 (2022-12-20)
--------------------------

### User changes:

- Now allow handling vertex-based fields when restarting from a different
  mesh. This can be sueful for ALE.

- Add `cs_function` API to refactor pre and user-defined location data
  evaluation. This defines function objects which can be evaluated
  on the fly, on the required mesh subset, and be otherwise used similarly
  to fields for postprocessing and logging.
  * This is a first step. Logging output is still missing, and most
    candidate functions are not migrated yet.

- Documentation: user guide fully migrated to HTML (based on Doxygen
  and Markdown).

- Preprocessor (mesh import): allow removal of cells with negative volume
  or incorrect topology.

- Groundwater flows: for the legacy module, now use global gravity
  definition rather than module-specific (`darcy_gravity`) variables;
  This may require updating case setups: using the GUI is recommended,
  but otherwise, the gravity may be defined in `cs_user_physical_model`.

- Add new optimized extended neighborhood reduction method;
  An additional option to keep a complete extended neighborhood on the
  boundary was improved and made available in the GUI.

- Change cs_user_boundary_conditions function arguments for cleaner
  icodcl/rcodcl access: pointers to the array sections associated to
  each field are now available in that field's bc_coeffs sub-structure,
  simplifying the index arithmetic.

- Simplify the use of soil-atmosphere model of the atmospheric module.

- Add experimental immersed boundary features. This can be computed from
  cloud of points.

- Groundwater flows: for the CDO module, modify the way the dispersive
  tensor is defined. Remove the theta scaling in front of the
  molecular diffusivity.

- Groundwater flows: Modify the default settings for a tracer equation
  (unsteady/reaction terms to voronoi Hodge, advection scheme - hybrid
  upwind/centered, linear solvers: jacobi preconditioner and tolerance
  to 10^-8). This follows from a first feedack on the VnV cases of the
  default options.

- Groundwater flows: Add a radioactive tracer. The default tracer has
  no first order decay coefficient from now. The first order decay
  coefficient is not anymore associated to a couple (tracer, soil) but
  only to a tracer.

- Groundwater flows: Easier handling of a radioactive decay chain. Add
  a function following the Bateman's solution to define the initial
  state of each radioactive tracer.

### Physical modeling:

- Add hybrid RANS/LES turbulence model (HTLES) from V. Duffal PhD.

- Add some atmospheric universal functions for large scale idealized wind
  profiles.

- Groundwater flows: Modify the precipitation model (add radioactive
  decay, work with mass of precipitate rather than concentration)

### Bug fixes:

- Restart from different mesh: fix values shift in variable exchange when
  some points are not located (updated PLE library).

- Fix performance issues with K-cycle and pairwise aggregation and MSR
  matrices in case of penalized diagonal entries or isolated diagonal
  values.

- Fix the stopping criterion for Uzawa-CG algorithm (used by CDO
  schemes for Stokes problems).

- Fix the precipitation model for radionuclide tracers in the GWF
  module (CDO). Add the missing radioactive decay phenomenon.

- Fix the treatment of Robin BCs in CDO schemes.

### Architectural changes:

- Completed migration of turbulence models source code from Fortran to C.

- Remove libtool dependency.

- FMI: modified FMU client and code_saturne controller communication
  protocol to regroup input/output variables exchange and reduce latency.

- Drop support for older EOS versions (prior to 1.8.0).

- Add the possibility to define the values at the boundary faces for a
  cs_property_t structure

### Studymanager:

- New option (--slurm-batch-size=N with N>0) to submit batches of cases using
  the SLURM resource manager on cluster. Cases are sorted by number of
  processors and level of dependency.

- New option (--report) to generate description report of studies based on latex
  file in the STUDY/REPORT folder.

### Numerics:

- Add MUMPS as a preconditioner of a Krylov method.

- Add SPD matrices as a new type of input matrices for MUMPS. This
  distinction corresponds to the key SYM in MUMPS. Now, "mumps" refers
  to SYM=0 (not symmetric), "mumps_sym" refers to SYM=2 (general
  symmetric) and "mumps_ldlt" refers to SYM=1 (SPD matrices). The same
  behavior applies to the single-precision version.

- Several optimizations/fixes for MUMPS (zero-value entries are not
  sent to MUMPS for instance).

- MUMPS user hook API has been changed

- Add a scaled pressure mass matrix as a Schur approximation for Uzawa
  CG algorithm (used by CDO schemes)

### CDO:

- Setup: Fix inconsistencies between cs_navsto_param_t and
  cs_equation_param_t for the momentum equation when settings a
  Navier-Stokes problem. This may induce a modification of user source
  files.


Release 7.2.0 (2022-06-30)
--------------------------

User changes:

- Fixes errors with Python >= 3.10 due to non-robust version tests.

- Add external matrix type for PETSc, handled in cs_matrix_petsc, so no
  additional copy of coefficients is needed in code_saturne.
  * Native matrices in the legacy code automatically switch to this mode,
    while matrices can still be provided in MSR or CSR form, as assembly is
    expected to be faster (though no copy is saved in this case).
  * Some PETSc structures are cas as (void *) pointers to avoid propagating
    compile warnings due to PETSc headers, which can now be included
    more sparingly.

- Subdivision of warped faces now only modifies boundary faces.

- Major changes to the studymanager GUI tool, which is now more complete
  and especially more robust.

- Replace `ussmag` Fortran user-defined subroutine with
  `cs_user_physical_properties_smagorinsky_c` C function.

- Remove legacy `cs_user_mass_source_terms` function.

Bug fixes:

- Balance by zone: fix crash in some configurations.

- Probes and profiles: fix possible location error on hybrid meshes
  depending on order of appearance of cell types.

- Preprocessor: automatically detect CGNS file type (ADF or HDF5).

- Preprocessor: fix crash when reading CGNS files with NGON (polyhedra)
  sections and vertex-based boundary condition definitions.

- Fix gradient reconstruction for anisotropic cases using
  Green-Gauss with least-squares gradient face values algorithm.

- Solidification/CDO: Several fixes to handle the melting process.

- CDO: Fix to handle properly the Neumann and full Neumann boundary conditions.

- Boundary layer insertion: Take properly into account vertices with a
  prescribed displacement.

Architectural changes:

- Debug wrapper and other tool wrappers are now configured in the
  `run.cfg` file, and can be directly defined through the
  `code_saturne run` command-line arguments.

- Drop support of AMD ACML (which is end-of-life) and IBM ESSL BLAS libraries,
  which were only used in unit tests comparisons.

Numerics:

- CDO: *Major reshape* Add a generic framework to manipulate linear systems
  arising from CDO schemes. This framework handles classical linear systems
  (scalar- or vector-valued) but also (2x2) saddle-point systems and linear
  systems of coupled equations (NxN). The notion of block is at the core of the
  framework and it can handle matrix in the standard way (cs_matrix_t) but also
  a collection of cs_matrix_t structures or a matrix which is unassembled in a
  part of the system or for the full system. The aim is to get a framework to
  manipulate block matrices and apply to it either external libraries or
  in-house algorithms to solve these systems.

- CDO: Add a resolution by increment in scalar-valued CDO-Vb schemes useful to
  handle non-linear algorithm (steady and unsteady cases are available)

- CDO: Add new ways to define source terms in CDO-Vb schemes (array defined on
  different locations)

- velocity-pressure algorithm: Add an option to solve the pressure correction
  with a CDO face-based scheme while still using the legacy FV scheme for the
  velocity prediction

Physical modeling:

- For VoF, add a surface tension force, accessible from the GUI. If null, one
  recovers previous behaviour.

- GWF/CDO: Add two new models for two-phase flows in porous media: an immiscible
  and a miscible model. Solved variables are the liquid and gas pressures. Only
  user-defined laws are available in this case to define the behavior of the
  soil. Several numerical algorithm are currently under tests (coupled or
  segregated solver; different ways of treating the non-linearities). This
  functionality is still under development and may evolve in depth.

Release 7.1.0 (2021-12-21)
--------------------------

User changes:

- UI: Add optional `--notebook-args [key1=val1] [key2=val]> ...`,
  `--parametric-args [arg1] [arg2] ...`, and `--kw-args [arg1] [arg2] ...`
  arguments to `code_saturne run` and `code_saturne submit` commands.
  * If provided, they are available as the `notebook` and `kw_args`
    attributes of each case domain object (as a Python dictionnary
    and list respectively). Otherwise, these attributes are `None`.
  * These options are applied when running the case, and allow modifying
    the setup "on the fly" at run initialization.
  * In the studymanager module, these options can be passed using
    the `<notebook>`, `<parametric>`, and `<kw_args>` tags respectively,
    providing a safe replacement to the `<prepro>` tag, which is now
    deprecated due to its possible side effects in the case directory.

- Turbulence: RSM models now default to isotropic (Shir) diffusivity model.

- Gas mix now accepts user definitions, using the
  `cs_gas_mix_add_species_with_properties` function.
  * The model currenty still assumes that of n associated species,
    n-1 will be solved variables, and the last one will be
    of property field, deduced from the others so the sum of
    mass fractions is 1.

Physical modeling:

- Solidification: Add a Stefan model for taking into account the phase
  change in a purely thermal model.
- Solidification: Add a non-linear update of the thermal source term
  with the Voller model relying on the variation of enthalpy.
- Solidification: Add a new strategy to update the thermal source term
  in the Voller (linear or non-linear) model relying on the
  solidification path.
- GWF: Major change in the way to use this module.

Numerics:

- Major change to the uncoupled RSM solver (`irijco = 0`):
  * The Rij components are handled as a single tensor, similarly to
    the coupled mode.
  * Extra-diagonal terms in local Rij component coupling matrices
    (such as boundary condition components) are set to zero
    in this uncoupled mode.
  * The numerical behavior should be basically unchanged, but switching
    from the coupled to the uncoupled solver does not require specific
    tests in user code any more.

- Add support for HYPRE solvers.
  * Most general purpose solvers and preconditioners supported by HYPRE
    for the IJ interface are supported. See
    [HYPRE documentation](https://hypre.readthedocs.io/en/latest/ch-solvers.html)
    for details.
  * GPU use may be enabled using the `cs_sles_hypre_set_host_device` function.
  * HYPRE uses its own matrices, handled in cs_matrix_hypre, so no additional
    copy of coefficients is needed in code_saturne.

- Add GCR (Generalized Conjugate Residual) which enables a flexible
  preconditioning for non-symmetric linear system. Equivalent to a flexible
  GMRES.

- Several improvements for the MUMPS interface
  * Add a MUMPS interface for native matrices (with or without symmetric
    storage).
  * Add a MUMPS interface to handle single-precision arithmetic. This
  functionnality allows one to get a better efficiency without degrading
  the accuracy when considering MUMPS as a preconditioner while reducing
  the memory consumption.

- CDO: Add block preconditioning for symmetric saddle-point problem and a
  MINRES algorithm (relying on a tuned storage of the saddle-point system).
  GCR algorithm for more complex block preconditionning such upper/lower Schur
  block preconditioning or Symmetric Gauss-Seidel block preconditioning.
  Several approximations of the Schur complement are available.

- CDO: Add an Uzawa algorithm accelerated with a CG (conjuguate gradient)
  strategy and preconditioned with a Cahouet-Chabard technique. This is
  another solver to solve saddle-point problem arising from the monolithic
  coupling of CDO face-based schemes for the (Navier-)Stokes equations.

- CDO: Add a transformation of saddle-point problem using the Notay's
  change of variable. Available only with a link to the PETSc library.
  This is a based on a common work with the ALGO team at CERFACS.

- CDO: The default settings for solving an equation is now a GCR
  algorithm using a 1st order Neumann polynomial preconditioner. Moreover,
  the tolerance is set to 1e-6 instead of 1e-8.

- CDO: Optimized pre-defined settings for Hypre/BoomerAMG and
  PETSc/GAMG.  This enables a significant improvement in terms of
  elapsed time when solving (at least) symmetric positive systems.

- CDO: Extend the notion of block preconditioning when solving
  vector-valued equations with PETSc. Now multiplicative and symmetric
  Gauss-Seidel block preconditioners can be easily defined. Previously,
  only additive block preconditioning was easily available.

- CDO: Add the Anderson acceleration algorithm to improve the
  convergence of the Picard algorithm when a non-linear problem has to
  be solved. Joint work with Qingqing FENG.

Architectural changes:

- Add optional MPI one-sided communications mode for ghost cells
  synchronization, using active target RMA "get" semantics.

- Simplify handling of timers, focusing on wall-time.

Release 7.0.0 (2021-06-15)
--------------------------

User changes:

- GUI: add the possibility to specify meteo profiles input directly in
  code_saturne rather than with a meteo file.

- Run script: `define_domain_parameters` from `cs_user_scripts.py` is
  now called after copying initial data to the case's execution directory,
  so as to allow case and run id-specific modifications to be done, avoiding
  modification of the base data and user sources.

- GUI: for gas combustion, allow creation of a thermochemistry data file
  based on user inputs.

- GUI: Remove legacy definitions of syrthes coupling when opening old XML file:
  * Starting from v7.0 a syrthes coupling is defined as a boundary
    condition and no longer as an independent model.
  * Pre v7.0 definitions are saved to a file
    'deprecated_syrthes_coupling_data.txt' and are then removed from the
    XML file.
  * A warning pop-up informs the user of the change.

- GUI: Add more options for probes/profiles:
  * User can now choose whether or not to snap a probe/profile point
    which is currently the default behavior. When snap is activated for
    a set of points and several are contained within the same cell, we
    only keep the one closest to the cell center.
  * User can now activate the interpolation for probes/profiles, hence
    providing an interpolated value at the chosen coordinate instead of
    the value of the closest cell center.
  * For better robustness for the interpolation, users can now only
    define fields as output for probes/profiles, and no longer only a
    component. For example, for the velocity, all 3 components will be
    written. For older xml files which contained components, the latter
    are kept for the sake of compatibility.

- GUI: User formulae for physical properties can now be defined
  over multiple zones.
  * This changes the handling of physical properties.
  * Physical properties are now defined as a volume condition, with
    the zone "all_cells" as the default one.
  * User-laws for physical properties can now be defined over multiple
    zones using the MEG mechanism. This requires the definition of a
    variable property over the main zone ("all_cells").
  * A zone can now be defined as solid, which then allows activating the
    internal coupling functions for a given set of scalars (see the
    "Coupling parameters" page).
  * The internal coupling function is now available in the GUI.

- Add a `cs_base_get_run_identity` utility function, allowing the
  user to query the run_id, case and study name from a running
  computation's user-defined functions.

- Major evolution concerning coupling with external codes
  * Create a generic "Coupling parameters" page after boundary
    conditions in the tree view. This page gathers the parameters for
    Syrthes, code_aster and cathare couplings.
  * The definition of a coupling with Syrthes is now done directly as a
    thermal boundary condition (like flux or temperature).
  * Display external coupling tab as deprecated.

- GUI: Change available options for time-stepping.
  * Previous steady option ("Steady (constant relaxation coefficient)"),
    `idtvar` = -1, is now only available if it is activated in the
    opened xml file. It is then named "Deprecated: Steady ..."
    to indicate its new status.
  * Previous "pseudo steady" option, `idtvar` = 2, is renamed as
    "Steady (local time step)" and is now the de-facto steady algorithm.

- Rename and reorganize `cs_stokes_model_t` and `cs_piso_t` structures.
  * Structures are renamed to `cs_velocity_coupling_model_t` and
    `cs_velocity_coupling_param_t`, and members are distributed in a more
    logical and consistent manner.

- When solving for enthalpy, multiplying a face's boundary condition
  code (icodcl) by -1 allows defining it by temperature instead.
  This option is also available for prescribed values using the GUI.

- Add a `code_saturne symbol2line` command.
  This command can be used with debug symbols from a stack trace in order
  to find the origin of the call.
  Usage example for the line:
  1: 0x55834c6102d7 <cs_function+0x17>             (cs_solver)
  Run:
  code_saturne symbol2line -s cs_function+0x17
  Additional info is provided using 'code_saturne symbol2line --help'

Bug fixes:

- Fix halo construction deadlock for some rare partitioning cases.

- Compressible: fix imposed inlet/outlet boundary condition.
  Boundary mass flux was not computed using the Rusanov scheme,
  and was not consistent with momentum and energy convective fluxes on
  the boundary faces with imposed inlet/outlet boundary condition.

- Fix overwrite of user-defined function settings for probe writer
  by GUI or default settings.

- GUI: fix incomplete handling of compressible boundary conditions
  preventing run unless completed in user subroutines.

- Fix off-by-one return values in cs_selector_get_family_list
  (advanced user selection).

- Gas combustion: do not ignore handle GUI-defined settings of model scalar
  numerical parameters.

- Avoid division by zero in parallel when trying to output empty profiles.

Numerics and physical modeling:

- Allow LES with Lagrangian module.

- Allow combination of radiative transfer and internal coupling.

- Change default precision for linear solvers
  from 1e-8 to 1e-5 except for the pressure field.

- Activate improved hydrostatic pressure interpolation (iphydr) by default.

- Add an explicit treatment of the advection term in CDO face-based equations.

- Atmo: Major modification the scheme for Solar radiation done the PhD of
  L. Asmar.
  In a previous version of the 1D solar radiation scheme of code_saturne, the
  effect of aerosols has been introduced by making the approximation that
  aerosols are purely diffusive. The same approximation was taken for
  water vapor absorption in the UV-visible band (O3 band) but with
  different optical properties.
  In this new version, the aerosols are taken into account as absorbing
  and diffusive particles characterized by their optical properties:
  aerosol optical depth, single scattering albedo and asymmetry factor. To
  do that, the adding method for multiple scattering has been added in the
  UV-Visible band.
  Other improvements have been made for:  the optical air mass, the
  Rayleigh diffusion and the direct radiation estimation. In addition,
  the absorption by minor gases has been introduced.

Architectural changes:

- Move salome_cfd extensions (CFDSTUDY module) to a separate repository.
  It should allow both simplifying installation in some cases
  and improving  maintainability.

Release 6.3.0 (2020-12-21)
--------------------------

User changes:

- GUI: volume and boundary zones are now defined under the "Mesh"
  section, so that their selection can be used in a consistent
  manner in all following sections.
  * some specific boundary conditions, such as those for the Lagrangian
    model, are now grouped with standard zone boundary conditions rather
    than appeating in a separate "Additional BC models" section.

- Add `--dest` option to `code_saturne run` and `code_saturne submit`
  commands. This allows specifying a top directory separate from the
  one in which a case is present.

- Change default options for the VOF model (convective scheme is now
  CICSAM by default, with a continuous limiter ensuring
  bounds).

- Make dynamic relaxation option for pure diffusion equation default option.
  (iswdyn = 2). This make the iterative process solving of the legacy
  solver more robust (especially for pressure).

- Add a `code_saturne parametric` command which currently allows the
  modification of several parameters of a case and the creation of a
  case  based on a reference case.
  Integration with studymanager and cfdstudy (for parametric studies)
  is to be added soon (before v7.0 release).

- Preprocessing: added mesh modification flag types to more easily
  activate re-partitioning.

- LES inflow synthetic turbulence configuration is now zone-based, with
  a new, more consistent user API.

- Add `cs_runaway_check` functions to try to detect diverging computations
  early enough for clean stop.
  * Default check is on velocity > 1e4 for incompressible computations,
    1e5 for compressible computations.

- Add possibility to generate a cartesian mesh on the fly
  * The GUI now allows to define a cartesian mesh which will be
    generated during runtime.
  * Needed parameters are minimal and maximal values for X,Y and Z
    coordinates, and number of cells per direction.
  * User can define the progression law for cells' size :
    - Constant step
    - Geometric step
    - Parabolic step (Geometric step with a symmetry w.r.t center-coordinate
      of the direction)
  * The newly created mesh contains 6 groups for boundary faces: X0, X1, Y0,
    Y1, Z0 and Z1 for the respectively the min and max coordinate of each
    direction (X,Y and Z)

- When creating a case with `code_saturne create`, the `--noref` option is
  enabled by default. It can be cancelled using the `--copy-ref` option.

- GUI: user source file editor can now also use reference files from the
  installation directory.

- code_saturne mesh format (.csm) is now treated in the same way as external
  formats.
  * Multiple meshes can be added in the "Mesh" page.
  * Preprocessor (read operation) is not done for these mesh, though mesh
    modification is.
  * For back compatibility matters, '.csm' files can still be used as a
    "mesh_input"

- When no boundary condition is provided for some boundary faces, a default
  (wall or symmetry) is used. This can be set using `cs_boundary_set_default`.

- Modify the way the code is handling input file.
  * If no xml file is provided by the 'run.cfg' file and a 'setup.xml' exists,
    the latter is used (previous behaviour). If it does not exist, the
    code prints a warning message indicating that no xml was found.
  * If an input file, other than 'setup.xml', is provided using the 'run.cfg'
    or '-p' option (for code_saturne run only), and a 'setup.xml' exists inside
    the DATA folder, the provided input file is used. Before this change
    'setup.xml' was the one used. A warning is printed to warn the user that
    having the two is against code_saturne BPG.

- Modify available C-API for usage of ParaMEDMEM coupling (MEDCoupling MPI).
  Specific user functions are added (cs_user_paramedmem_coupling.c) to allow
  an easier definitions of coupling. Examples are provided in the
  'cs_user_paramedmem_coupling-base.c' user_example file.
  Send/recieve operations still need to be done by the user, but are simplified
  by allowing send/recv based on a cs_field_t pointer.

- GUI: present volume and boundary conditions using sub-nodes in the
  left-hand tree for the appropriate zones. This is a first step in
  a change of zone setup presentation.

Bug fixes:

- GUI: Fix backward compatibility update function for NCFD v6.1 and later.

Numerics and physical modelling:

- CDO schemes: Add different strategies for the treatment of the advection term
  in Navier-Stokes equations (linearized implicit and explicit treatment). The
  Picard algorihtm (implicit) is the default as in the previous version. This is
  based on results obtained during R. MILANI's PhD.

- Amtospheric module: Major modification the scheme for Solar radiation based
  on the PhD of L. Asmar.

  In a previous version of the 1D solar radiation scheme of code_saturne, the
  effect of aerosols has been introduced by making the approximation that
  aerosols are purely diffusive. The same approximation was taken for
  water vapor absorption in the UV-visible band (O3 band) but with
  different optical properties.

  In this new version, the aerosols are taken into account as absorbing
  and diffusive particles characterized by their optical properties:
  aerosol optical depth, single scattering albedo and asymmetry factor. To
  do that, the adding method for multiple scattering has been added in the
  UV-Visible band.

  Other improvements have been made for:  the optical air mass, the
  Rayleigh diffusion and the direct radiation estimation. In addition,
  the absorption by minor gases has been introduced.

- Atmospheric module: add 3D radiative transfer for solar wave (shortwaves):
  direct and diffuse.
  The 1D atmospheric model is solved before to get boundary conditions at
  the top of the CFD domain. It is also used for the absorption
  coefficients.
  This work is part of M. Millez PhD, Y. Qu PhD and Y. Maanane internship.
  Also rename verbose mode for Radiative transfer (iimlum -> verbosity).

- Remove "extrag" option for the pressure gradient, as its use was
  limited to orthogonal meshes and was long superseded by the
  `iphydr=1` option.

Architectural changes:

- For parallel IO, allow ranks stepping at the `cs_file` level so as to
  allow using rank steps > 1 with all blocks, reducing the overhead
  of using fewer and larger blocks.

- In-situ postprocessing: add support for ParaView 5.9/Catalyst 2 pipelines.

- Add AmgX library support for linear system resolution on NVIDIA GPU's.
  * Associated matrices should be forced in CSR format.
  * In parallel, Mesh renumbering should be set
    to use `CS_RENUMBER_ADJACENT_LOW`.

- Extend `cs_xdef_t` structure with definitions by DoF (degrees of freedom) for
  advection fields, boundary and initial conditions (CDO schemes). This allow
  one to add more generic/complex definitions at the user level or in existing
  modules.

Release 6.2.0 (2020-08-27)
--------------------------

User changes:

- Volume mass injections can now be defined using zone and
  equation parameter based cs_equation_add_volume_mass_injection_*
  functions.
  * Previous definitions remain compatible but are deprecated.
  * Using the legacy mathod, if at least one zone is defined using
    the `CS_VOLUME_ZONE_MASS_SOURCE_TERM` flag, these zones will be used
    and the first 2 calls to cs_user_mass_source_terms are not needed.

- Add possibility to edit coupling parameters in the GUI using an editor.
  * Editor can be launched from the GUI of any Code_Saturne coupled case
    (CASE/DATA/code_saturne).
  * Editor can be launched from the top-level directory by launching the
    top-level 'code_saturne' launcher using the following command:
    code_saturne [gui|cplgui] run.cfg

- Thermal model:
  * The iscacp(iscal) array is replaced by the field key-word "is_temperature".
    This allows using it also from C code.

- Various improvments to CGNS output:
  * Added output of boundary condition sections (on by default)
  * Allow discarding steady data (to avoid isues with ParaView when both
    steady and unsteady data s present)
  * Use single precision by default for variables (leads to smaller files).

- Add example for computing porosity values from a CAD shape.

- When a non-transient probe set is clipped to cell centers and multiple
  probes located in a given cell, only the one closest to the cell center
  is retained.
  * This is similar to the behavior of the old GUI-defined profile system,
     but is cleaner and more controlled.
  * When the mash changes over time, all probes are kept as they may not
    be in the same cell on the next call.

- Allow activating P1 interpolation for probes and profiles using a simple
  function call: `cs_probe_set_option(pset, "interpolation", "1");`

- Add the possibility of defining fields to output with a given
  postprocessing mesh using the new `cs_post_mesh_attach_field`
  and `cs_probe_set_associate_field` functions.
  * This allows a finer control than the global 'auto_var' options
  * For profiles, the curvilinear and cartesian coordinates may also
    be output using the `cs_probe_set_auto_curvilinear_coords` and
    `cs_probe_set_auto_cartesian_coords` functions.

- Removed the `cs_user_turbulence_source_terms` user subroutine.
  Turbulent source terms can be added using the generic (C)
  `cs_user_source_terms` function.

- Add a code_saturne FMU example with full implementation of the methods
  allowing to operate with it from a master FMU. extras/FMI/README.md provides
  brief indications on how to create a code_saturne FMU. The example is based
  on the simulation of a stratified flow upstream and downstream a T-junction.

- Add new features to control file API aiming at using code_saturne as an FMU
  through socket communications:
  * retrieve a notebook parameter value ("get")
  * send back a message if an iteration has been done subsequently
    to an "advance" message to allow a finer control.

Numerics and physical modelling:

- Convection-diffusion multigrid: restore original aggregation criteria.
  * This criteria seemed to strangely limit the grid hierarchy depth
    on a provided test case, but seems to provide much better performance
    on Rij-epsilon-EBRSM test cases during the initial iterations (with
    no purely diffusive zone, but we assume a solver adapted to both
    convection and diffusion should behave well in convective zones).
  * This solver may now also be selected in the GUI for non-scalar variables.

- Turbulence: remove the vortex method for LES.

- CDO: Add the treatment of the non-linear advection term in
  Artificial Compressibility coupling algorithm with a Picard
  algorithm. This enables to treat incompressible Navier-Stokes
  equations with this coupling algorithm (based on the work performed
  during Riccardo Milani's PhD).

- CDO: Add the treatment of segregation of binary alloys. More
  validations are required but first results are promising.

- CDO: Add _solid cells_ for the treatment of permanent solid walls
  (by opposition to cells in the fluid domain which may become
  solid/mushy/fluid w.r.t. the thermophysical constraints) during the
  solidification or melting process.

- Add a first interface with the MUMPS library to get access to sparse
  direct solvers without the PETSc library (experimental).

Architectural changes:

- Add cs_array.c/cs_array.h for array utility functions.

- For coupled cases, replace `coupling_parameters.py` file by settings
  in the top-level `run.cfg` (see Doxygen documentation for details).
  Cases must be updated manually.

- Documentation: add optional MathJax support for Doxygen. Currently
  uses an external MathJax server, so works on-line only.

- [GUI]: enable Qt translation mechanism.
  Translation files are not maintained, but may be built using  pylupdate5,
  linguist-qt5, and lrelease, and files named code_saturne_<lang>.qm
  under <install_prefix>/share/translation will be handled.

- Use `cs_probe_t` type profiles for GUI-defined profiles.

- Fully replace MEI by generated (MEG) code. This removed the bootstrap
  dependency on Bison and Flex.

- Add pragmas in reference C user file to try to assign weak linkage
  to those files. This seems to work with gcc and clang, and is also
  mentioned in the IBM xl C documentation. This is a cleaner approach
  than allowing multiple definitions and should work on systems whose
  dynamic linkers cause issues with this.

- Documentation: move the installation documentation to a top-level
  INSTALL.md file. This replaces both the matching sections in the
  README.md and the LaTeX-based install.pdf file.

- Replace `SCRIPTS/runcase` with `DATA/run.cfg` logic.
  * DATA/SaturneGUI is replaced by DATA/code_saturne, which can both launch
    the GUI or run the code (using "./code_saturne" or "./code_saturne run").
  * The new logic allows keeping sections specific to various batch systems
    or named resource environments in the "run.cfg" file, as well as better
    align the "run" and "submit" commands, and is more easily extensible.

- Now require Python 3.4 or higher.
  Python 2 is still supported for the moment, but the configure checks
  require Python 3 to "strongly encourage" migration if not done yet.

- Add the possibility to create a desktop shortcut on Linux systems.
  Creation is done using the `make install-shortcut` command after the
  `make install` command. A `code_saturne_($version).desktop` file is
  created in the desktop folder and double-clicking on it launches
  the GUI.

Bug fixes:

- Fix convective outlet operator setting boundary coefficients for vectors
  (with or without anisotropic diffusion) and tensors.

- Major fix in hydrostatic pressure algorithm in the Rhie and Chow filter.
  It only impact cases where external force is not colinear to boundary
  faces where pressure Dirichlets are imposed (e.g. in atmospheric flows
  at the top of the domain with automatic BCs).

- Fix field output numbering on polyhedra or polygons when writer
  subdivision is activated in serial mode.

- Fix computation of Rhie & Chow first term in cases where diffusion in
  pressure correction step is anisotropic.
  This could have an impact on calculations with either head losses,
  tensorial porosity or the pseudo-coupled velocity-pressure solver.

Release 6.1.0 (2020-04-15)
--------------------------

User changes:

- Add lambda0 field for thermal conductivity to the
  cs_physical_properties structure. This is used by the GUI, but
  not yet used for the general setup, where visls0 is still used
  directly, and lambda0 left to 1. Usage of lambda0 instead of
  direct usage of visls0 should be generalized.

- Add the possibility of adding partial mesh modifications after the
  main mesh preprocessing stage (cs_user_mesh_modify_partial).
  This allows separating selected boundary faces (such as symmetry faces)
  from others, so as to ignore them in the main finite volume scheme.

- Add the possibility to keep more than just the last checkpoint file.

- Modify the GUI mathematical expressions syntax to allow only for C syntax for
  if loops.

- Add cs_mesh_remove_cells and cs_mesh_remove_cells_negative_volumes
  functions to allow removal of selected or degenrate cells in preprocessing.

- Add P1 interpolation function and associated example for probes output.

- Correctly handle mixed code_saturne/neptune_cfd couplings in run script.

- Add the possibility to compute a porosity from a file containing
  points (*.pts) coming from a 3D scan. To use it set the scan file
  with "cs_porosity_from_scan_set_file_name("my_scan.pts")" in
  cs_user_parameters.c.

- Multigrid: simplify plotting behavior.
  * Only convergence of cycles is now plotted, as many smoother
    operators do not compute the residual.
  * With a high enough verbosity, the residual after each level's
    smoother or solver application is logged.

- Add the possibility to compute each term of the balance equations of
  the Reynolds tensor and/or the turbulent flux of an additional transported
  variable when performing an LES simulation.

- Merge rough and smooth wall functions.
  User setting of roughness is now using the dedicated fields
  ("boundary_roughness" and "boundary_thermal_roughness").
  Buoyancy modification for atmospheric flows (Louis) are also available for
  dry and humid atmosphere for rough-smooth wall functions.

Numerics and physical modelling:

- Remove former slope test option for vector variable (component by component,
  isstpc = -1).

- Atmospheric model: add 3D radiative transfer for solar wave (shortwaves),
  with direct and diffuse components.
  The 1D atmospheric model is solved first to get boundary conditions at
  the top of the CFD domain. It is also used for the absorption
  coefficients.
  This is based on the PhD work of M. Millez and Y. Qu PhD and on
  the internship of Y. Maanane.

- The GMRES solver may now also be used for vector or tensor fields.
  The maximum number of restarts is set at 40 instead of 75
  (which is still quite high compared to most of the litterature).

- Remove hydrostatic pressure gradient computation (iphydr = 2).

- Wall pressure extrapolation is now considered advanced/deprecated,
  and is not set through the GUI anymore.

- LES: for dynamic Smagorinsky model, add a new filter for cases where
  the full extended cell neighborhood is not available.

- Add tensor gradient clipping.

- Add local (single cell) gradient reconstruction options for
  postprocessing.

- Gradient reconstruction and extended cell neighborhood:
  * allow new options for extended neighborhood reduction.

- Merge Atmospheric module and Cooling towers module air properties.

- Add the capability to solve CDO equations with an advection field
  defined by the mass flux arising from the FV legacy scheme.

- Add CDO edge-based schemes to solve vector-valued equations with curl terms

- Add a Maxwell module to solve electrostatic with CDO vertex-based schemes and
  magnetostatic equations with CDO edge-based schemes.

- Add the resolution of Oseen and Navier-Stokes equations with CDO
  face-based schemes (based on the work of Milani's PhD).

- Remove the files related to the Uzawa velocity-pressure coupling
  algorithm as well as the ones related to the Vector Penalty
  Projection (VPP) algorithm. The uzawa coupling is now merged to the
  monolithic approach.

- Add several new features for solving linear systems arising from a
  discretization of (Navier--)Stokes with CDO face-based schemes:
  Golub--Kahan biorthogonalization, Augmented Lagrangian--Uzawa
  algorithm (in incremental form), MUMPS with the PETSc interface,
  block preconditioning through PETSc

- Add a thermal module based on CDO schemes.

- Add a solidification module relying on CDO schemes (this module
  calls the thermal and the Navier-Stokes modules as well as the
  enforcement of internal degrees of freedom). A specific treatment of
  the phase changing is performed in this module which induces source
  terms and reaction terms in the related equations.

- Add new Hodge operators for vertex-based and face-based schemes based on a
  bubble stabilization or a sub-stabilization.

Architectural changes:

- Allow short path 'salome' for medcoupling, coolprop, metis and
  scotch in installation script.

- Move some parts of the user documentation to Doxygen.

- Remove support for atmospheric aerosol library SIREAM. Replaced by the
  aerosol library SSH-aerosol (EDF - ENPC - INRIA). The latter is released
  under GNU GPL v3, see https://sshaerosol.wordpress.com/ and
  https://github.com/sshaerosol/ssh-aerosol/

- Removed localization (French or English), as it was very incomplete
  and limited to a log file.

- Change the writing process of checkpoint files: when a checkpoint
  already exists, the previous one is renamed before writing the new one.
  The old one is then removed if the writing of the new file is
  successful, thus avoiding possible corrupted files.

- Update Melissa output driver to current Melissa API (c8608e1a6, January 2020)
  * Ensure components of multidimensional arrays are named distinctly, as
    Melissa ignores field dimension.
  * Add trace option to log main operations to file.
  * Add dry_run options to simulate output without Melissa (activates trace)
  * Add rank_step option to allow working with coarser communicators.
  * When a single MPI rank would be used, use the serial API instead.
  * If Melissa is build without MPI, use a ranks step so data is aggregated
    to the first ranks (rank_step = n_ranks) and use serial API.

- Addition of file extensions for both checkpoint (.csc) and mesh_input
  files (.csm). Code_Saturne still accepts old files
  (mesh_input/mesh_output/main/auxiliary/...) without the extensions.
  The GUI now also accepts "mesh_input" as a folder containing
  multiple mesh_input files, while waiting for a full integration of
  the Code_Saturne mesh file format into the GUI mesh handling page.

- All to all operation now all use the cs_all_to_all API, allowing better
  instrumentation and choice of underlying all to all algorithm.

- In field pointer API, CS_F_(e) and CS_FI_(e, 0) are now equivalent.

- code_aster coupling: remove the "milieu" module and move its functionnality
  into the main cs_ast_coupling.c code.

- Simplify handling of SpMV tuning, as various contraints on the type
  of matrix used for a given solver limited to use of the cross-type
  tuning, while that tuning added significant complexity.

- GUI: force SIP API version to 2 for PqQt4. This allows removing
  the "to_qvariant" wrapper function and makes code more readable;
  the GUI will fail if called from external code using PyQt4 which
  has loaded the default API, so compatibility with Salome
  versions 7 or older is dropped.

- Added optional support for CUDA in build system.
  At this stage, only device info is added to system information,
  no compute kernels are added yet.

Default option changes:

- Set clipping at zero for atmospheric chemical species

- Set k-epsilon turbulence models to uncoupled option by default
  (it was already uncoupled by default for all models except standard k-epsilon)

- Change a default setting (GUI) of the compressible model:
  for hydro. pressure treatment at walls (disabled by default now).

Bug fixes:

- Fix reconstruction gradient weighting in computation of face diffusion
  fluxes from pressure increment in pressure correction step. Impacts,
  for example, VoF calculation on non cartesian meshes.

- Fix launching of ncfd doxygen from GUI (help menu).

- Fix crash in ALE using internal structures coupling.

- Fix reading of GUI setting of turbulent flux model for additional model
  variables (specific physics).

- Fix in turbulent isotropic diffusion (Shir model) with EB-RSM. The turbulent
  viscosity did not take into account the effect of the wall.

- Fix to set number of iterations to 0 by studymanager command line (option -n).

Release 6.0.0 (2019-09-26)
--------------------------

User changes:

- Changed Catalyst writer behavior:
  * The name of the writer now only defines the name of the input referred by
    the coprocessing scripts.
  * All Catalyst coprocessing Python scripts present in the DATA (and execution)
    directory are used, irrespective of their name.
  This allows defining multiple pipelines with shared outputs.

- Add the possibility to define in the GUI postprocessing meshes based on the
  volume and boundary zones defined beforehand by the user.

- Lagrangian module: add particle event structure and associated statistics.
  * This allow handling boundary statistics in a manner more consistent
    with volume statistics.
  * Additional boundary statistics may be defined by the user.
  * The major boundary statistics are now handled through this system,
    though some need to be updated in the future.
  * This system could also be used to track volume events in the future.

- GUI: improve global workflow by adding two new functions (buttons):
  * A built-in text editor in the GUI, which allows editing of user routines
    and/or functions inside the SRC folder, as well as opening log files
    in a RESU folder.
  * A compiling test for possible user functions. In case of a compiling error,
    the message is transmitted to the GUI window for the user to analyse.

- GUI: allow disabling paralel IO for MED output writer.

- GUI: significant reorganization
  * Run or submit page moved to pop-up dialog callable from toolbar
  * Folders replaced by active pages (with new icons reflecting this)
  * Preprocessor/calculation modes replaced by run type in mesh page
  * Many minor changes

- GUI: add Volume of Fluid feature to GUI.

- Save mesh_input in restart by default. To avoid using excess disk space when
  meshes do not change, use hard links where appropriate (and move mesh_output
  to checkpoint/mesh_input upon checkpointing).

- Provide clean (advanced) option for deactivation of modified mesh output,
  usable from the GUI.
- Allow stopping criteria based on physical time and/or additional time.

- Add cs_restart_map_set_mesh_input function to allow mapping restarts
  from computation using a different mesh.

- Add option --create-xml to studymanager command allowing to generate a xml
  studymanager parameter file automatically from within a study directory.

- Add the possibility to visualize the turbulent production
  and buoyant terms for DRSM models
  (the user only has to create "rij_production" and/or "rij_buoyancy" field).

- Keep user reconstruction sweeps number in case dynamic sweeping (iswdyn) is
  enabled. Fix for issue #86.

Physical modeling:

- Add an option to disable dynamics in solid cells for internal coupling or for
  porous modelling.

- Add non-linear (quadratic) eddy viscosity model k-epsilon of Baglietto et al.
  * to enable it, set cs_glob_turb_model->iturb = 23
  * it is a Low Reynolds model, compatible with adaptative wall functions.

- Add eddy viscosity model k-epsilon Launder-Sharma
  * to enable it, set cs_glob_turb_model->iturb = 22
  * all y+ wall model is possible (cs_glob_wall_functions->iwallf = 7)
  * low Reynolds model behavior can be forced
    (cs_glob_wall_functions->iwallf = 0)
  * uses a segregated scheme to solve k-epsilon equations system
    (cs_glob_turb_rans_model->ikecou = 0).

- Add Boussinesq approximation as a variable density model. To activate it, set
  cs_glob_stokes_model->idilat to 0 in C or idilat = 0 in Fortran.

- Lagrangian module:
  * Add agglomeration and fragmentation algorithms.
  * A global bit-mask based CS_LAGR_P_FLAG attribute is now always present;
    It includes the features CS_LAGR_DEPOSITION_FLAG, and also includes a
    CS_LAGR_PART_FIXED bit replacing the use of a negative particle cell number.
  * Local zero-based particle cell id is now used instead of a
    signed 1-based cell-number
  * Advanced user functions are added for particle/face interaction, replacing
    the (never implemented) JBORD codes.

- Add a compressible two-phase homogeneous model
  * This model is solved using a fractional step method borrowing mass,
    momentum, energy balance steps from the compressible algorithm (single
    phase).
  * Convection and source terms (relaxation towards equilibrium) step for each
    fraction (volume, mass, energy) follow.
  * Thermodynamic of the mixture is generic, though initial values of some
    iterative processes (secant, dichotomy) used to solve implicit thermodynamic
    relations remain not generic.
  * Each phase thermodynamics follow a stiffened gas EOS which parameters are at
    hand for the user.
  * Relaxation time scale of return to equilibrium is also at hand for the user.
  * The model can be activated by setting to 2 the physical model parameter
    cs_glob_physical_model[CS_COMPRESSIBLE].
  This integrates the work of Olivier Hurisse.

- ALE module: use CDO vertex based numerical schemes (more robust) for ALE
  displacement. To activate it use "cs_glob_ale = 2" in C or "iale = 2" in
  Fortran.

- Make clipping of Reynolds stresses more robust for the coupled solver.
  When there is a negative eigen value, Rij is now clipped with a term of
  return to isotropy (hence preventing increase of turbulent kinetic energy).

Default option changes:

- Disable CS_FACE_RECONSTRUCTION_CLIP bad cells correction by default
  (clipping of face reconstruction distances |II'| and |JJ'|).
  This reduces precision (consistancy loss) and can impair space convergence
  on several verification test cases run on tetrahedral meshes
  (INTERNAL_COUPLING, PERMEABILITY_GRADIENT, PLANE_COUETTE_FLOW).

- Change a default setting (GUI) of the compressible model:
  for hydro. pressure treatment at walls (disabled by default now).

Bug fixes:

- Fix shifted index in imposed pressure outlet (GUI).

- Fix restart from LES to LES (do not add noise on velocity) and
  from RSM to LES.

- Fix for wall boundary conditions with the compressible model / algorithm.
  Always use a homogeneous Neumann on pressure for the mass balance step.
  This fixes mass conservation for cases:
  * with gravity and hydrostatic pressure treatment
  * without gravity and without hydrostatic pressure treatment.
  Note that the default behaviour (GUI) was to enable the hydro. pressure
  treatment at walls. Hence in all cases with this default setting and no
  gravity, this fix has no impact.

- Fix hydrostatic equilibrium option of compressible model / algorithm.
  Hydrostatic pressure was contributing twice to right hand side in mass
  balance step.

- Fix missing variable change for computed pressure with EVM in a
  postprocessing utility function.

- Fix for turbulent heat fluxes computation when Boussinesq hypothesis is used.
  Variance was not taken into account.

- Fixes for GMSH import.
  * Fix crashes for binary file import
  * Fix handling of physical groups for GMSH V4 format; names are also handled.

- Atmospheric model: fix issue combining meteo profiles and local injections.

- Lagrangian module: bug fixes in the deposition model.

- Compressible: fix density time scheme in transported passive scalar/vector
  balance to ensure conservativity with compressible algorithm.

- Fixes in time stepping:
  * Fix momentum conservation for the variable density deprecated algorithm
    with mass flux prediction (option idilat>=2, ipredfl=1). The density
    which is coherent with the mass balance needs to be updated right after
    mass flux prediction step.
  * Fix time stepping for buoyancy term. Compute it with last density value
    (EOS) when using first order time stepping, to avoid a time shift.
  * Fixes for time scheme of density in unsteady term of momentum, transported
    scalars / vectors balances, especially use density at time n in momentum
    balance unsteady term if mass flux is predicted.
  * Use EOS density in momentum balance unsteady term for the weakly variable
    density algorithm (idilat=1). Use it also in correction step to be coherent.
    This algorithm is now the same as pre v5.3 changes on time stepping.
  * Update density which is coherent with mass balance after velocity update for
    low Mach algorithm.
  * Use same time schemes for vector transport equation as for scalar transport
    equation.

- Fix time stepping of VoF / Cavitation algo. (density in unsteady term of
  momentum equation) and ensure variable density indicator (irovar) is set to 1
  when VoF / cavitation algorithm is enabled. Idem for variable viscosity
  indicator.

- Fix update of density in VoF / Cavitation algorithm in case the convection
  scheme for the void fraction is not upwind. Update was not ensuring proper
  mass conservation at solver precision.

- Fix time stepping in pressure gradient weighting for VoF / Cavitation algo.

- Fixes for parallel runs in sedimentation source term with humid atmosphere
  model.

Numerics:

- Change the way the dimensionless wall distance is computed (LES Smagorinsky
  model) i.e. convection flux is a face gradient and is now computed using a
  two point flux approximation (TPFA). Set solver precision to a higher value
  since there is no need for a high precision (1e-8 to 1e-5).

Architectural changes:

- Allow sourcing environment before launching with the
  --with-shell-env configure option.

- Scripts: add library dependency paths to LD_LIBRARY_PATH for better
  robustness when "rpath" does not have priority.

- Update CGNS support for CGNS 3.4 compatibility.

- Preprocessor: update Gmsh reader to handle GMSH v4.1 format.

- Handle mathematical expression in GUI by generating corresponding C code
  then inserted in a cs_meg_..._function.c file compiled with other run
  sources. This drastically improves performance when using MEI.

- Move Reynolds stress tensor transformation matrix (alpha in clca66)
  computation to C. C translation is taken from NEPTUNE_CFD.

- Move part of VoF module from Fortran 90 to C.

- Remove VOFI (VoF initialization) library detection as it is not used anymore.

Release 5.3.0 (2018-10-26)
--------------------------

User changes:

- GUI: when the NEPTUNE_CFD module is available, the GUI
  can now switch directly between Code_Saturne and NEPTUNE_CFD
  setups.

- Move velocity-pressure algorithm settings from global numerical
  parameters view to time step view in GUI. "Time step" view is
  subsequently renamed "Time settings".

- For studymanager, add default destination and repository directories
  if process is launched from a study directory. Destination directory
  is based on study directory name (RUN_study_name).

- Boundary layer insertion: added optional cell volume ratio limiter
  to reduce the extrusion near cells that would be excessively
  flattened or entangled.

- For new cases, the default postprocessing output writer uses the
  'separate_meshes' option, which avoids requiring the 'extract blocs'
  filter on ParaView.

- Add mesh refinement engine.
  * Coarsening not available at this stage.
  * Load balancing currently handled through complete repartitioning.

- Coal combustion:
  * The CO2 transport kinetic model may now be set with the GUI.
  * Fuel definitions must now always be done through
    the GUI; new cases do not include a reference dp_FCP.xml file
    anymore, though older such files may still be read.

- Allow short path 'salome' for hdf5, cgns, med in installation script.

- GUI: steady/unsteady computation type now handled with time step.

- Allow joining of meshes including isolated faces.

- Add postprocessing of temperature and flux at internal coupling interface.

Physical modelling:

- Make particle tracking compatible with transient turbomachinery model.

- Add a new continuous "all-y+" 2-scale wall model (iwallf=7) available
  with the EB-RSM and set by default for this model:
  * Ensure convergence towards standard EB-RSM when mesh is refined.
  * Degenerate in a SSG-like model on high Reynolds meshes.
  This was adapted from developpements done in J.F. Wald PhD.

- Major modification for K-omega SST (iturb=60) boundary condition.
  * Switch from a Neumann boundary condition to a Dirichlet
    boundary condition on omega.

Numerics:

- Use left anisotropic diffusion scheme (legacy FV) for mesh velocity solving
  in ALE framework.

- Major change in the time stepping to ensure 2nd time order for
  variable density flow if 2nd time order is activated.
   It impacts all the second time order time stepping and also all variable
   densities algorithms.
  * If you want to go back to the previous algorithm for variable
    density, you can specify ipredfl = 0 in the usipsu
    (cs_user_parameters.f90) subroutine.
  * The momentum equation is staggered in time, that is to say, when
    2nd order is activated, velocity is solved from time n-1/2 to n+1/2.
    A special care should be done for time averaged quantities.

- Porous modelling: adapte the numerics to discontinous porosity.
  * The velocity is interpolated at faces using mass conservation and the
    momentum is corrected so that the steady state of Euler equations is
    retrieved.
    This can be activated using iporos = 3, the improved hydrostatic treatment
    will then be activated.
    This was developped in the PhD of C. Colas.

- Improvements in mesh quantity computations.
  * Previous face center adjustment for volume removed. Adjustment
    used in versions 1.1 to 5.2 may be restored using
    the cs_mesh_quantities_face_cog_choice function.
  * A refinement of the face center computation (for warped faces)
    may be activated using the CS_FACE_CENTER_REFINE mesh quantities
    computation flag.
  * An option to compute cell centers based on the actual center of gravity is
    available. It can be enabled using the cs_mesh_quantities_cell_cen_choice
    function. The method based on face centers is kept as default.

- Added dispersion modeling option to DOM radiative model:
  * May be activated by setting
    cs_glob_rad_transfer_params->dispersion to true.
  * The associated dispersion coefficient may be set by modifying
    cs_glob_rad_transfer_params->dispersion_coeff.

- Added K-cycle multigrid type as an option.

- Added new multigrid coarsening algorithm compatible with matrices
  in MSR and native formats, usable with CDO/HHO matrices in addition
  to 2-point FV matrices.
  * Removed rarely used and inconclusive face traversal options from
    previous (default) coarsening algorithm.

- Add choice of inexact (flexible) preconditioned congugate gradient.
  * When using multiple threads and multigrid preconditioning with a
    Gauss-Seidel smoother, this choice is set by default over the
    standard preconditioned congugate gradient.

- Add vector-valued Laplacian for HHO schemes (case k=1). Based on
  Daniel Castenon's work.

- Add vector-valued Laplacian for CDO vertex-based schemes

- Add the steady Stokes equations with CDO Face-based schemes (based
  on the work of Riccardo Milani).
  * Velocity-pressure coupling is handled thanks to an Uzawa-Augmented
    Lagrangian algorithm

- Several new features for scalar-valued CDO Face-based scheme
  * Add an advection term: Upwind scheme in non-conservative
    formulation (Joint work with Hanz CHENG (Univ. Monash, Australia)
  * Add unsteady and reaction term for scalar-valued CDO Face-based
  scheme: Theta-scheme (including implicit and Crank-Nicolson)
  * Add a weak enforcement of Dirichlet boundary conditions (Nitsche
  technique) and also its symmetric version (based on Riccardo
  Milani's work)

- Add enforcement of internal degrees of freedom in scalar-valued CDO
  Vertex-based schemes

- Add Robin boundary conditions for scalar-valued CDO Vertex-based and
  Vertex+Cell-based schemes

Architectural changes:

- Move mesh velocity solving (ALE) from Fortran to C.

- Remove dependency to the libxml2 library.

- Preprocessor: update Gmsh reader to handle GMSH v4 format.

- Add "--disable-backend" configure option to build and install only
  front-end and documentation.

- Add CGNS writer "links" option to write mesh data in a separate file,
  mapped transparently in the main file through CGNS links.

- Add the possibility to create a CFD/Cathare2 coupled case using the
  "--cathare <name>" similarly to CFD/Syrthes coupling. This coupling
  is only functionnal for NEPTUNE_CFD at the moment.

- Add "--with-vofi=path_to_vofi_build" configure option to build and install
  with the VOFI (VoF initialization) library (only for some research uses;
  that library's licence is too restrictive for most uses).

- Add the possibility to conduction Uncertainty studies using OpenTURNS within
  the salome_cfd platform. The full workflow, Code_Saturne case creation to the
  UQ study, is available within the salome_cfd GUI. Code_Saturne unitary
  evaluations can be run on either the local workstation or distant machines
  such as computing clusters.

Default option changes:

- Change default options for bad meshes:
  CS_BAD_CELLS_WARPED_CORRECTION, CS_FACE_DISTANCE_CLIP, CS_FACE_RECONSTRUCTION
  are switched on by default.

Bug fixes:

- Fix initialisation of some turbulence constants. Constant/Models concerned:
  * sigmae (for all RSM models)
  * csrij (EB-RSM model).

- Fix CGNS reader so as to handle cases with unordered sections.

- Fix allocation size for values at injection in case of mass source terms
  with coupled Reynolds stress solver.

- Minor bug fix updates to Melissa writer.

- Fix face external force projection with tensorial diffusion and porous models 1, 2.
  This was impacting cases with head losses, improved pressure interpolation, and
  scalar or tensorial volume porosity models (iporos=1, 2).

- Fix in the Lagrangian particle tracking:
  * some minor inconsistencies were introduced on wraped faces
  * local normal (of the crossed sub-triangle of the crossed face) is stored for
    particle rebound at the boundary
  * if the maximum number of sweeps is reached, the put the particle at the last
    cell center.

Release 5.2.0 (2018-03-30)
--------------------------

Physical modelling:

- Atmospheric model: add new BCs for open-boundary flow such as atmospheric
  flows. It consists in changing the solved pressure. This change result
  in adding a (constant over space) momentum source terme. This momentum
  source term is designed to reach a target bulk momentum.
  The solved pressure has a zero Dirichlet as BC in case of constant
  density. In case of dry atmosphere, an hydrostatic pressure is
  pre-computed and imposed as BC.

- Turbomachinery module:
  * Improved legacy coupling-based turbomachinery option, frozen rotor
    is now also possible.
  * Allow time-varying rotor velocity for turbomachinery.
  * Account for Smirnov and Menter (ASME, 2009) modifications of Spalart-Shur eddy
    viscosity correction for rotation and curvature.

- Add internal coupling for vector quantities (e.g. for plasma - weldpool monolitic
  coupling).

- Make compatible Volume Of Fluid algorithm and option iphydr, improving face
  pressure interpolation by taking into account hydrostatic pressure gradient.

Numerics:

- Add the possibility to put buoyant scalars and density update in
  velocity-pressure loop.

- Integration of a vector Laplacian for CDO-Fb schemes.

User changes:

- Allow stopping at a given wall-clock time using the control_file.

- The GUI enforces the 'runcase' name for the script. This allows more robust
  handling of the runcase/XML file link, as one of the 2 names is fixed.

- Mesh-joing based turbomachinery model robustness is also improved by
  allowing multiple joining retries using a slightly perturbed position.

- Refactored extrusion handling to allow easier and finer control of extrusion
  vector generation.
  * a cs_mesh_extrude_face_info_t structure allows simple settings for each
    boundary face, including number of layers, thickness, and expansion factor,
    and optionally a start and end thickness. Utility functions allow
    "per zone" settings.
  * a cs_mesh_extrude_vectors_t structure defining local extrusion vectors,
    which is built automatically based on the face settings but can be
    modified by the user.

- Added boundary layer insertion at selected faces, using a mesh deformation
  (CDO vertex-based), followed by boundary layer extrusion.

- MEDCoupling: improve detection and add MEDLoader support. This allows e.g.
  to use MED for interpolating from one mesh to another.

Architectural changes:

- Add Support for reading CGNS files with NFACE_n elements.

- Partly unify boundary and volume zones, renaming the n_faces/n_cells
  members to n_elts and face_ids/cell_ids to elt_ids.

- Complete cs_all_to_all API to handle indexed cases.

- Use matrix assembler instead of matrix extension function for
  internal coupling, to allow for a broader range of linear system
  solvers and preconditioners.

- Add "Melissa" type writers for coupling with in-situ statistics
  Melissa server.

- Generate MEI parser/scanner upon bootstrap. This assumes Bison and Flex are
  available (but previous generated files are kept otherwise), and reduces
  dependencies for configure and build.

- CFD_Study: translate from Python2 to 3.

Bug fixes:

- Fix spurious mapping of previous mass flux with order 1 in time scheme.

Release 5.1.0 (2017-10-10)
--------------------------

Physical modelling:

- Turbulence: add DDES tubulence model
  * Activated using k-w SST model (iturb = 60) with keyword IDDES = 1
  * Hybrid Numerical scheme CD/SOLU can be activated with ISCHCV = 3

- Atmospheric model: add an atmostpheric infra-red absorption feature.
  Based on L. Makke PhD work.

Numerics:

- Atmospheric model: added partial radiative transfer properties model.

- Make supersonic outlet boundary condition implicit (compressible module).

- Added fallback mechanism to allow less stable but faster Krylov solvers
  (such as BiCGstab) to switch temporarily to GMRES in case of
  convergence error.

- Switch from Jacobi to processor-local symmetric Gauss-Seidel as default
  linear solver type for non-symmetric systems.

- Simplify selection of local symmetric Gauss-Seidel linear solver variant.

- Make vectorial and tensorial slope tests invariant by rotation and
  identical for all components (previous version available for the
  velocity with isstpc = -1 in var_cal_opt structure).

- Add an option to blend with upwind when the slope test is activated. The
  key word is blend_st in var_cal_opt structure.

Default option changes:

- Make coupled solving for DRSM turbulence models default option.

- Make prediction of mass fluxes before momentum solving optionnal. This
  impact variable density algorithm (idilat>=2).

- Make advanced turbulence inititalisation for k-omega, Bl-v2/k and
  EBRSM turbulence models default option.

- Slope test on vector (and tensor) is now invariant by rotation by default

User changes:

- Change low-level handling of mesh extrusion to allow easier fine
  control and prepare for boundary layer insertion.

- When importing CGNS files and requiring additional groups for mesh
  sections or zones, do not try to build boundary face groups based
  on vertex boundary conditions.

- Disable cases instead of exiting when a run directory check fails in
  studymanager.

- Separate compilation test and repository update features and add tags
  feature to studymanager:
  * --test-compilation or -t to test compilation for all cases
  * --update or -u to update GUI parameter files and set paths in scripts
    SaturneGui and runcase
  * --update-xml or -x to update GUI parameter files only
  * --with-tags to only process runs with all passed tags
  * --without-tags to exclude runs with one of the passed tags.

- Legacy coupling-based transient turbomachinery option is now defined in a
  manner similar to join-based coupling, simply replacing simply the joining
  by an internal coupling. Coupling is internal to a code instance, requiring
  a single data setup, and switching from a joining to a coupling-based
  scheme is now easy.

- Changes in Lagrangian module boundary and injection definitions.
  * renamed boundary particle classes to sets to avoid confusion with
    statistical classes.
  * some API changes in user functions; the cs_user_lagr_new_p_attr function
    is removed in favor of the more general cs_user_lagr_in.
  * boundary zones are now based on the general zones, as defined using the
    GUI or cs_user_zones.c.
  * definition of boundary injections is now more streamlined. Complex coal
    (non-raw) coal definitions, which were little tested, can be done using
    cs_user_lagr_in.
  * mesh-based profiles (injection rate) may be defined using user functions.
  * volume injections may now be done using the same mechanism.

- Removed handling of VTK-based scalar slice outputs in studymanager;
  graphical files may now be included, and using Catalyst to produce
  images is recommended.

- Removed "interp1d" (1d-D interpolation) option from MEI (Mathematical
  Expresion Interpreter) syntax, as it was poorly documented and
  added complexity; a more general interpolation framework, usable also
  from user functions, is planned for the future.

- Add option for postprocessing output of linear solver residuals.

Architectural changes:

- Add "histogram" type writers for easy output of variable histograms.
  Options include text, TixZ, or PNG (when Catalyst support is available).

- When coupling with Syrthes, source Syrthes library path changes only locally.

- Improve Catalyst environment handling for better automation with SALOME.

- Default to PyQt5 instead of PyQt4 when both are available.

- Remove support for restart files older than version 4.0.

Bug fixes:

- Fix broken path in update feature of studymanager.

- Fix possible crash or error when computing wall distance for 2D
  (single cell layer) meshes.

- Lagrangian module: fix characteristic time using temperature.

- Fix for issue #236 (crash on restart using periodicity).

- Fix for idilat>3 algo. and isstpp>1 (beta or NVD/TVD limiters) compatibility.

Release 5.0.0 (2017-06-02)
--------------------------

User changes:

- Allow user definitions of head losses by zones (see cs_user_head_losses.c).

- Do not relax k and epsilon by default when using Bl v2 k turbulence model.

- When inserting boundaries (shuch as thin walls), vertices are now separated
  on each side of those boundaries, for a true topological mesh modification.
  cs_user_mesh_thinwall is also renamed to cs_user_mesh_boundary.

Numerics:

- Improve time stepping of the coupled solver for the components of Reynolds
  stress tensor (irijco = 1).
  * No clipping on Rij should appear when no flux rectonstruction is needed and
    an upwind scheme is used (except those due to boundary conditions on R12).
  * Add eigen_max function to cs_math.c
  * EBRSM and SSG models are treated.
  * Add a clipping indicator which can be activated to check where clipping
    occurs in the output files

- Add and reorganize convection schemes (especially for VOF method)
  * reorganize Roe-Sweby limiters.
  * add Compressive convection schemes for accurate modeling of free surface
    with VOF method.

- Add a conservative gradient reconstructed with least squares gradient for
  vectors and scalars, in place of iterative gradient initialized by least
  squares gradient (imrgra=4,5,6).

- Add a cooling towers module to simulate exchange zones.

- Add a convection-diffusion equation solver for additional vector variables.

Physical modelling:

- Add a Volume Of Fluid algorithm:
  * partly merge cavitation algorithm with VOF algorithm.

- Add data assimilation feature to atmospheric module
  (optimal interpolation and nudging):
  * copy LU utilities to cs_math and keep static inline version of them
    in cs_sles_it
  * an optimal interpolation structure is created
  * interpol grid and measures set structures are used as well
  * multidimensional analysis are computed for multidimensional variables.

- Modification of the LES dynamic Smagorinsky clippings.

- Add internal coupling for scalars of two domains (for instance temperature
  between solid and liquid, or enthalpy for electric arcs between plasma and
  weldpool).

- Add sorption model treating non-equilibrium between solid and liquid phases
  to the ground water flow module.

- Lagrangian module:
  * add deposition and resuspension models on internal faces.
    The user can the impose the motion of deposited particles. If integral
    approach for porous modelling is set up (iporos=3), then the internal fluid
    section is reduced by particle deposition.
  * injection is now pseudo-continuous when injecting at every time step.
    To revert to the previous behavior, the CS_LAGR_RESIDENCE_TIME value
    must be set to 0 for newly injected particles.
  * 2 new attributes, CS_LAGR_TR_TRUNCATE and CS_LAGR_TR_REPOSITION, may
    be used to visualize particles with trajectory errors, rather than
    remove them. Particles are now only removed when "completely lost",
    which should never happen.

GUI changes:

- Add handling of multiple compute builds through the GUI.

- Add mesh extrusion capability to the GUI.

- Add a notebook to add global variables to be used in the GUI (such as variables
  in the physical laws).

- Add fans modelling (represented as momentum source terms) in the GUI.

- Add "mapped inlet" boundary condition (for recycled inlets) in the GUI.

- Add "imposed pressure" outlet boundary condition in the GUI.

- Add the verbosity mode for transported variables.

- Add "iterative process error estimators" in the GUI in "Volume solution contol".
  pannel.

- Add automatic transported scalar balance by zone and pressure drop by zone
  in the GUI.

- Physical properties files for the electric modules now use a fixed name,
  and are not handled by the GUI. If a local file is not present, the reference
  file will be used.

Architectural changes:

- Renamed Autovnv to studymanager.

- Remove dependency to Qt of GUI XMLengine.Case class, used by studymanager.
  This dependency was introduced in version 3.0.

- Paths needed for Catalyst are now sourced automatically in most cases;
  for movable builds, a CATALYST_ROOT_DIR environment variable may be used.

- Added general structures to handle boundary and volume zones.

- Extended run staging and submission to coupled cases. Temporary status file
  is now in execution ad not scripts directory.

- Add handling of lmod-based environment modules.

- Allow detection of HDF5 libraries named libhdf5-shared and libhdf5-static as well
  as libhdf5 (when built with CMake)

- Add functions for definition of mesh groups during mesh preprocessing.
  * moved most group handling functions to cs_mesh_group.*.

- Replaced icond keyword by icondb, icondv to allow to enable wall condensation
  and condensation on internals at the same time.

- Probes output activation is now based on field "post_vis" keyword, and does
  not allow fine-grained per-variable probes selection anymore (this being little-used,
  and feasible through use of additional probe sets).

- Prepared for cs_user_extra_op and cs_user_parameters translation to C in VnV base:
  * removed Fortran arrays corresponding to members of cs_var_cal_opt_t C structure
  * deployed use of get/set var_cal_opt everywhere
  * added cs_post_util functions (rij post-pro for EVM models on a subset of cells)
  * moved izfppp to C (cs_glob_face_zone) and map Fortran izfppp to it
  * added a cs_user_output function in cs_user_parameters.c (usipes equ.)
  * added a gas mix model structure in C

- Removed all the remaining mappings between fields and propce array:
  * removed propce arrays
  * removed ipproc array
  * removed nproce
  * removed iprpfl indirection array between field indices and properties
    numbers; iprpfl is kept for compatibility but is just an identity
    function (hence still known in user subroutines);
  * irom, iviscl, etc.. are now directly field indices (hence starting from 0)
  * test on variability of specific heat (icp), isochoric specific heat (icv)
    and volumetric viscosity (iviscv) are consequently shifted of 1
    (-1 : uniform field, >=0: non uniform)
  * renamed _owner utility subroutines
  * renamed most used add_property_field into add_property_field_1d
  * mesh viscosity now a 3-dimensional field when strictly orthotropic
  * removed iroma, ivisla, ivista, icpa variables; previous values of these
      fields now accessed by field_get_val_prev subroutine.

- Removed global_row number from matrix structures, as it was never used,
  and could cause confusion with global ids such as those used with
  range sets (which are realy used in some cases).

- Added range set structure, to ease operations related to handling of an
  owning rank for distributed entities.

- Added matrix assembler structure, to allow more general assembly of
  matrices, especially in parallel cases with rows matching elements
  shared on ranks (needed for parallelization of CDO operators).

- Refactored parallel numbering for space-filling curves and added
  a numbering generation based on a 1D series of real values.

Bug fixes:

- Fix missing reading and writing of mass fluxes in checkpoint/restart.

- Allow mesh joining when no vertices are modified.

- Allow use of Modak absorption using the GUI for gas combustion.

- Fix turbulent exchange coefficient of condensation in case of forced
  convection (the wall function computes Sch*y+/Y+ and not y+/Y+).

- Fix Neumann boundary condition for vectors with anisotropic diffusion.
  A priori, could be used only if the thermal turbulent flux was solved.

- Fix initialisation of flux and divergence boundary coefficients for coupled
  vector variables (field with dim=3).

Release 4.3.0 (2016-07-29)
--------------------------

User changes:

- Allow user activation or deactivation of parallel domain
  visualization output (cs_post_mesh_set_post_domain function).

- Automatic computation of head loss for a given volumic zone.

- When using boundary condition zones defined by the GUI, error reporting
  for overlapping zones is improved, with associated visualization.

- When extruding boundary faces, the group class id of interior faces
  previously on the boundary may by kept as an option.
  This is the default using the GUI.

Changes:

- Turbomachinery module: for better restart behavior, the joined mesh
  is now also handled using checkpoint/restart.

- GUI: remove titles from profile definitions, as these were only handled
  by the obsolete Salome VISU module.

- Added some turbomachinery post-processing utility functions.

- Allow coupling of radiative transfer with 1d wall thermal module.

- Extend automatic postprocessing output to fields defined at vertices.

- Add boundary mesh extrusion function callable from user subroutines.
  (in case of periodicity, rebuilding the periodicity in a later
  preprocessing stage may be necessary).

- When calling cs_selector_get_b_face_num_list or
  cs_selector_get_i_face_num_list, it is no longer necessary to build
  selection mechanisms during the optional mesh modification step.

Numerics:

- Add optional multigrid solver for scalars with convection and diffusion,
  based on the PhD work of Sana Khelifi.

- Lagragian module: new trajectory algorithm which does not loose particles even
  for warped faces (cs_lagr_tracking).

- Improve robustness of Boundary conditions for wall functions of scalars.
  Mainly impact Atmospheric flows.

- Add many convective flux limiters for convective schemes:
  * Fix in Convection scheme when the density is variable and the theta
  scheme is 1/2 so that the scheme is now conservative over time.
  * add inline functions for vectors
  * simplify non-relaxed inline functions
  * use inline functions for temperature explicit balance (bilsct)
  * Add an add hoc limiter, called Beta, which ensures that, for any
  convective scheme and any time stepping, that the variable remains
  between "min_scalar" and "max_scalar", to be given by the user.
  * Implementation of the Roe-Sweby Limiters for all the convective
schemes (Original SOLU, CENTER, Standard SOLU)
  * Creation of a new keyword: kst, for key source term;
  we use it to add source terms in the Beta limiter

- Handling of linear solver errors is improved to allow re-try with
  different solvers or parameters.

- Default to multigrid-preconditioned conjugate gradient for
  most symetric linear systems.

- Improve performance and possible precision of reduction operations.
  * improve OpenMP load balancing of SuperBlock reductions
  * add optional Kahan compensated sum-based BLAS operations;
  * add cs_blas_set_reduce_algorithm() function to choose algorithm

Architectural changes:

- Added parallel IO for MED output when underlying MED library is built
  with parallel (HDF5) IO support.

- Extend control file mechanism with Python module communicating through
  sockets (which may in turn be extended by an RPC mechanism); the
  commands are currently limited to connecting, advancing a time step,
  and disconnecting (in addition to control_file commands), but may
  be extended.

- Simplify MEDCoupling output writer so as to generate MEDCoupling meshes
  on the current partition rather than aggregating them to the root rank,
  as this should soon be directly usable by a VTK adaptor.

- Convert Lagrangian and radiative module implementations to C.
  This affects the associated user subroutines.

- Add basic (non-GUI) support for OAR resource manager.

- Add mesh adjacencies structure to allow direct access to a cell's
  neighbor cells and boundary faces. This allows looping on
  cells, as when computing a SpMV product using a CSR matrix,
  rather than looping successively on faces.

- Improvements to postprocessing writer handling.
  * add 'plot' and 'time_plot' writer types;
  * time plot flushing behavior is now determined by the
    cs_time_plot_set_flush_default() function, replacing the
    Fortran 'nthsav' and 'tplflw' values;
  * writers can now use the 'separate_meshes' option to create a
    separate format-specific writer per mesh.

- SALOME version required for MEDCoupling is now version 8.0
  in this version, MEDCoupling is a full standalone library and may
  be used separately from the MED module.

Release 4.2.0 (2015-12-23)
--------------------------

User changes:

- Automatic initialization of the Turbulence for EBRSM and k-omega models.
  From a reference velocity (uref), the turbulence profiles are reset
  after the first iteration. The velocity magnitude is also changed
  so that the Reichard profile is imposed next to walls.
  Activate it with reinit_turb=1 (in usipsu).

- Merge the bad cell and the mesh quality criterion for offsetting.

- Merge general boundary temperature handling with the radiative
  "wall temperature", for unified logging and post-processing.

- Added optional saving of scalar variable boundary values as fields
  (also done for temperature when a property).

- GUI: "check mesh" option is replaced by a new "preprocessor view": when
  building a new case, the GUI only shows sections relative to mesh
  selection and preprocessing, showing only the steps necessary up to
  preprocessing. "Tools"  menu entries and toolbar icons allow switching
  from the preprocessing mode to the computation mode.
  This new view encourages users to use the "preprocessing" run type,
  which was previously available at the end of the GUI options,
  instead of the "check mesh" option, which did not handle batch runs
  or user subroutines.

- Add boundary cell thickness computation to mesh quality criteria.

Physical modelling:

- Add gas mix thermo. and fix mass source term with the compressible module
  (internship of A. Menasria):
  * allow to select both gas mix (igmix) specific physics with compressible
    (icompf) spec. physics
  * add property field for deduced mass fraction (iddgas) and mixture molar
    mass (igmxml)
  * add a field pointer mapping for gas mix (mapping mix molar mass field)
  * add mapping of specific heat at constant volume field pointer
  * add one gas mix composed of helium, N2 and O2 (i.e. Helium+Air), O2 is the
    deduced species
  * add Sutherland behavior law for viscosity and thermal conductivity of
    gas mix (ivsuth option)
  * move setting of ieos from uscfx1 to usppmo
  * add thermo. relations for perfect gas mix in compressible thermo.
    (also for BCs),
  * igmix specific physics triggers switch of ieos to 3
  * treat conduction term in total energy equation correctly in gas mix
    (Cv depends on species)
  * fix mass source term in compressible (momentum source term coming from
    added mass has to be accounted for in simplified momentum equation during
    mass balance step).

- Add precipitation/dissolution modeling for particle tracking.

- Add the added-mass term in Lagrangian particle tracking (activate it
  with iadded_mass=1).

- Add FSCK radiative model for coal combustion.
  Activate it in cs_user_parameters.f90 with imfsck=1.

- Add a wall function for the velocity based on scalable wall function
  which is valid for both rough and smooth walls (activate it with
  iwallf = 6 in cs_user_parameters.f90).
  Moreover, the continuous wall function based on Van Driest (iwallf=5)
  is extended to Eddy Viscosity Models, and the use of roughness is allowed.
  For both wall functions, roughness must be specified in the field nammed
  "boundary_roughness" in cs_user_boundary_conditions.f90 for instance.

Numerics:

- Add new preconditioning management layer to allow for use of multigrid
  as a preconditioner.

- Extend weighted coefficients for lsq gradient calculation in case of
  anisotropic diffusion (this option can be enabled with the option
  iwgrec(ivar) = 1).

- Add a coupled component solver for symmetric tensors for the Reynolds
  stress model (irijco=1).
  * convection diffusion brick for tensors
  * add a gradient of a tensor field
  * only iterative methode up to now, least square gradients not yet
    available
  * add the advection/diffusion matrix for a tensor field
  * adapt matrix tuning fill and Jacobi and other linear solvers
  * update rotation functions
  * restart from an uncoupled Rij tensor is possible
  * Daly Harlow and Shir models are available

- CDO module:
  * Add reaction term to CDO module.
    Improve the computation of discrete Hodge op. when using
    conforming reconstruction function.
    Adapt the computation of the source term when using this kind of discrete
    Hodge operator.
  * Add time scheme and make some factorizations in vertex-based schemes.
    Add functionalities to cs_locmat_t structure.
  * Add convection for CDO vertex-based schemes.
    Add advection field along with its post-processing.
    Two advection operators are available related to the conservative or
    non-conservative form of the advection operator.
    Add post-processing of an equation (i.e. that of the related field).
  * Add different ways to enforce boundary conditions for CDO vertex-based
    schemes. Weak (Nitsche) enforcement, penalization and strong enforcement.
  * Add a domain and an equation structures
    Move cs_param_eq_t structure as a member of the cs_equation_t structure
  * Add the computation of the wall distance for CDO vertex-based schemes
  * Add the computation of the wall distance for face-based CDO schemes
  * Change the way user sets the paramters of a computation
  * Integrate PETSc solver calls
  * Add face-based CDO schemes for pure diffusion problem.

Architectural changes:

- Add --compute-build option to "code_saturne run" command to allow choosing
  one of several compute builds at runtime.

- Add "code_saturne submit" command to submit a batch job.
  The GUI now uses this to prepare data upon submission.

- Some renames and simplifications in matrix fill types:
  * to make things more generic, CS_MATRIX_33_BLOCK* is replaced by
    CS_MATRIX_BLOCK*
  * tuning for blocks assumes 3x3 blocks, but a specific CS_MATRIX_BLOCK_D_66
    type is added (as a subcase of CS_MATRIX_BLOCK_D) to allow for specific
    tuning for RSM coupled component matrixes.

- When using old (non-parallel) algorithm to compute wall distance, save
  distance in same array as for main algorithm.

- Implement handling of point and element tags so as to ignore location
  on elements sharing tag.

- Simplify Lagrangian particle tracking, removing a linked list.
  This is transparent for the user, except for a possible differences
  in particle numbering.

- Major changes to all_to_all API, so as to make its usage easier.
 * new API is similar to cs_part_to_block / cs_block_to_part, so
   as to be aligned more with sparse or neighborhod collectives
   then Crystal Router in the future.
 * place Crystal Router code in separate files from cs_all_to_all,
   to allow lower-level use. Also simplify its usage, and make its
   data access logic higher level.

- Added cs_debug_wrapper.py debugger wrapper, to assist in running the code
  under a gdb-based debugger, under Valgrind, or a combination of the two.
  * replaces the "valgrind" only debug in the GUI advanced options for execution.
  * may also be used in standalone mode.

- Do not build with BLAS by default, as only MKL is used outside of unit tests, and
  using it requires providing its path to "--with-blas" anyways.


Release 4.1.0 (2015-07-31)
--------------------------

User changes:

- Renamed 'efforts' to 'stress', which should be less confusing.

- Add simple fan effects modeling as explicit momentum source term in
  regions defined by fan characteristics.

- Add a code icodcl=11, allowing to easily impose a boundary value
  consisting in an affine function of the next cell value
  (used at wall faces in the compressible module).

- Add timer statistics and plots for different stages and operators.
  * Base timers for each time step (with initialization as time step 0)
    are now available in the timer_stats.csv file (or timer_stats.dat
    if the default format is changed).
  * Final performance data is also moved from listing to performance.log.

- Allow condensation model to use zone-based definitions. The examples
  are updated so as to recommend the zone-based setup, though previous
  single-zone setups definitions remain compatible.

- Add "--import-only" option to "code_saturne create" command so as to
  rebuild SaturneGUI and runcase scripts for a case which was copied
  from a different system.

- Add higher level functions for turbulent boundary condition settings.
  This allows moving tests on the current turbulence model inside
  the user-callable functions, for more concise and safer programming.

Physical modelling:

- Add 2-scales wall function (with V. Driest mixing length) and its consistant
  wall function on scalars (keyword iwallf in the doc.).

- Changes in Lagrangian Particle tracking.
  * modification of multi-layer model (iclogst = 1)
  * compute porosity in each cell based on mean deposition height only
    at boundary faces (porcel.f90)
  * influence of deposited layers on the flow (iflow = 1 option); addition of
    laghlo.f90 to compute the head loss due to deposition based on porosity)
  * remove lagbar.f90 (energy barrier now computed by cs_lagr_barrier for
    smooth wall/irough = 0, and cs_lagr_roughness_barrier for rough wall/
    irough = 1).

- Add stiffened gas thermodynamic law for the compressible module (ieos=2).
  Two new parameters are accessible in the user subroutine uscfx2:
  gammasg ("pseudo" specific heat ratio) and psginf (infinite pressure).

- Add ADF models for radiative transfers.

Numerics and linear solvers:

- Add vertex-based CDO (Compatible Discrete Operator) numerical scheme
  for diffusion problems. This is a first step, and does not handle
  parallelism yet. Integration with the rest of the code will be improved
  progressively.

- Change multigrid default settings based on recent performance benchmarks.

- Add 3-layer conjugate directions (residual) algorithm for non-symmetric
  linear systems (this algorithm was used in older EDF codes, and is
  different from the "standard" conjugate residual, which is limited to
  Hermitian matrixes).

- Add support for using PETSc in linear solvers.

- Convergence of linear solvers may now be plotted using CSV files.
  When solver verbosity is high, residues are not printed anymore
  (except in case of non-convergence), as plotting provides the same
  info in a more usable manner.

Architectural changes:

- Add --sig-defaults option to cs_solver executable to allow use of
  default signal handlers.

- Only trap floating-point exceptions by default for code compiled in
  debug mode, to avoid interfering with optimizer speculation.
  Exception trapping can be enabled or disabled locally using the
  cs_fp_exception_...() functions.

- Remove support for some older MPI libraries or process managers
  (which can still be used if needed using post-install settings).

- Complete Python 3 compatibility, except for builds with SALOME support
  (as the SALOME platform is not yet Python-3 compatible).

- Refactor EnSight, MED and CGNS output, replacing slice by block logic.
  The fvm_gather_...(*) API is now removed, and using parallel APIs
  for MED and CGNS should now be relatively straightforward.


Release 4.0.0 (2015-04-30)
--------------------------

Changes:

- Added boundary condition to flow mapping (feedback) functionnality
  (see boundary_conditions_map and boundary_conditions_mapped_set functions).

- By default, do not force use of an iterative gradient reconstruction method
  for pressure gradients, or other gradients deriving from a potential.
  To force it, a negative value of the IMRGRA keyword may be used.

- Coprocessing: if the main (default) writer uses Catalyst output, ignore
  joining output and use EnSight Gold format for joining and error writers
  if no matching Python script is found, to avoid uneeded abort.

- Turbomachinery modeling:
  * enable multiple rotors for rotor-stator model based on mesh joining.

- Generalize double backward implicit Euler time scheme for all variables.
  It can be activated with the keyword ibdtso(ivar) > 1.

- Add a "slope_test_upwind_id" field keyword, allowing postprocessing output
  of the contribution of slope tests to convected variables.

- Add a new dilatable (non conservative) algorithm for fire modelling.
  Activate it with idilat=4 (the formulation is in div(u) instead of
  div(rho u)). You can access to the previous dedicated algorithm by
  setting idilat to 5.

- Add L2 time normalized residual to the log.

- Major changes for solvers:
  * Unify handling of linear solvers, so as to allow finer user control, and
    enable future additions of solver options and user-defined or external
    solvers.
  * Add Block Gauss-Seidel linear solver, which may be accelerated for "upwind"
    type systems by a matrix line ordering. This is used (by default) for the
    DOM radiation module, using the ordering defined by the radiation direction,
    and leads to a factor of 2 to 4 improvement over the Jacobi solver for this
    case (tested on a small number of MPI ranks).
  * Single-reduction conjugate gradient is now an option rather than a separate
    solver. This allows switching automatically from one to the other based on
    computation vs. communication cost.
  * Add BiCGstab2 linear solver.
  * Change default solver for pure diffusion problems (from Conjugate gradient to
    multigrid solver).

- Add a new experimental module for Darcy flows. It can be
  activated with the keyword usppmo(idarcy) = 1.
  This path includes new developments to improve reconstruction gradient
  calculation with heterogeneous diffusion coefficients. This feature is only
  available for standard LSQ gradients and can be activated with
  the keyword iwgrec(ivar) = 1.

- Add drift modelling for coal combustion, and clean up the module:
  * Now, the enthalpy of the continuous phase (gas phase, h1) is transported
    rather than deduced (to be precises, X1.h1 is transported)
  * The convective flux for the gas phase is deduced from the convective fluxes
    of the particle classes and the convective flux of the bulk: therefore, the
    algorithm is fully conservative over time and space.
    (Note that the conservatism is obtained mainly due to having upwind schemes;
     transport equation on X1 with convective field V1 and
     transport equation on h1 with convective field X1.V1
     is equivalent to transport equation of X1.h1 with convective field V1...)
  * Some coal combustion fields are renamed
  * A model with a transported particle velocity per class is added (in fact,
    this velocity is handled as 3 scalars).
  * Correct cs_user_initialization-unified_combustion_coal.f90 user example.

- Add a condensation wall model.

- Added cell and face renumbering options to try to improve performance.

Architectural changes:

- New PLE (parallel location and exchange) library API, allowing new locator
  features and algorithm versioning. This requires an update on the Syrthes
  side for conjugate heat transfer.

- The preprocessor now uses integers of the same size as the kernel's
  global numbers (cs_gnum_t). This allows handling of meshes above
  200 million cells, at the cost of about 50% higher memory usage for
  this stage (unlesss the code is configured  with --disable-long-gnum).

- Remove rtp and rtpa arrays. Variables are stored with field
  structure owning its value arrays.

- Rename and merge user subroutines:
  * User scalars are not declared in usinsc (in cs_user_parameters.f90) anymore,
    but in cs_user_model (in cs_user_parameters.c). They are not named based
    on their number, but both their name and label is determined by the user.
    This should make it easier to combine data setups.
  * Remove the cs_user_field_parameters subroutine from cs_user_parameters.f90.
    user code in that subroutine may now be placed in usipsu or usipes.
    The usipsc subroutine is also removed, and scalar variable diffusivity
    behavior may be activated through usipsu.
  * merge (usebu1, uslwc1, usd3pt1, uscpl1, user_coal_ini1, user_fuel_ini1)
          in cs_user_combustion (in cs_user_paramters.f90)
  * split and rename usray1 in cs_user_radiative_transfer
    (in cs_user_parameters.f90)
          (Note that the declaration of the use of radiative transfer is in
                usppmo, as the other specific physics models)
  * rename usray2 in cs_user_radiative_transfer_bcs
  * merge usalin into usipph (in cs_user_parameters.f90)
  * rename ustsma into cs_user_mass_source_terms
  * rename uskpdc into cs_user_head_losses
  * remove the usipgl subroutine from cs_user_parameters.f90.
    Specific, model oriented options are defined in usipph, general options
    in usipsu.

- Lagrangian particle data is now shared between Fortran and C parts.
  Fortran arrays have been replaced by pointers, which map to the C data.
  arrays itepa, tepa, ettp, and ettpa are replaced respectively by pointers
  ipepa, pepa, eptp, eptpa, which use interleaved data; for example, the
  equivalent of itepa(np,jisor) would now be ipepa(jisor,np).
  Allocation is automatic, and npmabs is no more.

- When building with both lib CCM-IO and CGNS, do not try to use CGNS's ADF
  library for CCM-IO. This used to work in the past, but now leads to failure
  opening CCM files. The change might lead to issues for such builds when using
  ADF output for CGNS files, but HDF5 is now the default for CGNS (where
  available), and the hdf2adf and adf2hdf tools of the CGNS library suite allow
  conversion, so this should not be an issue. Builds with both libraries with
  linkers not allowing or configured to allow multiple definitions might be
  an issue.

- Restart files now use a new section naming scheme, at least for field data.
  this allows more automated handling of variables and properties in
  checkpoint/restart.

- Major indexing change: the mesh connectivity arrays (in cs_mesh_t
  structures, i_face_cell, b_face_cell, i_face_vtx_idx, i_face_vtx_lst,
  b_face_vtx_idx, and b_face_vtx_lst) now use 0-based indexing.
  This simplifies expressions in C code, but may break existing C code
  referencing those arrays. On the Fortran side, ifacel, ifabor, ipnfac,
  nodfac, ipnfbr, and nodfbr are no longer arrays, but elemental functions,
  which add 1 to the C array values, so as to make this transparent to the
  user. As most accesses are now in C, and the ifacel and ifabor functions are
  defined as "elemental pure" in the mesh module, the performance impact
  should not be significant.

- Some arrays for which there is a field keyword equivalent, such as ivarpr,
  are replaced by functions, which allow similar read-only access, but not
  modification, which must be done through the field API.

Bug fixes:

- Face mean secondary viscosity is now computed consistantly with face mean
  viscosity (namely harmonic if imvisf equal to 1).

- Fix bug introduced in rev. 4491 (before release of 3.0.0):

  - in rev 4491, the call to specific physics laws was moved from BEFORE
    the user call to AFTER the user call, which was only needed for the
    compressible module, whereas it impinges ALL the specifics physics.
  - therefore ppphyv is split into two.
  - WARNING: this change MAY impige all the specific physics
    (Atmospheric, Electric Arcs, Combustion, etc.)

- Many more bug changes (see ChangeLog for details)


Release 3.3.0 (2014-05-16)
--------------------------

Changes:

- Add 2 order backward Euler scheme in time for velocity prediction.
  Set ibdtso = 2 in usipsu to activate it.

- Add cavitation models. See the documentation (theory, user, Doxygen) for more details.
  You can activate it in cs_user_parameters.f90 with icavit=1.

- Add the Shir model for the turbulent diffusion of the second moment turbulence models
  (Rij and epsilon, isotropic diffusion). Available with the idirsm=0.

- Enable use of sparse matrix-vector tuning mechanism in regular calculations.
  Automatic tuning or variant selection may be done using the
  cs_user_matrix_tuning() function, in cs_user_performance_tuning.c.

- For non-batch systems, handling of the number of MPI ranks is
  based on a "code_saturne run" option, --nprocs, and is set in the
  runcase file, not in the XML file anymore.

- For couplings, runcase_coupling is replaced by a standard runcase,
  mechanism, using "code_saturne run --coupling coupling_parameters.py",
  where coupling_parameters.py contains the domains dictionnary previously
  in runcase_coupling. This allows unifying the handling of "code_saturne run"
  options.

- Rewrite of temporal moments handling. Moments handling is now more modular,
  and allows for variances in addition to means. Also, numerically stable recurrence
  relations are used to update moments, whose values are now directly usable at any
  given time. Weight accumulators are handled inside the module, and not seen as
  fields anymore. Also, support for user functions is added.
  Currently, this is mapped to the legacy data setup, and tested only in this
  context, but the added functionality will be exposed with future changes in
  case setup.

- Lagrangian module: Improvements in roughness and resuspension models
   * added a user keyword 'irough' for roughness surface (calculation of the
     energy barrier in the case of rough wall)
   * consideration of the electrostatic force in the adhesion force
     for the resuspension
   * mass flux update for particles rolling on the wall

- Ensure postprocessing of symmetric tensors for EnSight and Catalyst is consistent
  with {xx, yy, zz, xy, yz, xz} internal component representation.

- Update some external package versions for automatic installer.

Architectural changes:

- Remove rtp, rtpa and propce in user subroutines API. The access to the
  variables and the properties is carry out using field structure.

- Velocity field is now created as "owner field" in C (removed from rtp
  and rtpa). Velocity is now fully interleaved.

- Ensure only platform-independent files are installed to datarootdir
  (either as defined by ${prefix}/share or as ${datarootdir}. Files
  which contain install-related paths are now installed in sysconfdir.

- Fields representing variables are now declared in fldvar, previously call 1 of varpos.
  Fields representing properties are now declared in fldprp, previously call 2 of varpos.
  Fields for the radiation model are now declared in rayprp, previously call 4 of varpos.
  The only remaining call of varpos matches previous call 3.
  The nomvar array is replaced by field labels for fields, and nomva0 for temporaries.
  Field pointers for all variables are accessible from C using the macros in
  cs_field_pointer.h. The ichrvr, ilisvr, scamin and scamax arrays are also replaced by
  field keywords.

- Lagrangian module: for tracking, particle structure is replaced by raw data
  mapped at run-time, allowing for allocation only of useful data.

- Remove support for obsolete LAM and MPICH1 MPI libraries.

- Finalized cleanup of handling of boundary condition coefficients. Fields now allocate
  their coefficients, rather than map them from a common array.

Bug fixes:

 - Fix bug in the Generalized Gradient Diffusion Hypothesis (GGDH).
      - swich between two components (R13 and R23): this is not impacting the validation
        database because this model is only validated on 2D test cases
        where the R13 and R23 components are both 0.
      - from version 3.1, it is impacting all the Rij-epsilon calculations
        because from this version, the Daly Harlow model is passed by default on Rij
        (this model induces a GGDH on Rij and epsilon).


Release 3.2.0 (2013-12-04)
--------------------------

Changes:

- Documentation:
   * moved tutorials outside the codebase, and into a separate base.
     This allows looser synchronization with the code base, as tutorials may be
     updated somewhat less frequently.
   * complete the Doxygen documentation of
     - Fortran modules
     - user examples
     - Fortran routines
   * install Doxygen documentation from tarball (as built by "make dist").

- Remove deprecated uncoupled velocity solver (ivelco=0).

- Thermal model:
   * the thermal model is now defined by the "itherm" keyword/variable, which
     replaces iscsth(iscalt). In the case of temperature, the scale used is
     defined by a separate variable (itpscal). For additional user scalars,
     a new array iscacp is defined, such that iscacp(iscal) defines whether the
     scalar behaves like a temperature, so the possibility of modeling multiple
     passive "temperatures" is not lost.
     This change allows for better consistency between the standard and specific
     physics, as the thermal variable is now always a "model" scalar, and user
     scalars remain separate. It also allows querying the thermal model with one
     less indirection.
- Atmospheric module:
   * add gaseous chemistry models.
   * plug the SIze REsolved Aerosol Model (SIREAM).

- Turbomachinery modeling:
   * added a rotor-stator model based on mesh joining.

- Particle tracking module:
   * add a modeling of the drying phase of the coal particle combustion
   * add a new boundary condition to simulate coal fouling mechanism
   * bug fixes related to radiative transfer
   * implementation of a particle discretization in the coal combustion model:
     - backwards compatibility is ensured (set nlayer = 1)
     - computation of intra-particle thermal gradients
     - adaptation of chemical source terms to temperature discretization
     - reworked of the particle injection for coal (clear difference between standard
       and user-defined coal composition)
     - adapted the particles and trajectories export routines to be able to output variable
       information for a specific layer.

- Compressible module:
   * change the compressible algorithm from a density formulation to a pressure formulation
   * merge the compressible algorithm with the coupled velocity components algorithm
   * adapt standard bricks (codits, bilsc*) in order to make them compatible with the
     compressible algorithm
   * implement analytical flux boundary condition
     (plus a new total enthalpy / total pressure boundary condition with a fixed point algo,
      generalization of the subsonic outlet)
   * creation of a new set of BC coefficient for the convection operator for compressible flows
   * remove density in the variables array (rtp) and keep it only in the properties array.

- Coal combustion module:
  - added new NOx model for coal combustion;
  - introduction of the coal thermal conductivity
    (for the calculation of intra particle gradients in particle-tracking module)

- Boundary conditions:
   * Fix in the wall boundary conditions for the viscous boundary term
     (the viscous boundary term is not always parallel to the wall). This is
     mainly impactant for verification testcases.
   * add a new Boundary Condition type for free inlet,
     - this BC can be used for natural convective flows in free atmosphere for
       instance (plumes, flame, etc.).

- Turbulence:
   * Major change in Rij-epsilon models:
     - the Daly Harlow model on the diffusive term is now by default for SSG
     - the GGDH brick is used for all the models (LRR, SSG, EBRSM)
     - the "diffusivity_tensor" is added as a field key word
     - Rij-epsilon routines are cleaned up and doxygened.

Architectural changes:

- For CFDSTUDY, use PARAVIS instead of VISU.

- Only search for SCOTCH/PT-SCOTCH, Metis/ParMetis, Catalyst, and EOS
  if explicitely required to avoid detection of incorrect versions
  on clusters.

- Added cs_c_bindings.f90 module for general definitions of C bindings.
  For large modules, it is recommended to use separate files (see field.f90
  and post.f90 for example), but for smaller modules, this avoids requiring
  the definition of specific module files.

- Added cs_field_pointer API for quick access to main fields from C.

- Added experimental ParaView Catalyst co-processing output option.

- Added new all-to-all infrastructure, allowing alternative algorithms such as
  Crystal Router to MPI_Altoall/MPI_alltoallv as an option.

- Move the convection-diffusion balance (bilsc2.f90) to C.

Bug fixes:

- Fixes in the wall boundary conditions (validation cases re-run for safety,
  and were not impacted by these changes)
  - for the velocity, when the wall velocity condition is non-zero,
    the gradient BC was wrong.
  - for the scalars, a term (proportional to the wall value) was missing
    in the gradient BC
  - for scalars with small molecular Prandtl number, the value in the viscous
    sub-layer was not consistent between the gradient BC and the flux BC.

- Many bug fixes (see ChangeLog for details).


Release 3.1.0 (2013-06-11)
--------------------------

Changes:

- Many documentation updates

- Documentation: correction of the wall functions provided by D. Laurence,
  description of the BL-v2-k model provided by F. Billard.

- Improve k-omega robustness when yplus is low.

- Use field structure for the scalar with drift and add a drift
  model for coal combustion.

- New key words are added to the field structure:
  * scalar_id (inverse of the array isca(iscal))
  * scalar_class (class of the scalar which belongs to, in particular sharing
    the convective mass flux)
  * inner_mass_flux_id (index of the convective field at inner faces)
  * boundary_mass_flux_id (index of the convective field at boundary faces)
  * drift_scalar_model (type of model for drift scalar:
    0 no drift, > 0 with drift)
  The different models may be activated in cs_user_parameters
  (user_field_parameters)), and one example is available in src/user_examples
  and documented with Doxygen.

- Lagrangian module: implement a zero-flux particle boundary condition to be
  applied with eulerian symmetries.

- Lagrangian module with combustion: use of a formulation of the coal density
  local to a particle and improve the numerics.

- Lagrangian module: add a particle resuspension model.

- Lagrangian module: implementation of a wall law for fluid velocity, k and
  epsilon for the deposition sub-model

- Radiative model: fix existing T2 and T4 quadratures, and add new
  S4 S6 S8 and Tn quadratures .

- Add numbering options for threads.

- Use relative paths when possible in execution script.

- Lagrangian deposition sub-model: improvement of the initialization of the
  particle state vector.

- Increase max coal number to 5.

- Add coke composition for coal combustion.

- Remove old coal combustion model

- Implementation of a Lagrangian boundary condition based on the DLVO theory

GUI changes:

- If the GUI is launched through SALOME, update the object browser in order to
  display the results.

- Possible activation of the Lagrangian particle deposition sub-model
  through the GUI.

- Add new case creation capability by GUI.

- CFD_STUDY: display the listing in the standard dialog window.

User and pre/post processing changes:

- Lagrangian model: full rewrite of the postprocessing output,
  which is now usable in parallel.

- Add handling of white-space in paths.

- GUI and CFD_STUDY: display the monitoring points on the SALOME VTK viewer.

AutoVnV changes:

- Allow to call Code_Saturne GUI functions in external preprocessing scripts
 (automatically set the PYTHONPATH variable).

- Add file input feature to the postprocessing environment in autovnv.

- Add use Agg and suppress usetex to run autovnv on cluster with no X server.

- Merge the preprocessing and case running steps (update documentation example).

- Add the possibility to apply a global post-processing to a study composed
  of several cases.

- Add the possibility to run the same case several times and to prescribe
  impose the name of results directory.

Architectural changes:

- Automatic installer changes. The installer is now in the top-level directory,
  and does not download Code_Saturne anymore. The setup file template is
  generated by a first call to install_saturne.py. MPI should now be installed
  upstream, but PT-SCOTCH and ParMetis are now handled.

- Improve robustness of SALOME detection for YAMM builds. Sourcing the SALOME
  environment should not be required anymore when the --with-salome configure
  option is used with a YAMM salome build (which incluides its own Python
  interpreter by default).

- To simplify current and future tests, a build-aux/cs_config_test.py
  file was added to allow writing of configuration tests in Python.

- Improve MPI variant detection.

- Add EOS detection and support.

- Add MPMD support for Blue Gene/Q.

- Remove support for legacy (version 1-1 to 1.3) restart files, simplifying code.

- Remove support for obsolete METIS 4/ParMETIS 3 versions.

- Add "flush" stage times to postprocessing writer.

- Add CS_LOG_DEFAULT type to use logging API to listing file.

Bug fixes:

- Many bug fixes (see ChangeLog for details).


Release 3.0.0 (2013-03-22)
--------------------------

Changes:

- Many documentation updates.

- Lagrangian deposition submodel: improved parallel handling of the
  deposition submodel.

- Always perform the first iteration of the iterative processes (to prevent
  problems if the user set a too large value for the solver precision)

- LES: restore the wall shear stress in the velocity components coupled
  solver (ivelco = 1) when the first cell is in the viscous sub-layer.

- Compressible module is now usable again.

- Lagrangian module: forbid the use of the broken trajectory and displacement
  post-processing in parallel mode.

- Allow user selection of single reduction congugate gradient algorithm.

- Add a Boundary type to vectors (generalized symmetry) where a Dirichlet
  is imposed on the normal component and a flux is imposed on the tangential
  ones (icodcl=14). Needed for Marangoni effect.

- MAJOR fix in Head Losses:
  - itrgpv.f90, itrmav.f90 and projtv.f90 are now adapted to anisotropic
    diffusion (same algorithm as diften.f90 for scalars) for ivelco=1
  - fix in BCs on the rpessure field for ivelco=0 with head losses.
  - the symmetric tensor K for head losses is stored as follows
    [k11, k22, k33, k12, k23, k13] instead of [k11, k22, k33, k12, k13, k23]
    to be consistent with the anisotropic diffusion operator.
  - pseudo pressure velocity coupling is also corrected.

- Replace symmetry/block flags for matrix benchmarking by fill type.
  This allows finer control and easier addition of fill types,
  such as diagonal/extradiagonal block/no block combinations,
  and possible future special block types.

- Add finer control for computation and visualization  of bad cells.

- Add parallel checkpoint/restart for particle tracking.
  (this breaks the contents compatibility for the restart content).

- Unset default "relaxation" on k-omega.

- MAJOR fix in k-omega turbulence model. Buoyancy term and -2/3divu term are
  implicit in the equation of k if they are negative.
  It prevents from clippings on k (if an upwind scheme is used and if ircflu=0).

- MAJOR change in "k-epsilon" turbulence models (iturb = 20, 21, 50, 51):
  - the time stepping on k-epsilon Linear Production (iturb=21) is revisited
    so that any negative source term is taken implicit (it could be the case of
    "-2/3*k* trace(Grad u)" or of the buoyant term). No Clipping on k or epsilon
    should now appear if the convective scheme is upwind and if there is no
    reconstruction on the diffusive term (ircflu=0).
  - The limitation done by the k-epsilon Linear Production (iturb=21) model when
    the Strain rate is big if performed on the deviatoric strain rate and not
    the strain rate itself, as proposed by the model. That was an issue because
    trace(Grad u) is not exactly 0 even for incompressible flows.
  - The subroutine turbke is partially rewritten using additional arrays making
    the code clearer. The order of the steps is then more logical.
  - The relaxation which was performed by default on k and epsilon in tridim is
    deactivated for iturb=21.
  - The negative source terms ("-2/3*k* trace(Grad u)" or the buoyant term) are
    also implicit for models iturb=50 and iturb=51.
  - The buoyant term for atmospheric flows has been adapted.

- Lagrangian module: seed the random number generator with (num_rank+1)
  to avoid statistical bias.

- Lagrangian module: reimplementation of the counter of depositing particles.

- New Fortran bindings for postprocessing, and update user output of variables
  (this breaks the previous usvpst routine, but should be clearer).

- Move definition of variables output of all physics to usipes and separate
  detailed examples.

- Add flags for condensation source terms and exchange coefficient correlation
  (icond and iwallt flags).

- Always use an iterative gradient reconstruction method for pressure gradients.
  In fact, this is done for all gradients deriving from a potential. To
  allow the previous behavior, a negative value of the IMRGRA keyword may be
  used (forcing the method matching its aboslute value).
  Iterative gradients may now also be initialized by least squares using an
  extended neighborhood.

- Set iclsyr=1 as default option (improved symmetry boundary conditions on
  Reynolds stress tensor).

- Adapt diffusive modelization of turbulent flux (T'u') to any scalar.
  This is performed using the new field structure.
  For each scalar, the user can choose 3 options:
  - iturt(iscal) = 0 SGDH (default) muT * GradT
  - iturt(iscal) = 10 GGDH, T'u' = ctheta*k/eps/Rij * GradT
  - iturt(iscal) = 20 AFM, (Alegbraic model)
  - iturt(iscal) = 30 DFM, (Transport equation on u'T').
  Note that the user BCs on u'T' are set at the end of icocl and rcodcl
  (if iturt=30).

- Add a hydrostatic pressure gradient computation (iphydr = 2).
  This allows handling the imbalance between the pressure gradient and gravity
  source term. This pressure gradient currently requires an orthogonal mesh.

- Remove limit on number of classes for Lagrangian inlets.

- Syrthes now only allows nearest-neighbor search for non-matching faces
  as an advanced option. By default, only matching faces are detected, and
  non-matching faces provoke an error. This avoids silently switching to a much
  slower algorithm, while non-matching faces are most often due to inconstitent
  coordinates between domains, rather than different mesh feature detail levels.

- Add Fortran derived types to manage arrays of pointers.

- Add several boundary postprocessing options (notably wall T+, temperature, and Nusselt).

- Add GGDH-AFM and DFM models for the thermal scalar. This options is available
  with the key word iturbt (0, 10, 20, 30).
  Warning: the transport equation on turbulent fluxes (DFM) requires changing
           the gradient boundary conditions on the thermal scalar to have the
           proper production term. It is now consistent with the gardient
           boundary conditions on the velocity.

- Improve user scripts, to add data preparation and results copying user hooks.
  An example is provided, allowing simplified successive restart handling.

- All files in case's DATA and not in subdirectories are now copied
  automatically to execution directory.

- Move initialization of IO logging to main initialization block.

- Add a soot model for gas combustion.

- Add computation of combustion source terms for the Low Mach algorithm (idilat=4).

- Change the default value for the precision of the reconstruction iterative
  process (codits, coditv): it is now by default 10 times the precision of the
  linear solver, to be consitent with resopv.f90.

- Add an extra block size to matrix (used for solving vectorial field with
  tensorial diffusion). Add bricks for tensorial diffusion for vectorial fields
  (matrvv, vistnv, diftnv).

- Modify iterative process to print the proper residual in the listing.
  That means we always compute the residual (whereas we didn't before), so
  it might be a little bit more costly.

- Adapt the dynamic relaxation to vectorial field (from codits to coditv).
  The log messages are also changed so that the normed residual in the listing
  is now the one of the iterative process (codits and coditv).
  The nuber of iterations for the solver is now the sum of iterations for the
  solver over the sweeps in codits or coditv, as it was already done in resopv.

- Add a new algorithm to solve diffusion term with anisotropic & heterogenous
  viscosity for scalars (for use with GGDH);
  to activate it set idften(ivar) = 6 (for symmetric tensor diffusion)
                                  = 1 (for standard scalar diffusion)
                                  = 3 (be for orthotropic diffusion, diagonal tensor)

- Add a new dynamic relaxation for solving the Poisson equation on the pressure;
  to activate it set iswdyn(ivar) = 1 for a relaxation with the last increment
                                  = 2 for a relaxation with the two last increments.
  should be used only for transport equation without advection.

- Replace depecrated activated() SIGNALs by triggered() ones as advised by
  the Qt documentation.

- Add rotation/curvature correction (Spalart-Shur and Cazalbou) for Eddy viscosity
  turbulence models.

- Add (exact) Coriolis terms in Rij-epsilon models (LRR, SSG and EBRSM).

- Use vector coefa/coefb for postprocessing output when applicable.

- Major Fix on Boundary conditions with the coupled velocity components algo (ivelco=1).
  - The Gradient boundary term for the velocity in k-epsilon, 2 scales of velocity was
    not computed correctly.
  - The wall shear stress in Rij-epsilon was counted twice.
    Note: The previous versions of CS for Rij didn't ensure that rho uk*uet was
          the wall shear stress.
    Note 2: The 2 scales of velocity for Rij does not work properly in a channel,
            because the scale uk is underestimated.

GUI changes:

- Integration of UNDO/REDO function.

- Allow setting of RHS reconstruction sweeps for pressure.

- Add electric models.

- Change partitioning page into tab of a more complete performance tuning page.

- Add a popup for print XML function

- suppress Current species class.

User and pre/post processing changes:

- Change default v2f model from phi model to BL-v2k model.

- Output probe set coordinates in CSV format.

- Regroup usatsoil and usatdv user subroutines in cs_user_atmospheric_model.f90.

- Replace ficstp by control_file, with extended options. This file may now
  also be used to force postprocessing or checkpoint output.

- Introduce physical time for mesh rotation start.

- Allow using SYRTHES 4 even if matching environment is not sourced.

- Allow user definition of additional compile/link options to simplify usage of
  external libraries.

- Move some post-processing functionnality to utility functions, for easier usage.

- Add user example for visualization of space-filling curves.

- Add a starting time for the computation of moments.

- Read JANAF file from reference directory directly instead of using a local copy.

- Use reconstructed value for postprocessing of heat flux.

- Remove coefa/coefb arguments from user subroutines.
  We now encourage the use of the field API. Legacy coefa/coefb
  arguments are nonetheless kept in cs_user_extra_operations.f90,
  as we may still need to compute a gradient of a single velocity
  component even when ivelco = 1 (at least until examples are improved).

- Add CCMIO output for main volume and boundary meshes.

- Move user control of mesh warping handling to cs_user_mesh.c

- Improve defaults for cs_user_parameters and make examples more consistent.
  The "--nogui" option of "code_saturne create" is finally removed.
  Tests on "iutile" are replaced by .true. or .false. for clarity.

- Remove obsolete igghexa to MED 2.3 converter tool.

- Lagrangian module: removal of the now-useless routine uslabo.

Architectural changes:

- Move field properties from cstphy to optcal.

- Update detection of PT-SCOTCH/SCOTCH to handle new version 6.0.

- Allow definition of a user or site shell rcfile to modify the execution environment.

- Add finer control for parallel block I/O.
  Global data is now never read in a collective manner, but broadcast,
  while additional modes and hints for MPI-IO may be chosen by the user
  (control is moved from the command line to the GUI and user files).
  It is now possible to use a specific I/O communicator of smaller
  size than the main calculation communicator.
  These additions are part of a "Performance Tuning" group of user options.

- First version of the Windows port of code_saturne.
  This brings up the following adaptations: the installation can now be relocated
  (also available on Unix systems), the scripts (SaturneGUI, runcase,
  run_solver.sh) are now generated on the fly and depend on the underlying
  system and shell, and the interaction with NEPTUNE_CFD package is also changed
  for a correct behavior.

- Add C API for access to time step information.

- Use ISO_C_BINDINGS to improve Fortran API for fields

- Python cleanup and partial Python 3 compatibility.
  Except for Unicode and Qt4 aspects of the GUI, Python code
  should now be compatible both with Python 2 and Python 3.

- Move some base Fortran API functions from cs_base.c to cs_base_fortran.c,
  and replace CPU times by elapsed time in Fortran timings.

- Default logging mechanism may now switch between multiple modes.
  This allows combining a C stdio-based API an a Fortran IO based API,
  switching to Fortran where it is dominant (in which case bft_printf()
  uses Fortran IO, as before), and back to C where preferred.
  This removes the need for flushing of Fortran output, as bft_printf_flush()
  is used mainly by the C API, and a switch back to the C API is forced
  (with the Fortran output being closed) before error logging.

- Allow query of locale directory and package data directory from main executable.

- Rewrite gradient quality tests in C.

Bug fixes:

- Many bug fixes (see ChangeLog for details).


Release 2.3.0 (2012-07-23)
--------------------------

Changes:

- Set the coupled velocity component solver as default option (ivelco=1).

- Default turbulence model is now linear production variant of k-epsilon.

- MAJOR atmospheric/meteo addition of new features (humid atmosphere and
  soil module).
  Note: humididy, soil and 1D radiation models are experimental.

- MAJOR change: the Temperature transport equation is multiplied by
  Cp (specific Heat) so that the equation has now the dimension of energy.
  The change is expected to have NO influence when Cp is constant.
  When Cp is variable in space, this change fixes the error done on this
  equation.

- Change default values for NSWRSM.

- Add synthetic turbulence inflow methods for LES:
  - random method (Gaussian noise),
  - Batten method (based on Fourier decomposition of turbulent fluctuations),
  - Synthetic Eddy Method (SEM).

- Adapt the Multigrid algorithm to vectorial Poisson equation
  (such as mesh velocity for ALE).
  - The aggregation criterion is based on the trace of the diagonal block
    DA (3x3), but can be a changed (We could test n.DA.n),
  - The multigrid algorithm for scalars (such as the pressure field) is
    rigourously unchanged.

- Add a low Mach algorithm (semi-analitical: idilat=4).

- The subroutines codits and inimas/inimav are transparently changed:
  - codits save and return the last increment,
  - inimas/inimav can compute a velocity flux OR a mass flux regarding
    the value of itypfl.

- First stage of the implementation of the diffusion-inertia model of
  aerosol deposition. Still under development in this revision.

- Add an multi-species algorithm for low-Mach number algorithm (idilat=3).

- MAJOR changes to the formulation of boundary conditions for diffusive part.
  The changes impact the subroutines where the BCs are computed (condli,
  clptur, clptrg, clsyvt) which are rewritten (and doxygened).
  The changes should not impact the results (to the truncature error precision).
  A new boundary condition is added: convective/radiative outlet (icodcl=2).
  The radiative transfer module has been widely modified:
  - a new solved variable (in rtp) call ilum has been created,
  - the boundary coefficient are now in the same array as for the
    other variables,
  - the new Boundary condition (radiative) is used and user set BCs with
    rcodcl and icodcl
  This may have introduced bugs in the compressible module.

- Lagrangian module: removal of the 'snap_to_grid' (legacy) method of
  particle localization

- Lagrangian module: first stage of the implementation of the parallelism.

- Add bad cells detection and post-processing.

- Lagrangian module:
  - English translation of the main messages
  - Suppression of the display of the ambiguous mean values of the stats

- Clean the coupled velocity component version or the correction step
  of pressure. The loop over non orthogonalities is performed in a
  clearer manner. The updating of the mass flux is always performed so
  that the continuity equation is fullfilled exactly (at the pressure
  precision) even when the iterative process have converged.

- Pass the mass aggregation term (-div(rho u)T) directly into the linear
  system (matrix and bilsc) to be coherent with NCFD.

- Add a Low Mach compressible algorithm conservative in time for the
  momentum equation and the transport equation of any scalar.
  It adds a prediction step of the mass flux. Available with the
  key word idilat (2 or 3).

- Add a low-Mach algorithm to account for the mass equation for dilatable
  flows. It is only available for mono-species flows at the moment.
  This can be activated with the idilat keyword set to 3 (1 being the
  current default).

GUI changes:

- Add handling of gas combustion model.

- Add handling of solid fuel combustion.

- Add compressible algorithm in GUI (not fully operational in this version
  yet: version 3.0 recommended).

- New presentation and use MEI to initialize scalar boundary conditions.

- Volumic initialization for meteo variables.

- Add choice for NSWRSM and IRESOL.

- Use splines to define profiles.

- Modification of min, max and initialization for scalars and variance.

- Langrangian module: modification of the oundary conditions:
  - Rename of the classical boundary conditions

- Lagrangian module: modification of the volume and boundary
  statistics management:
  - Rename of the default names
  - Names non-modifiable in the GUI
  - Post-processing or not or the default variables
  - Move of the names from uslag1 to lagopt.f90
  - Rename in lagopt to be consistent with the GUI

- Pulverized-coal model not activatable in the GUI (deprecated).

- Use MEI to initialize turbulence, velocity and thermal variables by zone.

- Add use of MEI for turbulence boundary conditions.

- Added piso and modification of control time step with
  velocity-pressure algorithm choice.

- Move hydrostatic pressure option from body forces to numerical parameters.

- Use MEI only for deformable mesh and control access view.

- Add support for Rij-epsilon EBRSM and Spalart-Allmaras turbulence models.

User and pre/post processing changes:

- Separate usproj.f90 into multiple examples.

- Add a user example for LES inflow.

- SYRTHES coupling: in case of unlocated Code_Saturne elements in SYRTHES mesh,
  these elements are now post-processed.

- In case of a single SYRTHES coupling, adjust name to automatic match.

- Split usray5 ("User" subroutine for radiative transfert) into 2 parts:
  - raycll, which is not a User subroutines, which sets BCs on luminance
  - usray5 (iappel=2) where a net flux is computed.
  This should allow to remove luminance from rtp (ilum).

- Add variable/port name to SALOME Kernel Calcium API messages.

- Postprocess deformation in ALE mode.

- Add user function to disable or force mesh_output.

- Replace usphyv by cs_user_physical_properties.

- Force use atphyv for density property for meteo physics.

- Update of the max number of particles to visualize from 500 to 100000

- EnSight Gold variable description limit is now 49 characters, not 19.

- Autovnv: add new option --update in order to upgrade a repository
  of test cases by reload files of parameters (i.e. run the
  backwardCompatibility method) and changes pathes for SaturneGUI and
  runcase. Add capability to mix Code_Saturne and NEPTUNE_CFD test
  cases in the same Study.

- New Autovnv. functionality: 2D view of scalar (needs pyvtk).

- Remove Fortran utility subroutines used to obtain global mesh element numbers.

Architectural changes:

- Use 16 characters instead of 8 for variable and property names.
  Improve formatting, and make default variable names witout GUI more
  consistent.

- Add C parallel operation wrappers.

- Add new structures for measures set -> global mesh interpolation
  (cs_measures_set_t) and global mesh -> point cloud interpolation
  (cs_interpol_grid_t).
  These features will be used soon in atmospheric module.

- Dependencies on libraries are now handled in the standard
  Automake manner, to avoid issues with the gold linker
  (this means building a dynamic libsaturne on top of static dependency
  may no longer be possible).

- Detection of dynamic versions of SCOTCH, METIS, and HDF5 is
  now made possible.

- On Mac OS X, build is static by default, and compiling of user
  subroutines for a static build requires unarchiving the library and
  overwriting selected object files, to avoid issues with multiple
  definitions not being handled by the Mac OS X linker.

- Move partitioning to main solver executable.
  This simplifies the toolchain, as a separate partitioner is no
  longer required. Additional options for finer-grained control
  are provided, and parallel partitioning is encouraged.

- libPLE: use same FLAGS as parent Code_Saturne build.
  This helps ensure subconfigure is consistent. Bootstrapping also
  inherits prior cleaning stage from that of parent.
  Also, dependencies on MPI are handled in the standard Automake
  manner when building libPLE (this means building a dynamic
  libPLE on top of a static MPI may no longer be possible).
  This implies changing at least the patch release number for PLE.

- Port to Blue Gene/Q.

- Moved block distribution functions from fvm to base.

- Allow user configuration of compilation flags for performance-critical files.

- Merge fvm_parall.* and cs_parall.*.

- Remove support for coupling with (obsolete) version 3.4 of SYRTHES.

- Reduce use of external BLAS functions based on recent comparisons
  with internal functions.
  Reduce usage of external BLAS to benchmarking (to allow for occasional
  comparisons).

- Replace external BLAS dot product with superblock variant for
  better expected precision.

- Move definition of parallel rank and thread status from cs_base.*
  to cs_defs.*.

Bug fixes:

- Many bug fixes (see ChangeLog for details).


Release 2.2.0 (2012-03-30)
--------------------------

Changes:

- Add a porosity formulation.
  The transport equations (espacially in turbulence) have to be checked.
  We also have to check if the formulation is all right in presence of
  Coriolis forces.

- Add Least square method for gradients of a vector.
  Add a clipping for gradients of a vector. Only available for
  coupled velocity components algo (ivelco=1).

- Multigrid: iagmax variable in autmgr.f90/_automatic_aggregation should
  be reset to for each coarsening, and should thus be an internal variable,
  not an argument.
  This leads to more regular aggregation patterns, though often 1 or 2
  more grid levels.

- Make the Cocg matrix for the vectorial iterative gradient DIMENSIONLESS.

- Update handling of periodicity of rotation for more consistent usage
  of halo synchronization. The parcom and percom routines are now
  fully replaced by the halo syn*** series of routines.
  At this stage, the effective operations are unchanged, although
  the API should be clearer.

- Lagrangian module: Implementation of a simpler way to calculate the
  determinant in the trajectography sub-module (default choice from now on)

- Make the Coriolis source term partially implicit with the coupled solver.

- Make free-surface flow independant of the axis direction but dependant of the
  gravity direction.  A checking is added in vericl.

- Add initial version of renumbering for hypbrid parallelism using OpenMP.
  The renumbering algorithm for interior faces is based upon one of the simpler
  algorithms provided by IBM, and is mainly destined for tests, as better
  performing algorithms will be added as a second step.
  The renumbering for boundary faces is simpler and is not based on the
  IBM library.

- New option added: dynamic relaxp in resopv (swpdyn = 1).
  Only available with ivelco = 1.

- Add the Rij EBRSM model (iturb = 32)

- CS-CS coupling avalaible in ivelco=1.

- Remove the multigrid algorithm by the default for eletric variables.

- Use reserved name for temperature or enthalpy field.

- Interleaved ALE displacement array for ivelco=1.

- Log SYRTHES 4 overheads timing information.

- Improve the temperature calculation  in the Libby-Williams model when
  the model is not adiabatic by using the transported enthalpy.

- Improve SALOME CFDSTUDY module

- Ensure Laplacian used to compute wall distance is positive.
  A second pass without reconstruction is run if this is not the case.

- Add matrix dump to file function.

- Add a free surface boundary conditions for the mesh velocity (only
  available with the coupled framework ivelco = 1). Also add a utilitarian
  function to find the closest node of a given point.

- Add a coupled framework for the ALE module.

- Make implicit the handling of the transpose gradient of the velocity
  as well as the secondary viscosity within the diffusive flux
  computation routine (bilsc4).

- Add mesh smoother with initial unwarping algorithm.

- Many minor changes (see ChangeLog for details).

User and pre/post processing changes:

- Allow choice of ADF or HDF5 format for CGNS.

- Add postprocessing options for Syrthes in runcase_coupling.

- Move all user files in a single dedicated directory named src/user.

- Empty reference boundary conditions and initialization user subroutines.
  Examples are now given separately in a SRC/EXAMPLES case subdirectory.
  - Replace all the specific boundary conditions routines by a single
    cs_user_boundary_conditions one.
  - Replace all the specific initialization routines by a single
    cs_user_initialization one.

- Merge all the option initializations for the specific physics
  into the existing usini1.f90, including usppmo.f90

- Merge mesh-related user files in a single cs_user_mesh.c

- Merge nearly all source terms definition in a single cs_user_source_terms.f90

- The mechanism for advanced selection or modification of postprocessing
  meshes is now based on user-defined selection functions.

- Remove never-used user function for multigrid coarsening.

Architectural changes:

- Remove support for CGNS versions < 3.1

- Remove support of PROSTAR/ngeom input format.

- Remove MED 2.3 support.

- Added function attributes to bft_printf() and bft_error() with GNU or Intel
  compilers so as to check format arguments, and fixed all errors and
  warnings subsequently reported.

- Rewrite multigrid autmgr.f90 and crstgr.f90 in C (respectively
  _automatic_aggregation and _build_coarse_lvl in cs_grid.c).

- Rewrite gradient reconstruction in C, with interleaving and OpenMP loops.

- Add typedefs for multidimensional arrays. This allows for much
  clearer syntax for interleaved multidimensional arrays with
  fixed "local" dimensions. Examples of their use may be found
  in gradient computation functions.

- Force link with C++ when using MEDCoupling or PARAMEDMEM.

- Remove obsolete/unused parallel call counters in Fortran wrappers.

- Move additional quantities such as COCG, II', and JJ' from Fortran to
  cs_mesh_quantities structure.

- Added field (and mesh location) API for both C and Fortran
  Field maintain their own metadata, including keywords and boundary condition
  coefficients information.

- Add initial MEDCoupling output plugin writer.

- Add plugin mechanism when dlopen/dlsym/dlclose are available.

- Add a variable related to the Fortran modules directory so as to enable
  NEPTUNE_CFD to compile Code_Saturne Fortran files.

- ChangeLog is now now automatically generated either by running
  'make changelog' on a Subversion checkout, or when
  distributing an archive by running 'make dist'.

- Separation of postprocessing into a common part (usable by other codes)
  and specific Code_Saturne additional default outputs.

- Separate extended neighborhood management from LES filter.

- Maintain vertices interface set in mesh structure for use by mesh
  modification or ALE.

- Replace fvm_interface_t with cs_interface_t structure.
  The API is slightly modified, as cs_interface_t uses 0-based ids for
  local elements, uses a send order to fix bugs with the previous
  API when periodic local and matching elements could not both be sorted by
  increasing global number, and includes some utility functions,
  for copies and sums.

- Use cs_lnum_t instead of cs_int_t for mesh structure.

- Merge system info from bft_sys_info.* and cs_base.* to cs_sytem_info.*.

- Rename cs_perio_* to cs_halo_perio_* (prepares move of fvm_perio_*
  to cs_perio_*).

- Migrate some fvm_ functions and types to the cs_ name prefix.
  fvm_order_* is also moved to cs_order_*.

- New LaTeX classes are added.
  A new plan is added in the theory guide.
  A note on basic rules to write the theory guide is added in the
  developer guide.

- Improve documentation tools detection and add a doc rule at the top level
  directory, and rename the doc directory to docs to avoid Makefile
  incoherency). Also add an update-po rule at the top level directory.

- Add a --with-salome configure option that takes care of SALOME modules
  detection avoiding the use of the different --with-salome-xxx.

- Remove tests using diagonal matrix BLAS 2 from benchmark mode,
  as initial tests show they are slower than specific code, while
  they would require more costly storage than a blocked diagonal.

- Improve handling of external BLAS.
  Either wrappers with the "cs_" rather than "cblas_d" prefix are used
  (for daxpy and ddot with strides of 1), or those functions are
  provided if no external BLAS is available.
  Support for Fortran BLAS is dropped, except for IBM ESSL and AMD ACML
  (which is added), as most other modern BLAS provide C API's.

- Add logging API, which is used for improved multigrid performance information.

Bug fixes:

- Many bug fixes (see ChangeLog for details).


Release 2.1.0 (2011-10-20)
--------------------------

Changes:

- Add thermochemistry reference files for the unified combustion modeling.

- Add coupling with SYRTHES 4.0
  This includes boundary and volume coupling features.
  A conservativity flag ay also be used to force energy conservation.

- Set the relaxation coefficient to the classical value of 0.7 instead of
  0.9 for the steady algorithm.

- Make the thermochemistry file name coherent with the scripts.

- New turbulence model in the framework of the v2-f models.
  This is a "blended" v2f (iturb = 51).

- Add a prototype of solving the velocity components in a coupled way.
  It can be tested by setting the ivelco variable to 1 but is not yet
  compatible with every feature (use with caution...).

- Add an experimental combustion model for coal and heavy fuel oil
  with a unified gas-combustion modelling.
  This new modelling replaces the previous fuel modelling (which is
  still available as a backup).

- Add thin wall insertion.

- Implement the correct behavior to handle thermal wall-function in
  the framework of the scalable wall-functions.

- First implementation of a V&V automatization tool.

- Non-symmetric matrix coefficients are now interleaved.
  (the Fortran API still allows non-interleaved matrixes as an option).

- Remove synchronization of postprocessing with code_aster.

- Add a particle deposition model to the Lagrangian module.

- Prefer initializing the random-number generator with a seed equal
  to 1 instead of 0 (as advised in the zufall.f90 comments).

- Delay discarding of isolated faces to just prior calculation phase.
  This allows post-processing of those faces, as well as their usage
  in mesh joining (not updated for this yet) or user modification.

- Add possibility of defining a default destination rank for
  block to part distributor based on strided adjacency, so as to
  also distribute elements with no adjacency.

- Added choice of Hilbert SFC for domain partitioning.

- Added GMRES to linear solvers.

- Enable build of Doxygen documentation (currently minimal).

- Finalize Code_Saturne/Code_Aster coupling s now known to work
  with SALOME 5.1.5 and Code_Aster NEW10 recent enough on Calibre 5.

- Add a one-equation turbulence  model: the Spalart-Allmaras model.

- Added single-reduction variant of preconditioned conjugate gradient solver.

- Optional merge of coarse grids across processors so as to allow coarsening
  beyond domain boundaries.

- Remove the MATISSE module.
  The MATISSE module can be found in its last version by the Subversion
  tag pre_removing_matisse.

- Fix a mesh-numbering dependancy (and number-of-processors dependancy) in the
  choice of the face at which one takes the reference presure, when the user
  defined a free outlet boundary condition.
  Results may differ from previous calculations but should be not much wrong
  than before and will be less dependant on other parameters.

- Log number of cycles instead of equivalent iterations for multigrid.

- Handle multiple batch systems through the GUI depending on the configuration file.

- Add SALOME module generation, for future Code_Saturne/Code_Aster coupling.

- Move the mass-flux update in case of a rotating mesh, only relevant when
  disabling the pressure reconstruction (no impact on current simulations).

- Remove rarely used interior faces selection for mesh checking.

- Many minor changes (see ChangeLog for details).

User and pre/post processing changes:

- Many updates and additions to GUI.

- Switch to new case directory hierarchy, allowing a new coupling
  directory structure using one subdirectory per domain.
  A domain's name matches its subdirectory, and a new
  runcase_coupling may be generated by cs_create when multiple domains
  are present.
  File and directories are not renamed by the script anymore, so
  copying of results from the execution directory is better automated
  and the latter may be purged when no error has occured. Better,
  the case may be run directly in the results directory.
  The user does not need to specify files to retrieve, but may specify
  scratch files not to retrieve (which should be much rarer).

- Add a master script Python module to call the different modules.

- Scripts and GUI overhauled, with user Python functions and script XML reader.
  Restart behavior is now based on the presence of a restart sub-directory in
  the execution directory. A "--preprocess" option has also been added to the
  Kernel so as to handle purely preprocessing runs.

- Remove inconsistant 'check_mesh' command from the main script (superseded
  by both the graphical interface and the 'run' command). Also remove
  the 'check_consistency' command which was a pale copy of the verification
  stage (--quality option of the solver).

- Allow selection of an alternate build for compute tasks.
  This is useful mainly for supercomputers with separate front-ends, such
  as IBM Blue Gene machines.
  Also allow disabling build of front-end tools.

- Try to use the same shell as the user's current shell in generated
  scripts used by the "code_saturne run" command.

- Add a SALOME command to the main Code_Saturne script so as to launch the
  SALOME platform with the CFDSTUDY module enabled.

- Remove the now useless "code_saturne plot_probes" command since the probes
  now have only one first column for either the time-step number of the time
  value (it was only a wrapper aroun xmgrace -nxy $file).

- Replace the drag'n drop feature in the graphical interface by a more
  intuitive add/remove mechanism for the profile and time average pages.

- Allow reading of CGNS files with BC's referring non-existant edges, as a
  workaround for a bug in ICEM CFD 13 CGNS export.

- Use SCOTCH instead of METIS by default when partitioning.

- Update ParMETIS support to version 4.0.

- Add the possibility for the user to define advanced post-processing,
  writers and meshes, from the GUI.

- ficstp.MOD renamed to ficstp_updated to avoid confusion with Fortran modules.

- Use the "executable" mode of YACS instead of the "loader" one
  for both Code_Saturne and Code_Aster.

- Separate mesh joining verbosity and visualization level. Above verbosity
  level 2, output is directed to a processor-local log file.

- Move activation of Cooling Towers postprocessing to that module.

- Make activation of SYRTHES postprocessing coupling-instance specific.

- MED output: a field must have different names on different MED
  meshes within a single file, so postpend mesh name if necessary.

- Add the possibility of setting up a frequency output in seconds in the GUI.

- Remove all the occurences of the user and developer arrays. User should now
  use the user_modules.f90 suggested in the user files, and developers should
  use allocatable arrays in one of the modules (see src/Makefile.am).

- Transform user files for joining, periodicity, mesh modification,
  SYRTHES and Code_Saturne coupling from Fortran to C.

- Preprocessor: merging of coincident vertices or multiply-referenced
  vertices in face definitions moved inside generation of descending
  connectivity so as to handle cases where quadrangle faces transformed
  into triangles by vertex merging may be identical with faces already
  defined as triangles.

- Preprocessor: extracted meshes should only contain connected vertices.

- Preprocessor: removed postprocessing format options
  and additional minor source code cleanup.

- Preprocessor only creates visualization files if necessary.

- Removed preprocessor meta-file and command-line input file support.

- Preprocessor now converts all colors to groups on mesh import.

- Add output of default volume and boundary mesh groups to EnSight format.
  This requires the usage of additional cases, so as to avoid time dependency
  issues and group name clashes between surface and volume parts.
  Removed group output by Preprocessor for EnSight format now it is in Kernel.

- Improve the cs_solver options management so that the error message is clearer
  in case the code is not compiled with X feature support.

- For code coupling, remove all occurences of the application number and
  replace it by the, application name, which is based on the matching
  computation domain case's direcory name.

- The monitoring files are now created in a specific directory.

- The checkpoint files are now created in a specific directory. Likewise, the
  restart files are read from a specific directory.
  Moreover, the filenames cannot be changed by the user anymore (it was
  already breaking the scripts anyway).

- Add bash completions for Code_Saturne executables.

- GUI: handle batch system based on batch template and code_saturne.cfg.

- Add optional code_saturne.cfg configuration file:
    - select batch type and default batch card.
    - handle the temporary directory from the config file.
    - handle the notion of mesh database directory.

- Change the default behavior for copying the result files (without suffix).

- Make the particles output compliant with EnSight 6 format and improve the
  robustness of the output for visualizing with ParaView.

- Add diff capabilities to the I/O dump utility.

- Move the warped-face-cutting setup from the command line option to a Fortran
  user routine and XML file.

- Switch the joining and periodicity support from the runcase script to the
  XML reader when using the graphical interface.

- Switch to the Python runcase script; the GUI now handles this script,
  and the old (shell) script is removed.

- Add timing information for mesh input I/O.

Architectural changes:

- Add a developers guide with coding standards.

- Add an installation documentation.

- Enable build with MED 3.0 (based on 2.9) as well as with MED 2.3.

- Enable build with CGNS 3.1 as well as CGNS 2.5.

- Print additional MPI environment info for Blue Gene/P.

- Use specific intracommunicator instead of MPI_COMM_WORLD for coupling with
  SYRTHES 3, to avoid bug due to non-contiguous MPMD rank-numbering on
  BLue Gene/P or similar machines.

- Many changes to matrix API:
  - Move matrix type definitions to cs_matrix_priv.h,
    so as to allow implementing matrix operations over several files
    (adding cs_matrix_util.* here).
  - Assign coefficients to matrices outside linear solvers.
    This and allows for simpler and more consitent calls to linear solvers.
    An added benefit is the fact that the overhead of assigning coefficients
    is applied once per time step instead of once per multigrid cycle.
  - Remove alpha.A.x + beta.y operation type for matrixes.

- Added Matrix tuning infrastructure and updated benchmark mode.
  Benchmark mode now uses a fixed minimum time rather than number of
  passes, so as to be less dependant on mesh size.
  Matrixes now include an MSR structure, (modified CSR, with separate diagonal),
  and new operators are added.

- Add OpenMP directives for linear solvers.

- Express Jacobi solvers using y <- A.x rather than y <- alpha.A.x - beta.y
  formulations, and take advantage of the new (A-D).x operator to
  remove extra matrix for polynomial preconditioning.

- Remove sockets support for SYRTHES3 coupling.

- Remove support for IRIX and Tru64 Unix, as these systems are now obsolete.

- Merge bft_file.* and fvm_file.* into cs_file.*, replacing bft_file_printf()
  by simple fprintf(), and making endian-swapping code local where needed.

- Replace cs_bool_t with bool.

- Merge bft_config_defs.h, fvm_config_defs.h, and some definitions from
  cs_base.h into a single cs_defs.h header. This also includes preparations
  for replacing fvm_*num_t by cs_*num_t.

- Un-version TODO file.

- Add --with-modules configure option to allow override of automatic detection.

- Update METIS support for METIS 5.0.

- Upgrade post-processing management API for better mapping with the
  GUI, and for better clarity of user options. All related settings
  are now local to cs_post.c, and are removed from entsor.f90.

- Add a relocatable option to the configure script to handle relocatable
  installation (of executables). It is based on the $ORIGIN variable for
  RPATHs that should be understood by at least GNU and SOLARIS linkers.

- Split the base directory into 4 directories with 3 new ones dedicated
  to the mesh handling, the turbulence modelling, and the Finite Volume
  schemes and resolution.

- Major memory management update: a dynamic allocation when
  necessary, remove useless work arrays and move work arrays where
  there are really used.
  Remove the allocation of main arrays in the C part and the
  corresponding user interface (both graphical and Fortran).
  The only remaining bits of the macro arrays are the ones related to
  the post-processing management (aka ipp2ra): rtp, rtpa, propce, dt and
  tpucou.

- Split the gradient API into subroutines to highlight the pressure
  gradient computation and facilitate future simplifications and tests.
  Only the new 'grdpot' subroutine is now able to account for body forces.

- Remove explicit support of n phases (support was never used, and is still
  possible using adequate variables and defining relations between them).
  Phase index is removed from all arrays, improving code readability and
  avoid bugs in which iphas was not re-set to 1 after a previous loop on
  phases.

- Add a developers guide with coding standards.

- Remove the vectorization directives (when forced) so as to simplify
  the code maintainance and since we do not have  access to a vector
  computer at the moment. Thus, the ivecti/b indicators are now useless.

- Remove matrix structure symmetry flag, replacing if by using a specific
  symmetric matrix type where applicable (added symmetric CSR type).

- Move the calculation of geometric quantities (useful for handling the
  non-orthogonalities) from the Fortran part to the C part.
  The arrays are now accessed directly instead of through a pointer to ra.

- Make use of the "package" concept in the graphical interface to have an
  easier management of both Code_Saturne and NEPTUNE_CFD codes.

- Added detection and handling of environment modules. The configuration
  detected at configure/build time are saved and applied when compiling
  user subroutines or when running a calculation script. A --with-modules
  configure option allows overriding of automatic detection.

- Add the CFD Proxy library to Code_Saturne repository.
  This library is used when Code_Saturne must be wrapped by YACS (SALOME
  component) as a shared library, providing specific Calcium wrappers, for
  example in the Code_Aster / Code_Saturne coupling case. This library
  should be removed when we switch to the executable mode of YACS.

- Use C++ wrapper to link with MED library if necessary.
  This avoids requiring the definition of MED dependency
  libraries, which could interfere with libtool in cases
  mixing static and dynamic libraries and cross-compilation.

- Added C API to handle probes.
  This allows not writing to a temporary file anymore and enabling both
  the legacy .dat (text block data/XmGrace) format and CSV files.

- Add members to the cs_vars_t structure to account for the number of phases,
  in order for NEPTUNE_CFD to be able to link against Code_Saturne source code.

- Introduce a package-specific Python module to avoid too many overloading of
  Python classes by NEPTUNE_CFD. It replaces the cs_config.dirs class and adds
  some useful specific information like the package name, version, ...
  Some work remains to be done in the graphical user interface code, because
  the cs_package module is directly used instead of being passed on.

- C++ compiler detection for a possible link with this compiler when linking
  with MED support and to ease NEPTUNE_CFD configury.

- Move the bison/yacc compilation from the bootstrap stage to the
  compilation stage. Also generate the corresponding C files at dist stage.

- Generate the documenation at dist stage.
  Install the documentation with install target if present.

- Added selection functions for faces at cell criteria boundaries
  and for pre-selection of families.

- Attributes are removed from group classes.:
  Attributes read in the input mesh are converted to group names so as
  to be renamable by the user and listed in the meshe's groups, but
  additional mesh family items may be built when group names are convertible
  to integers so that user subroutines looping on families and family
  properties still work as before.

  Also, selectors still contain internal attributes so as to ensure a
  range[] of integers works as expected, but these are now built directly
  from group names which are convertible to integers.

- EnSight writer now enables parallel IO in binary mode.

- Output information on groups and periodic faces in solver log file,
  and add selection functions to list periodic faces.

- Enable saving of meshes in preprocessor data output whenever
  the input mesh has been modified or meshes have been concatenated.

- Add preprocessor_output files concatenation to the kernel (works in parallel).
  Appending multiple meshes is now done through a user function:
  coordinate transformations and group renames are now possible.
  Mesh concatenation support is removed from Preprocessor.

- Transform all Fortran common blocks to modules.
  Fold dimfbr into dimens and vector into parall.

- Move all Fortran include files from include/* to src/* and rename
  their extension to .f90 in preparation for switch to modules.
  Rename vortex.h to vorinc.f90 to avoid conflict with existing vortex.f90.

- Instrumentation of the gradient computation. Also add a first
  implementation of a pure C interface fo the gradients.

- Add a new API for halo synchronization. Four new functions are available for
  scalars, vectors, diagonal tensor and 9x9 tensor. The old API
  "parcom/percom" remains available.

- Make the CCM files reader able to read both 2.6.1 and later versions.

- Move the preprocessor into the Code_Saturne tree.
  Major code refactoring to ease the maintenance.

- Rename the package from ncs to code-saturne instead.

- Change the behavior of the configure option for finding a Python interpreter.

- Migrate the BFT (Basic Functions and Tools) library to the Code_Saturne tree.

- Migrate the MEI (Mathematical Expression Interpreter) library to the
  Code_Saturne tree.

- Migrate the FVM (Finite Volume Mesh) library to the Code_Saturne tree.

- Replace the FVM coupling library by the new PLE one, which may be installed
  either as part of code_saturne, or separately (Code_Saturne may use either
  its internal PLE library or an external build).

- Replace the SWIG dependency for the syntax checking in the graphical
  interface by a small executable called at each check.

- Add the periodicity joining handling to the solver's parallel joining
  algorithm.

- Move the partitioning and I/O dump tools from the preprocessor to the kernel
  package. Add ParMETIS and PT-SCOTCH/SCOTCH support. By default, serial
  partitionning is done, but parallel partitioning may be selected with
  the user script.

Bug fixes:

- Many bug fixes (see ChangeLog for details).


Release 2.0.0 (2010-08-20)
--------------------------

Code_Saturne 2.0 is a major (fully validated) production release, replacing
version 1.3 (which is to be maintained until the release of version 3.0).

- Correct value for the diffusion coefficient in k-omega SST.

- Move the mass-flux update in case of a rotating mesh, only relevant when
  disabling the pressure reconstruction (no impact on current simulations).

- Translate variable and property names in the example Xml files.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 2.0.0-rc2 (2010-06-22)
------------------------------

- By-pass a possible bug with the flush Fortran 2003 statement on BG systems.

- Make the particles output compliant with EnSight 6 format and improve the
  robustness of the output for visualizing with ParaView.

- Revert a changeset regarding the turbulence relaxation (k/epsilon/omega) to
  avoid stability issues in some calculations. This still needs to be
  investigated in order to have a correct unsteadyness when needed.

- Under HP-UX, handle both hppa and Itanium systems.

- Miscellaneous improvements for the parallel joining algorithm,
  added the periodicity "joining" handling in parallel by the solver.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 2.0.0-rc1 (2010-02-19)
------------------------------

- English translation of more Fortran subroutines headers.

- Move upwards the definition of additional scalars in the GUI because the
  physical properties can depend on them.

- Change the behavior of the xml file relatively to the time step min/max
  factor so that it corrresponds to what is given by the user in the gui.

- Add rotor/stator interaction (specific interpolation) for relative frame
  calculation.

- Add new interpolation schemes to code/code coupling for
  rotor/stator interaction modelling.

- Improve robustness for wall boundary conditions for atmospheric modelling.

- Use new discovery mechanism for coupling applications.

- Remove hyphens and dots from variable names defined by the GUI.

- Rename the main script in code_saturne to avoid a conflict
  with the cs executable of csound package.

- Set the relaxation factor between P0/P1 interpolations to 0.95 in multigrid
  algorithm in order to stabilize the algorithm in some calculations.

- Set the aggregation limit to 3 instead of 8 in the multigrid solver
  for better robustness (but slightly lower performance).

- Upgrade to FVM 0.15

- The space-filling-curve algorithm for backup domain splitting is activated
  by default.

- Add a mpi_io option to the command line so that one can use FVM compiled
  with MPI/IO support and decide at runtime which type of I/O to use.

- Merge the preprocessor user guide into Code_Saturne user guide.

- Add a new Python runcase script with more code coupling capabilities
  and better handling of MPI and batch environments. This scipt is enabled by
  the --new-runcase option to the "code_saturne create" command.
  The GUI still only handles the legacy (shell) script.

- Disable by default the multigrid algorithm for the potential vector
  in MHD as it does not seem to work correctly.

- Many GUI improvements and fixes.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 2.0.0-beta1 (2009-07-29)
--------------------------------

- Move loop on velocity/pressure system so that it includes
  the boundary conditions calculation.

- Use selection mechanism for exchange zone definition.

- Improve cooling tower module.

- Information on SYRTHES installation and compilers can now be
  given after Code_Saturne installation.

- Fix a wrong assumed behavior on x86_64 and IA64 computers
  in the particles tracking algorithm. One assumed that the internal
  FP precision was always 80 bits (true for FPU x87 coprocessor on
  x86 processors but wrong with SSE optimizations)...
  This is a short term fix (assembler + compiler-dependent option).

- Add Coriolis source terms in standard calculations.

- Do not relax turbulent variables (k/epsilon/omega) in unsteady
  simulations.

- Adapt 1d-profiles file to SALOME format for title and labels.

- Add radiative transfer support for fuel combustion.

- Add a space-filling-curve algorithm for parallel partitioning.
  This algorithm is based on a Morton curve.

- On Linux systems, use fenv library to trap floating point exceptions.

- Add parall mesh joining feature in the Kernel (by default, the legacy
  joining feature from the Preprocessor is still used.

- Add non-neutral atmosphere modelling (both dry and humid, though
  humid atmosphere modelling is not yet fully functional).

- When detecting a divergence in the linear system solver stage, abort
  and write graphical data of the given matrix (rhs, diagonal, ...).

- Simplify running a SYRTHES 3 / Code_Saturne coupling through sockets.

- Add a mesh category for an easier management of user-defined
  post-processing.

- Add a new turbulence model for LES: the WALE model.

- Fix when testing for radiative transfer.

- Enable the multigrid algorithm for the diffusion equation in the
  electric arcs, Lagrangian, radiative transfer (P1 model) and ALE modules.

- Replace all the different user scripts by a single one.

- The GUI is now based on QT4 rather than Tk, and is now part of the
  kernel source tree to simplify version control and installation.

- The Syrthes coupling library is now included in the kernel source
  tree for easier versioning and installation.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 2.0.0-beta1 (2009-05-26)
--------------------------------

- Add a quick reference card documentation.

- Add man pages for the main scripts to the documentation.

- Do not compile the code if there is no user file.

- Create a graphical post-processing view of the boundary

- New build system, based on the GNU Autotools.

- Conversion of all the Fortran files (headers included)
  form fixed F77 format to free f90 format. A small-case
  extension is used so that one can, in the future, skip the
  preprocessing stage.

- Change the mesh quality verification behavior. The different
  gradient calculation modes are now done in one pass. Thus, only
  one option remains for the command line (-q, with no sub-option).

- The kernel is now able to post-process edges for a given
  mesh via the PSTEDG subroutine. This functionnality is
  thus removed from the preprocessor.

- English translation of many C source files.

- Add cell renumbering capability and re-organize face renumbering.

- SYRTHES coupling upgraded to syr_cs 2.3.0.

- Upgrade checkpoint/restart to allow MPI-IO.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 1.4.0 (2008-11-28)
--------------------------

Code_Saturne 1.4 is an intermediate development release, very similar to
version 1.3.3 except for the items below:

- Add atmospheric flow modelling for neutral atmosphere flows.

- Always set the head loss tensor to be a full symmetric tensor
  in order to ease the setup within the Graphical User Interface
  (NCKPDC has been removed and replaced by the constant 6).

- Make the test on the convergence of the gradient reconstruction
  method dimensionless if the maximum control volume is greater
  than 1 (useful for some atmospheric simulations).

- Remove the forcing to ASCII mode for post-processing when
  the calculation is in verification mode.

- Remove the pipe communication mode for SYRTHES coupling.

- Add oxycombustion and fuel classes management.

- Separate preprocessing from partitioning. There is now a single
  preprocessor_output file, so that mesh joining does not need to be run again
  for a different number of processors.

- Preprocessor output, partitioning, and restart files a now use a common
  binary format (independent of the number of processors), which may be
  written or read using MPI-IO (not enabled for restart files yet).

- Scripts translation and improvements.

- The reconstruction of the right-hand side has now its own precision
  parameter (EPSRSM) instead of using EPSILO.

- Re-write of the algebric multigrid algorithm for the resolution of
  the of pressure equation. The algorithm works in serial and parallel,
  and is now the default, leading to run-times reduced by a factor or
  2 to 3 on average.

- Add wet cooling tower heat exchange specific physics model.

- The size of the work arrays IA and RA is now defined with a Fortran
  user file and no more via the scripts, making them independent of
  the number of sub-domains in case of a parallel run.

- Simplified the command line option for a parallel runs.

- Several bug fixes and minor improvements (see ChangeLog for details).


Release 1.3.3 (2008-11-27)
--------------------------

- Upgrade dependencies to FVM 0.12.0 and BFT 1.0.8.

- Upgrade SYRTHES to 3.4.2 version (syr_cs 2.1.0).  Restarting
  from SYRTHES 3.3 is not possible due to incompatible format.

- Add a QUALITY_ASSURANCE file which precises if the current
  version of code_saturne is validated under EDF quality
  assurance.

- Portability updates.

- Add a patch directory in which the user can find
  some untested patches fixing tricky issues.

- Updating the halos at the beginning of an outer iteration is now
  done on extended halos (reducing error with rotational periodicity).

- Improved mesh and coherency tests.

- Added localization (French or English) to the Code_Saturne Kernel
  and translated all scripts to English.

- Several bug fixes (see ChangeLog for details).


Release 1.3.2 (1008-04-16)
--------------------------

- Port to BlueGene/P and Cray XT

- Update FVM API to take into account some particular cases
  where two periodicities are not commutative.

- Add a warning when changing the mesh vertices coordinates along
  with periodicity, which can break the periodicity parameters.

- Correct distance use at the first time-step when using the k-w SST
  model and the old algorithm to compute the wall distance.

- Change the behavior when trying to get a list of cells or faces
  with the selector. Revert to the old behavior when the selection
  was done through the mesh properties. A warning is issued instead
  of an error when the criteria returns an empty selection.

- Move unmaintained macros files to aux/macros_old.

- Make Code_Saturne stop when a Lagrangian calculation is run with the
  steady-state algorithm.

- Change the name Syrthes to SYRTHES following the trademark registration.

- Add documentation on faces/cells selection with GETxxx via fvm_selector

- Add comments on steady-state algorithm IDTVAR=-1.

- Add a temporary "tmp_Saturne/$ETUDE.$CAS.$DATE" directory when running
  a case in $TMPDIR (when defined by the system) to avoid an issue
  with the linux distribution CAELinux

- Complete reorganisation of the theory documentation to follow the template
  of the user and tutorial documents

- Variant for Intel compiler on Itanium only, optimizations by BULL.

- Rename the Preprocessor listing from listenv to listpre.

- Portability improvements in scripts.

- Several bug fixes (see ChangeLog for details).


Release 1.3.1 (2007-11-28)
--------------------------

- Storage on preprocessor files in a subdirectory of
  the temporary execution directory.

- Added code_Saturne tutorial.

- synchronization with Syrthes only if ITRALE>0 in caltri.F

- creation of erreur_n*** files only for processor 0 and the
  processors receiving a segmentation fault or floating point
  exception message.

- change of kernel options "couleur" and "groupe" to "color"
  and "group" for Syrthes options.

- modification of residue test to detect when conjugate
  gradient diverges.

- Several bug fixes (see ChangeLog for details).


Release 1.3.0 (2007-08-02)
--------------------------

Code_Saturne 1.3 is a major (fully validated) production release.
Release 1.3.0 is the feature-freeze release, and versions 1.3.1
and 1.3.2 are validation bug-fix versions. 1.3.3 is the first production
version of code_saturne 1.3.

This version of the code is also the first version of code_saturne distributed
under a free software (GPL + LGPL) licence.
