Release 2.0.4 (unreleased)
==========================

Bug fixes:
----------

- Fix detection of GCC compiler with some Cray versions.

Changes:
--------

- Remove some definitions for some ancient compilers and MPI libraries.

- Add `ple_locator_exchange_point_var_all` function to handle
  exchanges with both located and unlocated points.

- Extend `ple_coupling_mpi_set` features:
  * Add `ple_coupling_mpi_set_compute_timestep` function to compute
    a recommended time step for the current application based on
    provided flags and values of applications in a set.
  * Add `PLE_COUPLING_TS_FOLLOWER` and `PLE_COUPLING_TS_FOLLOWER`
    synchronization flag defintitions for additional time step
    update schemes.

- Remove libtool from build system.
  This requires slightly more commplex rules in the Makefile.am's, but
  should avoid issues on some systems where lbitool's absolute refusal
  to ignore `.la` files even when they contain incorrect values
  could cause more problems than it solves.

- Add key/value-based settings functions:
  * `ple_locator_set_default_option`
  * `ple_locator_set_options`

Release 2.0.3 (November 8, 2021)
================================

Changes:
--------

- Add communication ordering to reduce serialization in parallel.
  Algorithm versioning ensures this is not used when combined
  with an older PLE library version.

Bug fixes:
----------

- Fix PLE locator variable exchange in asynchronous reverse mode.

- Really fix PLE locator for local cases when points fall outside domain.

- Avoid crash in `ple_locator_shift_location` for empty locator.

Release 2.0.2 (March 16, 2018)
==============================

Bug fixes:
----------

- Fix in `ple_locator_extend_search` for non-local cases.

Changes:
--------

- Use asyncronous locator variable exchange by default when
  number of communicating ranks is not too large (< 128).

- Implement handling of point tags so as to ignore location
  on elements sharing tag (also requires matching features
  in `ple_mesh_elements_locate_t` function used.

Architectural changes
---------------------

- Add support for internationalization.

Release 2.0.1 (May 21, 2015)
============================

Bug fixes:
----------

- Fix in PLE locator indexing for `ple_locator_extend_search`.

Changes:
--------

- Improve Cray support.

Release 2.0.0 (April 30, 2015)
==============================

Changes:
--------

- Upgrade PLE API so as to allow multiple searches. This replaces the use of
  functions to locate closest elements, as using successive searches with
  increasing tolerance requires no specific function and is expected to have
  more stable behavior.

Release 1.0.4 (August 18, 2014)
===============================

Bug fixes:
----------

- Fix asynchronous exchange for `ple_locator_exchange_point_var`.

Changes:
--------

- Add Clang compiler detection.

- Updates in build scripts and MPI detection.

Release 1.0.3 (April 03, 2013)
==============================

Bug fixes:
----------

- Fix a test for "locate on closest" and another for empty realloc.

Release 1.0.2 (October 29, 2012)
================================

Bug fixes:
----------

- Fix extra count (should be 3, not 4) in PLE flags exchange in `ple_locator.c`.

- Fixes for mixed "locate on closest / do not locate on closest" cases.

Changes:
--------

- Avoid deprecated MPI datatypes to prepare for MPI 3.

- Improve robustness of MPI detection tests.

Release 1.0.1 (June 26, 2012)
=============================

Bug fixes:
----------

- Fix link problem in PLE unit tests.

- Fix MPI autodetection on FEDORA systems.

Changes:
--------

- Use same FLAGS as parent Code_Saturne build.

- Remove support for IRIX and Tru64 Unix, as these systems are
  now obsolete.

Release 1.0.0 (July 28, 2011)
=============================

- Initial release of the "Parallel Location end Exchange" library,
  based on the separation of the related functionality from the
  "Finite Volume Mesh" library.
